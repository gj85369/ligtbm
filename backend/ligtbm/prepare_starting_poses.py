import os
import click
import json
import multiprocessing as mp
import shutil

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
CHEMID_TO_LIGINFO_FILE = os.path.join(FILE_DIR, 'data', 'Ligand_Expo', 'chemid_to_liginfo.json')
PDBID_TO_CHEMID = os.path.join(FILE_DIR, 'data', 'Ligand_Expo', 'pdbid_to_chemids.json')

import logging
import logging.config
logging.config.fileConfig(os.path.join(FILE_DIR, 'logging.conf'), disable_existing_loggers=False)
user_logger = logging.getLogger('user')
dev_logger = logging.getLogger('dev')

import pymol2
from rdkit import Chem
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.Data.SCOPData import protein_letters_3to1

import tools
from BLAST_tools import run_blastp, BLASTError
from gen_opt_confs import gen_opt_confs
from collect_ligs_prepare_mcs import collect_ligs_prepare_mcs
from fetch_complex import fetch_complex, PyMolError
from CS_confs_align import min_cs_rmsd_conf_align

try:
    from modeller_tools import build_model, ModellerError
except ImportError:
    have_modeller = False
else:
    have_modeller = True


def get_seq(structure_file, chain):
    max_ratio_of_bad_residues_warning = 0.1  # 10%, chosen arbitrary
    min_num_of_good_residues_warning = 15  # chosen arbitrary
    with pymol2.PyMOL() as p:
        p.cmd.load(structure_file, 'orig')
        num_receptor_atoms = p.cmd.select('rec', 'orig and polymer.protein and chain {}'.format(chain))
        if num_receptor_atoms == 0:
            user_logger.error('Cannot find protein chain {} in the receptor'.format(chain))
            return None
        seq = ''.join(p.cmd.get_fastastr('rec').splitlines()[1:])
        num_nonstandard = seq.count('?')
        if num_nonstandard > max_ratio_of_bad_residues_warning * len(seq):
            user_logger.warning('More than {:.1f}% of the receptor residues are non-standard',
                                max_ratio_of_bad_residues_warning * 100)
        if '?' in seq:  # Try to use extended alphabet containing non-standard residues from SCOP
            space = {'residues': []}
            p.cmd.iterate('rec', 'residues.append((resi, resn))', space=space)
            residues = []
            for resi_resn in space['residues']:
                if resi_resn not in residues:
                    residues.append(resi_resn)
            for resi, resn in residues:
                if resn not in protein_letters_3to1.keys():
                    user_logger.warning('Unknown residue {} in the receptor'.format(resn))
            seq = ''.join([protein_letters_3to1.get(resn, 'X') for resi, resn in residues])
    seq_known = seq.replace('X', '')
    if len(seq_known) == 0:
        user_logger.error('Unable to find any known residue in the receptor')
        return None
    if len(seq_known) < (1-max_ratio_of_bad_residues_warning) * len(seq):
        user_logger.warning('More than {:.1f}% of the receptor residues are unknown'.format(
            max_ratio_of_bad_residues_warning * 100))
    if len(seq_known) < min_num_of_good_residues_warning:  # 15 residue threshold chosen arbitrary
        user_logger.warning('Less than {} known residues in the receptor'.format(min_num_of_good_residues_warning))
    return seq


def mp_func(ref_rec_pdbid, cases, args):
    rec_pdb_file = args['rec_pdb_file']
    rec_chain = args['rec_chain']
    rec_seq = args['rec_seq']
    alt_loc = args['alt_loc']
    contact_dist = args['contact_dist']
    lig_mol = args['lig_mol']
    no_remodeling = args['no_remodeling']

    pdbid_msg_exc_list = []
    pdbid_outputs = {}

    with pymol2.PyMOL() as p:
        local_pdb_path = os.getenv('LOCAL_PDB')
        p.cmd.set('fetch_type_default', 'cif')
        if local_pdb_path:
            p.cmd.set('fetch_host', 'file:///')
            p.importing.hostPaths['cif'] = local_pdb_path
        else:
            dev_logger.warning('Local PDB database not provided, will download structures from the web')

        p.cmd.load(rec_pdb_file, 'target_orig')
        p.cmd.select('target_rec', 'target_orig and polymer.protein and chain {}'.format(rec_chain))
        for case in cases:
            ref_rec_chain = case['blast_alignment']['sbjct_chains'][0]  # take only one
            ref_rec_seq = case['blast_alignment']['sbjct']
            ref_lig_chemid = case['refinfo']['ref_chemid']

            try:
                # TODO sele_complex, replace "within dist" exception with None output
                fetch_complex(ref_rec_pdbid, ref_rec_chain, ref_lig_chemid, 'ref', seq=None, alt_loc=alt_loc,
                              contact_dist=contact_dist, pymol_instance=p)
            except PyMolError as exc:
                msg = ''
                pdbid_msg_exc_list.append((msg, exc))
                continue

            p.cmd.create('complex', 'ref.rec or ref.lig')
            p.cmd.delete('ref.rec ref.lig ref.env')
            cealign_result = p.cmd.cealign('target_rec', 'complex')
            cealign_len = cealign_result['alignment_length']
            cealign_rmsd = cealign_result['RMSD']

            # TODO potential problem with several hsps in BLAST
            case_dir = '{}_{}_{}'.format(ref_rec_pdbid, ref_rec_chain, ref_lig_chemid)
            os.makedirs(case_dir, exist_ok=True)

            ref_rec_pdb_file = os.path.join(case_dir, 'ref_rec.pdb')
            ref_lig_mol_file = os.path.join(case_dir, 'ref_lig.mol')
            rec_pdb_file = os.path.join(case_dir, 'rec.pdb')
            p.cmd.save(ref_rec_pdb_file, 'complex and polymer.protein')
            p.cmd.save(ref_lig_mol_file, 'complex and organic')
            p.cmd.save(rec_pdb_file, 'target_rec')
            p.cmd.delete('complex')

            ref_liginfo = case['refinfo']['ref_liginfo']
            try:
                ref_lig_mol = tools.mol_from_mol_file(ref_lig_mol_file, ref_liginfo)
            except tools.RDKitError as exc:
                shutil.rmtree(case_dir)
                msg = ''
                pdbid_msg_exc_list.append((msg, exc))
                continue

            mcs_smarts = case['refinfo']['mcs']['smarts']
            try:
                aligned_lig_mol, mcs_rmsd, mcs_mapping = min_cs_rmsd_conf_align(lig_mol, ref_lig_mol, mcs_smarts)
            except tools.RDKitError as exc:
                shutil.rmtree(case_dir)
                msg = ''
                pdbid_msg_exc_list.append((msg, exc))
                continue
            aligned_lig_mol_file = os.path.join(case_dir, 'lig.mol')
            Chem.MolToMolFile(aligned_lig_mol, aligned_lig_mol_file)

            if no_remodeling:
                rec_pdb_file_to_use = rec_pdb_file
            else:
                modeller_rec_pdb_file = os.path.join(case_dir, 'modeller_rec.pdb')
                try:
                    build_model(modeller_rec_pdb_file, rec_seq, ref_rec_pdb_file, ref_rec_chain)
                except ModellerError as exc:
                    msg = 'MODELLER failed, using user-submitted structure'
                    pdbid_msg_exc_list.append((msg, exc))
                    rec_pdb_file_to_use = rec_pdb_file
                else:
                    rec_pdb_file_to_use = modeller_rec_pdb_file

            pdbid_outputs[case_dir] = {'rec': rec_pdb_file_to_use,
                                       'lig': aligned_lig_mol_file,
                                       'ref_parts': [{'ref_rec': ref_rec_pdb_file,
                                                      'ref_lig': ref_lig_mol_file,
                                                      'blast_alignment': case['blast_alignment'],
                                                      'refinfo': case['refinfo'],
                                                      'rmsd': mcs_rmsd,
                                                      'mapping': [{'target_atom_id': target_atom_id,
                                                                   'ref_atom_id': ref_atom_id
                                                                   } for target_atom_id, ref_atom_id in mcs_mapping]
                                                      }]
                                       }
    return pdbid_outputs, pdbid_msg_exc_list


def prepare_starting_poses(lig_smiles, rec_pdb_file, rec_chain,
                           mcs_coverage_threshold=0.5, fp_tanimoto_threshold=0.4, num_confs=1000,
                           evalue_threshold=1e-20, identity_threshold=0.3, database='pdb',
                           chemid_to_liginfo_file=CHEMID_TO_LIGINFO_FILE, pdbid_to_chemids_file=PDBID_TO_CHEMID,
                           contact_dist=8, alt_loc='A', excl_pdbids=None, top_n_outputs=20, no_remodeling=False,
                           num_threads=1):
    rec_seq = get_seq(rec_pdb_file, rec_chain)
    if rec_seq is None:
        user_logger.error('Unable to extract sequence from the provided PDB')
        return None

    rec_fasta_file = 'rec.fasta'
    seq_record = SeqRecord(Seq(rec_seq, generic_protein))
    # TODO remove tmp fasta file
    SeqIO.write(seq_record, rec_fasta_file, "fasta")

    if not no_remodeling and not have_modeller:
        dev_logger.error('Unable to initialize MODELLER, even though initial tests passed')
        user_logger.error('MODELLER requested, but failed to initialize')
        raise RuntimeError('Unable to initialize MODELLER')

    try:
        pdbid_to_blast_alignments, unparsed_titles = run_blastp(rec_fasta_file, evalue_threshold, identity_threshold, database, num_threads)
    except BLASTError as exc:
        user_logger.error('Unable to run BLAST on the receptor: {}'.format(str(exc)))
        dev_logger.exception('Unable to run BLAST on the receptor')
        return None  # or maybe raise an exception

    if len(unparsed_titles) != 0:
        dev_logger.warning('Unparsed titles: {}'.format(unparsed_titles))
    # TODO temporary checkpoint for debugging
    with open('pdbid_to_blast_alignments.json', 'w') as f:
        json.dump(pdbid_to_blast_alignments, f, indent=4)
    pdbids = list(pdbid_to_blast_alignments.keys())

    if excl_pdbids is not None:
        excl_pdbids = [pdbid.upper() for pdbid in excl_pdbids]
        pdbids = list(set(pdbids) - set(excl_pdbids))

    # Ligand Expo
    with open(chemid_to_liginfo_file) as f:
        chemid_to_liginfo = json.load(f)
    with open(pdbid_to_chemids_file) as f:
        pdbid_to_chemids = json.load(f)

    try:
        lig_mol = tools.safe_mol_init(lig_smiles)
    except tools.RDKitError as exc:
        user_logger.error('Incorrect input SMILES: {}'.format(str(exc)))
        return None  # or maybe raise an exception

    try:
        lig_mol = gen_opt_confs(lig_mol, num_confs, num_threads)
    except tools.RDKitError as exc:
        user_logger.error('Failed generation/optimization of conformers (ligand too big?): {}'.format(str(exc)))
        dev_logger.exception('Failed generation/optimization of conformers')
        return None  # or maybe raise an exception

    pdbid_to_refinfos = collect_ligs_prepare_mcs(lig_mol, chemid_to_liginfo, pdbid_to_chemids, pdbids,
                                                 mcs_coverage_threshold, fp_tanimoto_threshold, num_threads)
    # TODO temporary checkpoint for debugging
    with open('pdbid_to_refinfos.json', 'w') as f:
        json.dump(pdbid_to_refinfos, f, indent=4)
    pdbids = list(pdbid_to_refinfos.keys())

    pdbid_to_cases = {}
    for pdbid in pdbids:
        cases = []
        for refinfo in pdbid_to_refinfos[pdbid]:
            for blast_alignment in pdbid_to_blast_alignments[pdbid]:
                cases.append({'refinfo': refinfo,
                              'blast_alignment': blast_alignment})
        pdbid_to_cases[pdbid] = cases
    # TODO temporary checkpoint for debugging
    with open('pdbid_to_cases.json', 'w') as f:
        json.dump(pdbid_to_cases, f, indent=4)

    identity = lambda case: case['blast_alignment']['identity']
    mcs_coverage = lambda case: case['refinfo']['mcs']['coverage']
    case_score = lambda case: (identity(case) - identity_threshold)**2 +\
                              (mcs_coverage(case) - mcs_coverage_threshold)**2
    pdbid_case_list = [(pdbid, case) for pdbid, cases in pdbid_to_cases.items() for case in cases]
    pdbid_case_list.sort(key=lambda pdbid_case: case_score(pdbid_case[1]), reverse=True)

    pdbid_to_cases = {}
    for pdbid, case in pdbid_case_list[:2 * top_n_outputs]:
        if pdbid not in pdbid_to_cases:
            pdbid_to_cases[pdbid] = []
        pdbid_to_cases[pdbid].append(case)

    args = {'rec_pdb_file': rec_pdb_file,
            'rec_chain': rec_chain,
            'rec_seq': rec_seq,
            'alt_loc': alt_loc,
            'contact_dist': contact_dist,
            'lig_mol': lig_mol,
            'no_remodeling': no_remodeling}

    mp_inputs = [pdbid_cases + (args,) for pdbid_cases in pdbid_to_cases.items()]
    with mp.Pool(num_threads) as pool:
        mp_outputs = pool.starmap(mp_func, mp_inputs)

    case_name_to_output = {}
    for pdbid_outputs, pdbid_msg_exc_list in mp_outputs:
        for msg, exc in pdbid_msg_exc_list:
            dev_logger.warning(msg, exc_info=exc)
        case_name_to_output.update(pdbid_outputs)

    case_name_output_list = list(case_name_to_output.items())
    output_score = lambda output: case_score(output['ref_parts'][0])
    case_name_output_list.sort(key=lambda case_name_output: output_score(case_name_output[1]), reverse=True)
    for case_dir, output in case_name_output_list[top_n_outputs:]:
        shutil.rmtree(case_dir)
    case_name_to_output = {case_name: output for case_name, output in case_name_output_list[:top_n_outputs]}

    with open('case.json', 'w') as f:
        json.dump(case_name_to_output, f, indent=4)


@click.command()
@click.argument('param_json', type=click.Path(exists=True))
@click.option('--num_threads', type=int, default=1, help='Number of threads for parallel processing.')
@click.option('--chemid_to_liginfo_file', type=click.Path(exists=True), default=CHEMID_TO_LIGINFO_FILE,
              help='Ligand Expo json file with chemid-to-ligandinfo dict.')
@click.option('--pdbid_to_chemids_file', type=click.Path(exists=True), default=PDBID_TO_CHEMID,
              help='Ligand Expo json file with pdbid-to-chemids dict.')
@click.option('--database', type=str, default='pdb', help='Database for BLAST search.')
def _prepare_starting_poses_click(param_json, num_threads, chemid_to_liginfo_file, pdbid_to_chemids_file, database):
    with open(param_json) as f:
        param_dict = json.load(f)
    lig_smiles = param_dict['lig_smiles']
    rec_pdb_file = param_dict['rec_pdb_file']
    rec_chain = param_dict['rec_chain']
    no_remodeling = param_dict.get('no_remodeling', False)
    mcs_coverage_threshold = param_dict.get('mcs_coverage_threshold', 0.5)
    fp_tanimoto_threshold = param_dict.get('fp_tanimoto_threshold', 0.4)
    num_confs = param_dict.get('num_confs', 1000)
    evalue_threshold = param_dict.get('evalue_threshold', 1e-20)
    identity_threshold = param_dict.get('identity_threshold', 0.3)
    contact_dist = param_dict.get('contact_dist', 8)
    alt_loc = 'A'  # tmp
    excl_pdbids = param_dict.get('pdb_exc', [])
    top_n_outputs = param_dict.get('top_n_outputs', 20)

    prepare_starting_poses(lig_smiles, rec_pdb_file, rec_chain, mcs_coverage_threshold, fp_tanimoto_threshold,
                           num_confs, evalue_threshold, identity_threshold, database, chemid_to_liginfo_file,
                           pdbid_to_chemids_file, contact_dist, alt_loc, excl_pdbids, top_n_outputs, no_remodeling,
                           num_threads)


if __name__ == '__main__':
    _prepare_starting_poses_click()
