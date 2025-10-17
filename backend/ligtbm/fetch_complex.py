import click
import pymol
import pymol2
import Bio.SeqIO


class PyMolError(Exception):
    pass


def fetch_complex(pdbid, chain, chemid, name, seq=None, alt_loc='A', contact_dist=8, pymol_instance=None):
    if pymol_instance:
        cmd = pymol_instance.cmd
    else:
        # With PyMOL 2.1, calling any pymol.cmd function will automatically start a backend process without the GUI in the main thread.
        cmd = pymol.cmd

    cmd.fetch(pdbid, 'orig')
    # TODO Split all states and alternative locations into separate pymol objects
    cmd.remove('orig and not (alt ''+{})'.format(alt_loc))
    cmd.alter('orig', "alt=''")

    cmd.show_as('sticks', 'orig and organic')
    cmd.show_as('cartoon', 'orig and polymer.protein')
    cmd.show_as('nb_spheres', 'orig and (solvent or inorganic)')

    if seq:
        rec_sele = 'orig and polymer.protein and chain {} and pepseq {} and present'.format(chain, seq)
        rec_check = cmd.select('rec', rec_sele)
        if rec_check == 0:
            cmd.delete('orig rec')
            raise PyMolError('cannot select {} sequence in {}_{}'.format(seq, pdbid, chain))
            # cmd.select('rec', 'orig and polymer.protein and chain {}'.format(chain))
            # chain_seq = ''.join(cmd.get_fastastr('rec').splitlines()[1:])  # There are still some problems with the HEADER
            # match = difflib.SequenceMatcher(None, seq, chain_seq).find_longest_match(0, len(seq), 0, len(chain_seq))
            # match_seq = seq[match.a: match.a + match.size]
            # rec_sele = 'orig and polymer.protein and chain {} and pepseq {} and present'.format(chain, match_seq)
            # rec_check = cmd.select('rec', rec_sele)
        else:
            match_seq = ''.join(cmd.get_fastastr('rec').splitlines()[1:])  # There are still some problems with the HEADER
    else:
        match_seq = None  # There are still some problems with the HEADER
        rec_check = cmd.select('rec', 'orig and polymer.protein and chain {} and present'.format(chain))
        if rec_check == 0:
            cmd.delete('orig rec')
            raise PyMolError('cannot select {}_{}'.format(pdbid, chain))

    lig_check = cmd.select('lig', '(orig and organic and resn {}) within {} of rec'.format(chemid, contact_dist))
    if lig_check == 0:
        cmd.delete('orig rec lig')
        if match_seq:
            raise PyMolError('Cannot find {} ligand within {} of {} in {}_{}'.format(chemid, contact_dist, match_seq, pdbid, chain))
        else:
            raise PyMolError('Cannot find {} ligand within {} of {}_{}'.format(chemid, contact_dist, pdbid, chain))
    space = {'lig_atoms': []}
    cmd.iterate('lig', 'lig_atoms.append((chain, resi, resn))', space=space)
    lig_residues = list(set(space['lig_atoms']))
    if len(lig_residues) == 1:
        lig_res = lig_residues[0]
    else:
        pocket_sele = 'rec within {} of (orig and organic and chain {} and resi {} and resn {})'
        options = [(cmd.select('tmp', pocket_sele.format(contact_dist, *lig_res)), lig_res) for lig_res in lig_residues]
        cmd.delete('tmp')
        lig_res = max(options, key=lambda x: x[0])[1]
    cmd.select('lig', 'orig and organic and chain {} and resi {} and resn {}'.format(*lig_res))

    env_check = cmd.select('env', 'byres ((orig and (solvent or inorganic or organic)) near_to {} of lig)'.format(contact_dist))

    cmd.create('{}.rec'.format(name), 'rec')
    cmd.create('{}.lig'.format(name), 'lig')
    cmd.create('{}.env'.format(name), 'env')

    cmd.delete('orig rec lig env')
    return match_seq


@click.command()
@click.argument('pdbid', type=str)
@click.argument('chain', type=str)
@click.argument('chemid', type=str)
@click.argument('rec_file', type=click.Path())
@click.argument('lig_file', type=click.Path())
@click.argument('env_file', type=click.Path())
@click.option('--seq_fasta_file', type=click.Path(exists=True), default=None, help='FASTA file with receptor sequence.')
@click.option('--contact_dist', type=float, default=8, help='Receptor-ligand contact cut-off distance in Angstroms.')
@click.option('--alt_loc', type=str, default='A', help='Alternate location.')
def _fetch_complex_click(pdbid, chain, chemid, rec_file, lig_file, env_file, seq_fasta_file, contact_dist, alt_loc):
    '''
        PyMol-based script for fetching protein-ligand complexes.\n
        Receptor part of the protein is described by PDBID, CHAIN and AA subsequence (optionally) of the CHAIN written
        in SEQ_FASTA_FILE. Receptor will be saved in REC_FILE.\n
        Ligand part is described by CHEMID. It will be saved in LIG_FILE.\n
        Additionally, the environment of the pocket (solvent, non-ligand organic, and inorganic molecules) will be
        saved in ENV_FILE.\n
        REC_FILE, LIG_FILE, and ENV_FILE can be written in any PyMol-available format (e.g. pdb and mol).
    '''
    if seq_fasta_file:
        seq_record = next(Bio.SeqIO.parse(seq_fasta_file, 'fasta'))
        seq = str(seq_record.seq)
    else:
        seq = None
    with pymol2.PyMOL() as pm:
        fetch_complex(pdbid, chain, chemid, 'abc', seq, alt_loc, contact_dist, pm)
        pm.cmd.save(rec_file, 'abc.rec')
        pm.cmd.save(lig_file, 'abc.lig')
        pm.cmd.save(env_file, 'abc.env')


if __name__ == '__main__':
    _fetch_complex_click()
