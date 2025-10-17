import multiprocessing as mp
import logging
user_logger = logging.getLogger('user')
dev_logger = logging.getLogger('dev')

from rdkit import DataStructs
from rdkit.Chem import rdFMCS
from rdkit.Chem.Fingerprints import FingerprintMols

import tools


def calc_mcs(mol1, mol2, mcs_flags=[]):
    if 'aa' in mcs_flags:
        atomcompare = rdFMCS.AtomCompare.CompareAny
    elif 'ai' in mcs_flags:
        # CompareIsotopes matches based on the isotope label
        # isotope labels can be used to implement user-defined atom types
        atomcompare = rdFMCS.AtomCompare.CompareIsotopes
    else:
        atomcompare = rdFMCS.AtomCompare.CompareElements

    if 'ba' in mcs_flags:
        bondcompare = rdFMCS.BondCompare.CompareAny
    elif 'be' in mcs_flags:
        bondcompare = rdFMCS.BondCompare.CompareOrderExact
    else:
        bondcompare = rdFMCS.BondCompare.CompareOrder

    if 'v' in mcs_flags:
        matchvalences = True
    else:
        matchvalences = False

    if 'chiral' in mcs_flags:
        matchchiraltag = True
    else:
        matchchiraltag = False

    if 'r' in mcs_flags:
        ringmatchesringonly = True
    else:
        ringmatchesringonly = False

    if 'cr' in mcs_flags:
        completeringsonly = True
    else:
        completeringsonly = False

    maximizebonds = True

    mols = [mol1, mol2]
    try:
        mcs_result = rdFMCS.FindMCS(mols,
                                    atomCompare=atomcompare,
                                    bondCompare=bondcompare,
                                    matchValences=matchvalences,
                                    ringMatchesRingOnly=ringmatchesringonly,
                                    completeRingsOnly=completeringsonly,
                                    matchChiralTag=matchchiraltag,
                                    maximizeBonds=maximizebonds)
    except:
        # sometimes Boost (RDKit uses it) errors occur
        raise tools.RDKitError('MCS calculation was failed')
    if mcs_result.canceled:
        raise tools.RDKitError('MCS calculation ran out of time')
    return mcs_result.smartsString, mcs_result.numAtoms, mcs_result.numBonds


def calc_fp_tanimoto(mol1, mol2):
    fp1 = FingerprintMols.FingerprintMol(mol1)
    fp2 = FingerprintMols.FingerprintMol(mol2)
    fp_tanimoto = DataStructs.FingerprintSimilarity(fp1, fp2)
    return fp_tanimoto


def mp_func(ref_chemid, ref_liginfo, args):
    target_mol = args['target_mol']
    mcs_coverage_threshold = args['mcs_coverage_threshold']
    fp_tanimoto_threshold = args['fp_tanimoto_threshold']

    output = {'ref_chemid': ref_chemid,
              'ref_liginfo': ref_liginfo}

    try:
        ref_mol = tools.safe_mol_init(ref_liginfo)
    except tools.RDKitError as exc:
        output.update({'exc': exc})
        return output

    try:
        mcs_smarts, mcs_num_atoms, mcs_num_bonds = calc_mcs(target_mol, ref_mol, mcs_flags=['v'])
        mcs_strict_smarts, mcs_strict_num_atoms, mcs_strict_num_bonds = calc_mcs(target_mol, ref_mol, mcs_flags=['v', 'cr'])
    except tools.RDKitError as exc:
        output.update({'exc': exc})
        return output

    mcs_coverage = float(mcs_num_atoms) / target_mol.GetNumAtoms()
    mcs_strict_coverage = float(mcs_strict_num_atoms) / target_mol.GetNumAtoms()
    mcs_tanimoto = float(mcs_num_atoms) / (target_mol.GetNumAtoms() + ref_mol.GetNumAtoms() - mcs_num_atoms)
    # abs(float(target_mol.GetNumAtoms() - ref_mol.GetNumAtoms())) / mcs_num_atoms

    fp_tanimoto = calc_fp_tanimoto(target_mol, ref_mol)
    if (mcs_coverage >= mcs_coverage_threshold) and (fp_tanimoto >= fp_tanimoto_threshold):
        output.update({'fp_tanimoto': fp_tanimoto,
                       'mcs': {'coverage': mcs_coverage,
                               'strict_coverage': mcs_strict_coverage,
                               'tanimoto': mcs_tanimoto,
                               'smarts': mcs_smarts,
                               'num_atoms': mcs_num_atoms,
                               'num_bonds': mcs_num_bonds}})
        return output
    return None


def collect_ligs_prepare_mcs(target_mol, chemid_to_liginfo, pdbid_to_chemids, pdbids=None,
                             mcs_coverage_threshold=0.7, fp_tanimoto_threshold=0.4, num_threads=1):
    if pdbids is None:
        pdbids = list(pdbid_to_chemids.keys())  # Ligand Expo might be out-of-date

    chemids = set()
    for pdbid in pdbids:
        chemids.update(pdbid_to_chemids.get(pdbid, []))  # Ligand Expo might be out-of-date
    chemids = list(chemids)

    chemid_liginfo_list = [(chemid, chemid_to_liginfo.get(chemid, None)) for chemid in chemids]
    absent_chemids = [chemid_liginfo[0] for chemid_liginfo in chemid_liginfo_list if chemid_liginfo[1] is None]
    if len(absent_chemids) != 0:
        dev_logger.warning('the following chemids are absent in Ligand Expo Database: {}'.format(' '.join(absent_chemids)))

    args = {'target_mol': target_mol,
            'mcs_coverage_threshold': mcs_coverage_threshold,
            'fp_tanimoto_threshold': fp_tanimoto_threshold}
    mp_inputs = [chemid_liginfo + (args,) for chemid_liginfo in chemid_liginfo_list if chemid_liginfo[1]]
    with mp.Pool(num_threads) as pool:
        mp_outputs = pool.starmap(mp_func, mp_inputs)

    chemid_to_refinfo = {output['ref_chemid']: output for output in mp_outputs if output and 'exc' not in output}
    for output in mp_outputs:
        if output and 'exc' in output:
            dev_logger.warning('reference chemid: {}'.format(output['ref_chemid']), exc_info=output['exc'])

    pdbid_to_refinfos = {}
    for pdbid in pdbids:
        chemids = pdbid_to_chemids.get(pdbid, [])
        refinfos = [chemid_to_refinfo[chemid] for chemid in chemids if chemid in chemid_to_refinfo]
        if len(refinfos) != 0:
            pdbid_to_refinfos[pdbid] = refinfos
    return pdbid_to_refinfos
