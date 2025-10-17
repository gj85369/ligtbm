import json
import click
import tools
from rdkit import Chem
from rdkit.Chem import AllChem


def cs_sym_mappings(mol1, mol2, cs_smarts):
    cs_patt = Chem.MolFromSmarts(cs_smarts)
    matches1 = mol1.GetSubstructMatches(cs_patt, uniquify=False, useChirality=True)
    matches2 = mol2.GetSubstructMatches(cs_patt, uniquify=False, useChirality=True)
    # Maybe redundant to use uniquify=True second time

    mappings = []
    for match1 in matches1:
        for match2 in matches2:
            mapping = set(zip(match1, match2))
            if mapping not in mappings:
                mappings.append(mapping)
    mappings = [list(mapping) for mapping in mappings]
    return mappings


def best_rmsd_mapping_conf(mol, ref_mol, mappings):
    variants = []
    for cid in range(mol.GetNumConformers()):
        for mapping in mappings:
            rmsd = AllChem.AlignMol(mol, ref_mol, prbCid=cid, atomMap=mapping)  # aligned RMSD
            variants.append((rmsd, cid, mapping))
    min_rmsd, min_rmsd_cid, min_rmsd_mapping = min(variants, key=lambda x: x[0])
    return min_rmsd, min_rmsd_cid, min_rmsd_mapping


def min_cs_rmsd_conf_align(mol, ref_mol, cs_smarts):
    mappings = cs_sym_mappings(mol, ref_mol, cs_smarts)
    if len(mappings) == 0:
        raise tools.RDKitError('no mappings have been found between target mol and reference mol')
    min_rmsd, min_rmsd_cid, min_rmsd_mapping = best_rmsd_mapping_conf(mol, ref_mol, mappings)
    min_rmsd = AllChem.AlignMol(mol, ref_mol, prbCid=min_rmsd_cid, atomMap=min_rmsd_mapping)
    mol = tools.get_conf_mol(mol, min_rmsd_cid)
    return mol, min_rmsd, min_rmsd_mapping


@click.command()
@click.argument('target_sdf_file', type=click.Path(exists=True))
@click.argument('target_smiles', type=str)
@click.argument('ref_mol_file', type=click.Path(exists=True))
@click.argument('ref_smiles', type=str)
@click.argument('cs_smarts', type=str)
@click.argument('aligned_mol_file', type=click.Path())
@click.argument('mapping_rmsd_json_file', type=click.Path())
def _cs_confs_align_click(target_sdf_file, target_smiles, ref_mol_file, ref_smiles, cs_smarts,
                          aligned_mol_file, mapping_rmsd_json_file):
    mol = tools.mol_confs_from_sdf_file(target_sdf_file, target_smiles)
    ref_mol = tools.mol_from_mol_file(ref_mol_file, ref_smiles)
    mol, rmsd, mapping = min_cs_rmsd_conf_align(mol, ref_mol, cs_smarts)
    Chem.MolToMolFile(mol, aligned_mol_file)
    with open(mapping_rmsd_json_file, 'w') as f:
        mapping_rmsd_dict = {'CS_mapping': mapping, 'CS_RMSD': rmsd}
        json.dump(mapping_rmsd_dict, f, indent=4)


if __name__ == '__main__':
    _cs_confs_align_click()
