from rdkit import Chem
from rdkit.Chem import AllChem


class RDKitError(Exception):
    pass


def mol_confs_to_sdf_file(mol, sdf_file):
    writer = Chem.SDWriter(sdf_file)
    for cid in range(mol.GetNumConformers()):
        writer.write(mol, confId=cid)
    writer.close()


def get_conf_mol(mol, cid):
    mol = Chem.Mol(mol)
    for i in range(mol.GetNumConformers()):
        if i != cid:
            mol.RemoveConformer(i)
    mol.GetConformer(cid).SetId(0)
    return mol


def safe_mol_init(arg):
    if type(arg) == dict:
        if 'SMILES_CACTVS' not in arg:  # for now we take into consideration only CACTVS SMILES from Ligand Expo
            raise RDKitError('CACTVS SMILES string is not available')
        smiles = arg['SMILES_CACTVS']
    elif type(arg) == str:
        smiles = arg
    else:
        raise RDKitError('the argument is not a liginfo dictionary nor a SMILES string')

    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return mol
    else:
        raise RDKitError('SMILES {} is not RDKit-correct'.format(smiles))


def mol_from_mol_file(mol_file, arg):
    mol = safe_mol_init(arg)
    mol_unsanitized = Chem.MolFromMolFile(mol_file, sanitize=False)
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(mol, mol_unsanitized)
        # TODO stereochemistry assignment
        Chem.SanitizeMol(mol)
        return mol
    except:
        raise RDKitError('failed attempt to assign bond orders and stereochemistry')