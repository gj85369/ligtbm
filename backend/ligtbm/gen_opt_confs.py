import click
import tools
import os
from rdkit import Chem
from rdkit.Chem import AllChem

import logging.config
FILE_DIR = os.path.dirname(os.path.abspath(__file__))
logging.config.fileConfig(os.path.join(FILE_DIR, 'logging.conf'), disable_existing_loggers=False)
dev_logger = logging.getLogger('dev')


def gen_etkdg_mmff(mol, num_confs=1000, num_threads=1):
    num_confs_acceptable = int(num_confs * 0.25)  # 25% threshold chosen arbitrary
    mol = Chem.AddHs(mol)
    # Since the 2018.09 release of the RDKit, ETKDG is the default conformation generation method.
    AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, numThreads=num_threads)
    num_confs_generated = mol.GetNumConformers()
    if num_confs_generated != num_confs:
        dev_logger.warning('Generated only {} out of {} conformers'.format(num_confs_generated, num_confs))
    if num_confs_generated < num_confs_acceptable:
        raise tools.RDKitError('Cannot generate enough (at least {}) conformers for this molecule'.format(
            num_confs_acceptable))
    AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=num_threads)
    mol = Chem.RemoveHs(mol)
    return mol


def gen_opt_confs(mol, num_confs=1000, num_threads=1):
    patt = Chem.MolFromSmarts('[r;!r3;!r4;!r5;!r6;!r7]')
    if False:  # mol.HasSubstructMatch(patt)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        Chem.MolToMolFile(mol, 'BRIKARD_input.mol')
        # TODO call BRIKARD
    else:
        mol = gen_etkdg_mmff(mol, num_confs, num_threads)
    return mol


@click.command()
@click.argument('smiles', type=str)
@click.argument('sdf_file', type=click.Path())
@click.option('--num_confs', type=int, default=1000, help='Number of conformers.')
@click.option('--num_threads', type=int, default=1, help='Number of threads for parallel generation.')
def _gen_opt_confs_click(smiles, sdf_file, num_confs, num_threads):
    '''
        RDKit & BRIKARD(temporarily unavailable)-based script for conformers generation.
    '''
    mol = tools.safe_mol_init(smiles)
    mol = gen_opt_confs(mol, num_confs, num_threads)
    tools.mol_confs_to_sdf_file(mol, sdf_file)


if __name__ == '__main__':
    _gen_opt_confs_click()
