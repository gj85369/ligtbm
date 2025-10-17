#!/usr/bin/env python3

import numpy as np
import tempfile
import subprocess
import json
import prody
import sys
import shutil

ref_pdb = prody.parsePDB('data/old_rmin_out.pdb')

with tempfile.NamedTemporaryFile() as fp_pdb, tempfile.NamedTemporaryFile() as fp_dat:
    subprocess.check_call(['../../rmin/build/rmin',
        'data/0_123.8X1.5V59_A.rec.pdb', 'data/0_123.8X1.5V59_A.rec.psf', 'data/0_123.8X1.5V59_A.lig.json', 'data/0_123.8X1.5V59_A.restraints.pdb',
        'data/parm_new.prm', 'data/pdbamino_new.rtf', 'data/libmol.params',
        fp_pdb.name, fp_dat.name, '-s', '10', '-c', '5', '-p', '1'], stdout=subprocess.DEVNULL)
    test_pdb = prody.parsePDB(fp_pdb.name)
    fp_dat.seek(0)
    energies = json.load(fp_dat)
    shutil.copy(fp_pdb.name, 'output.pdb')

en_total = energies['total']
rmsd = prody.calcRMSD(ref_pdb.noh, test_pdb.noh)
lig_sel = 'heavy and chain X'
lig_rmsd = prody.calcRMSD(ref_pdb.select(lig_sel), test_pdb.select(lig_sel))
print('HEAVY ATOM RMSD', rmsd)
print('HEAVY ATOM LIGAND RMSD', lig_rmsd)
print('TOTAL ENERGY', en_total)

if rmsd > 0.5:
    sys.exit(-1)

