#!/bin/bash
python ../../modeller_tools.py 2w8t/2w8t_A_model.pdb CAPRI_T152_CASP_T1003.fasta 2w8t/2w8t.pdb A &
python ../../modeller_tools.py 2x8u/2x8u_B_model.pdb CAPRI_T152_CASP_T1003.fasta 2x8u/2x8u.pdb B &
python ../../modeller_tools.py 3tqx/3tqx_B_model.pdb CAPRI_T152_CASP_T1003.fasta 3tqx/3tqx.pdb B &
python ../../modeller_tools.py 4iw7/4iw7_A_model.pdb CAPRI_T152_CASP_T1003.fasta 4iw7/4iw7.pdb A &
python ../../modeller_tools.py 5txr/5txr_B_model.pdb CAPRI_T152_CASP_T1003.fasta 5txr/5txr.pdb B &
