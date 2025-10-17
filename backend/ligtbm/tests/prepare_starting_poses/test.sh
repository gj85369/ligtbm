#!/bin/bash -l

#$ -pe omp 8
#$ -l h_rt=12:00:00
#$ -N prepare_starting_poses
#$ -m ea
#$ -j y

module load blast+
export BLASTDB=/projectnb/k3domics/ligtbm/blast/
export PATH="/projectnb/k3domics/sergeik/miniconda3/envs/ligtbm/bin:$PATH"
python ../../prepare_starting_poses.py --num_threads $NSLOTS --database pdbaa param.json
