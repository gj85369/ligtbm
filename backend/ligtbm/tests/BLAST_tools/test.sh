#!/bin/bash -l

#$ -pe omp 8
#$ -l h_rt=12:00:00
#$ -N blast_example
#$ -m ea
#$ -j y

module load blast+
export BLASTDB=/projectnb/k3domics/ligtbm/blast/
export PATH="/projectnb/k3domics/sergeik/miniconda3/bin:$PATH"
python BLAST_tools.py --database pdbaa --num_threads $NSLOTS --identity_threshold 0.7 example.fasta example.json
