#!/bin/bash -l
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 200:00:00
#SBATCH --mem=8G
#SBATCH -o out_msd.out
#SBATCH -e out_msd.err
#SBATCH --job-name=submit_20m_msd

cd $SLURM_SUBMIT_DIR
conda --version
#source activate /home/siepmann/singh891/.conda/envs/halogen37
conda activate litfsi37
date
python msd.py
date

