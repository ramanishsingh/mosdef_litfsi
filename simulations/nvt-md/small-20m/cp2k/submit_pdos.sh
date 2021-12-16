#!/bin/bash -l
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 32:00:00
#SBATCH --mem=20G
#SBATCH -o out_pdos.out
#SBATCH -e err_pdos.err
#SBATCH --job-name=submit_pdos

cd $SLURM_SUBMIT_DIR
module purge
module load conda
conda --version
#source activate /home/siepmann/singh891/.conda/envs/halogen37
conda activate litfsi37
date
python pdos_submitter.py
date

