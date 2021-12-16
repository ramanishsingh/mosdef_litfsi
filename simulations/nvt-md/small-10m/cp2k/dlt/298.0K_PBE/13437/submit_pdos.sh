#!/bin/bash -l
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 30:00:00
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --job-name=pdos_idx
#SBATCH --partition=small
#SBATCH --mem=20gb

cd $SLURM_SUBMIT_DIR
module purge
module load intel/cluster/2018
module load mkl
module load fftw
date
cp2k.popt -i litfsi_pdos_input.inp -o litfsi_pdos_input.out
date

