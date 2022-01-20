#!/bin/sh -l
#PBS -j oe
#PBS -l nodes=1:ppn=16
#PBS -l walltime=96:00:00
#PBS -q standard
#PBS -V
#PBS -m abe
#PBS -e ./error/
#PBS -o ./output/
#PBS -M xiaobo.lin@vanderbilt.edu
#PBS -N small_analy

cd $PBS_O_WORKDIR/
conda activate myenv
source activate myenv
date
python analysis.py

date

