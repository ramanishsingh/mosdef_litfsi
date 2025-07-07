{% extends "base_script.sh" %}
{% block header %}
#!/bin/sh -l
#PBS -j oe
#PBS -l nodes=1:ppn=16
#PBS -l walltime={{ walltime|format_timedelta }}
#PBS -q standard
#PBS -V
#PBS -m abe
#PBS -e ./error/
#PBS -o ./output/
#PBS -M xiaobo.lin@vanderbilt.edu
#PBS -N small_litfsi20

conda activate myenv
module load gromacs/2020

{% endblock %}