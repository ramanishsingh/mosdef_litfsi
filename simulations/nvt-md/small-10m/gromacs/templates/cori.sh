{% extends "base_script.sh" %}
{% block header %}
#!/bin/bash
#SBATCH --job-name=sample
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time={{ walltime|format_timedelta }}
#SBATCH -C haswell
#SBATCH --error=/global/project/projectdirs/m1046/Xiaobo/project/cosolvent_project/output/error/error_%j.err
#SBATCH --output=/global/project/projectdirs/m1046/Xiaobo/project/cosolvent_project/output/run_out_%j.log 
#SBATCH -q regular

echo working > test.txt

module load python
source activate myenv
module load gromacs/2020.1.hsw

{% endblock %}

