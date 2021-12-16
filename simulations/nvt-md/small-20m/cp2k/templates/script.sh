{% extends base_script %}
{% block project_header %}
#SBATCH -t 96:00:00
#SBATCH --mem-per-cpu=3g
#SBATCH --partition=small
module purge
module load intel/cluster/2018
module load mkl
module load fftw
module load conda
source activate /home/siepmann/singh891/.conda/envs/litfsi37
date >> execution.log
{{ super() }}
{% endblock %}

