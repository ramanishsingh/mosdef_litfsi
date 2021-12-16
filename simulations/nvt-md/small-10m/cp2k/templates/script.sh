{% extends base_script %}
{% block project_header %}
#SBATCH -t 96:00:00
#SBATCH --mem-per-cpu=2500MB
#SBATCH --partition=small
module purge
module load intel/cluster/2018
module load mkl
module load fftw
module load conda
conda --version
source activate /home/siepmann/singh891/.conda/envs/litfsi37
#source activate litfsi37
date >> execution.log
{{ super() }}
{% endblock %}

