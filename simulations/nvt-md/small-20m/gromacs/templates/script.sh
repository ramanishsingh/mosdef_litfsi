{% extends base_script %}
{% block project_header %}


#SBATCH -p k40
#SBATCH --gres=gpu:k40:1
#SBATCH -t 23:30:00   #Walltime in hh:mm:ss or d-hh:mm:ss
#SBATCH --ntasks 1   #Number of cores

module purge
module load intel

module load gromacs/2019.1-ompi-gpu

. /home/rs/anaconda3/etc/profile.d/conda.sh
conda activate

conda --version
#source activate /home/siepmann/singh891/.conda/envs/litfsi37
conda activate litfsi38
date >> execution.log
{{ super() }}
{% endblock %}

