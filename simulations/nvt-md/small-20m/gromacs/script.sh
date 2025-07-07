{% extends base_script %}
{% block project_header %}
#SBATCH -t 192:00:00
#SBATCH --mem-per-cpu=2500MB
#SBATCH --partition=silver
#SBATCH --ntasks=24                  # Number of MPI tasks (i.e. processes)
#SBATCH --cpus-per-task=1            # Number of cores per MPI task 

module purge
module load gromacs
module load cuda/10.2 
module load gcc
export GMXLIB=$GMXLIB:/soft/gromacs/5.0.0-impi-single/share/gromacs/top

. /home/rs/anaconda3/etc/profile.d/conda.sh
conda activate

conda --version
#source activate /home/siepmann/singh891/.conda/envs/litfsi37
conda activate litfsi37
date >> execution.log
{{ super() }}
{% endblock %}

