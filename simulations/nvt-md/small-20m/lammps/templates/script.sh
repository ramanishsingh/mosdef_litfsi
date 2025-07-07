{% extends base_script %}
{% block project_header %}
#SBATCH -t 192:00:00
#SBATCH --mem-per-cpu=2500MB
#SBATCH --partition=silver
module purge
module load intel
#module load mkl
#module load fftw
#module load conda
module load ompi
module load lammps
. /home/rs/anaconda3/etc/profile.d/conda.sh
conda activate

conda --version
#source activate /home/siepmann/singh891/.conda/envs/litfsi37
conda activate litfsi37
date >> execution.log
{{ super() }}
{% endblock %}

