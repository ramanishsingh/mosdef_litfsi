
# Project Setup and Simulation Instructions

## Conda Environment Creation

Setting up, running, and analyzing the simulations require various Python packages which can be downloaded and installed from conda channels. Creating a conda environment helps minimize interference with other installed packages.

You can create the conda environment using either the `environment.yml` file or the set of commands below.

### Method 1: Using `environment.yml`

```bash
# Clone the repository
git clone https://github.com/ramanishsingh/mosdef_litfsi.git
cd mosdef_litfsi

# Create the environment
conda env create -f environment.yml --yes
conda activate mosdef_litfsi311
cd ..
```

### Method 2: Using Commands

```bash
mamba create -n mosdef_litfsi311 python=3.11
mamba activate mosdef_litfsi311
mamba install -c conda-forge mbuild
mamba install -c conda-forge foyer
pip install signac-flow==0.18.1 signac==1.7.0 mdtraj freud-analysis unyt ele==0.2.0 progressbar periodictable
conda install -c conda-forge openbabel --yes

# Software from source
mkdir software
cd software

# Scattering
git clone https://github.com/ramanishsingh/scattering.git
cd scattering/
git checkout partial_sq_atom_name
pip install -e .
cd ..

# ilforcefields
git clone https://github.com/mattwthompson/ilforcefields.git
cd ilforcefields
pip install -e .
cd ..

# mtools
git clone https://github.com/XiaoboLinlin/mtools
cd mtools
pip install -e .
cd ..

cd ..
```

## CP2K Installation

### Using a Package Manager

#### On Linux:

```bash
sudo apt install cp2k
```

#### On macOS:

```bash
brew install cp2k
```

For installation from source on clusters, visit: https://www.cp2k.org/howto

## mosdef_cp2k_writer and cp2kmdpy Installation

These tools generate input files and run CP2K simulations.

Activate the `mosdef_litfsi311` conda environment and run:

```bash
cd software

# cp2kmdpy
git clone https://github.com/ramanishsingh/cp2kmdpy.git
cd cp2kmdpy
# Modify runners.py or runners_mpi.py to replace "cp2k.popt" with your preferred command
pip install -e .
cd ..

# mosdef_cp2k_writer
git clone https://github.com/ramanishsingh/mosdef_cp2k_writer.git
cd mosdef_cp2k_writer
pip install -e .
cd ..
```

## Running and Analyzing Simulations

### CP2K Simulations

```bash
cd mosdef_litfsi/simulations/nvt-md/small-10m/cp2k
python init.py
python project.py run -o copy_setter
python project.py run -o copy_structure
python project.py run -o md_files
python project.py run -o run_md

# To submit jobs to cluster:
# python project.py submit -o run_md

# To restart simulations:
python project.py submit -o restart_md
```

**Cluster-specific settings (queue name, walltime, etc.) should be edited in** `templates/script.sh`.

### Analysis

```bash
python sfac.py              # Computes structure factor
python analysis.py          # Computes RDFs
python li_o_neighbor_analysis.py  # Neighbor analysis
```

Repeat the same steps for the 20m system.

### LAMMPS Simulations

```bash
cd mosdef_litfsi/simulations/nvt-md/small-10m/lammps
# Update LAMMPS run command in project.py accordingly

python init.py
python project.py run       # Run operations, may need to be rerun until completion
python project.py submit    # To submit jobs to cluster
```

### GROMACS Simulations

```bash
cd mosdef_litfsi/simulations/nvt-md/small-10m/gromacs
# Update GROMACS run command in project.py accordingly

python init.py
python project.py submit
```

Analysis commands are the same as those used for CP2K.
