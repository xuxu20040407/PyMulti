#!/bin/bash

#SBATCH --job-name=AI4ICF_Simu1D       # Name
#SBATCH --partition=64c512g          # cpu
#SBATCH -n 50                # total node
#SBATCH --ntasks-per-node=50     # core per node
#SBATCH --output=output/%j.out
#SBATCH --error=output/%j.err

module load miniconda3/4.10.3
source activate multinet
which python
nproc
python 1D_database_final.py
