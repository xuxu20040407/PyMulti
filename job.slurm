#!/bin/bash

#SBATCH --job-name=AI4ICF_Simu1D       # Name
#SBATCH --partition=64c512g          # cpu 
#SBATCH -n 100                 # total cpu 
#SBATCH --ntasks-per-node=64     # cpu/node
#SBATCH --output=out/%j.out
#SBATCH --error=err/%j.err 

python 1D_database.py
