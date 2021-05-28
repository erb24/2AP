#!/bin/bash
#SBATCH --job-name="CDcalc-python"
#SBATCH --output="CDcalc-python.log"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=ALL
#SBATCH --mem 25G 
#SBATCH -t 48:00:00

BD="../"
python ${BD}/gen_unformatted_traj.py
python ${BD}/CDcalc_extended_dipole.py

exit