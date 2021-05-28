#!/bin/bash
#SBATCH --job-name="CDcalc-python"
#SBATCH --output="CDcalc-python.log"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=ALL
#SBATCH --mem 25G 
#SBATCH -t 48:00:00

module load gromacs
echo "[ State 1 ]" > frames.ndx
echo "1" >> frames.ndx

traj="after_rot.xtc"
top="NPT_prod.tpr"
gmx_mpi trjconv -f $traj -s $top -o Atoms-anly.pdb -n Atoms.ndx
echo "1" | gmx_mpi trjconv -f Atoms-anly.pdb -s Atoms-anly.pdb -o top.pdb -fr frames.ndx
echo "1" | gmx_mpi trjconv -f Atoms-anly.pdb -s Atoms-anly.pdb -o Atoms-anly.g96

python gen_unformatted_traj.py
python CDcalc_extended_dipole.py 

exit
