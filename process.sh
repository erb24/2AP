#!/bin/bash -l
#SBATCH --job-name="process"
#SBATCH --output="./process.log"
#SBATCH -A uoo104
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=ALL
#SBATCH -t 0:30:00

if [ "$#" -ne 3 ]; then
        echo "Script needs two command line arguments: trajectory (.xtc) file, topology (.tpr) file, and name of the system (string)."
        echo "Example: sh process.sh 1UBQ.xtc 1UBQ.tpr 1UBQ"
        exit 1
fi

#Code to process the trajectory for analysis
traj=$1 
top=$2 
protname=$3

module load gromacs

gmx=`which gmx_mpi`
echo $gmx

#Dump centered topology file
echo "1" | $gmx trjconv -f $traj -s $top -dump 0 -o top.pdb

$gmx editconf -f top.pdb -c -center 0 0 0 -o top.pdb

#top='top.pdb'

#Correct PBCs
echo "1 1 1" | $gmx trjconv -quiet -f $traj -s $top -o after_pbc.xtc -pbc cluster -center -boxcenter zero

echo "1 1" | $gmx trjconv -quiet -f after_pbc.xtc -s $top -fit rot -o after_rot.xtc

rm -v ./#*#
