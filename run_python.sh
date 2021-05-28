#!/bin/bash
#SBATCH --job-name="CDcalc-python"
#SBATCH --output="CDcalc-python.log"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=ALL 
#SBATCH -t 24:00:00

module load gromacs
echo "[ State 1 ]" > frames.ndx
echo "1" >> frames.ndx

traj="after_rot.xtc"
top="NPT_prod.tpr"
BD=$PWD
for i in `seq 10 10 1000`
do
  mkdir -v ${i}ns/
  dir=${i}ns
  gmx_mpi trjconv -f $traj -s $top -o ${dir}/Atoms-anly.pdb -n Atoms.ndx
  echo "1" | gmx_mpi trjconv -f ${dir}/Atoms-anly.pdb -s ${dir}/Atoms-anly.pdb -o ${dir}/top.pdb -fr frames.ndx
  echo "1" | gmx_mpi trjconv -f ${dir}/Atoms-anly.pdb -s ${dir}/Atoms-anly.pdb -o ${dir}/Atoms-anly.g96
  cd $dir/
  cp -v ${BD}/run_CD.sh ./
  sbatch run_CD.sh
  cd $BD
done
exit
