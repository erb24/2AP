#!/bin/bash -l

BD=$PWD
cd $BD
sys="2AP"

#################################################################################
#       SOLVATION AND ENERGY MINIMIZATION
#################################################################################

echo $sys > sys.txt
echo "SOL" > SOL
echo "1" > pdb2gmx.inp # selects local force-field directory
echo "1" >> pdb2gmx.inp # selects water model
gmx_mpi pdb2gmx -f ${sys}.pdb -p ${sys}.top -o ${sys}.gro -ignh  < pdb2gmx.inp
gmx_mpi editconf -f ${sys}.gro -o out_${sys}.gro -d 1 -bt cubic
gmx_mpi solvate -cp out_${sys}.gro -cs -p ${sys}.top -o psolv_${sys}.gro
touch ions.mdp
gmx_mpi grompp -f ions.mdp -p ${sys}.top -c psolv_${sys}.gro -o ions.tpr -maxwarn 5
echo "SOL" | gmx_mpi genion -s ions.tpr -o ionized.pdb -p ${sys}.top -neutral -conc 0.1
gmx_mpi grompp -f em.mdp -p ${sys}.top -c ionized.pdb -o em.tpr -maxwarn 5
gmx_mpi mdrun -v -s em.tpr -deffnm em -c em.pdb


exit

