#!/bin/bash
#SBATCH --job-name=ALK
#SBATCH --output=output.out
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=4GB
#SBATCH --time=02-00:00:00
#SBATCH --partition=main
#SBATCH --mail-user=hendrykn123@physik.fu-berlin.de
#SBATCH --mail-type=fail

module load gromacs/single/2021.3
#module load gromacs/single/2019

nLoop=20
nCore=1

bash cleaner.sh
python Run_Topology_Writer.py

gmx solvate -maxsol 1 -box 4.100 4.100 4.100 -cs pentane.gro -o pentane_box.gro -p topol.top
gmx solvate -maxsol 1000 -cp pentane_box.gro -cs spc216.gro -o pentane_solvate.gro -p topol.top

gmx grompp -f em.mdp -c pentane_solvate.gro -p topol.top -o en_min.tpr -maxwarn 1
gmx mdrun -nt $nCore -deffnm en_min

gmx grompp -f npt.mdp -c en_min.gro -p topol.top -o npt.tpr -maxwarn 1
gmx mdrun -nt $nCore -pin on -deffnm npt

# Run the first nvt sequence outside of loop, with npt.gro as input
gmx grompp -f nvt.mdp -c npt.gro -p topol.top -o nvt_1.tpr -maxwarn 1
gmx mdrun -nt $nCore -deffnm nvt_1 

# run the remaining nvt sequences inside loop
for j in $(eval echo "{1..$[$nLoop-1]}")
do
    gmx grompp -f nvt.mdp -c nvt_$j.gro -p topol.top -o nvt_$[$j + 1].tpr -maxwarn 1 
    gmx mdrun -v -nt $nCore -deffnm nvt_$[$j + 1]
done
