#!/bin/bash
#SBATCH --job-name=XYZ
#SBATCH --output=output.out
#SBATCH --nodes=1
#SBATCH --ntasks=1 
#SBATCH --mem=1GB
#SBATCH --time=4:00:00
#SBATCH --partition=main

mkdir dihedrals
module load gromacs/single/2021.3           

#left dihedral
for j in {1..2}
do
    echo "2" | gmx angle -type dihedral -n index.ndx -f nvt_$j.xtc -ov Dihed.xvg
    sed '1,17d' Dihed.xvg > dihedrals/Dihedral_$j.txt
    rm *.xvg

done

#right dihedral
for j in {1..2}
do
    echo "3" | gmx angle -type dihedral -n index.ndx -f nvt_$j.xtc -ov Dihed.xvg
    sed '1,17d' Dihed.xvg > dihedrals/Dihedral2_$j.txt
    rm *.xvg

done

