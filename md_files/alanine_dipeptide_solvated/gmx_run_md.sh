#!/bin/bash
#set -eux

#SBATCH --job-name=alanine2
#SBATCH --output=output.out
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=6GB
#SBATCH --time=5-00:00:00
#SBATCH --partition=main
#SBATCH --mail-user=hendrykn123@physik.fu-berlin.de
#SBATCH --mail-type=fail

module load gromacs/single/2021.3

NATIVE=alanine-dipeptide.pdb
EM_NATIVE=em.pdb
FORCEFIELD=amber99
WATERMODEL=tip3p
DISTANCE=1.2
IGNORE_H=true

nLoop=50
nCore=16 

bash cleaner.sh
mkdir dihedral

# convert initial PDB:
gmx pdb2gmx -f $NATIVE -o start.gro $($IGNORE_H && echo -ignh) -ff $FORCEFIELD -water $WATERMODEL

# define simulation box:
gmx editconf -f start.gro -o box.gro -bt cubic -d $DISTANCE

# add solvent:
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top

# minimize energy:
gmx grompp -f em.mdp -c solvated.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# save the energy-minimized state:
echo 1 | gmx trjconv -f em.gro -s em.tpr -o $EM_NATIVE

# equilibrate with NVT:
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt

# NPT equilibrate:
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt

# Run the first md sequence outside of loop, with nvt.gro as input
gmx grompp -f nvt.mdp -c nvt.gro -p topol.top -o nvt_1.tpr -maxwarn 1
gmx mdrun -nt $nCore -deffnm nvt_1


# correct first trajectory (for plotting):
echo 1 | gmx trjconv -f nvt_1.xtc -s nvt_1.tpr -pbc nojump -o nvt_nojump_1.xtc
echo 1 | gmx trjconv -f nvt_1.xtc -s nvt_1.tpr -pbc mol -o nvt_center_1.xtc
echo 1 1 | gmx trjconv -f nvt_center_1.xtc -s nvt_1.tpr -fit rot+trans -o nvt_fitted_1.xtc

#compute dihedrals of first trajectory
gmx rama -f nvt_1.xtc -s nvt_1.tpr -xvg none -o Dihed.xvg
sed '1,17d' Dihed.xvg > dihedral/phi_psi_1.txt
rm *.xvg

# run the remaining nvt sequences inside loop
for j in $(eval echo "{1..$[$nLoop-1]}")
do
    gmx grompp -f md.mdp -c nvt_$j.gro -p topol.top -o nvt_$[$j + 1].tpr -maxwarn 1
    gmx mdrun -v -nt $nCore -deffnm nvt_$[$j + 1]
    gmx rama -f nvt_$[$j + 1].xtc -s nvt_$[$j + 1].tpr -xvg none -o Dihed.xvg
    sed '1,17d' Dihed.xvg > dihedral/phi_psi_$[$j + 1].txt
    rm *.xvg
done

