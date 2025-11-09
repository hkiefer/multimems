#!/bin/bash
#SBATCH --job-name=XYZ
#SBATCH --output=output.out
#SBATCH --nodes=1
#SBATCH --ntasks=1 
#SBATCH --mem=1GB
#SBATCH --time=4:00:00
#SBATCH --partition=main

module load gromacs/single/2021.3

mkdir pentane_xyz

for j in {1..10}
do
    gmx trjconv -f nvt_$j.xtc -pbc nojump -ndec 3 -o trj.xtc
    echo "2" | gmx traj -f trj.xtc -n index.ndx -s npt.gro -ox trjxyz.xvg
    sed '1,35d' trjxyz.xvg > pentane_xyz/xyz_traj_$(printf "%02d" $j).txt
    rm trjxyz.xvg trj.xtc
done

