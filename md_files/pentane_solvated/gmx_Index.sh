
module load GROMACS/2020-fosscuda-2019b

#gmx make_ndx -f npt.gro -o index_freeze.ndx
gmx make_ndx -f npt.gro -o index.ndx
