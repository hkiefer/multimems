# md_files

This folder contains input files for molecular dynamics (MD) simulations with GROMACS, used to generate dihedral angle trajectories for the axiv paper (https://arxiv.org/abs/2506.05966). The provided systems are solvated pentane and dialanine dipeptide.

## Folder Structure

- `alanine_dipeptide_solvated/`
    - `alanine-dipeptide.pdb`: Structure file for dialanine dipeptide.
    - `start.gro`: Starting configuration.
    - `em.mdp`, `npt.mdp`, `nvt.mdp`: GROMACS parameter files for energy minimization, NPT equilibration, and MD runs.
    - `gmx_run_md.sh`: Bash script to run the MD simulation.
    - `cleaner.sh`: Utility script for cleaning up files.


- `pentane_solvated/`
    - `pentane.gro`, `pentane_box.gro`, `pentane_solvate.gro`: Structure files for pentane and solvated system.
    - `pentane.itp`, `spce_custom.itp`: Topology files for pentane and solvent.
    - `em.mdp`, `npt.mdp`, `nvt.mdp`: GROMACS parameter files for energy minimization, NPT equilibration, and MD runs.
    - `index.ndx`: Index file for atom selection.
    - `gmx_run_md.sh`: Bash script to run the MD simulation.
    - `gmx_dihedral.sh`: Script to calculate dihedral angles.
    - `cleaner.sh`: Utility script for cleaning up files.
    - Additional utilities: `gmx_Index`, `gmx_XTC_to_XYZ`, `Run_Topology_Writer.py`.

## Running MD Simulations

To run the MD simulations, navigate to the respective folder and execute:

    bash gmx_run_md.sh

This will perform the MD workflow using GROMACS.

## Dihedral Angle Calculation

After the MD simulations, e.g. in `pentane_solvated/`, calculate dihedral angles using:

      bash gmx_dihedral.sh

This will extract the relevant dihedral angle trajectories for further analysis. For `alanine_dipeptide_solvated/`, this is already done running `gmx_run_md.sh`.