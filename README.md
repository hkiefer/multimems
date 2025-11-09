# Multimems

Python 3 tool suite for the computation of memory kernel matrices in a multi-dimensional generalized Langevin equation (GLE) from a set of time-series data. A publication about the multi-dimensional GLE is in preparation; a preprint version is available at: https://arxiv.org/abs/2506.05966.

To get started, run

    pip install .

to install. See `example/example_extraction.ipynb` for an introduction of the extraction scheme, and `example/example_embedding.ipynb` for an introduction of the Markovian embedding scheme. The notebook `example/fitting_multi_dim_kernels_exp_matrices.ipynb` includes an instruction how to fit multi-dimensional memory kernels with multi-exponential matrices.

In `md_files/`, input files for the molecular dynamics simulations with GROMACS, to obtain the dihedral angle trajectories used in https://arxiv.org/abs/2506.05966, are provided.

