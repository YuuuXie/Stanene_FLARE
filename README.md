# FLARE Force Field for Stanene

- This repository provides the code for developing FLARE force field of stanene, which is used to generate results in this [paper [1]](https://arxiv.org/abs/2008.11796).
- The FLARE package is installed from [MIR-FLARE](https://github.com/mir-group/flare).
  ```
  pip install mir-flare
  ```
- The generated data reported in the paper has been uploaded to Materials Cloud Archive. You can see it [here [2]](https://archive.materialscloud.org/record/2020.99), including LAMMPS coefficient file for the force field, trajectories and videos.

## About the versions

- Part of the data is NOT generated from the latest version of FLARE. The original version of code for the data is released as [v1.0](https://github.com/YuuuuXie/Stanene_FLARE/releases/tag/v1.0), which does not predict atomic or global energy, and only predicts atomic forces. 
- But we will try catching up with the latest FLARE repository, and we recommend you to use the latest code that we provide.

## The repo consists of

- OTF: on-the-fly (OTF) training scripts of stanene, and postprocessing
- Hierarchical: hierarchical training of stanene, and postprocessing
- Convergence: force & uncertainty convergence of MGP w.r.t. GP in accuracy
- Bulk: investigation of alpha & beta phases including phonon and bulk modulus, and molecular dynamics of stanene-bcc transformation
- Liquid: LAMMPS molecular dynamics of stanene melting process, and ab-initio molecular dynamics of melted tin
- Timing: speed comparison of different methods
- Plot.ipynb: a notebook for data analysis, and postprocessing

## The repo doesn't consist of

- The dumped GP file (`.pickle` or `.json`) with training data of Sn
- The LAMMPS coefficient file of MGP
- The (LAMMPS & ab-initio) molecular dynamics trajectories

because they are too large. The above files are available in [Materials Cloud Archive doi: 10.24435/materialscloud:cs-tf](https://archive.materialscloud.org/record/2020.99) [2]

## How to cite?

- Paper:

  [1] Xie, Yu, et al. "Fast Bayesian Force Fields from Active Learning: Study of Inter-Dimensional Transformation of Stanene." arXiv preprint [arXiv:2008.11796](https://arxiv.org/abs/2008.11796) (2020).

- Data:

  [2] Yu Xie, Jonathan Vandermause, Lixin Sun, Andrea Cepellotti, Boris Kozinsky, Fast Bayesian force fields from active learning: study of inter-dimensional transformation of stanene, Materials Cloud Archive 2020.99 (2020), [doi: 10.24435/materialscloud:cs-tf](https://archive.materialscloud.org/record/2020.99).
