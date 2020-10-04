# FLARE Force Field for Stanene

- This repository provides the code for developing [FLARE](https://github.com/mir-group/flare) force field of stanene, which is used to generate results in this [paper [1]](https://arxiv.org/abs/2008.11796).
- The generated data reported in the paper has been uploaded to Materials Cloud Archive. You can see it [here [2]](https://archive.materialscloud.org/record/2020.99), including LAMMPS coefficient file for the force field, trajectories and videos.

## About the versions

- The data is NOT generated from the latest version of FLARE. The original version of code for the data is released as [v1.0](). 
- But we will try catching up with the latest FLARE repository, and we recommend you to use the latest code that we provide.

## The code consists of

1. On-the-fly (OTF) training of stanene
2. Postprocessing of OTF training
3. Hierarchical training of stanene
4. Postprocessing of Hierarchical training

## How to cite?

- Paper:

  [1] Xie, Yu, et al. "Fast Bayesian Force Fields from Active Learning: Study of Inter-Dimensional Transformation of Stanene." arXiv preprint [arXiv:2008.11796](https://arxiv.org/abs/2008.11796) (2020).

- Data:

  [2] Yu Xie, Jonathan Vandermause, Lixin Sun, Andrea Cepellotti, Boris Kozinsky, Fast Bayesian force fields from active learning: study of inter-dimensional transformation of stanene, Materials Cloud Archive 2020.99 (2020), [doi: 10.24435/materialscloud:cs-tf](https://archive.materialscloud.org/record/2020.99).
