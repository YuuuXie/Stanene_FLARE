# FLARE Force Field for Stanene

- This repository provides the code for developing [FLARE](https://github.com/mir-group/flare) force field of stanene, which is used to generate results in this [paper [1]](https://arxiv.org/abs/2008.11796).
- The generated data reported in the paper has been uploaded to Materials Cloud Archive. You can see it [here [2]](https://archive.materialscloud.org/record/2020.99), including LAMMPS coefficient file for the force field, trajectories and videos.

## About the versions

- The data is NOT generated from the latest version of FLARE. The original version of code for the data is released as [v1.0](https://github.com/YuuuuXie/Stanene_FLARE/releases/tag/v1.0), which does not predict atomic or global energy, and only predicts atomic forces. 
- But we will try catching up with the latest FLARE repository, and we recommend you to use the latest code that we provide.

## The code consists of

1. On-the-fly (OTF) training of stanene
2. Postprocessing of OTF training
3. Hierarchical training of stanene
4. Postprocessing of Hierarchical training

## How to cite?

- Paper:

  [1] Xie, Y., Vandermause, J., Sun, L. et al. Bayesian force fields from active learning for simulation of inter-dimensional transformation of stanene. npj Comput Mater 7, 40 (2021). https://doi.org/10.1038/s41524-021-00510-y

- Data:

  [2] Yu Xie, Jonathan Vandermause, Lixin Sun, Andrea Cepellotti, Boris Kozinsky, Fast Bayesian force fields from active learning: study of inter-dimensional transformation of stanene, Materials Cloud Archive 2020.99 (2020), [doi: 10.24435/materialscloud:cs-tf](https://archive.materialscloud.org/record/2020.99).
