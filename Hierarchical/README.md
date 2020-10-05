# Hierarchical training of stanene force field

After the OTF training is complete, we use MGP force field in LAMMPS for long time-scale simulations.
Snapshots are picked up from the LAMMPS trajectory and calculated in DFT for validation and re-training.

1. `add_snap` in mgp-"n": Add snapshots from DFT calculations in mgp-"n-1"
2. `train_hyps`: Train hyperparameters of GP
3. `build_mgp`: Use optimized hyps to build MGP
4. `lmp`: Run LAMMPS with MGP
5. `uncertainties`: Predict uncertainties 
6. `dft`: Pick up three snapshots for DFT
