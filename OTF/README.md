# On-the-fly training of the stanene force field

1. `atom_setup.py` constructs the supercell of stanene
2. `flare_setup.py` initializes parameters for FLARE force field
3. `qe_setup.py` sets input parameters of Quantum Espresso (QE). Set to use 288 CPUs in total.
4. `otf_run.py` sets parameters for OTF training and molecular dynamics. Set to run 100 ps in total.
5. `postprocess.py` analyzes the OTF trajectory and generate a figure of mean square displacement

To launch the training, use
```
python otf_run.py
```
