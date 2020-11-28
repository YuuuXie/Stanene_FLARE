import os
from copy import deepcopy
import numpy as np
from ase.io import read, write

from ase.calculators.lammpsrun import LAMMPS
from ase import Atoms

species = ["Sn"]
pot_file = "abs_f.mgp"
parameters = {
    "command": os.environ.get("lmp"),  # set up executable for ASE
    "newton": "off",
    "pair_style": "mgp",
    "pair_coeff": [f"* * {pot_file} Sn yes yes"],
}

files = [pot_file]


def compute_mae(atoms, lmp_calc):
    dft_f = atoms.get_forces()
    atoms.pbc = True
    atoms.calc = lmp_calc
    lmp_f = atoms.get_forces()

    mae = np.mean(np.abs(dft_f - lmp_f))
    maf = np.mean(np.abs(dft_f))
    print("MAE MAF", mae, maf)
    return mae

atoms_list = read("../OTF/AL-3/4-DFT/dft_frames.xyz", index=":")
total_mae = 0.0
for a, atoms in enumerate(atoms_list): 
    # create ASE calc
    lmp_calc = LAMMPS(
        label=f"tmp",
        keep_tmp_files=True,
        tmp_dir="./tmp/",
        parameters=parameters,
        files=files,
        specorder=species,
    )
    
    mae = compute_mae(atoms, lmp_calc)
    total_mae += mae

total_mae /= len(atoms_list)
print("Total MAE", total_mae)
