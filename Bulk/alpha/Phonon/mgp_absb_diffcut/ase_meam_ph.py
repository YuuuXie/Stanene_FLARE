import os
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read, write
from ase.phonons import Phonons
from ase.calculators.lammpsrun import LAMMPS

species = ["Sn"]
parameters = {
    "command": os.environ.get("lmp"),  # set up executable for ASE
    "pair_style": "meam/c",
    "pair_coeff": [f"* * library.meam Sn Sn.meam Sn"],
}

files = ["library.meam", "Sn.meam"]

# create ASE calc
calc = LAMMPS(
    label=f"tmp",
    keep_tmp_files=True,
    tmp_dir="./tmp/",
    parameters=parameters,
    files=files,
    specorder=species,
)

atoms = read("../scf/scf.in", format="espresso-in")
N = 8

# first relax
....

ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.05)
ph.run()
ph.read(acoustic=True)
ph.clean()

path = atoms.cell.bandpath('GLXG', npoints=154)
bs = ph.get_band_structure(path)
xcoords, label_xcoords, orig_labels = bs.get_labels()
label_xcoords = list(label_xcoords)
np.save("xcoords", xcoords)
np.save("label_xcoords", label_xcoords)
np.save("orig_labels", orig_labels)
np.save("energies", bs.energies)

plt.switch_backend("agg")
fig, axes = plt.subplots(figsize=(7, 4))
emax=0.03
bs.plot(ax=axes, emin=0.0, emax=emax)
plt.savefig("mgp_ph.png", dpi=200)
