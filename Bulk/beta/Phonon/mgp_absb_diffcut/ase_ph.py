import os
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read, write
from ase import Atoms
from ase.phonons import Phonons
from ase.calculators.lammpsrun import LAMMPS
from ase.optimize import BFGS

species = ["Sn"]
pot_file = "absb_diffcut.mgp"
parameters = {
    "command": os.environ.get("lmp"),  # set up executable for ASE
    "pair_style": "mgp",
    "pair_coeff": [f"* * {pot_file} Sn yes yes"],
}

files = [pot_file]

# create ASE calc
calc = LAMMPS(
    label=f"tmp",
    keep_tmp_files=True,
    tmp_dir="./tmp/",
    parameters=parameters,
    files=files,
    specorder=species,
)

#atoms = read("/n/kozinsky_lab/Users/xiey/Stanene/beta-0/scf/scf.in", format="espresso-in")
#atoms = read("../scf/scf.in", format="espresso-in")
a = 2.9 
c = 1.61
scaled_positions = np.array([[0.0, 0.0, 0.0], [0.25, 0.75, 0.5]])
cell = np.array([[-a, a, c], [a, -a, c], [a, a, -c]])
atoms = Atoms(
    symbols="Sn2",
    scaled_positions=scaled_positions,
    cell=cell,
    pbc=True,
)

N = 8
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.1)
ph.run()
ph.read(acoustic=True)
ph.clean()

path = atoms.cell.bandpath('GXNG', npoints=154)
bs = ph.get_band_structure(path)
xcoords, label_xcoords, orig_labels = bs.get_labels()
label_xcoords = list(label_xcoords)
np.save("xcoords", xcoords)
np.save("label_xcoords", label_xcoords)
np.save("orig_labels", orig_labels)
np.save("energies", bs.energies)

plt.switch_backend("agg")
fig, axes = plt.subplots(figsize=(7, 4))
emax=0.02
bs.plot(ax=axes, emin=0.0, emax=emax)
plt.savefig("mgp_ph.png", dpi=200)
