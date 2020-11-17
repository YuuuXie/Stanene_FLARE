import os
from copy import deepcopy
import numpy as np
from ase.io import read, write

from ase.calculators.lammpsrun import LAMMPS
from ase import Atoms


species = ["Sn"]
pot_file = "absb_diffcut_3frames.mgp"
parameters = {
    "command": os.environ.get("lmp"),  # set up executable for ASE
    "pair_style": "mgp",
    "pair_coeff": [f"* * {pot_file} Sn yes yes"],
}

files = [pot_file]

# create ASE calc
lmp_calc = LAMMPS(
    label=f"tmp",
    keep_tmp_files=True,
    tmp_dir="./tmp/",
    parameters=parameters,
    files=files,
    specorder=species,
)


def compute_energy(atoms, v):
    en_list = []

    cell = np.array(atoms.get_cell())
    volume = atoms.get_volume()
    new_volume = v * len(atoms)
    new_cell = cell * (new_volume / volume) ** (1 / 3)

    atoms.set_cell(new_cell, scale_atoms=True)

    atoms.calc = deepcopy(lmp_calc)
    energy = atoms.get_potential_energy()

    E_per_atom = energy / len(atoms)
    en_list.append(E_per_atom)
    print(atoms.get_volume(), E_per_atom, flush=True)
    return E_per_atom

a = 2.91 
c = 1.61
#a = 3.028
#c = 1.65
#scaled_positions = np.array([[0.0, 0.0, 0.0], [0.25, 0.75, 0.5]])
#cell = np.array([[-a, a, c], [a, -a, c], [a, a, -c]])
#atoms = Atoms(
#    symbols="Sn2",
#    scaled_positions=scaled_positions,
#    cell=cell,
#    pbc=True,
#)
atoms = read("scf.pwi", format="espresso-in")
new_cell = np.diag([2 * a, 2 * a, 4 * c])
atoms.set_cell(new_cell, scale_atoms=True)

#v = 207.40020436537114 / len(atoms)
v = atoms.get_volume() / len(atoms)
#compute_energy(atoms, v)

V_list = []
E_list = []

v_range = np.linspace(v * 0.8, v * 1.2, 11)
for v in v_range:
    scaled_atoms = deepcopy(atoms)
    E_per_atom = compute_energy(scaled_atoms, v)
    E_list.append(E_per_atom)
    V_list.append(v)



from ase.units import kJ
from ase.eos import EquationOfState

V = V_list
E = E_list

eos = EquationOfState(V, E)
v0, e0, B = eos.fit()
print(B / kJ * 1.0e24, 'GPa')
