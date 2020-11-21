import os
from copy import deepcopy
import numpy as np
from ase.io import read, write
from ase import Atoms

from ase.calculators.lammpsrun import LAMMPS

species = ["Sn"]
pot_file = "absb_diffcut_6frames.mgp"
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

a_list = np.linspace(3.2, 3.3, 11)
scaled_positions = np.array([[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]])

E_list = np.zeros(len(a_list))
V_list = np.zeros(len(a_list))
for ind_a, a in enumerate(a_list):
    cell = np.array([[0, a, a], [a, 0, a], [a, a, 0]])
    atoms = Atoms(
        symbols="Sn2",
        scaled_positions=scaled_positions,
        cell=cell,
        pbc=True,
        calculator=deepcopy(calc),
    )
    E_list[ind_a] = atoms.get_potential_energy()
    V_list[ind_a] = atoms.get_volume()

print(atoms.positions)
e_min = np.min(E_list)
e_argmin = np.argmin(E_list)
print(e_min, e_argmin)
print("a", a_list[e_argmin])


from ase.units import kJ
from ase.eos import EquationOfState

V = V_list
E = E_list

eos = EquationOfState(V, E)
v0, e0, B = eos.fit()
print(B / kJ * 1.0e24, 'GPa')
