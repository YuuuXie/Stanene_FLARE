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

a_list = np.linspace(2.85, 2.95, 11)
c_list = np.linspace(1.55, 1.65, 11)
# a, c, c/a 3.028 1.65 0.5449141347424042
scaled_positions = np.array([[0.0, 0.0, 0.0], [0.25, 0.75, 0.5]])

E_list = np.zeros((len(a_list), len(c_list)))
V_list = np.zeros((len(a_list), len(c_list)))
for ind_a, a in enumerate(a_list):
    for ind_c, c in enumerate(c_list):
        cell = np.array([[-a, a, c], [a, -a, c], [a, a, -c]])
        atoms = Atoms(
            symbols="Sn2",
            scaled_positions=scaled_positions,
            cell=cell,
            pbc=True,
            calculator=deepcopy(calc),
        )
        E_list[ind_a, ind_c] = atoms.get_potential_energy()
        V_list[ind_a, ind_c] = atoms.get_volume()

e_min = np.min(E_list)
e_argmin = np.unravel_index(np.argmin(E_list, axis=None), E_list.shape)
print(e_min, e_argmin)
a_min = a_list[e_argmin[0]]
c_min = c_list[e_argmin[1]]
print("a, c, c/a", a_min, c_min, c_min / a_min)


from ase.units import kJ
from ase.eos import EquationOfState

V = V_list.ravel()
E = E_list.ravel()

eos = EquationOfState(V, E)
v0, e0, B = eos.fit()
print(B / kJ * 1.0e24, 'GPa')
