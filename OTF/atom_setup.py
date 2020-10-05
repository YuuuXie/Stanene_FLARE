import numpy as np
from ase import Atoms
from ase.build import bulk, make_supercell, add_adsorbate
from ase.calculators.espresso import Espresso
from ase.lattice.hexagonal import Graphene

symbol = "Sn"
a = 4.675238246
c = 20.0  # vaccum
d0 = 0.425705071
unit_cell = Graphene(symbol, latticeconstant={"a": a, "c": c})
unit_cell.positions[0][2] -= d0
unit_cell.positions[1][2] += d0

# make super cell
multiplier = np.array([[4, 0, 0], [0, 4, 0], [0, 0, 1]])
super_cell = make_supercell(unit_cell, multiplier)
