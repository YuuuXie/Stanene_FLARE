import numpy as np
from ase import Atoms
from ase.build import bulk, make_supercell, add_adsorbate
from ase.lattice.hexagonal import Graphene

symbol = 'Sn'
a = 4.581579042
c = 20.0 # vaccum
d0 = 0.415267584 # buckling
unit_cell = Graphene(symbol, latticeconstant={'a':a,'c':c})
unit_cell.positions[1][2] += 2 * d0
multiplier = np.array([[4,0,0],[0,4,0],[0,0,1]])
super_cell = make_supercell(unit_cell, multiplier)
super_cell.positions[:, 2] -= d0 

positions = []
for _ in range(10):
    pos = super_cell.positions + 0.02 * (2 * np.random.rand(*super_cell.positions.shape) - 1)
    positions.append(pos)

np.save('test_frames_32atoms', positions)
