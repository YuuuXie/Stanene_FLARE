import numpy as np
from ase import Atoms
from ase.build import bulk, make_supercell, add_adsorbate
from ase.calculators.espresso import Espresso
from ase.lattice.hexagonal import Graphene
from ase import units

def get_atoms(temperature, multiplier):
    symbol = "Sn"
    a = 4.581579042
    c = 20.0  # vaccum
    d0 = 0.415267584
    unit_cell = Graphene(symbol, latticeconstant={"a": a, "c": c})
    unit_cell.positions[0][2] -= d0
    unit_cell.positions[1][2] += d0
    
    # make super cell
    multiplier = np.diag(multiplier)
    super_cell = make_supercell(unit_cell, multiplier)
    
    # jitter atomic positions
    nat = len(super_cell)
    dt = 1 * units.fs
    xi = np.random.standard_normal((nat, 3))
    kT = temperature * units.kB
    init_v = np.sqrt(kT / super_cell.get_masses())[:, np.newaxis] * xi
    delta_pos = init_v * dt
    super_cell.positions += delta_pos
    return super_cell


if __name__ == "__main__":
    temp = 200
    multiplier = [4, 4, 1]
    atoms = get_atoms(temp, multiplier)
