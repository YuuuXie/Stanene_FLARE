import numpy as np
import os
from copy import deepcopy

from ase import units
from ase.md.velocitydistribution import (
    MaxwellBoltzmannDistribution,
    Stationary,
    ZeroRotation,
)

from flare.ase.otf import ASE_OTF
import atoms_setup, flare_setup, dft_setup

def init_md(super_cell, temperature):
    # set up MD engine
    MaxwellBoltzmannDistribution(super_cell, temperature * units.kB)
    Stationary(super_cell)  # zero linear momentum
    ZeroRotation(super_cell)  # zero angular momentum



print("modules loaded")

temperature = 200
super_cell = atoms_setup.get_atoms(temperature, multiplier=[2,2,1])

flare_calc = flare_setup.get_flare_calc()
super_cell.set_calculator(flare_calc)

#n_cpus = 256
#npool = 16
#kpts = (4, 4, 4)
n_cpus = 32
npool = 1
kpts = [1, 1, 1]
dft_calc = dft_setup.get_dft_calc(n_cpus, npool, kpts)

init_md(super_cell, temperature)
md_engine = "VelocityVerlet"
md_kwargs = {}

# set up OTF engine
N_steps = 3
atoms_added = 8
otf_params = {
    "init_atoms": [i for i in range(atoms_added)],
    "output_name": "otf",
    "std_tolerance_factor": 2.0,
    "max_atoms_added": atoms_added,
    "freeze_hyps": 100,
    "write_model": 3,
    "rescale_temps": [600, 800, 1300],
    "rescale_steps": [10000, 20000, 30000],
    "store_dft_output": [["scf.pwi", "scf.pwo"], "dft"],
}

test_otf = ASE_OTF(
    super_cell,
    timestep=1 * units.fs,
    number_of_steps=N_steps,
    dft_calc=dft_calc,
    md_engine=md_engine,
    md_kwargs=md_kwargs,
    **otf_params
)

print("otf run!")
test_otf.run()
