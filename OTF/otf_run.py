import numpy as np
from copy import deepcopy

from ase import units
from ase.md.velocitydistribution import (
    MaxwellBoltzmannDistribution,
    Stationary,
    ZeroRotation,
)

from flare.ase.otf_md import otf_md
from flare.ase.logger import OTFLogger

import atom_setup, flare_setup, qe_setup

np.random.seed(12345)

super_cell = atom_setup.super_cell
nat = len(super_cell.positions)

flare_calc = deepcopy(flare_setup.flare_calc)
super_cell.set_calculator(flare_calc)

dft_calc = qe_setup.dft_calc

# set up OTF MD engine
md_engine = "VelocityVerlet"
dt = 1 * units.fs
md_params = {"timestep": dt, "trajectory": None, "dt": None}

atoms_added = 1
otf_params = {
    "dft_calc": dft_calc,
    "init_atoms": list(range(atoms_added)),
    "std_tolerance_factor": 1.0,
    "max_atoms_added": atoms_added,
    "freeze_hyps": 0,
    "restart_from": None,
    "use_mapping": super_cell.calc.use_mapping,
    "non_mapping_steps": [i for i in range(100)]
    + [i + 20000 for i in range(200)]
    + [i + 40000 for i in range(500)]
    + [i + 70000 for i in range(500)],
}

temperature = 200

# jitter atomic positions
xi = np.random.standard_normal((nat, 3))
kT = temperature * units.kB
init_v = np.sqrt(kT / super_cell.get_masses())[:, np.newaxis] * xi
delta_pos = init_v * dt
super_cell.positions += delta_pos

# intialize velocity
MaxwellBoltzmannDistribution(super_cell, kT)
Stationary(super_cell)  # zero linear momentum
ZeroRotation(super_cell)  # zero angular momentum

test_otf = otf_md(md_engine, super_cell, md_params, otf_params)

# set up logger
otf_logger = OTFLogger(
    test_otf, super_cell, logfile="otf_train.log", mode="w", data_in_logfile=True
)
test_otf.attach(otf_logger, interval=1)

# run otf
number_of_steps = 100000
rescale_temp = [600, 800, 1300]
rescale_steps = [20000, 40000, 70000]
test_otf.otf_run(number_of_steps, rescale_temp, rescale_steps)
