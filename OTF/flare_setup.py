import numpy as np

from flare.gp import GaussianProcess
from flare.mgp.mgp import MappedGaussianProcess
from flare.ase.calculator import FLARE_Calculator

# ---------- create gaussian process model -------------------
hyps = np.array([0.1, 1.0, 0.001, 1, 0.05])
two_cut = 7.2
three_cut = 7.2
cutoffs = np.array([two_cut, three_cut])
hyp_labels = ["sig2", "ls2", "sig3", "ls3", "noise"]
opt_algorithm = "BFGS"

gp_model = GaussianProcess(
    kernel_name="2+3",
    hyps=hyps,
    hyp_labels=hyp_labels,
    cutoffs=cutoffs,
    opt_algorithm=opt_algorithm,
    parallel=True,
    n_cpus=32,
)

# ----------- create mapped gaussian process ------------------
struc_params = {
    "species": [50],
    "cube_lat": np.eye(3) * 100,
    "mass_dict": {"0": 118.71},
}

# grid parameters
lower_cut = 2.35
grid_num_2 = 128
grid_num_3 = 96
grid_params = {
    "bounds_2": [[lower_cut], [two_cut]],
    "bounds_3": [[lower_cut, lower_cut, -1], [three_cut, three_cut, 1]],
    "grid_num_2": grid_num_2,
    "grid_num_3": [grid_num_3, grid_num_3, grid_num_3],
    "svd_rank_2": 0,
    "svd_rank_3": 0,
    "bodies": [2, 3],
    "load_grid": None,
    "update": False,
}

mgp_model = MappedGaussianProcess(
    gp_model.hyps,
    gp_model.cutoffs,
    grid_params,
    struc_params,
    mean_only=False,
    container_only=False,
    GP=gp_model,
    lmp_file_name="lmp.mgp",
)

# ------------ create ASE's flare calculator -----------------------
flare_calc = FLARE_Calculator(gp_model, mgp_model, par=True, use_mapping=True)
