import os
from flare.gp import GaussianProcess
from flare.mgp import MappedGaussianProcess
from flare.utils.parameter_helper import ParameterHelper
from flare.ase.calculator import FLARE_Calculator
from flare.output import Output, set_logger


def get_flare_calc():

    # set up GP hyperparameters
    kernels = ["twobody", "threebody"]
    parameters = {"cutoff_twobody": 7.2, "cutoff_threebody": 6.0}
    pm = ParameterHelper(kernels=kernels, random=True, parameters=parameters)

    hm = pm.as_dict()
    hyps = hm["hyps"]
    cut = hm["cutoffs"]
    print("hyps", hyps)

    gp_model = GaussianProcess(
        kernels=kernels,
        component="mc",
        hyps=hyps,
        cutoffs=cut,
        hyp_labels=["sig2", "ls2", "sig3", "ls3", "noise"],
        opt_algorithm="L-BFGS-B",
        parallel=True,
        per_atom_par=True,
        n_cpus=32,
    )

    grid_params = {
        "twobody": {"grid_num": [128]},
        "threebody": {"grid_num": [96, 96, 96]},
    }
    mgp_model = MappedGaussianProcess(
        grid_params, unique_species=[50], n_cpus=32, var_map="pca",
    )

    # set up flare calculator
    flare_calc = FLARE_Calculator(
        gp_model, mgp_model=mgp_model, par=True, use_mapping=True
    )
    return flare_calc
