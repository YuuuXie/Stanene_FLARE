import numpy as np
from time import time
import os
import multiprocessing as mp
import sys

sys.path.append('../build_mgp')
from build_mgp import build_mgp

sys.path.append('../lmp')
from lmp_2_xyz import parse_lmp

from flare.mgp.mgp import MappedGaussianProcess
from flare.mgp.utils import get_l_bound
from flare.env import AtomicEnvironment
from flare.struc import Structure
from flare.otf_parser import OtfAnalysis
from flare import kernels
from flare.dft_interface.qe_util import parse_dft_forces_and_energy


def test_uncertainties(traj_file, gp_model, mgp_model):
    t0 = time()
    positions_list, cell = parse_lmp(traj_file, skip=10)
    f_mgp = []
    v_mgp = []
    for frame in positions_list:
        nat = len(frame)
        species = [50 for i in range(nat)]
        struc_curr = Structure(cell, species, frame)
        forces = np.zeros((nat, 3))
        stds = np.zeros((nat, 3))
        for n in range(nat):
            chemenv = AtomicEnvironment(struc_curr, n, gp_model.cutoffs)
            f, v = mgp_model.predict(chemenv, mean_only=False)
            forces[n] = f
            stds[n] = np.sqrt(np.absolute(v))
    
        f_mgp.append(forces)
        v_mgp.append(stds)
    
    print('chemenv species', chemenv.species)
    np.save('mgp_f', np.array(f_mgp))
    np.save('mgp_v', np.array(v_mgp))
    print('\nprediction time:', time()-t0)


# ------------------------- test structure -----------------------------------
if __name__ == "__main__":

    traj_file = '../../mgp-1/lmp/traj.lammps'
    mgp_model, gp_model = build_mgp('../train_hyps/stn_2_trained.gp.pickle', 
                                    '../build_mgp/')
    test_uncertainties(traj_file, gp_model, mgp_model)
