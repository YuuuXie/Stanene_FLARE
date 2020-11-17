import numpy as np
from time import time
import os
import multiprocessing as mp
import sys

from build_mgp import build_mgp, build_gp

from flare.mgp.mgp import MappedGaussianProcess
from flare.mgp.utils import get_l_bound
from flare.env import AtomicEnvironment
from flare.struc import Structure
from flare.otf_parser import OtfAnalysis
from flare import kernels
from flare.dft_interface.qe_util import parse_dft_forces_and_energy
from flare.predict import predict_on_structure_par

def test_uncertainties(traj_file, gp_model, mgp_model):
    positions_list = np.load(traj_file)
    t_mgp = []
    t_mgp_v = []
    t_gp = []
    t_env = []

    train_size = len(gp_model.training_data)
    test_size = nat = len(positions_list[0])

    if test_size == 32:
        # cell for 32 atoms
        cell = np.array([[18.326316168, 0.0, 0.0],
                         [-9.163158084, 15.87105535927349, 0.0],
                         [0.0, 0.0, 20.0]])
    elif test_size == 200:
        # cell for 200 atoms
        cell = np.array([[39.6776,  0.0   ,  0],
                         [ 0.0   , 45.8158,  0],
                         [ 0.0   ,  0.0   ,  50]])
    n = 0
    t_frame = time()
    for frame in positions_list:
        nat = len(frame)
        species = [50 for i in range(nat)]
        struc_curr = Structure(cell, species, frame)
        for n in range(nat):
            t0 = time()
            chemenv = AtomicEnvironment(struc_curr, n, gp_model.cutoffs)
            t_env.append(time()-t0)

            t0 = time()
            mgp_model.predict(chemenv, mean_only=True)
            t_mgp.append(time()-t0)

            t0 = time()
            mgp_model.predict(chemenv, mean_only=False)
            t_mgp_v.append(time()-t0)
  
        t0 = time()
        predict_on_structure_par(struc_curr, gp_model)
        t_gp.append(time()-t0)
        n += 1
        print(train_size, test_size, n, 'predicted:', time()-t_frame)
 
    np.savetxt(f't_mgp-{train_size}-{test_size}', np.array(t_mgp))
    np.savetxt(f't_mgp_v-{train_size}-{test_size}', np.array(t_mgp_v))
    np.savetxt(f't_gp-{train_size}-{test_size}', np.array(t_gp))
    np.savetxt(f't_env-{train_size}-{test_size}', np.array(t_env))


# ------------------------- test structure -----------------------------------
if __name__ == "__main__":
    train_size = int(sys.argv[1])
    test_sizes = [32, 200]
    gp_model = build_gp('stn_3_trained.gp.pickle', train_size)
    mgp_model, gp_model = build_mgp(f'gp_train{train_size}.pickle', None) 
    for test_size in test_sizes:
        traj_file = f'test_frames_{test_size}atoms.npy'
        test_uncertainties(traj_file, gp_model, mgp_model)
