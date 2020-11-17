import numpy as np
from time import time
import os
import multiprocessing as mp
import pickle

from flare.mgp.mgp import MappedGaussianProcess
from flare.mgp.utils import get_l_bound
from flare.env import AtomicEnvironment
from flare.struc import Structure
from flare.otf_parser import OtfAnalysis
from flare import kernels
from flare.dft_interface.qe_util import parse_dft_forces_and_energy

def build_gp(filename, crop_train_set):
    gp = pickle.load(open(filename, 'rb'))
    gp.training_data = gp.training_data[:crop_train_set]
    gp.training_labels = gp.training_labels[:crop_train_set]
    gp.training_labels_np = np.hstack(gp.training_labels[:crop_train_set])
    gp.n_cpus = 32
    gp.par = True
    gp.per_atom_par = True
    gp.set_L_alpha()
    gp.write_model(f'gp_train{crop_train_set}', format='pickle')

def build_mgp(filename, load_grid=None):
    # load trained gp model
    t0 = time()
    gp_model = pickle.load(open(filename, 'rb'))
    cutoffs = gp_model.cutoffs
    print('gp loaded:', time()-t0)

    # mgp params
    train_size = len(gp_model.training_data)
    grid_num_2 = 144
    grid_num_3 = 128
    grid_params = {'bounds_2': np.array([[2.35], [cutoffs[0]]]),
                   'bounds_3': np.array([[2.35, 2.35, -1], 
                                         [cutoffs[1], cutoffs[1], 1]]), 
                   'grid_num_2': grid_num_2, 
                   'grid_num_3': [grid_num_3 for i in range(3)], 
                   'svd_rank_2': np.min([grid_num_2, train_size*3]), 
                   'svd_rank_3': np.min([grid_num_3**3, train_size*3]),
                   'bodies': [2, 3], 
                   'load_grid': load_grid, 
                   'update': False}
    
    struc_params = {'species': [50], 
                    'cube_lat': 100 * np.eye(3), 
                    'mass_dict': {'0': 12}}
    
    # build mgp
    t0 = time()
    mgp_model = MappedGaussianProcess(gp_model.hyps, gp_model.cutoffs,
                grid_params, struc_params, mean_only=False, container_only=False,
                GP=gp_model, lmp_file_name='tmp.mgp')
    print('map built:', time()-t0)
    
    return mgp_model, gp_model


if __name__ == "__main__":
    build_mgp('stn_3_trained.gp.pickle')
