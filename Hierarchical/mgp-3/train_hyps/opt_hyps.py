import numpy as np
from time import time
import os
import multiprocessing as mp
import pickle

from flare.mgp.mgp import MappedGaussianProcess
from flare.env import AtomicEnvironment
from flare.gp import GaussianProcess
from flare.struc import Structure
from flare.output import Output
from flare.otf_parser import OtfAnalysis
from flare import kernels
from flare.dft_interface.qe_util import parse_dft_forces_and_energy
from flare.output import Output


def parse_gp():
    kernel = kernels.two_plus_three_body
    kernel_grad = kernels.two_plus_three_body_grad
    hyps = np.array([0.2, 1., 0.0001, 1, 0.05]) 
    cutoffs = np.array([7.2, 7.2])
    hyp_labels = ['sig2', 'ls2', 'sig3', 'ls3', 'noise']
     
    # parse traj file
    t0 = time()
    old_otf = OtfAnalysis('../../../../../../../OTF_train/otf_run_100ps.out')
    print('old otf ready:', time()-t0)
    
    # construct gp
    t0 = time()
    call_no = len(old_otf.gp_position_list)
    gp_model = old_otf.make_gp(kernel=kernel,
                               kernel_grad=kernel_grad,
                               call_no=call_no, 
                               cutoffs=cutoffs, 
                               hyps=hyps)
    train_size = len(gp_model.training_data)
    print('hyps:', gp_model.hyps, 'likelihood:', gp_model.likelihood)
    print('gp ready:', time()-t0)
    print('training size from old: {}, min dist: {}'.format(train_size, 'unknown'))
    print('training data species:', gp_model.training_data[0].species)

    # save gp model
    gp_model.par = True
    gp_model.no_cpus = 32
    gp_model.ncpus = 32
    gp_model.maxiter = 40
    gp_model.write_model('stn_1_original.gp', format='pickle')
    
    return gp_model


def train_gp(gp_model):
    # train gp model
    t0 = time()
    logfile = Output('train', always_flush=True)
    gp_model.train(output=logfile)
    gp_model.write_model('stn_1_trained.gp', format='pickle')
    print('training time:', time()-t0)
    print('gp hyps:', gp_model.hyps)
    print('gp likelihood:', gp_model.likelihood)

    return gp_model


def load_gp(filename):
    # load trained gp model
    t0 = time()
    gp_model = pickle.load(open(filename, 'rb'))
    print('gp loaded:', time() - t0)
    return gp_model


if __name__ == "__main__":
#    parse_gp()
#    gp_model = load_gp('../add_snap/stn_3_original.gp.pickle')
    gp_model = load_gp('stn_3_trained.gp.pickle')
    train_gp(gp_model)
