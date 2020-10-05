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
from flare.util import is_std_in_bound
from flare.predict import predict_on_structure_par
from flare.dft_interface.qe_util import parse_dft_forces_and_energy,\
        dft_input_to_structure


def parse_xyz(filename):
    f = open(filename)
    lines = f.readlines()
    nat = int(lines[0].strip())
    labels = np.zeros(nat)
    positions_list = []
    for l,line in enumerate(lines):
        line = line.split()
        if len(line) <= 1:
            positions = np.zeros((nat, 3))
        else:
            atom = l % (nat+2) - 2
            positions[atom] = np.array([float(line[i]) for i in range(1, 4)])
            if atom == nat-1:
                positions_list.append(positions)
    return positions_list

def add_atoms_to_train(xyz_file, pwo_file, add_list, test_var, cell, gp_model,
                       max_atoms_added):

    positions = parse_xyz(xyz_file)[0]
    nat = len(positions)
    species = [50 for i in range(nat)]
    struc_curr = Structure(cell, species, positions)
    dft_forces, _ = parse_dft_forces_and_energy(pwo_file)
    dft_forces = dft_forces[:struc_curr.nat]
    gp_model.l_bound = get_l_bound(10, struc_curr, True)
    
    # add atoms to gp
    t0 = time()
    added_atoms = add_list 
    if test_var:
        # test uncertainty
        _, stds = predict_on_structure_par(struc_curr, gp_model)
        in_bound, added_atoms = is_std_in_bound(0.1, 
                                                gp_model.hyps[-1], 
                                                struc_curr,
                                                max_atoms_added)

        print('Add atoms', added_atoms)
    
    else:
        # add to database if uncertainty is high
        print('Add atom', atom) 


    for atom in added_atoms:
        print('Add atom', atom) 
        atom = int(atom)
        env_curr = AtomicEnvironment(struc_curr, atom, gp_model.cutoffs)
        forces_curr = np.array(dft_forces[atom])
        gp_model.training_data.append(env_curr)
        gp_model.training_labels.append(forces_curr)
        gp_model.training_labels_np = gp_model.force_list_to_np(gp_model.training_labels)
        gp_model.update_L_alpha() 
    
    print('gp from snap {} ready: {}'.format(snap, time()-t0))
    print('training size from old: {}, min dist: {}'.format(len(gp_model.training_data), gp_model.l_bound))


if __name__ == '__main__':
    # load gp model
    gp_model = pickle.load(open('../../mgp-1/train_hyps/stn_1_trained.gp.pickle', 'rb'))

    # training dft data
    path = '../../mgp-1/dft/'
    scf_pwi = [path+'1000/1000ps.xyz', path+'1980/1980ps.xyz', path+'2900/2900ps.xyz']
    scf_pwo = [path+'1000/scf.pwo', path+'1980/scf.pwo', path+'2900/scf.pwo']
    cell = np.array([[ 39.6776,  0.00000,  0.0000],
                     [  0.0000, 45.81580,  0.0000],
                     [  0.0000,  0.00000, 50.0000]])

    for snap in range(len(scf_pwi)):
        add_list = np.arange(200)
        test_var = True
        add_atoms_to_train(scf_pwi[snap], scf_pwo[snap], 
                add_list, test_var, cell, gp_model, 30)
    
    # save gp model 
    gp_model.write_model('stn_2_original.gp', format='pickle')
