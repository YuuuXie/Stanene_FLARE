import sys
import numpy as np
import datetime
import time
from typing import List
import copy
import multiprocessing as mp
import subprocess
import concurrent.futures
from flare import struc, gp, env, qe_util, md, output
from flare.mff.utils import get_l_bound

class OTF(object):
    def __init__(self, qe_input: str, dt: float, number_of_steps: int,
                 gp: gp.GaussianProcess, pw_loc: str,
                 std_tolerance_factor: float = 1,
                 prev_pos_init: np.ndarray=None, par: bool=False,
                 skip: int=0, init_atoms: List[int]=None,
                 calculate_energy=False, output_name='otf_run.out',
                 backup_name='otf_run_backup.out',
                 max_atoms_added=None, freeze_hyps=False,
                 rescale_steps=[], rescale_temps=[], add_all=False,
                 no_cpus=1, 
                 use_mapping: bool=False, # below are added for mff
                 non_mapping_steps: list=[], 
                 l_bound: float=None, two_d: bool=False):

        self.qe_input = qe_input
        self.dt = dt
        self.number_of_steps = number_of_steps
        self.gp = gp
        self.pw_loc = pw_loc
        self.std_tolerance = std_tolerance_factor
        self.skip = skip
        self.dft_step = True
        self.freeze_hyps = freeze_hyps

        # parse input file
        positions, species, cell, masses = \
            qe_util.parse_qe_input(self.qe_input)

        self.structure = struc.Structure(cell=cell, species=species,
                                         positions=positions,
                                         mass_dict=masses,
                                         prev_positions=prev_pos_init)

        self.noa = self.structure.positions.shape[0]
        self.atom_list = list(range(self.noa))
        self.curr_step = 0

        # --------- params for mapped force field --------------
        # whether and when to begin using mapped force field
        self.use_mapping = use_mapping
        self.non_mapping_steps = non_mapping_steps
        self.two_d = two_d
        struc_l_bound = get_l_bound(100, self.structure, self.two_d)
        if l_bound:
            self.l_bound = np.min((l_bound, struc_l_bound))
        else:
            self.l_bound = struc_l_bound
        self.is_mff_built = False
        # -------------------------------------------------------

        if max_atoms_added is None:
            self.max_atoms_added = self.noa
        else:
            self.max_atoms_added = max_atoms_added

        # initialize local energies
        if calculate_energy:
            self.local_energies = np.zeros(self.noa)
        else:
            self.local_energies = None

        # set atom list for initial dft run
        if init_atoms is None:
            self.init_atoms = [int(n) for n in range(self.noa)]
        else:
            self.init_atoms = init_atoms

        self.dft_count = 0

        # set pred function
        if not par and not calculate_energy:
            self.pred_func = self.predict_on_structure
        elif par and not calculate_energy:
            self.pred_func = self.predict_on_structure_par
        elif not par and calculate_energy:
            self.pred_func = self.predict_on_structure_en
        elif par and calculate_energy:
            self.pred_func = self.predict_on_structure_par_en
        if self.use_mapping:
            self.mff = None
            if par:
                self.pred_func = self.predict_on_structure
                self.pred_func_gp = self.predict_on_structure
            else:
                self.pred_func = self.predict_on_structure
                self.pred_func_gp = self.predict_on_structure
        self.par = par

        # set rescale attributes
        self.rescale_steps = rescale_steps
        self.rescale_temps = rescale_temps

        self.output_name = output_name
        self.backup_name = backup_name
        self.last_backup_time = time.time()
        self.add_all = add_all

        # set number of cpus for qe runs
        self.no_cpus = no_cpus
        print('otf initialized, use mapping:', self.use_mapping)

    def run(self):
        output.write_header(self.gp.cutoffs, self.gp.kernel_name, self.gp.hyps,
                            self.gp.algo, self.dt, self.number_of_steps,
                            self.structure, self.output_name,
                            self.std_tolerance)
        counter = 0
        self.start_time = time.time()

        # relaunch mode
#        if self.curr_step > 0 and self.use_mapping:
#            self.train_mff(skip=False)

        while self.curr_step < self.number_of_steps:
            print('curr_step:', self.curr_step)
            # run DFT and train initial model if first step and DFT is on
            if self.curr_step == 0 and self.std_tolerance != 0:
                self.run_dft()
                dft_frcs = copy.deepcopy(self.structure.forces)
                new_pos = md.update_positions(self.dt, self.noa,
                                              self.structure)
                self.update_temperature(new_pos)
                self.record_state()

                # make initial gp model
                self.update_gp(self.init_atoms, dft_frcs)

                if not self.freeze_hyps:
                    self.train_gp()

                # build mapped force field
                if self.use_mapping and (not self.is_mff_built):
                    self.train_mff(skip=False)

                # check if remaining atoms are above uncertainty threshold
                if self.add_all:
                    std_in_bound = True
                else:
                    self.pred_func()
                    atom_list = self.init_atoms
                    std_in_bound, target_atom = self.is_std_in_bound(atom_list)

            # otherwise, try predicting with GP model
            else:
                if self.use_mapping:
                    if not self.mff: # self.mff == None
                        self.train_mff(skip=False)
                    if (self.curr_step-1) not in self.non_mapping_steps:
                        self.pred_func = self.predict_on_structure_mff
                    else:
                        if self.par: 
                            self.pred_func = self.predict_on_structure
                        else:
                            self.pred_func = self.predict_on_structure

                self.pred_func()
                self.dft_step = False
                new_pos = md.update_positions(self.dt, self.noa,
                                              self.structure)

                std_in_bound, target_atom = self.is_std_in_bound([])

                if std_in_bound: 
                    # check if the lower bound changes and mff needs re-build
                    if self.use_mapping:
                        if self.l_bound < self.mff.bounds_2[0,0]:
                            self.train_mff()
                        std_curr = np.max(self.structure.stds)
                        std_threshold = self.std_tolerance * np.abs(self.gp.hyps[-1])
                        if std_curr <= 0.8 * std_threshold:
                            if (not self.is_mff_built) and \
                                    (self.curr_step not in self.non_mapping_steps):
                                self.train_mff(skip=False)
                        else:
                            if not self.is_mff_built:
                                self.non_mapping_steps.append(self.curr_step)
                else:
                    atom_list = [target_atom]

                    # record GP forces
                    self.update_temperature(new_pos)
                    self.record_state()

                    # record DFT forces
                    self.dft_step = True
                    self.run_dft()
                    dft_frcs = copy.deepcopy(self.structure.forces)
                    new_pos = md.update_positions(self.dt, self.noa,
                                                  self.structure)
                    self.update_temperature(new_pos)
                    self.record_state()

                    # add atoms to training set until max error is below threshold
                    if self.add_all:
                        if not std_in_bound:
                            self.update_gp(self.atom_list, dft_frcs)
                            if not self.freeze_hyps:
                                self.train_gp()
 
                            self.is_mff_built = False
#                            if self.use_mapping:
#                                self.train_mff()
                            self.pred_func()
                    else:
                        atom_count = 0
                        while (not std_in_bound and atom_count <
                               self.max_atoms_added):
                            self.update_gp([target_atom], dft_frcs)
                            atom_list.append(target_atom)
                            if self.use_mapping:
                                self.pred_func_gp() # if use_mapping, then just use GP to predict here
                            else:
                                self.pred_func()
                            std_in_bound, target_atom = \
                                self.is_std_in_bound(atom_list)
                            atom_count += 1
    
                        if not self.freeze_hyps:
                            self.train_gp()
 
                        self.is_mff_built = False
                        self.non_mapping_steps.append(self.curr_step)
#                        if self.use_mapping:
#                            self.train_mff()
          
            # write gp forces only when counter equals skip
            if counter >= self.skip and not self.dft_step:
                self.update_temperature(new_pos)
                self.record_state()
                counter = 0

            counter += 1
            self.update_positions(new_pos)
            self.curr_step += 1

            # back up output every 1 hour
            if time.time() - self.last_backup_time >= 3600:
                subprocess.run(["cp", self.output_name, self.backup_name])
                self.last_backup_time = time.time()

        output.conclude_run(self.output_name)

    def predict_on_atom(self, atom):
        chemenv = env.AtomicEnvironment(self.structure, atom, self.gp.cutoffs)
        comps = []
        stds = []
        # predict force components and standard deviations
        for i in range(3):
            force, var = self.gp.predict(chemenv, i+1)
            comps.append(float(force))
            stds.append(np.sqrt(np.abs(var)))

        return comps, stds

    def predict_on_atom_en(self, atom):
        chemenv = env.AtomicEnvironment(self.structure, atom, self.gp.cutoffs)
        comps = []
        stds = []
        # predict force components and standard deviations
        for i in range(3):
            force, var = self.gp.predict(chemenv, i+1)
            comps.append(float(force))
            stds.append(np.sqrt(np.abs(var)))

        # predict local energy
        local_energy = self.gp.predict_local_energy(chemenv)
        return comps, stds, local_energy

    def predict_on_structure_par(self):
        n = 0
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for res in executor.map(self.predict_on_atom, self.atom_list):
                for i in range(3):
                    self.structure.forces[n][i] = res[0][i]
                    self.structure.stds[n][i] = res[1][i]
                n += 1

    def predict_on_structure_par_en(self):
        n = 0
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for res in executor.map(self.predict_on_atom_en, self.atom_list):
                for i in range(3):
                    self.structure.forces[n][i] = res[0][i]
                    self.structure.stds[n][i] = res[1][i]
                self.local_energies[n] = res[2]
                n += 1

    def predict_on_structure(self):
        for n in range(self.structure.nat):
            chemenv = env.AtomicEnvironment(self.structure, n, self.gp.cutoffs)
            for i in range(3):
                force, var = self.gp.predict(chemenv, i + 1)
                self.structure.forces[n][i] = float(force)
                self.structure.stds[n][i] = np.sqrt(np.abs(var))

    def predict_on_structure_en(self):
        for n in range(self.structure.nat):
            chemenv = env.AtomicEnvironment(self.structure, n, self.gp.cutoffs)
            for i in range(3):
                force, var = self.gp.predict(chemenv, i + 1)
                self.structure.forces[n][i] = float(force)
                self.structure.stds[n][i] = np.sqrt(np.abs(var))
            self.local_energies[n] = self.gp.predict_local_energy(chemenv)

    def run_dft(self):
        output.write_to_output('\nCalling Quantum Espresso...\n',
                               self.output_name)

        # calculate DFT forces
        forces = qe_util.run_espresso_par(self.qe_input, self.structure,
                                          self.pw_loc, self.no_cpus)
        self.structure.forces = forces

        # write wall time of DFT calculation
        self.dft_count += 1
        output.write_to_output('QE run complete.\n', self.output_name)
        time_curr = time.time() - self.start_time
        output.write_to_output('number of DFT calls: %i \n' % self.dft_count,
                               self.output_name)
        output.write_to_output('wall time from start: %.2f s \n' % time_curr,
                               self.output_name)

    def update_gp(self, train_atoms, dft_frcs):
        output.write_to_output('\nAdding atom {} to the training set.\n'
                               .format(train_atoms),
                               self.output_name)
        output.write_to_output('Uncertainty: {}.\n'
                               .format(self.structure.stds[train_atoms[0]]),
                               self.output_name)

        # update gp model
        self.gp.update_db(self.structure, dft_frcs,
                          custom_range=train_atoms)

        if self.curr_step == 0:
            self.gp.set_L_alpha()
        else:
            self.gp.update_L_alpha()

    def train_gp(self):
        self.gp.train()
        output.write_hyps(self.gp.hyp_labels, self.gp.hyps,
                          self.start_time, self.output_name,
                          self.gp.like, self.gp.like_grad)

    def is_std_in_bound(self, atom_list):
        # set uncertainty threshold
        if self.std_tolerance == 0:
            return True, -1
        elif self.std_tolerance > 0:
            threshold = self.std_tolerance * np.abs(self.gp.hyps[-1])
        else:
            threshold = np.abs(self.std_tolerance)

        # find max std
        max_std = 0
        for atom, std in enumerate(self.structure.stds):
            std_curr = np.max(std)

            if std_curr > max_std and atom not in atom_list:
                max_std = std_curr
                target_atom = atom

        # if above threshold, return atom
        if max_std > threshold:
            return False, target_atom
        else:
            return True, -1

    def update_positions(self, new_pos):
        if self.curr_step in self.rescale_steps:
            rescale_ind = self.rescale_steps.index(self.curr_step)
            temp_fac = self.rescale_temps[rescale_ind] / self.temperature
            vel_fac = np.sqrt(temp_fac)
            self.structure.prev_positions = \
                new_pos - self.velocities * self.dt * vel_fac
        else:
            self.structure.prev_positions = self.structure.positions
        self.structure.positions = new_pos
        self.structure.wrap_positions()
        
        # update mff lower bound
        self.l_bound = get_l_bound(self.l_bound, self.structure, self.two_d)

    def update_temperature(self, new_pos):
        KE, temperature, velocities = \
                md.calculate_temperature(new_pos, self.structure, self.dt,
                                         self.noa)
        self.KE = KE
        self.temperature = temperature
        self.velocities = velocities

    def record_state(self):
        output.write_md_config(self.dt, self.curr_step, self.structure,
                               self.temperature, self.KE,
                               self.local_energies, self.start_time,
                               self.output_name, self.dft_step,
                               self.velocities)
