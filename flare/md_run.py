from flare.gp import GaussianProcess
from flare.struc import Structure
from flare.env import AtomicEnvironment
import numpy as np
import time
import datetime
import concurrent.futures
from flare import md
from flare import output


class MD:
    """Generates NVE dynamics from a GP model."""

    def __init__(self, dt: float, number_of_steps: int, gp: GaussianProcess,
                 pos_init: np.ndarray, species, cell, masses,
                 prev_pos_init: np.ndarray=None, par: bool=False, skip: int=0,
                 output_name='otf_run.out'):

        self.dt = dt
        self.Nsteps = number_of_steps
        self.gp = gp

        self.structure = Structure(cell=cell, species=species,
                                   positions=pos_init,
                                   mass_dict=masses,
                                   prev_positions=prev_pos_init)

        self.noa = self.structure.positions.shape[0]
        self.atom_list = list(range(self.noa))
        self.curr_step = 0

        # choose prediction function
        if par is True:
            self.pred_func = self.predict_on_structure_par_en
        else:
            self.pred_func = self.predict_on_structure_en

        # initialize local energies
        self.local_energies = np.zeros(self.noa)

        self.pes = []
        self.kes = []

        self.output_name = output_name

    def run(self):
        output.write_header(self.gp.cutoffs, self.gp.kernel_name, self.gp.hyps,
                            self.gp.algo, self.dt, self.Nsteps, self.structure,
                            self.output_name)
        self.start_time = time.time()

        while self.curr_step < self.Nsteps:
            # verlet algorithm follows Frenkel p. 70
            self.pred_func()
            new_pos = md.update_positions(self.dt, self.noa, self.structure)
            self.update_temperature(new_pos)
            self.record_state()
            self.update_positions(new_pos)
            self.curr_step += 1

        output.conclude_run(self.output_name)

    def predict_on_structure_en(self):
        for n in range(self.structure.nat):
            chemenv = AtomicEnvironment(self.structure, n, self.gp.cutoffs)
            for i in range(3):
                force, var = self.gp.predict(chemenv, i + 1)
                self.structure.forces[n][i] = float(force)
                self.structure.stds[n][i] = np.sqrt(np.absolute(var))
            self.local_energies[n] = self.gp.predict_local_energy(chemenv)

    def predict_on_structure_par_en(self):
        n = 0
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for res in executor.map(self.predict_on_atom_en, self.atom_list):
                for i in range(3):
                    self.structure.forces[n][i] = res[0][i]
                    self.structure.stds[n][i] = res[1][i]
                self.local_energies[n] = res[2]
                n += 1
        self.structure.dft_forces = False

    def predict_on_atom_en(self, atom):
        chemenv = AtomicEnvironment(self.structure, atom, self.gp.cutoffs)
        comps = []
        stds = []
        # predict force components and standard deviations
        for i in range(3):
            force, var = self.gp.predict(chemenv, i+1)
            comps.append(float(force))
            stds.append(np.sqrt(np.absolute(var)))

        # predict local energy
        local_energy = self.gp.predict_local_energy(chemenv)
        return comps, stds, local_energy

    def update_positions(self, new_pos):
        self.structure.prev_positions = self.structure.positions
        self.structure.positions = new_pos
        self.structure.wrap_positions()

    def update_temperature(self, new_pos):
        KE, temperature = \
                md.calculate_temperature(new_pos, self.structure, self.dt,
                                         self.noa)
        self.KE = KE
        self.temperature = temperature

    def record_state(self):
        self.pes.append(np.sum(self.local_energies))
        self.kes.append(self.KE)
        output.write_md_config(self.dt, self.curr_step, self.structure,
                               self.temperature, self.KE, self.local_energies,
                               self.start_time, self.output_name)
