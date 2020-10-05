import sys
sys.path.append('..')
from flare.struc import Structure

import numpy as np
from ase.calculators.espresso import Espresso
from ase.calculators.eam import EAM
from ase.md.npt import NPT

class UQ_NPT(NPT):

    def __init__(self, atoms, timestep, temperature, 
            externalstress, ttime, pfactor, mask=None, 
            trajectory=None, logfile=None, loginterval=1,
            # on-the-fly parameters
            std_tolerance_factor=None, max_atoms_added=1,
            freeze_hyps=False, dft_input={}):

        super().__init__(atoms, timestep, temperature, 
                externalstress, ttime, pfactor, mask=None, 
                trajectory=None, logfile=None, loginterval=1)

        self.std_tolerance = std_tolerance_factor
        self.max_atoms_added = max_atoms_added
        self.freeze_hyps = freeze_hyps
        self.dft_input = dft_input

        # initialize gp
        self.std_in_bound = False
        self.target_atom = 0
        self.stds = []
        dft_forces = self.call_DFT()
        self.update_GP(dft_forces)

#    def forcecalculator(self):
#        forces, stds = self.atoms.get_forces()
#        self.stds.append(stds)
#        return forces 

    def is_std_in_bound(self, atom_list):
        # set uncertainty threshold
        if self.std_tolerance == 0:
            self.std_in_bound = True
            self.target_atom = -1
        elif self.std_tolerance > 0:
            noise = np.abs(self.atoms.calc.gp_model.hyps[-1])
            threshold = self.std_tolerance * noise
        else:
            threshold = np.abs(self.std_tolerance)

        # find max std
        max_std = 0
        for atom, std in enumerate(self.stds):
            std_curr = np.max(std)

            if std_curr > max_std:
                max_std = std_curr
                target_atom = atom

        # if above threshold, return atom
        if max_std > threshold and atom not in atom_list:
            self.std_in_bound = False
            self.target_atom = target_atom
        else:
            self.std_in_bound = True
            self.target_atom = -1

    def run(self, steps):
        """Perform a number of time steps."""
        if not self.initialized:
            self.initialize()
        else:
            if self.have_the_atoms_been_changed():
                raise NotImplementedError(
                    "You have modified the atoms since the last timestep.")

        for i in range(steps):
 #           self.stds = []
            print('step:', i)
            self.stds = np.zeros((32,3))
            self.step()
            self.nsteps += 1

            # figure out if std above the threshold
            #self.attach(self.is_std_in_bound, 1, [[]])
            self.call_observers() 

            # call dft/eam
            print('calling dft')
            dft_forces = self.call_DFT()

            # update gp
            print('updating gp')
            self.update_GP(dft_forces)
    
    def call_DFT(self):
        prev_calc = self.atoms.calc
        pseudopotentials = self.dft_input['pseudopotentials']
#        label = self.dft_input['label']
#        input_data = self.dft_input['input_data']
#        kpts = self.dft_input['kpts']
#        calc = Espresso(pseudopotentials=pseudopotentials, label=label, 
#                        tstress=True, tprnfor=True, nosym=True, 
#                        input_data=input_data, kpts=kpts) 
        calc = EAM(potential=pseudopotentials)        
        self.atoms.set_calculator(calc)
        forces = self.atoms.get_forces()
        self.atoms.set_calculator(prev_calc)
        return forces

    def update_GP(self, dft_forces):
        atom_count = 0
        atom_list = []
        calc = self.atoms.calc.gp_model
        while (not self.std_in_bound and atom_count <
               self.max_atoms_added):
            # build gp structure from atoms
            atom_struc = Structure(self.atoms.cell, 
                    ['A']*len(self.atoms.positions), 
                    self.atoms.positions)

            # update gp model
            calc.update_db(atom_struc, dft_forces,
                           custom_range=[self.target_atom])
    
            if not calc.alpha:
                calc.set_L_alpha()
            else:
                calc.update_L_alpha()

            atom_list.append(self.target_atom)
            self.forcecalculator()
            self.is_std_in_bound(atom_list)
            atom_count += 1

        if not self.freeze_hyps:
            calc.train()

