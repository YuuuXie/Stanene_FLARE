import numpy as np
import subprocess
import os


class EAM_Force_Calculator:
    def __init__(self, struc, style_string, coeff_string, lammps_folder,
                 lammps_executable):
        self.struc = struc
        self.style_string = style_string
        self.coeff_string = coeff_string
        self.lammps_folder = lammps_folder
        self.lammps_executable = lammps_executable

        self.input_file = lammps_folder + '/tmp.in'
        self.output_file = lammps_folder + '/tmp.out'
        self.dat_file = lammps_folder + '/tmp.data'
        self.dump_file = lammps_folder + '/tmp.dump'

        self.input_text = self.lammps_input()
        self.dat_text = self.lammps_dat()

    def get_forces(self):
        self.lammps_generator()
        self.run_ewald()
        forces = self.lammps_parser()
        return forces

    def run_ewald(self):
        # create input and data files
        self.lammps_generator()

        # run lammps
        lammps_command = '%s < %s > %s' % (self.lammps_executable,
                                           self.input_file,
                                           self.output_file)
        os.system(lammps_command)

    def lammps_generator(self):
        self.write_file(self.input_file, self.input_text)
        self.write_file(self.dat_file, self.dat_text)

    def lammps_parser(self):
        forces = []

        with open(self.dump_file, 'r') as outf:
            lines = outf.readlines()

        for count, line in enumerate(lines):
            if line.startswith('ITEM: ATOMS'):
                force_start = count

        for line in lines[force_start+1:]:
            fline = line.split()
            forces.append([float(fline[-3]),
                          float(fline[-2]),
                          float(fline[-1])])

        return np.array(forces)

    def lammps_input(self):
        input_text = """# lammps input file created with eam.py.
units metal
atom_style atomic
dimension  3
boundary   p p p
read_data %s

pair_style %s
pair_coeff %s

thermo_style one
dump 1 all custom 1 %s id type x y z fx fy fz
dump_modify 1 sort id
run 0
""" % (self.dat_file, self.style_string, self.coeff_string, self.dump_file)

        return input_text

    def lammps_dat(self):
        dat_text = """Header of the LAMMPS data file

%i atoms
1 atom types
""" % (self.struc.nat)

        dat_text += self.lammps_cell_text()

        dat_text += """
Masses

1  1

Atoms
"""

        dat_text += self.lammps_pos_text()

        return dat_text

    def lammps_cell_text(self):
        """ Write cell from structure object. Assumes orthorombic periodic
        cell."""

        cell_text = """
0.0 %f  xlo xhi
0.0 %f  ylo yhi
0.0 %f  zlo zhi
%f %f %f  xy xz yz
""" % (self.struc.cell[0, 0],
            self.struc.cell[1, 1],
            self.struc.cell[2, 2],
            self.struc.cell[1, 0],
            self.struc.cell[2, 0],
            self.struc.cell[2, 1])

        return cell_text

    def lammps_pos_text(self):
        pos_text = '\n'
        for count, pos in enumerate(self.struc.positions):
            pos_text += '%i 1 %f %f %f \n' % \
                (count+1, pos[0], pos[1], pos[2])
        return pos_text

    @staticmethod
    def write_file(file, text):
        with open(file, 'w') as fin:
            fin.write(text)
