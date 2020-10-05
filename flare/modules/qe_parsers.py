import numpy as np

# unit conversion from http://greif.geo.berkeley.edu/~driver/conversions.html
en_conv = 13.6056925330  # bohr to eV
force_conv = 25.71104309541616  # bohr/Ry to eV/A
stress_conv = 25.71104309541616 / (0.529177208**2)


# get energy from scf run
def parse_scf_energy(outfile):
    with open(outfile, 'r') as outf:
        lines = outf.readlines()

    for line in lines:
        if '!' in line:
            energy = float(line.split()[-2])

    return energy


def parse_scf(outfile: str):
    """
    Get forces from a pwscf file in eV/A

    :param outfile: str, Path to pwscf output file
    :return: list[nparray] , List of forces acting on atoms
    """
    forces = []
    total_energy = np.nan

    with open(outfile, 'r') as outf:
        for line in outf:
            if line.lower().startswith('!    total energy'):
                total_energy = float(line.split()[-2])

            if line.find('force') != -1 and line.find('atom') != -1:
                line = line.split('force =')[-1]
                line = line.strip()
                line = line.split(' ')
                line = [x for x in line if x != '']
                temp_forces = []
                for x in line:
                    temp_forces.append(float(x))
                forces.append(np.array(list(temp_forces)))

    assert total_energy != np.nan, "Quantum ESPRESSO parser failed to read " \
                                   "the file {}. Run failed.".format(outfile)

    total_energy = total_energy * en_conv
    forces = [force_conv * force for force in forces]

    return total_energy, forces


def parse_md_output(outfile):

    steps = {}

    with open(outfile, 'r') as outf:
        lines = outf.readlines()

    split_indexes = [N for N in range(len(lines)) if '!' == lines[N][0]]

    step_chunks = []
    for n in range(len(split_indexes)):
        step_chunks.append(lines[split_indexes[n]:split_indexes[n+1] if
                           n != len(split_indexes)-1 else len(lines)])

    for current_chunk in step_chunks:

        force_start_line = [line for line in current_chunk if
                            'Forces acting on atoms' in line][0]
        force_end_line = [line for line in current_chunk if
                          'Total force' in line][0]
        force_start_index = current_chunk.index(force_start_line)+2
        force_end_index = current_chunk.index(force_end_line)-2

        atoms_start_line = [line for line in current_chunk if
                            'ATOMIC_POSITIONS' in line][0]
        atoms_end_line = [line for line in current_chunk if
                          'kinetic energy' in line][0]
        atoms_start_index = current_chunk.index(atoms_start_line)+1
        atoms_end_index = current_chunk.index(atoms_end_line)-3

        stress_start_line = [line for line in current_chunk if
                             'Computing stress' in line][0]
        stress_start_index = current_chunk.index(stress_start_line) + 3

        temperature_line = [line for line in current_chunk if
                            'temperature' in line][0]
        dyn_line = [line for line in current_chunk if
                    'Entering Dynamics' in line][0]
        dyn_index = current_chunk.index(dyn_line)
        time_index = dyn_index+1

        forces = []
        for line in current_chunk[force_start_index:force_end_index+1]:
            forceline = line.split('=')[-1].split()
            forces.append(force_conv * np.array([float(forceline[0]),
                                                 float(forceline[1]),
                                                 float(forceline[2])]))
        forces = np.array(forces)
        total_force = float(force_end_line.split('=')[1].strip().split()[0])
        SCF_corr = float(force_end_line.split('=')[2].strip()[0])

        stress = []
        for line in current_chunk[stress_start_index:stress_start_index+3]:
            stress_line = line.split()
            stress.append(np.array([float(stress_line[0]),
                                    float(stress_line[1]),
                                    float(stress_line[2])]))
        stress = np.array(stress)
        stress = stress * stress_conv  # convert to eV/A^3

        positions = []
        elements = []
        for line in current_chunk[atoms_start_index:atoms_end_index+1]:
            atomline = line.split()
            elements.append(atomline[0])
            positions.append(np.array([float(atomline[1]), float(atomline[2]),
                             float(atomline[3])]))
        positions = np.array(positions)

        toten = \
            float(current_chunk[0].split('=')[-1].strip().split()[0]) * en_conv
        temperature_line = temperature_line.split('=')[-1]
        temperature = float(temperature_line.split()[0])
        iteration = int(dyn_line.split('=')[-1])
        timeline = current_chunk[time_index].split('=')[-1].strip().split()[0]
        time = float(timeline)
        Ekin = float(atoms_end_line.split('=')[1].strip().split()[0])

        steps[iteration] = {'iteration': iteration,
                            'forces': forces,
                            'stress': stress,
                            'positions': positions,
                            'elements': elements,
                            'temperature': temperature,
                            'time': time,
                            'energy': toten,
                            'ekin': Ekin,
                            'kinetic energy': Ekin,
                            'total energy': toten,
                            'total force': total_force,
                            'SCF correction': SCF_corr}

    return(steps)
