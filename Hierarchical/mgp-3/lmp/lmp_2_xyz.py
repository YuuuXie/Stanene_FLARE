import numpy as np
import sys

def parse_lmp(filename, skip):
    f = open(filename)
    
    lines = f.readlines()
    nat = int(lines[3].strip())
    cell = np.zeros((3,3))
    for i in range(3):
        line = lines[i+5].split()
        cell[i,i] = float(line[1])

    traj_len = len(lines) // (nat+9)
    xyz_list = []
    for t in range(0, traj_len, skip):
        head = t * (nat + 9) + 9
        tail = head + nat
        block = lines[head:tail]
        xyz = np.zeros((nat, 3))

        for l, line in enumerate(block):
            line = line.split()
            atom = int(line[0])-1
            xyz[atom] = np.array([float(l) for l in line[2:5]])

        xyz_list.append(xyz)

    return xyz_list, cell


def list_2_xyz(xyz_list, xyz_file):
    f = open(xyz_file, 'w')

    for frame in xyz_list:
        nat = len(frame)
        f.write('{} \n\n'.format(nat))
        for atom in frame:
            f.write('Sn {} {} {}\n'.format(atom[0], atom[1], atom[2]))

    f.close()

if __name__ == "__main__":
    snap = int(sys.argv[1])
    xyz_list, cell = parse_lmp('{}ps.lmp'.format(snap), snap)
    list_2_xyz(xyz_list, '{}ps.xyz'.format(snap))
