import numpy as np
from ase.io import read

from flare.gp import GaussianProcess
from flare.ase.atoms import FLARE_Atoms

hyps = np.array(
    [3.27429828e-02, 5.18480876e-01, 3.75272358e-05, 6.90918291e-01, 1.40038968e-01]
)
cutoffs = [7.2, 7.2]
gp = GaussianProcess(kernels=["twobody", "threebody"], hyps=hyps, cutoffs=cutoffs)

training_frames = read("../data/train_frames.xyz", index=":")
add_atoms = open("../data/add_atoms.dat").readlines()[1:]

for ind, struc in enumerate(training_frames):
    add = add_atoms[ind].split()
    add_list = [int(i) for i in add[2:]]

    struc = FLARE_Atoms.from_ase_atoms(struc)
    gp.update_db(struc, struc.get_forces(), add_list)
    if len(gp) >= 100:
        break

gp.set_L_alpha()

print('training_set_size', len(gp.training_data))
gp.write_model('stn_small.gp', format='json')
