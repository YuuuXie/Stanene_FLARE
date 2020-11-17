import numpy as np
import os
from flare.gp import GaussianProcess
from mgp_code_modified import MappedGaussianProcess
from flare.predict import predict_on_structure_efs_par, predict_on_structure_mgp
from ase.io import read
from flare.ase.atoms import FLARE_Atoms
from flare.env import AtomicEnvironment

os.environ["OMP_NUM_THREADS"] = "64"
gp = GaussianProcess.from_file("stn_small.gp.json")
print("training size", len(gp))

struc_xyz = read("../Data/test_frames.xyz", index="0")
struc = FLARE_Atoms.from_ase_atoms(struc_xyz)
_, gp_f, _, gp_v, _, _ = predict_on_structure_efs_par(struc, gp, 32, write_to_structure=False)
np.save("gp_f", gp_f)
np.save("gp_v", gp_v)
print("gp_f mean, max, min", np.mean(gp_f), np.max(gp_f), np.min(gp_f))
print("gp_v mean, max, min", np.mean(gp_v), np.max(gp_v), np.min(gp_v))

os.environ["OMP_NUM_THREADS"] = "2"
rank_list = [5 * (i+1) for i in range(60)]
#print("rank_list", rank_list, flush=True)
for n3_pow in range(1, 8):
    n3 = 2 ** n3_pow
    grid_params = {
        "twobody": {"grid_num": [144]},
        "threebody": {"grid_num": [n3, n3, n3]},
    }

    data = gp.training_statistics

    mgp_model = MappedGaussianProcess(
        grid_params=grid_params, unique_species=data["species"], n_cpus=32, #var_map="pca",
    )

    mgp_model.build_map(gp)

    mgp_f, mgp_v = predict_on_structure_mgp(struc, mgp_model, write_to_structure=False)
    mae = np.mean(np.abs(mgp_f - gp_f))
    print("force MAE", n3_pow, mae, np.mean(mgp_f), flush=True)
    continue

    for rank in rank_list:
        if rank > n3 ** 3:
            break

        mgp_v = []
        for a in range(len(struc)):
            atom_env = AtomicEnvironment(
                struc, a, mgp_model.cutoffs, cutoffs_mask=mgp_model.hyps_mask
            )
            _, v, _, _ = mgp_model.predict(atom_env, rank=rank)
            mgp_v.append(np.sqrt(np.abs(v)))

        mgp_v = np.array(mgp_v)[:,0]
        #assert np.all(mgp_v.shape == gp_v.shape), f"{mgp_v.shape}, {gp_v.shape}"
        mae = np.mean(np.abs(mgp_v - gp_v))
        print("MAE", n3_pow, rank, mae, np.mean(mgp_v), flush=True)
