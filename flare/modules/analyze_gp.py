import numpy as np
import sys
from flare.modules import crystals, qe_parsers
import copy
from flare import gp, struc, env


class MDAnalysis:
    def __init__(self, MD_output_file: str, cell: np.ndarray):
        self.MD_output_file = MD_output_file
        self.cell = cell
        self.MD_data = self.get_data_from_file()

    def get_data_from_file(self):
        data = qe_parsers.parse_md_output(self.MD_output_file)
        return data

    def get_structure_from_snap(self, snap):
        positions = self.MD_data[snap]['positions']
        species = self.MD_data[snap]['elements']
        structure = struc.Structure(self.cell, species, positions)
        return structure

    def get_forces_from_snap(self, snap):
        forces = self.MD_data[snap+1]['forces']
        return forces


# given a gp model, predict vacancy diffusion activation profile
def vac_diff_fcc(gp_model, gp_cell, fcc_cell, species, cutoff, nop):

    # create 2x2x2 fcc supercell with atom 0 removed
    alat = fcc_cell[0, 0]
    fcc_unit = crystals.fcc_positions(alat)
    fcc_super = crystals.get_supercell_positions(2, fcc_cell, fcc_unit)
    vac_super = copy.deepcopy(fcc_super)
    vac_super.pop(0)
    vac_super = np.array(vac_super)

    # create list of positions for the migrating atom
    start_pos = vac_super[0]
    end_pos = np.array([0, 0, 0])
    diff_vec = end_pos - start_pos
    test_list = []
    step = diff_vec / (nop - 1)
    for n in range(nop):
        test_list.append(start_pos + n*step)

    # predict force on migrating atom
    vac_copy = copy.deepcopy(vac_super)
    store_res = np.zeros((6, nop))

    cutoffs = np.array([cutoff, cutoff])

    for count, pos_test in enumerate(test_list):
        vac_copy[0] = pos_test

        struc_curr = struc.Structure(gp_cell, species, vac_copy)
        env_curr = env.AtomicEnvironment(struc_curr, 0, cutoffs)

        fc_comp_x, vr_x = gp_model.predict(env_curr, 1)
        fc_comp_y, vr_y = gp_model.predict(env_curr, 2)
        fc_comp_z, vr_z = gp_model.predict(env_curr, 3)

        store_res[0, count] = fc_comp_x
        store_res[1, count] = fc_comp_y
        store_res[2, count] = fc_comp_z
        store_res[3, count] = vr_x
        store_res[4, count] = vr_y
        store_res[5, count] = vr_z

    return store_res


# return untrained gp model
def get_gp_from_snaps(md_trajectory, training_snaps, kernel,
                      kernel_grad, hyps, cutoffs, algorithm='BFGS',
                      maxiter=1000, par=False):

    gp_model = gp.GaussianProcess(kernel, kernel_grad, hyps,
                                  cutoffs, opt_algorithm=algorithm,
                                  maxiter=maxiter, par=par)

    for snap in training_snaps:
        structure = md_trajectory.get_structure_from_snap(snap)
        forces = md_trajectory.get_forces_from_snap(snap)
        gp_model.update_db(structure, forces)

    return gp_model


# return untrained gp model
def custom_range_gp(md_trajectory, training_snaps, training_range, kernel,
                    kernel_grad, hyps, cutoffs, algorithm='BFGS',
                    maxiter=1000):

    gp_model = gp.GaussianProcess(kernel, kernel_grad, hyps,
                                  cutoffs, opt_algorithm=algorithm,
                                  maxiter=maxiter)

    for snap_count, snap in enumerate(training_snaps):
        structure = md_trajectory.get_structure_from_snap(snap)
        forces = md_trajectory.get_forces_from_snap(snap)
        gp_model.update_db(structure, forces, training_range[snap_count])

    return gp_model


def predict_forces_on_structure(gp_model, md_trajectory, snap, cutoffs):
    structure = md_trajectory.get_structure_from_snap(snap)
    forces = np.array(md_trajectory.get_forces_from_snap(snap))
    noa = len(structure.positions)

    predictions = np.zeros([noa, 3])
    variances = np.zeros([noa, 3])

    for atom in range(noa):
        env_curr = env.AtomicEnvironment(structure, atom, cutoffs)

        for n in range(3):
            d = n+1
            comp_pred, comp_var = gp_model.predict(env_curr, d)
            predictions[atom, n] = comp_pred
            variances[atom, n] = comp_var

    # mean absolute error
    MAE = np.mean(np.abs(forces-predictions))

    # mean absolute std
    stds = np.sqrt(variances)
    MAS = np.mean(stds)

    return predictions, variances, forces, MAE, MAS


def predict_forces_on_test_set(gp_model, md_trajectory, snaps, cutoff):
    all_predictions = []
    all_variances = []
    all_forces = []
    for snap in snaps:
        predictions, variances, forces, _, _ = \
            predict_forces_on_structure(gp_model, md_trajectory, snap, cutoff)
        predictions = predictions.reshape(predictions.size)
        variances = variances.reshape(variances.size)
        forces = forces.reshape(forces.size)

        for m in range(len(predictions)):
            all_predictions.append(predictions[m])
            all_variances.append(variances[m])
            all_forces.append(forces[m])

    all_predictions = np.array(all_predictions)
    all_variances = np.array(all_variances)
    all_forces = np.array(all_forces)

    return all_predictions, all_variances, all_forces


if __name__ == '__main__':
    pass
