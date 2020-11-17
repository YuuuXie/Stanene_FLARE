import os, time

from ase import io
from ase.calculators.espresso import Espresso


def get_dft_calc(n_cpus, npool, kpts):
    # set up qe file
    calculation = "scf"
    label = calculation
    qe_dir = "/n/home08/xiey/q-e"
    pseudo_dir = f"{qe_dir}/pseudo/"

    os.environ[
        "ASE_ESPRESSO_COMMAND"
    ] = f"srun -n {n_cpus} --mpi=pmi2 {qe_dir}/bin/pw.x -npool {npool} -in {label}.pwi > {label}.pwo"

    # set up .pwi input data
    input_data = {
        "control": {
            "prefix": label,
            "pseudo_dir": pseudo_dir,
            "outdir": "./out",
            "verbosity": "high",
            "calculation": calculation,
        },
        "system": {
            "ibrav": 0,
            "ecutwfc": 50,
            "ecutrho": 400,
            "smearing": "mv",
            "degauss": 0.01,
            "occupations": "smearing",
        },
        "electrons": {
            "conv_thr": 1.0e-08,
            #'startingwfc': 'file',
            "electron_maxstep": 200,
            "mixing_beta": 0.5,
        },
        "cell": {"cell_dofree": "2Dxy"},
    }

    # since there could be disordered phase, increase the k pts in z-direction
    pseudopotentials = {"Sn": "sn_pbesol_v1.4.uspp.F.UPF"}

    dft_calc = Espresso(
        pseudopotentials=pseudopotentials,
        label=label,
        directory="dft",
        tstress=True,
        tprnfor=True,
        nosym=True,
        input_data=input_data,
        kpts=kpts,
    )

    return dft_calc
