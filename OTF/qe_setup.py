import os, time

from ase import io
from ase.calculators.espresso import Espresso


# set up qe file
calculation = "scf"
label = calculation
input_file_name = label + ".pwi"
output_file_name = label + ".pwo"
qe_dir = "/n/home08/xiey/q-e"
pseudo_dir = f"{qe_dir}/pseudo/"

n_cpus = 288
npool = 36
os.environ[
    "ASE_ESPRESSO_COMMAND"
] = f"srun -n {n_cpus} --mpi=pmi2 {qe_dir}/bin/pw.x -npool {npool} -in scf.pwi > scf.pwo"

# set up .pwi input data
input_data = {}
input_data["control"] = {
    "prefix": label,
    "pseudo_dir": pseudo_dir,
    "outdir": "./out",
    "verbosity": "high",
    "calculation": calculation,
}

input_data["system"] = {
    "ibrav": 0,
    "ecutwfc": 50,
    "ecutrho": 400,
    "smearing": "mv",
    "degauss": 0.01,
    "occupations": "smearing",
}

input_data["electrons"] = {
    "conv_thr": 1.0e-08,
    #'startingwfc': 'file',
    "electron_maxstep": 200,
    "mixing_beta": 0.5,
}

input_data["cell"] = {"cell_dofree": "2Dxy"}

# since there could be disordered phase, increase the k pts in z-direction
kpts = (4, 4, 4)
pseudopotentials = {'Sn': 'sn_pbesol_v1.4.uspp.F.UPF'}

dft_calc = Espresso(
    pseudopotentials=pseudopotentials,
    label=label,
    tstress=True,
    tprnfor=True,
    nosym=True,
    input_data=input_data,
    kpts=kpts,
)
