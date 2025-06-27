import os
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.vasp.inputs import Incar, Poscar, Kpoints, Potcar
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.core import Structure

# Define paths
base_dir = os.getcwd()
structure_file = os.path.join(base_dir, "POSCAR")
incar_file = os.path.join(base_dir, "INCAR")
potcar_file = os.path.join(base_dir, "POTCAR")

# Load initial structure and INCAR
structure = Structure.from_file(structure_file)
incar = Incar.from_file(incar_file)
potcar = Potcar.from_file(potcar_file)

# Define parameter ranges
kpoint_grids = [(3, 3, 3), (4, 4, 4), (5, 5, 5)]
kpar_values = [1, 2, 4]
encut_values = [300, 400, 500]

# Function to run VASP
def write_slurm_script(job_name="vasp_job", nodes=1, ntasks=48, time="400:00:00"):
    with open("run.slurm", "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={ntasks}
#SBATCH --time={time}
#SBATCH --output=vasp.out
#SBATCH --error=vasp.err

module load vasp
srun vasp_std
""")

def run_vasp():
    write_slurm_script()
    subprocess.run(["sbatch", "run.slurm"])

# Function to extract total energy from OUTCAR
def get_total_energy(outcar_path):
    outcar = Outcar(outcar_path)
    return outcar.final_energy

# Optimization loops
def optimize_parameter(param_name, values, fixed_params):
    energies = []
    for val in values:
        dir_name = f"{param_name}_{val}"
        os.makedirs(dir_name, exist_ok=True)
        # Write POSCAR
        poscar = Poscar(structure)
        poscar.write_file(os.path.join(dir_name, "POSCAR"))
        # Write INCAR
        incar_copy = incar.copy()
        incar_copy.update(fixed_params)
        if param_name == "ENCUT":
            incar_copy["ENCUT"] = val
        elif param_name == "KPAR":
            incar_copy["KPAR"] = val
        incar_copy.write_file(os.path.join(dir_name, "INCAR"))
        # Write POTCAR
        shutil.copy(potcar_file, os.path.join(dir_name, "POTCAR"))
        # Write KPOINTS
        if param_name == "KPOINTS":
            kpts = Kpoints.gamma_automatic(kpts=val)
        else:
            kpts = Kpoints.gamma_automatic(kpts=fixed_params.get("KPOINTS", (3, 3, 3)))
        kpts.write_file(os.path.join(dir_name, "KPOINTS"))
        # Run VASP
        run_vasp(dir_name)
        # Extract energy
        energy = get_total_energy(os.path.join(dir_name, "OUTCAR"))
        energies.append((val, energy))
    return energies






# Optimize KPOINTS
print("Optimizing KPOINTS...")
kpoints_energies = optimize_parameter("KPOINTS", kpoint_grids, {"ENCUT": 400, "KPAR": 1})
# Select optimal KPOINTS
optimal_kpoints = min(kpoints_energies, key=lambda x: x[1])[0]

# Optimize KPAR
print("Optimizing KPAR...")
kpar_energies = optimize_parameter("KPAR", kpar_values, {"ENCUT": 400, "KPOINTS": optimal_kpoints})
# Select optimal KPAR
optimal_kpar = min(kpar_energies, key=lambda x: x[1])[0]

# Optimize ENCUT
print("Optimizing ENCUT...")
encut_energies = optimize_parameter("ENCUT", encut_values, {"KPAR": optimal_kpar, "KPOINTS": optimal_kpoints})
# Select optimal ENCUT
optimal_encut = min(encut_energies, key=lambda x: x[1])[0]

# Final relaxation
print("Running final relaxation...")
relax_dir = "final_relaxation"
os.makedirs(relax_dir, exist_ok=True)
# Write POSCAR
poscar = Poscar(structure)
poscar.write_file(os.path.join(relax_dir, "POSCAR"))
# Write INCAR
final_incar = incar.copy()
final_incar.update({
    "ENCUT": optimal_encut,
    "KPAR": optimal_kpar,
    "IBRION": 2,
    "ISIF": 3,
    "NSW": 100,
    "EDIFF": 1e-5,
    "EDIFFG": -0.02
})
final_incar.write_file(os.path.join(relax_dir, "INCAR"))
# Write POTCAR
shutil.copy(potcar_file, os.path.join(relax_dir, "POTCAR"))
# Write KPOINTS
kpts = Kpoints.gamma_automatic(kpts=optimal_kpoints)
kpts.write_file(os.path.join(relax_dir, "KPOINTS"))
# Run VASP
run_vasp(relax_dir)

print("Optimization and relaxation complete.")
