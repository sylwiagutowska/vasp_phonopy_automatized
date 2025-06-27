import os
from common import *
  
 


def encut_convergence(structure,incar, kpoints, kpar=1, energy_tol=1e-3):
    encut=int(incar["ENCUT"])
    encut_list = np.arange(encut, 2*encut, 20)
    previous_energy = None

    dir='encut_conv'
    os.system('mkdir '+dir+'; cp POSCAR POTCAR '+dir)
    kpoints.write_file(dir+"/KPOINTS")

    incar["KPAR"] = kpar

    for ne,encut in enumerate(encut_list):
        incar["ENCUT"] = encut
        incar.write_file(dir+"/INCAR")

        print(f"Running with ENCUT = {encut}")
        run_vasp(dir)
        energy = get_energy(dir)
        os.system('mv '+dir+'/OUTCAR '+dir+'/OUTCAR.'+str(encut))

        print(f"Energy = {energy} eV")

        if previous_energy is not None:
            diff = abs(energy - previous_energy)/len(structure)
            print(f"Energy difference per atom = {diff} eV")
            if diff < energy_tol:
                print(f"ENCUT converged at {encut}")
                break

        previous_energy = energy
    incar["ENCUT"] = encut_list[ne-1]
    incar.write_file("INCAR")
    clean_vasp_files(dir)