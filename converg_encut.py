import os
from common import *
  
 


def encut_convergence(structure,incar, kpoints, kpar=1, energy_tol=1e-3):
    my_print('ENCUT CONVERGENCE STARTED',0)
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

        my_print(f"Running with ENCUT = {encut}",1)
        run_vasp(dir)
        energy = get_energy(dir)
        os.system('mv '+dir+'/OUTCAR '+dir+'/OUTCAR.'+str(encut))

        my_print(f"Energy = {energy} eV",2)

        if previous_energy is not None:
            diff = abs(energy - previous_energy)/len(structure)
            my_print(f"Energy difference per atom = {diff} eV",2)
            if diff < energy_tol:
                my_print(f"ENCUT converged at {encut}",1)
                my_print('ENCUT CONVERGENCE FINISHED',0)
                break

        previous_energy = energy
    incar["ENCUT"] = encut_list[ne-1]
    incar.write_file("INCAR")
    clean_vasp_files(dir)
