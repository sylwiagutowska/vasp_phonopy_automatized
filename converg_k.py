import os
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar
from pymatgen.io.vasp.outputs import Outcar
from common import *

 


oldene,ene=2000,0
oldk=1


def kpoints_convergence(structure, incar, kpar=1, energy_tol=1e-3):
    k_density_list = [100, 200, 400, 600, 800, 1000, 1200]  #angstrom-1 or per atom
    previous_energy = None

    dir='kpoints_conv'
    os.system('mkdir '+dir+'; cp POSCAR POTCAR '+dir)

    incar["ENCUT"] = incar["ENCUT"]*1.3
    incar["KPAR"] = kpar
    incar.write_file(dir+"/INCAR")

    for nk,k_density in enumerate(k_density_list):
        print(f"Testing k-points density: {k_density}")
        kpoints = Kpoints.automatic_density(structure, k_density)
        kpoints.write_file(dir+"/KPOINTS")
        
        run_vasp(dir)
        energy = get_energy(dir)
        os.system('mv '+dir+'/OUTCAR '+dir+'/OUTCAR.'+str(k_density))
        print(kpoints) 
        print(f" Energy = {energy} eV")

        if previous_energy is not None:
            diff = abs(energy - previous_energy)/len(structure)
            print(f" energy difference per atom = {diff} eV")
            if diff < energy_tol:
                print(f"KPOINTS converged at density {k_density}")
                break

        previous_energy = energy
    kpoints=Kpoints.automatic_density(structure, k_density_list[nk-1])
    kpoints.write_file("KPOINTS")
    incar["ENCUT"] = incar["ENCUT"]/1.3
    clean_vasp_files(dir)