import os
import numpy as np
from common import *
 

def get_runtime(dir):
    outcar = Outcar(dir+"/OUTCAR")
    return outcar.run_stats["Total CPU time used (sec)"]


def get_nkpts(dir):
    with open(dir+"/OUTCAR") as f:
        for line in f:
            if "NKPTS =" in line:
                nkpts = int(line.split()[3])
                print(f'number of irreducible points:{nkpts}')
                return nkpts


def kpar_convergence(incar, kpoints,  npr):
    best_kpar = None
    best_time = float("inf")

    dir='kpar_conv'
    os.system('mkdir '+dir+'; cp POSCAR POTCAR '+dir)
    kpoints.write_file(dir+"/KPOINTS")

    nkpts=1
    for kpar in range(1,npr):
        if kpar>nkpts: break
        if npr%kpar!=0: continue
        incar["KPAR"] = kpar
        incar.write_file(dir+"/INCAR")

        print(f"Running with KPAR = {kpar}")
        run_vasp(dir)
        runtime = get_runtime(dir)

        print(f" Runtime = {runtime} sec")

        if kpar==1: nkpts=get_nkpts(dir)

        os.system('mv '+dir+'/OUTCAR '+dir+'/OUTCAR.'+str(kpar))

        if runtime < best_time:
            best_time = runtime
            best_kpar = kpar

    print(f" Optimal KPAR = {best_kpar} with runtime = {best_time} sec")
    incar["KPAR"] = best_kpar
    incar.write_file("INCAR")
    clean_vasp_files(dir)


 