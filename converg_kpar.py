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
    my_print('CHOICE OF KPAR STARTED',0)
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
        incar["ISTART"]=0
        incar.write_file(dir+"/INCAR")

        my_print(f"Running with KPAR = {kpar}",1)
        run_vasp(dir)
        runtime = get_runtime(dir)

        my_print(f"Runtime = {runtime} sec",2)

        if kpar==1: nkpts=get_nkpts(dir)

        os.system('mv '+dir+'/OUTCAR '+dir+'/OUTCAR.'+str(kpar))

        if runtime < best_time:
            best_time = runtime
            best_kpar = kpar

    my_print(f"Optimal KPAR = {best_kpar} with runtime = {best_time} sec",1)
    incar["KPAR"] = best_kpar
    incar["ISTART"]=1
    incar.write_file("INCAR")
    clean_vasp_files(dir)
    my_print('CHOICE OF KPAR FINISHED',0)


 
