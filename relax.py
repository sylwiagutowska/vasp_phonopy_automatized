import os 
from common import *
def relax_status(dir):
    with open(dir+"/OUTCAR") as f:
        for line in f:
            if "reached required accuracy - stopping structural energy minimisation =" in line:
                return 1
    return 0



def relaxx(incar,kpoints):
 print('relax starts')
 dir='relax'
 os.system('mkdir '+dir)
 incar_relax=incar.copy()
 incar_relax["ISIF"]=3
 incar_relax["NSW"]=100
 incar_relax["EDIFF"]=1e-7
 incar_relax["EDIFG"]=1e-5
 incar_relax["LCHARG"]=True
 incar_relax["LWAVE"]=True
 incar_relax.write_file(dir+'/INCAR')
 kpoints.write_file(dir+'/KPOINTS')
 os.system('cp POTCAR POSCAR '+dir)
 for m in range(20):
  clean_vasp_files(dir)
  os.system('cp POTCAR  '+dir)
  run_vasp(dir)
  if  relax_status(dir):    break
  os.system("cp {dir}/CONTCAR {dir} POSCAR")
  os.system(f"mv {dir}/OUTCAR {dir}/OUTCAR.{m}")
 os.system(f"mv POSCAR POSCAR_init; cp {dir}/CONTCAR POSCAR")
 print('relax ends')





