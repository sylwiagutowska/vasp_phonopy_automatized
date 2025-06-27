import sys
from common import *
from prepare_files import *
from converg_k import *
from converg_kpar import *
from converg_encut import *
from bands_dos import *
from relax import *
from run_phonopy import *
[name,vaspcomm,npr]=sys.argv[-3:]
npr=int(npr)

'''
print(name)
structure,energy_above_hull=download_poscar(name)
if structure==-1: 
  raise Warning(name+' not found in materialsproject') 
  structure=Structure.from_file('POSCAR')
else:
 print(name,'poscar downloaded successfully')

potcar,incar=prepare_potcar_incar(name,structure,energy_above_hull)  
set_vasp(vaspcomm,incar)


structure,potcar,incar= Structure.from_file('POSCAR'), Potcar.from_file('POTCAR'), Incar.from_file('INCAR')
kpoints_convergence(structure, incar, kpar=1, energy_tol=1e-4)

structure,incar,kpoints= Structure.from_file('POSCAR'), Incar.from_file('INCAR'), Kpoints.from_file('KPOINTS')
encut_convergence(structure,incar, kpoints, kpar=1, energy_tol=1e-3)

structure,potcar,incar,kpoints= Structure.from_file('POSCAR'), Potcar.from_file('POTCAR'), Incar.from_file('INCAR'), Kpoints.from_file('KPOINTS')
kpar_convergence(incar, kpoints,  npr)

'''

structure,potcar,incar,kpoints= Structure.from_file('POSCAR'), Potcar.from_file('POTCAR'), Incar.from_file('INCAR'), Kpoints.from_file('KPOINTS')
set_vasp(vaspcomm,incar)
#scf_bands_dos(structure,incar,kpoints)
#plot_bands_dos('bands','dos')

run_phonopy(structure, incar, kpoints, kpar=1,max_q=3)


#structure,potcar,incar,kpoints= Structure.from_file('POSCAR'), Potcar.from_file('POTCAR'), Incar.from_file('INCAR'), Kpoints.from_file('KPOINTS')
#relaxx(incar,kpoints)
