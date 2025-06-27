import numpy as np
import os
import glob
from pymatgen.io.vasp.inputs import Incar, Poscar, Kpoints, Potcar
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.core import Structure
import sys

def set_vasp(vaspcomm,incar):
 try:
  if incar["LNONCOLLINEAR"]: vaspcomm=vaspcomm.replace('std','ncl')
  print('i detected lnoncollinear so i changed vasp command to ncl:')
  print(vaspcomm)
 except: 0 
 global vasp_command
 vasp_command=vaspcomm



def run_vasp(dir):
    
    os.system('cd '+dir+'; '+vasp_command+';cd ..')
    h=open(dir+'/OUTCAR')
    tmp=h.readlines()
    h.close()
    for i in tmp:
        if 'The electronic self-consistency was not achieved in the given' in i or 'EEEEEEE  RRRRRR   RRRRRR   OOOOOOO  RRRRRR' in i:
         return 0 
    return 1

def read_energy(dir='./'):
    h=open(dir+'/OUTCAR')
    tmp=h.readlines()
    h.close()
    ene=None
    for i in tmp: 
        if 'free  energy   TOTEN' in i: ene=float(i.split()[4])
    return ene


def get_energy(dir='./'):
    outcar = Outcar(dir+"/OUTCAR")
    return outcar.final_energy


def clean_vasp_files(dir):
    files_to_remove = [
        "WAVECAR",  # large wavefunction file
        "CHGCAR",   # charge density file
        "CHG",      # charge density
        "IBZKPT",
        "EIGENVAL",
        "DOSCAR",
        "vasprun.xml",
        "OUTCAR",
        "OSZICAR",
        "vaspout.h5",
        "XDATCAR",
        "PCDAT",
        "POTCAR",
        "REPORT"
    ]

    for f in files_to_remove:
        if os.path.exists(dir+'/'+f):
            os.remove(dir+'/'+f)
