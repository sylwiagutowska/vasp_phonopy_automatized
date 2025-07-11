import os
import glob
import numpy as np
import sys
from common import *
from phonopy import Phonopy
from phonopy.interface.calculator import read_crystal_structure, write_crystal_structure
from phonopy.units import VaspToTHz

import xml.etree.ElementTree as ET
import numpy as np

def parse_forces_from_vasprun(xml_file):
    def strip_ns(tag):
        return tag.split('}', 1)[-1] if '}' in tag else tag

    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Find all varray elements with name="forces"
    force_arrays = []
    for varray in root.iter():
        if strip_ns(varray.tag) == 'varray' and varray.attrib.get('name') == 'forces':
            force_arrays.append(varray)

    if not force_arrays:
        raise ValueError("No forces found in vasprun.xml.")

    final_forces_raw = force_arrays[-1]
    forces = []
    for v in final_forces_raw:
        if strip_ns(v.tag) == 'v': 
            forces.append([float(x) for x in v.text.strip().split()])
    return forces

def parse_set_of_forces(filenames):
   forces=[]
   for filename in filenames:
     forces.append(parse_forces_from_vasprun(filename))
   return np.array(forces)*2 #.reshape(-1, 3) #i have to use factor 2 to have the same results as phonopy -fc -f ..., phonopy mesh.conf

#@profile
def run_phonopy(structure, incar, kpoints, kpar=1,max_q=8):
   my_print("PHONON CALCULATIONS",0)
   for q in range(1,max_q+1):
      if q==1 and len(structure)==1: continue #phonons for 1 atom dont exist
      #directory
      dir='q'+str(q)
      my_print('supercell No.:'+dir,1)
      if len(glob.glob(dir+'/total_dos.dat'))!=0: 
         my_print('total dos already calculated; i continue with next q if required',2)
        # olddos=np.loadtxt(dir+'/total_dos.dat',skiprows=1)[:,1]
         continue
      #os.system('mv '+dir+' '+dir+'_old')
      os.system('mkdir '+dir+'; cp POSCAR POTCAR '+dir)
      #kpoints
      my_print('Preparing kpoints',2)
      kpold=np.array(kpoints.kpts[0])
      kp=np.array([(int(m)) for m in (np.ceil(np.array(kpold)/q))])
      kpointsq=Kpoints.monkhorst_automatic(kp)
      kpointsq.write_file(dir+"/KPOINTS")
      #phonopy
      dim=(np.round(q*kpold/kpold[0]).astype(int))
      my_print('supercell size',dim,'kpoints size',kp,2)
      unitcell, _ = read_crystal_structure("POSCAR")
      my_print('Ready for phononpy with supercell of dim',dim,2)
      phonon = Phonopy(unitcell,  supercell_matrix=[[dim[0], 0, 0], [0, dim[1], 0], [0, 0, dim[2]]], factor=VaspToTHz)
      phonon.generate_displacements(distance=0.01)
      supercells = phonon.supercells_with_displacements 
      incar_phonopy=incar.copy()
      incar_phonopy["LWAVE"]=True
      incar_phonopy["LCHARG"]=True
      incar_phonopy["LWAP"]=True
      incar_phonopy["IBRION"]=-1
      incar_phonopy["EDIFF"]=1e-8
      incar_phonopy["NSW"]=1
      incar_phonopy["KPAR"]=int(np.ceil(incar_phonopy["KPAR"]/q))
      if "MAGMOM" in incar and len(incar["MAGMOM"])!=0: incar_phonopy["MAGMOM"]= len(Structure.from_file('SPOSCAR').sites)*"0 "
      incar_phonopy.write_file(dir+"/INCAR")
      my_print("Files ready. I calculate perfect cell...",2)
      #perfect run
      dir_perfect=dir+'/PERFECT/'
      if len(glob.glob(dir_perfect))==0:
       os.system('mkdir '+dir_perfect)
       os.system('cp '+dir+"/{POTCAR,INCAR,KPOINTS} "+dir_perfect)
       write_crystal_structure(dir_perfect+'/POSCAR',phonon.supercell)
       run_vasp(dir_perfect)
      my_print('Perfect calculated',2)

      dirs=[dir+'/disp-'+str(i+1) for i in range(len(supercells))]
      for ni,i in enumerate(dirs):
         my_print(i,'is being calculated...',2)
         if len(glob.glob(i))==0:
            os.system('cp -r '+dir_perfect+' '+i)
            write_crystal_structure(i+'/POSCAR',supercells[ni])
            run_vasp(i)
            os.system('rm '+i+'/{CHG,CONTCAR,EIGENVAL,OSZICAR,PCDAT,POTCAR,WAVECAR,vaspout.h5,CHGCAR,DOSCAR,IBZKPT,REPORT,XDATCAR}')
      my_print('All displacements calculated',2)
      nms=([str(i)+'/vasprun.xml' for i in  dirs])
      my_print('Forces will be read from following files', nms,2)
      forces =  parse_set_of_forces( filenames=nms)
      phonon.set_forces(forces)
      phonon.produce_force_constants()
      phonon.symmetrize_force_constants(level=1) #level :    Application of translational and permulation symmetries is repeated by this number. Default is 1.
      my_print("Forces finished. I calculate DOS now...",2)
   #   phonon.run_mesh([40, 40, 40],with_eigenvectors=True,is_mesh_symmetry=False) 
   #   phonon.run_total_dos(freq_min=-2.,freq_max=20.-0.005,freq_pitch=0.05)
      phonon.auto_total_dos(mesh=40,is_mesh_symmetry=False,plot=False,write_dat=True,filename=dir+"/total_dos.dat")
   #   phonon.write_total_DOS(filename=dir+"/total_dos.dat")
      my_print('DOS  calculated',2)
      phonon.save(filename=dir+"/phonopy_params.yaml",settings={'force_constants': True})
      
      #old method of termination: when DOS doesnt change much
      '''
      try: dos=np.loadtxt(dir+'/total_dos.dat',skiprows=1)[:,1]
      except: 
       if q==1: dos=np.zeros((1001))
      if q==1: olddos=dos
      if q>1:
         dif=(np.sum(np.abs(dos-olddos))/np.sum(np.abs(dos)))
         print('dif',dif)
         olddos=dos
         if dif<0.3: break
      ''' 

      #new method of termination: when forces constant drop by 3 orders of magnitude
      force_constants=phonon.get_force_constants()
      my_print('Dimension of force_constant matrix',force_constants.shape,2)
      #print('force constants:\n',force_constants)
      atoms=phonon._primitive._symbols
      primitive_matrix=structure.lattice.matrix 
      sc_len=len(force_constants)//len(atoms) 
      force_constants2=np.abs(np.array([[np.matmul(np.linalg.inv(primitive_matrix),np.matmul(m, primitive_matrix )) for m in i] for i in force_constants]))
      primitive_to_supercell_map=np.array([[m for m in range(i,i+sc_len)] for i in phonon._primitive._p2s_map])
      my_print( 'Map atoms from primitive to supercell',phonon._primitive._p2s_map,'reshaped to',primitive_to_supercell_map,2) 
      if_forces_reduced_enough=np.zeros((len(atoms),len(atoms),3))
      my_print('Force constants analysis',2)
      for nato1,ato1 in enumerate(primitive_to_supercell_map): #over atoms in primitive cell. ato1 is a list of atoms equivalent to given atom in primitive cell
        for nato2,ato2 in enumerate(primitive_to_supercell_map): #over elements
            my_print(atoms[nato1],atoms[nato2],':',2)
            for ndire,dire in enumerate(['x','y','z']):
               if ndire==0: select_atoms=[ato2[0]+m for m in range(dim[0])]
               elif ndire==1: select_atoms=[ato2[0]+m*dim[0] for m in range(dim[1])]
               elif ndire==2: select_atoms=[ato2[0]+m*dim[1]*dim[2] for m in range(dim[2])]
               my_print('atoms in the direction',dire,':',select_atoms,2)
               maxf,minf=np.max(force_constants2[ato1[0],ato2[select_atoms],ndire,ndire]),np.min(force_constants2[ato1[0],ato2[select_atoms],ndire,ndire])
               my_print(dire,':','max force constant: ',maxf,'min force constant:',minf,3)
               if maxf>minf*900: if_forces_reduced_enough[nato1,nato2,ndire]=1
            for ncell1,cell1 in enumerate(ato1):
              if ncell1!=0: break
              for ncell2,cell2 in enumerate(ato2):
               my_print(f'{ncell1}:{ncell2}: \n',force_constants2[cell1][cell2],3)
      if_converg=np.all(if_forces_reduced_enough)
      if if_converg:
        my_print('Convergence with respect to supercell size achieved. Now I run band and DOS calculations',1)
      if (not if_converg) and q==max_q:
        my_print('Convergence with respect to supercell size NOT achieved. Now I run band and DOS calculations',1)
      if np.all(if_forces_reduced_enough) or q==max_q:
        phonon.auto_band_structure(with_eigenvectors=True,write_yaml=True,plot=False)
        my_print('Band structure calculated',1)
        phonon.run_mesh([int(100)/sc_size for sc_size in dim],with_eigenvectors=True,is_mesh_symmetry=False) 
        phonon.run_total_dos(freq_pitch=0.01)
        phonon.write_total_DOS(filename=dir+"/total_dos.dat")
        my_print('DOS  calculated',1)
        phonon.run_projected_dos(freq_pitch=0.01)
        plot=phonon.plot_band_structure_and_dos() 
        plot.savefig("phonon_bands_and_dos.png", dpi=600, bbox_inches='tight')
        phonon.save(filename=dir+"/phonopy_params.yaml",settings={'force_constants': True})
        my_print('Data saved. PHONONS FINISHED',0)
        #other commands to deal with dos and bands
        # frequencies, dos = phonon.get_total_DOS()
        #  phonon.plot_total_DOS().show() 
        #  phonon.plot_projected_dos().show()
        break
      
