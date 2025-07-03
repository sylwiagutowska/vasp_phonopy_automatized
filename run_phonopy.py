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
   print("######PHONON CALCULATIONS#####")
   for q in range(1,max_q+1):
      if q==1 and len(structure)==1: continue #phonons for 1 atom dont exist
      #directory
      dir='q'+str(q)
      print(dir)
      if len(glob.glob(dir+'/total_dos.dat'))!=0: 
         print('total dos already calculated; i continue with next q if required')
         olddos=np.loadtxt(dir+'/total_dos.dat',skiprows=1)[:,1]
         continue
      #os.system('mv '+dir+' '+dir+'_old')
      os.system('mkdir '+dir+'; cp POSCAR POTCAR '+dir)
      #kpoints
      print(' Preparing kpoints')
      kpold=np.array(kpoints.kpts[0])
      kp=np.array([(int(m)) for m in (np.ceil(np.array(kpold)/q))])
      kpointsq=Kpoints.monkhorst_automatic(kp)
      kpointsq.write_file(dir+"/KPOINTS")
      #phonopy
      qs=(np.round(q*kpold/kpold[0]).astype(int))
      print(qs,kpold,kp)
      unitcell, _ = read_crystal_structure("POSCAR")
      print(' Ready for phononpy with supercell of dim',qs)
      phonon = Phonopy(unitcell,  supercell_matrix=[[qs[0], 0, 0], [0, qs[1], 0], [0, 0, qs[2]]], factor=VaspToTHz)
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
      print(" Files ready. I calculate perfect cell...")
      #perfect run
      dir_perfect=dir+'/PERFECT/'
      if len(glob.glob(dir_perfect))==0:
       os.system('mkdir '+dir_perfect)
       os.system('cp '+dir+"/{POTCAR,INCAR,KPOINTS} "+dir_perfect)
       write_crystal_structure(dir_perfect+'/POSCAR',phonon.supercell)
       run_vasp(dir_perfect)
      print(' Perfect calculated')

      dirs=[dir+'/disp-'+str(i+1) for i in range(len(supercells))]
      for ni,i in enumerate(dirs):
         print(' ',i,'is being calculated...')
         if len(glob.glob(i))==0:
            os.system('cp -r '+dir_perfect+' '+i)
            write_crystal_structure(i+'/POSCAR',supercells[ni])
            run_vasp(i)
            os.system('rm '+i+'/{CHG,CONTCAR,EIGENVAL,OSZICAR,PCDAT,POTCAR,WAVECAR,vaspout.h5,CHGCAR,DOSCAR,IBZKPT,REPORT,XDATCAR}')
      print(' All displacements calculated')
      nms=([str(i)+'/vasprun.xml' for i in  dirs])
      print(' Forces will be read from following files', nms)
      forces =  parse_set_of_forces( filenames=nms)
      phonon.set_forces(forces)
      phonon.produce_force_constants()
      phonon.symmetrize_force_constants(level=1) #level :    Application of translational and permulation symmetries is repeated by this number. Default is 1.
      print(" I calculate DOS now...")
   #   phonon.run_mesh([40, 40, 40],with_eigenvectors=True,is_mesh_symmetry=False) 
   #   phonon.run_total_dos(freq_min=-2.,freq_max=20.-0.005,freq_pitch=0.05)
      phonon.auto_total_dos(mesh=40,is_mesh_symmetry=False,plot=False,write_dat=True,filename=dir+"/total_dos.dat")
   #   phonon.write_total_DOS(filename=dir+"/total_dos.dat")
      print(' DOS  calculated')
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
      #print('force constants:\n',force_constants)
      atoms=phonon._primitive._symbols
      primitive_matrix=structure.lattice.matrix 
      sc_len=len(force_constants)//len(atoms) 
      force_constants2=np.abs(np.array([[np.matmul(np.linalg.inv(primitive_matrix),np.matmul(m,np.linalg.inv(primitive_matrix).transpose())) for m in i] for i in force_constants]))
      print( ' Map atoms from primitive to supercell',phonon._primitive._p2s_map) 
      primitive_to_supercell_map=[[m for m in range(i,i+sc_len)] for i in phonon._primitive._p2s_map]
      if_forces_reduced_enough=np.zeros((len(atoms),len(atoms),3))
      print(' Force constants analysis')
      for nato1,ato1 in enumerate(primitive_to_supercell_map):
        for nato2,ato2 in enumerate(primitive_to_supercell_map):
            print(atoms[nato1],atoms[nato2],':')
            for ndire,dire in enumerate(['x','y','z']):
               maxf,minf=np.max(force_constants2[ato1[0],ato2,ndire,ndire]),np.min(force_constants2[ato1[0],ato2,ndire,ndire])
               print('  ',dire,':','max force constant: ',maxf,'min force constant:',minf)
               if maxf>minf*900: if_forces_reduced_enough[nato1,nato2,dire]=1
            for ncell1,cell1 in enumerate(ato1):
              if ncell1!=0: break
              for ncell2,cell2 in enumerate(ato2):
               print(f'  {ncell1}:{ncell2}: \n',force_constants2[cell1][cell2])
      if_converg=np.all(if_forces_reduced_enough)
      if if_converg:
        print('Convergence with respect to supercell size achieved. Now I run band and DOS calculations')
      if not if_converg and q==max_q:
        print('Convergence with respect to supercell size NOT achieved. Now I run band and DOS calculations')
      if np.all(if_forces_reduced_enough) or q==max_q:
        phonon.auto_band_structure(with_eigenvectors=True,write_yaml=True,plot=False)
        print('Band structure calculated')
        phonon.run_mesh([int(100)/sc_size for sc_size in qs],with_eigenvectors=True,is_mesh_symmetry=False) 
        phonon.run_total_dos(freq_pitch=0.01)
        phonon.write_total_DOS(filename=dir+"/total_dos.dat")
        print('DOS  calculated')
        phonon.run_projected_dos(freq_pitch=0.01)
        plot=phonon.plot_band_structure_and_dos()
        tick_positions, tick_labels = plot.get_ticks()
        tick_positions = [round(pos, 2) for pos in tick_positions]
        plot.plot(ticks=[tick_positions, tick_labels])
        plot.savefig("phonon_bands_and_dos.png", dpi=600, bbox_inches='tight')
        phonon.save(filename=dir+"/phonopy_params.yaml",settings={'force_constants': True})
        print('Data saved. PHONONS FINISHED')
        #other commands to deal with dos and bands
        # frequencies, dos = phonon.get_total_DOS()
        #  phonon.plot_total_DOS().show() 
        #  phonon.plot_projected_dos().show()
        break
      