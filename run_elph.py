import os
import glob
import numpy as np
import sys
from common import *  



def read_fftmesh(dir='./'):
    h=open(dir+'/OUTCAR')
    tmp=h.readlines()
    h.close()
    ene=None
    for i in tmp: 
        if 'dimension x,y,z NGX ' in i: 
            j=i.split()
            return [ j[m] for m in [4,7,10]]
    return 0

def read_yaml(dir='./'):
    import yaml
    with open(dir+'/phonopy_disp.yaml', 'r') as file:
       data = yaml.safe_load(file)
    return data['phonopy']['configuration']['dim']


def read_a2f(dir='./'):
    import h5py
    h5 = h5py.File(dir+'/vaspout.h5','r')
    a2f = h5['results/electron_phonon/phonons/self_energy_1/a2F'][:]
    w= h5['results/electron_phonon/phonons/self_energy_1/frequency_grid'][:]

    a2f=a2f[0,0,:]
  #  plt.plot(w,a2f)
  #  plt.show()

    h=open('a2f.dat-vasp','w')
    for ni,i in enumerate(w):
     h.write(str(i)+' '+str(a2f[ni])+'\n')
    h.close()

def allen_dynes_tc(mus):
 h=open('a2f.dat-vasp','r')
 tmp=h.readlines ()[1:]
 tmp2=np.array([i.split() for i in tmp if float(i.split()[0])>0],dtype=float)
 W=tmp2[:,0]*47.9924/(2*np.pi)
 A2F=tmp2[:,1]
 a=np.trapz(A2F*np.log(W)/W,x=W)
 l=2*np.trapz(A2F/W,x=W)
 omlog=np.exp(2./l*a)
 for mu in mus:
  tc=omlog/1.2*np.exp(-1.04*(1+l)/(l-mu*(1+0.62*l))) # rownowazne wyrazeniu z ph$
  my_print('for mu*='+str(mu)+': lambda='+str(l)+' Tc='+str(tc)+' omlog='+str(omlog),2)

 


def run_phonopy(structure, incar, kpoints, kpar=1,max_q=8):
    my_print('ELPH CALC STARTED',0)
    fftmesh=read_fftmesh('relax')
    my_print('fft mesh',fftmesh,1)
    my_print('Preparing kpoints and incar',1)

    dir='ELPH'
    os.system('mkdir '+dir)
    os.system('cp POSCAR POTCAR '+dir)
    kp=np.array(kpoints.kpts[0],dtype=int)
    kp_dense=np.array([(int(m*2)) for m in (kp)])
    kp.write_file(dir+"/KPOINTS")
    kp_dense.write_file(dir+"/KPOINTS_DENSE") 
    incar_elph=incar.copy()
    incar_elph["ELPH_DRIVER"]='PH'
    incar_elph["ELPH_MODE"]='ELIASHBERG_ISO'
    incar_elph["ELPH_SELFEN_CARRIER_DEN"]=0
    incar_elph["ELPH_SELFEN_TEMPS"]=10
    incar.write_file(dir+"/INCAR")


    for q in range(1,max_q):
      dir_ph='q'+str(q)
      if len(glob.glob(dir_ph+'/total_dos.dat'))!=0: nq=q
    #nq=2
    dir_ph='q'+str(nq)
    my_print( "I run ELPH for phonons from ",dir_ph,1)
    dim=read_yaml(dir_ph)
    dispdirs=['../'+m for m in glob.glob(dir+'/disp*')]
    my_print('supercell size=',dim,'. Data will be read from ',dispdirs,1)
#os.system('cp '+dir+'/phonopy_disp.yaml ELPH/phelel_disp.yaml')
    phelel='cd '+dir+'; /fs/home/sylwia/src/miniconda3/bin/phelel -d --dim '+''.join([str(m) for m in dim])+' --amplitude 0.01; rm POSCAR*-*'
    os.system(phelel)
    phelel='cd '+dir+';/fs/home/sylwia/src/miniconda3/bin/phelel -c POSCAR --dim '+''.join([str(m) for m in dim])+' --create-derivatives ../'+dir+'/PERFECT '+' '.join(dispdirs)+' --fft-mesh '+' '.join(fftmesh)
    my_print(phelel,1)
    os.system(phelel)
    run_vasp('./ELPH/',24)
    read_a2f('./ELPH/')
    allen_dynes_tc([0.1,0.13])
    my_print("ELPH CALC FINISHED",0)