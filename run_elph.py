import os
import glob
import numpy as np
vasp_command='mpirun -np 24 /fs/home/sylwia/src/vasp-tmp1-gmatrix/IPA-AMD-OMP/std/vasp'
import sys
if len(sys.argv)>1: vasp_command=' '.join(sys.argv[1:])
print('vasp command:', vasp_command)




def run_vasp(dir,npc):
    vc=vasp_command #.replace('24',str(npc))
    os.system('a=$(pwd); echo $a; cd '+dir+'; '+vc+';cd $a')
    h=open(dir+'/OUTCAR')
    tmp=h.readlines()
    h.close()
    for m in tmp:
        if 'The electronic self-consistency was not achieved in the given' in m or 'EEEEEEE  RRRRRR   RRRRRR   OOOOOOO  RRRRRR' in m:
         return 0 
    return 1

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
  print('for mu*='+str(mu)+': lambda='+str(l)+' Tc='+str(tc)+' omlog='+str(omlog))

 



fftmesh=read_fftmesh('relax')
print(fftmesh)

os.system('mkdir ELPH; cp ../ins/INCAR.elph ELPH/INCAR; cp relax/CONTCAR ELPH/POSCAR; cp POTCAR ELPH/; cp relax/KPOINTS ELPH/')
os.system('grep ENCUT relax/INCAR >>ELPH/INCAR')
os.system("sed -i '3s/.*/Gamma/' ELPH/KPOINTS") #ELPH gives error with Monkhorst grid mthod
with open('ELPH/KPOINTS') as h:
    tmp=h.readlines()
    kp=[int(m) for m in tmp[3].split()]
kpd=[str(2*m) for m in kp]
h=open('ELPH/KPOINTS_DENSE','w')
for i in tmp[:3]:
    h.write(i)
h.write(' '.join(kpd)+'\n0 0 0')
h.close()


for q in range(1,9):
   dir='q'+str(q)
   if len(glob.glob(dir+'/total_dos.dat'))!=0: nq=q
#nq=3
print(nq)
dir='q'+str(nq)
dim=read_yaml(dir)
dispdirs=['../'+m for m in glob.glob(dir+'/disp*')]
print(dim)
#os.system('cp '+dir+'/phonopy_disp.yaml ELPH/phelel_disp.yaml')
phelel='cd ELPH; /fs/home/sylwia/src/miniconda3/bin/phelel -d --dim '+''.join([str(m) for m in dim])+' --amplitude 0.01; rm POSCAR*-*'
os.system(phelel)
phelel='cd ELPH;/fs/home/sylwia/src/miniconda3/bin/phelel -c POSCAR --dim '+''.join([str(m) for m in dim])+' --create-derivatives ../'+dir+'/PERFECT '+' '.join(dispdirs)+' --fft-mesh '+' '.join(fftmesh)
print(phelel)
os.system(phelel)
run_vasp('./ELPH/',24)
read_a2f('./ELPH/')
allen_dynes_tc([0.1,0.13])