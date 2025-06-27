import h5py
import matplotlib.pyplot as plt
import numpy as np

def print_attrs(name, obj):
    print (name)
    for key, val in obj.attrs.items():
        print( "    %s: %s" % (key, val))

h5 = h5py.File('vaspout.h5','r')
h5.visititems(print_attrs)
a2f = h5['results/electron_phonon/phonons/self_energy_1/a2F'][:]
w= h5['results/electron_phonon/phonons/self_energy_1/frequency_grid'][:]

a2f=a2f[0,0,:]
plt.plot(w,a2f)
plt.show()

h=open('a2f.dat-vasp','w')
for ni,i in enumerate(w):
 h.write(str(i)+' '+str(a2f[ni])+'\n')
h.close()





