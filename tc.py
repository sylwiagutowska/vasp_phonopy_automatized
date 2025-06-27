import os
#import math
import numpy as np

mu=float(input('Enter mu* value: \n'))
#read a2f's
 



#calc. of lambda and Tc for all 10 a2f functions
if 1==1:
 h=open('a2f.dat-vasp','r')
 tmp=h.readlines ()[1:]
 tmp2=np.array([i.split() for i in tmp if float(i.split()[0])>0],dtype=float)
 W=tmp2[:,0]*47.9924/(2*np.pi)
 A2F=tmp2[:,1]

# W=
 a=np.trapz(A2F*np.log(W)/W,x=W)
 l=2*np.trapz(A2F/W,x=W)
 omlog=np.exp(2./l*a)
 tc=omlog/1.2*np.exp(-1.04*(1+l)/(l-mu*(1+0.62*l))) # rownowazne wyrazeniu z ph$
 print('lambda='+str(l)+' Tc='+str(tc)+' omlog='+str(omlog))

 
