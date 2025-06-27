from mp_api.client import MPRester
import numpy as np
import os
elements=["Ac","Ag","Al","Am","Ar","As","At","Au","B","Ba","Be","Bh","Bi","Bk","Br","C","Ca","Cd","Ce","Cf","Cl","Cm","Cn","Co","Cr","Cs","Cu","Db","Ds","Dy","Er","Es","Eu","F","Fe","Fl","Fm","Fr","Ga","Gd","Ge","H","He","Hf","Hg","Ho","Hs","I","In","Ir","K","Kr","La","Li","Lr","Lu","Lv","Mc","Md","Mg","Mn","Mo","Mt","N","Na","Nb","Nd","Ne","Nh","Ni","No","Np","O","Og","Os","P","Pa","Pb","Pd","Pm","Po","Pr","Pt","Pu","Ra","Rb","Re","Rf","Rg","Rh","Rn","Ru","S","Sb","Sc","Se","Sg","Si","Sm","Sn","Sr","Ta","Tb","Tc","Te","Th","Ti","Tl","Tm","Ts","U","V","W","Xe","Y","Yb","Zn","Zr"]

 

def download_poscar(formula,name=''):
 with MPRester("rsnyxATXEleWPypUdNDxszAd2BiKhqVu") as mpr:
  docs = mpr.summary.search(formula = formula)
  #print(docs)
  if len(docs)==0: return -1
  no=np.argsort([i.energy_above_hull for i in docs])[0]
  structure = docs[no].structure
  print('energy_above_hull',docs[no].energy_above_hull)
#  print(structure)
#  structure.to(fmt='poscar', filename='POSCAR'+name)
  return structure
 
 
h=open('from_pub.dat')
#data=np.loadtxt('from_pub.dat')
data=[i.split() for i in h.readlines()]
head=data[0]
data=data[1:]
h.close()
names=[i[0] for i in data]
data=np.array([i[1:] for i in data])
 

for i in names[1:]:
 structure=download_poscar(i)
 if structure==-1: 
  print(i+' not found in materialsproject')
  continue
# print(structure.elements)
# prepare_files(i,structure) 





