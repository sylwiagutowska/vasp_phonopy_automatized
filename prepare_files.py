from mp_api.client import MPRester
from common import *

elements=["Ac","Ag","Al","Am","Ar","As","At","Au","B","Ba","Be","Bh","Bi","Bk","Br","C","Ca","Cd","Ce","Cf","Cl","Cm","Cn","Co","Cr","Cs","Cu","Db","Ds","Dy","Er","Es","Eu","F","Fe","Fl","Fm","Fr","Ga","Gd","Ge","H","He","Hf","Hg","Ho","Hs","I","In","Ir","K","Kr","La","Li","Lr","Lu","Lv","Mc","Md","Mg","Mn","Mo","Mt","N","Na","Nb","Nd","Ne","Nh","Ni","No","Np","O","Og","Os","P","Pa","Pb","Pd","Pm","Po","Pr","Pt","Pu","Ra","Rb","Re","Rf","Rg","Rh","Rn","Ru","S","Sb","Sc","Se","Sg","Si","Sm","Sn","Sr","Ta","Tb","Tc","Te","Th","Ti","Tl","Tm","Ts","U","V","W","Xe","Y","Yb","Zn","Zr"]

 


def read_enmax(file):
   with open(file) as h:
     for l in h.readlines():
      if 'ENMAX' in l: return float(l.split()[2].replace(';',''))
       
def choose_potcar(i):
  if glob.glob('/fs/home/sylwia/src/potpaw_PBE/'+str(i)): 
   fil= '/fs/home/sylwia/src/potpaw_PBE/'+str(i)+'/POTCAR'
   return fil, read_enmax(fil)
  else: 
   files=glob.glob('/fs/home/sylwia/src/potpaw_PBE/'+str(i)+'*/POTCAR') 
   en=[read_enmax(m) for m in files]
   num=np.argsort(en)[0]
   return files[num],en[num]
    
def prepare_incar(structure,para):
 incar_dict = {
    "ISTART": 1,
     "ISPIN": 1,
    "ENCUT": 520,
    "EDIFF": 1e-7,
    "EDIFFG": -2e-5,
    "ISMEAR": 1,
    "SIGMA": 0.15,
    "IBRION": -1, 
    "ISIF": 3,
    "NSW": 1,
    "PREC": "Accurate",
    "LREAL": "FALSE",
    "ALGO": "Normal", 
   
 } 
 for i in para.keys(): incar_dict[i]=para[i]
 incar = Incar(incar_dict)
 incar.write_file("INCAR")
 return incar


def download_poscar(formula,name=''):
 my_print('PREPARING THE POSCAR',0)
 with MPRester("rsnyxATXEleWPypUdNDxszAd2BiKhqVu") as mpr:
  docs = mpr.summary.search(formula = formula)
 # print(docs.material_id)
  if len(docs)==0: return -1,-1
  no=np.argsort([i.energy_above_hull for i in docs])[0]
  structure = docs[no].structure
  my_print('energy_above_hull',docs[no].energy_above_hull,1)
  print(structure)
  structure.to(fmt='poscar', filename='POSCAR'+name)
  my_print('POSCAR PREPARED',0)
  return structure,docs[no].energy_above_hull
 
def prepare_potcar_incar(formula,structure,energy_above_hull):
 my_print('PREPARING POTCAR AND INCAR',0)
 #read poscar
 structure_file = os.path.join("POSCAR")
 structure = Structure.from_file(structure_file)
 
 #prepare potcar
 enmaxs=[]
 potcar_file='POTCAR'
 os.system('rm '+potcar_file)
 for i in structure.elements: 
  fil,en=choose_potcar(i)
  os.system('cat '+fil+' >> '+potcar_file)
  enmaxs.append(en)
  my_print(f"chosen potcar {fil}",1)
 encut=np.ceil(np.min(enmaxs))
 potcar = Potcar.from_file(potcar_file)

 my_print(f"chosen ENCUT {encut}",1) 
 #prepare incar
 dic={"ENCUT": encut }
 if np.any([site.specie.Z>73 for site in structure.sites]):  #higher than Ta
  my_print('heave element (Z>73) detected, i add spin-orbit coupling',1)
  dic["LSORBIT"]= ".TRUE."; dic["LNONCOLLINEAR"]= ".TRUE."; dic["SAXIS"]=[0, 0, 1]
  dic["MAGMOM"]=(int(len(structure.sites)))*"0 "
 incar=prepare_incar(structure,dic)
 my_print('FILES PREPARTED',0)
 
  
 #os.system('mv {POSCAR,POTCAR,INCAR} '+formula+'/.')

 

 return potcar,incar






