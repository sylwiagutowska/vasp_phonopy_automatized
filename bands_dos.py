import os
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar
from pymatgen.io.vasp import Vasprun, BSVasprun
from pymatgen.electronic_structure.plotter import BSDOSPlotter
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.core import Spin,Orbital
import matplotlib.pyplot as plt
from common import *
import xml.etree.ElementTree as ET


def write_incar_static(incar,dir):
    incar["LWAVE"]= True
    incar["LCHARG"]= True
    incar.write_file(dir+'/INCAR')

def write_incar_bands(incar,dir):
    incar_bands=incar.copy()
    incar_bands["LWAVE"]= False
    incar_bands["LCHARG"]= False
    incar_bands["ICHARG"]=11
    incar_bands["ISMEAR"]= -5
    incar_bands.write_file(dir+'/INCAR')

def write_incar_dos(incar,dir):
    incar_dos=incar.copy()
    incar_dos["LWAVE"]= False
    incar_dos["LCHARG"]= False
    incar_dos["ICHARG"]=11
    incar_dos["ISMEAR"]= -5
    incar_dos["NEDOS"]= 2000
    incar_dos["LORBIT"]= 11
    incar_dos.write_file(dir+'/INCAR')

def find_enrange(dir):
  tree = ET.parse(dir+'/vasprun.xml')
  root = tree.getroot()
  ENE_v=[]
  a=root.find('calculation/eigenvalues/array/set/set')
  for i in a: 
    if i.tag=='set' and 'kpoint' in i.attrib['comment']: 
      ENE_v.append([])
      for j in i.findall('r'):
        ENE_v[-1].append(j.text.split()[0])

  ENE_v=np.array(ENE_v,dtype=float) #.transpose()

  a=root.find('calculation/dos/i')
  ef=float(a.text.split()[0])

  ENE_v-=ef
  
  erange=[-4,4]

  for e in np.arange(0,20,0.2):
     if not np.any(np.round(ENE_v-e) == 0): 
        erange[1]=np.round(e,1)
        break
  for e in np.arange(0,20,0.2):
     if not np.any(np.round(ENE_v+e) == 0): 
        erange[0]=np.round(e,1)
        break
   
  return erange

def save_bands(dir='bands'):
  # Load vasprun.xml and get band structure
  vasprun = BSVasprun(dir+"/vasprun.xml", parse_projected_eigen=True)
  bs = vasprun.get_band_structure(kpoints_filename="KPOINTS", line_mode=True)
  fermi = bs.efermi
  # Get band energies (assuming spin-unpolarized)
  bands = bs.bands
  spin = list(bands.keys())[0]  # Spin.up or Spin.down
  band_energies = bands[spin] - fermi  # Shape: (n_bands, n_kpoints)
  # Transpose to shape (n_kpoints, n_bands)
  band_energies = band_energies.T
  # Get k-point coordinates
  kpoints = np.array([kp.frac_coords for kp in bs.kpoints])
  n_kpoints = len(kpoints)
  # Build header
  header = "k_index  kx  ky  kz  " + "  ".join([f"Band{i+1}" for i in range(band_energies.shape[1])])
  # Prepare data
  k_indices = np.arange(n_kpoints).reshape(-1, 1)
  data = np.hstack((k_indices, kpoints, band_energies))
  # Save to file
  np.savetxt(dir+"/bands.dat", data, header=header, fmt="%.6f")


def save_dos(dir='dos'):
  vasprun = Vasprun(dir+"/vasprun.xml")
  dos = vasprun.complete_dos  # Or vasprun.tdos for total DOS only
  energies = dos.energies - vasprun.efermi  # Align with Fermi level
  print(dos.densities.keys())
  total_dos = dos.densities[Spin.up]        
  if Spin.down in dos.densities:
    total_dos += dos.densities[Spin.down]
  # Save energies and total DOS to file
  np.savetxt(dir+"/dos.dat", np.column_stack((energies, total_dos)), header="E-E_F (eV)    DOS (1/eV) integral", fmt="%.6f")

def save_pdos(dir='dos'):
# Load vasprun.xml
  vasprun = Vasprun(dir+"/vasprun.xml", parse_projected_eigen=True)
  cdos = vasprun.complete_dos
  energies = cdos.energies - vasprun.efermi  # Shift energies to Fermi level
  efermi_index=np.argsort(np.abs(energies))[0] #ene closest to EF
  # Initialize data: first column is energy
  data = [energies]
  headers = ["Energy"]
  print("PDOS available for sites:", list(cdos.pdos.keys()))
  # Loop over all atoms and orbitals
  # Iteruj po atomach, dla których mamy PDOS
  for site, spin_orb_dict in cdos.pdos.items():
    print(f"Atom {site.specie} at {site.frac_coords}")

    for spin in spin_orb_dict:
        for orb, dos in spin_orb_dict[spin].items():
            print(f"  Orbital {orb.name}, Spin: {spin.name},  DOS(EF): {dos[efermi_index]:.3f}")
            data.append(dos)
            headers.append(f"{site.specie}_{orb.name}") 
  # Transpose data to get one row per energy point
  data = np.array(data).T
  # Save to file
  np.savetxt(dir+"/projected_dos.dat", data, header="  ".join(headers), fmt="%.6f")


def plot_bands_dos(bands_vasprun_path='bands', dos_vasprun_path='dos', output="bands_dos.png"):
    # Load Band structure
    bs_vrun = BSVasprun(bands_vasprun_path+'/vasprun.xml', parse_projected_eigen=True)
    bs = bs_vrun.get_band_structure(line_mode=True)

    # Load DOS
    dos_vrun = Vasprun(dos_vasprun_path+'/vasprun.xml')
    dos = dos_vrun.complete_dos
    
    [emin,emax]=find_enrange(dos_vasprun_path)
    # Plot Bands + DOS + PDOS
    plotter = BSDOSPlotter(bs_projection="elements", dos_projection="elements", vb_energy_range=abs(emin), cb_energy_range=(emax))
    plot = plotter.get_plot(bs, dos)
    plt.tight_layout()
    plt.savefig(output, dpi=300)
 #   plt.show()


def scf_bands_dos(structure,incar,kpoints):
    # SCF run
    dir='SCF'
    os.system('mkdir '+dir)
    os.system('cp POSCAR POTCAR '+dir)
    write_incar_static(incar,dir)
    kpoints.write_file(f"{dir}/KPOINTS")
    incar.write_file(f"{dir}/INCAR")
    run_vasp(dir)

    # Band structure run
    dirb='bands'
    os.system(f'rm -r {dirb}; cp -r SCF {dirb}')
    write_incar_bands(incar,dirb)
    kpath = HighSymmKpath(structure)
    # Create KPOINTS file for band structure along that path
    kpoints_bands = Kpoints.automatic_linemode(divisions=50, ibz=kpath)
    kpoints_bands.write_file(f"{dirb}/KPOINTS")
    run_vasp(dirb)
    save_bands(dirb)

    # DOS run (dense mesh)
    dird='dos'
    os.system(f'rm -r {dird}; cp -r SCF {dird}')
    write_incar_dos(incar,dird)
    kpoints_dos = Kpoints.automatic(' '.join([str(4*m) for m in kpoints.kpts[0]]))
    run_vasp(dird)
    save_dos(dird)
    save_pdos(dird)

    
 
