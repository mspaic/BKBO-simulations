#imports
import numpy as np
from javelin.structure import Structure
from javelin.unitcell  import UnitCell
from javelin.modifier  import SetDisplacementNormalXYZ
from javelin.modifier  import SetDisplacementNormal
from custom.energies  import MorseEnergy
from custom.energies  import YukawaEnergy
from javelin.mc import MC
import multiprocessing as mp 


""" 
   Contains a Javelin implementation of Yukawa based model for 
   observed size-effect distortion involving Bi and Ba/K sublattices. Program performs 
   MC relaxation on a specified number of independent samples in parallel.
   Output: - Formatted .txt files which can be read using read_structure function 
             from custom.inout for further analysis and simulation 
           - Formatted Scatty input files for diffuse scattering calculations 
           - Plots the distribution of Bi-Ba and Bi-K pair distances 

"""


def lattice_constant(x): 
    """
    Empirical doping dependence of the 
    pseudocubic lattice constant. 
    """
    return 4.3548 - 0.1743*x


def CreateStructure(ncells = 10, doping = 0): 
    """
    Creates a javelin.structure supercell with 'ncells' unitcells 
    in each spatial direction. Potassium dopant is randomly
    distributed throughout the supercell with composition ratio 
    given by 'doping' argument. 
    """
    pos = [(0,0,0),(0.5,0,0),(0,0.5,0),(0,0,0.5),(0.5,0.5,0.5)]
    cubic_cell = UnitCell(lattice_constant(doping))
    BaK = np.random.choice([19,56], size =ncells**3 , p=[doping, 1-doping])
    structure = Structure(numbers=[83,8,8,8,1]*ncells**3,
    positions = pos*ncells**3, unitcell=cubic_cell)
    structure.reindex([ncells,ncells,ncells,5])
    structure.get_atomic_numbers()[4::5] = BaK
    structure.update_atom_symbols()
    return structure



def Yukawa_simulation(structure = None):
    """ Performs a Monte Carlo simulation of the model"""
    dop = structure.get_average_site(4)['K']['occ']

    Range = 1
    sum_distance = 3*Range

    # Defining neighbors 
    BaBi = structure.get_neighbors(4,0, maxD = sum_distance)
    BiBa = structure.get_neighbors(0,4, maxD = sum_distance)
    BaBa = structure.get_neighbors(4,4, maxD = sum_distance)
    BiBi = structure.get_neighbors(0,0, maxD = sum_distance)

    BiO =  structure.get_neighbors(0,1) + structure.get_neighbors(0,2) + structure.get_neighbors(0,3)  
    BaO =  structure.get_neighbors(4,1) + structure.get_neighbors(4,2) + structure.get_neighbors(4,3)

    E_BiO =  MorseEnergy(D = 10,  b = 0.318/lattice_constant(dop), desired = 1/2 ,         atom_type1 = 83, atom_type2 = 8)
    E_ko =   MorseEnergy(D = 2 , b = 0.398/lattice_constant(dop), desired = np.sqrt(2)/2, atom_type1 = 19, atom_type2 = 8)
    E_BaO =  MorseEnergy(D = 2, b = 0.406/lattice_constant(dop), desired = np.sqrt(2)/2, atom_type1 = 56, atom_type2 = 8)


    # Yukawa interactions -> Effective charges taken to be proportional to average oxidation states 
    k = 30
    E_BiBi =  YukawaEnergy(length = Range, g= k*(4+dop)*(4+dop), atom_type1=83 , atom_type2=83 )
    E_BiBa =  YukawaEnergy(length = Range, g= k*(4+dop)*2      , atom_type1=83 , atom_type2=56 )
    E_BiK  =  YukawaEnergy(length = Range, g= k*(4+dop)*1      , atom_type1=83 , atom_type2=19 )
    E_BaK  =  YukawaEnergy(length = Range, g= k*2              , atom_type1=56 , atom_type2=19 )
    E_BaBa  = YukawaEnergy(length = Range, g= k*4              , atom_type1=56 , atom_type2=56 )
    E_kk    = YukawaEnergy(length = Range, g= k*1              , atom_type1=19 , atom_type2=19 )


    Size_effect = MC() 
    Size_effect.temperature = .1
    Size_effect.cycles = 300

    # Defining modifiers used to generate MC moves 
    Sigma_Bi = 0.02
    Sigma_Ba = 0.01

    Size_effect.add_modifier(SetDisplacementNormal(0,0.0,Sigma_Bi))
    Size_effect.add_target(BiO,   E_BiO)

    Size_effect.add_target(BiBi,  E_BiBi)
    Size_effect.add_target(BiBa,  E_BiBa)
    Size_effect.add_target(BiBa,  E_BiK)



    Size_effect.add_modifier(SetDisplacementNormalXYZ(4,0.5,Sigma_Ba, 0.5, Sigma_Ba, 0.5, Sigma_Ba)) 
    Size_effect.add_target(BaO ,  E_BaO)
    Size_effect.add_target(BaO ,  E_ko)


    Size_effect.add_target(BaBa,  E_BaBa)
    Size_effect.add_target(BaBa,  E_BaK)
    Size_effect.add_target(BaBa,  E_kk)
    Size_effect.add_target(BaBi,  E_BiBa)
    Size_effect.add_target(BaBi,  E_BiK)
    structure = Size_effect.run(structure)
    return structure


   
""" Running the simualtion """
n_samples = 10
samples = [CreateStructure(15, 0.35) for i in range(n_samples)]

results = []
count = 0
with mp.Pool() as pool:
    # call the function for each item in parallel
    for result in pool.map(Yukawa_simulation, samples):
        results.append(result)
        count = count + 1

from custom.inout import write_structure


"""Writing structures to a formatted .txt file"""
for n in range(n_samples): 
    write_structure(results[n],"Samples/yukawa_model"+str(n)+".txt")


"""Calculating standard deviations of displacements and pair distances. """
U_Bi = 0.0 
U_Ba = 0.0
U_K = 0.0
BiK_distances = np.array([])
BiBa_distances = np.array([])

n = results[0].get_neighbors(0,4)
import custom.analysis as analysis
for sample in results: 
    U_Bi += analysis.std_dev(sample, site = 0, atom = 83)
    U_Ba += analysis.std_dev(sample, site = 4, atom = 56)
    U_K  += analysis.std_dev(sample, site = 4, atom = 19)

    BiK_distances = np.concatenate([BiK_distances, analysis.get_distances_bw_atoms(sample,n, 83, 19)], axis = 0)
    BiBa_distances= np.concatenate([BiBa_distances, analysis.get_distances_bw_atoms(sample,n, 83, 56)], axis = 0)


U_Bi = U_Bi/n_samples 
U_Ba = U_Ba/n_samples
U_k  = U_K/n_samples

print("Bi rms displacement: ", U_Bi )
print("Ba rms displacement: ", U_Ba )
print("K rms displacement: ",  U_k )

print("Mean BiBa distance: ",  np.mean(BiBa_distances))
print("Mean BiK  distance: ",  np.mean(BiK_distances))



""" Writing scatty input files for diffuse scattering calculation. """
from custom.inout import Scatty as SC 
sc = SC() 
sc.basis = [(0.0,0.0,0.0),(0.5,0.0,0.0),(0.0,0.5,0.0),(0.0,0.0,0.5),(0.5,0.5,0.5)]
sc.centre = (0,0,0)
sc.x_axis = (6,0,0,600)
sc.y_axis = (0,6,0,600)
sc.z_axis = (0,0,0,1)
sc.remove_bragg = True
sc.cell_centering = "P"
sc.expansion_error = 0.001
sc.expansion_order = 10
sc.window = 2
sc.cutoff = 2
sc.ppm_output = True
sc.ppm_colourmap = "viridis"
sc.radiation = "X"
sc.symmetry = "m3m"
sc.path = "Samples/"
sc.title = "BKBO"
sc.name = "HK"
sc.write_input(results)
from matplotlib import pyplot as plt 


"""
Plotting pair distances for Bi-Ba and Bi-K pairs 
in the form of a histogram. 
"""
plt.figure(figsize=[8,11])
plt.hist(BiBa_distances, bins = 90 , density = True, color = "green", fill = False, 
label="Bi-Ba pairs", histtype = "step", lw = 3)
plt.hist(BiK_distances , bins = 90 , density = True, color = "orange",fill = False, 
label="Bi-K pairs",  histtype = "step", lw = 3)
plt.axvline(np.sqrt(3)/2 , 0, 1, label='', ls = ':', color = "firebrick", lw = 2.5 )
plt.legend(edgecolor = "b", fancybox = True ,fontsize = 15)
plt.tick_params(axis = 'both', labelsize=15)
plt.xlabel("Distance (l.u.)", fontsize = 15)
plt.ylabel("Pair distribution (arb. units)", fontsize = 15)
plt.ylim((0,50))
plt.show() 


