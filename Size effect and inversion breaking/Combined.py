#imports
import numpy as np
from javelin.structure import Structure
from javelin.unitcell  import UnitCell
from javelin.modifier  import SetDisplacementNormalXYZ
from javelin.modifier  import SetDisplacementNormal
from custom.energies  import MorseEnergy
from custom.energies  import YukawaEnergy
from javelin.mc import MC
import BvsOptimizer as bvo
import multiprocessing as mp
from time import sleep

"""
Combines SizeEffect.py and InversionBreaking.py into a single program 
which produces the supercells with both size-effect and inversion 
breaking disorder present. 
Output:   - Formatted .txt files which can be read using read_structure function 
               from custom.inout for further analysis and simulation 
          - Formatted Scatty input files for diffuse scattering calculations
          - Plots the distribution of BiBa, BiK, BaO, KO and OO pair distances.

"""


def lc(x): 
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
    cubic_cell = UnitCell(lc(doping))
    BaK = np.random.choice([19,56], size =ncells**3 , p=[doping, 1-doping])
    structure = Structure(numbers=[83,8,8,8,1]*ncells**3,
    positions = pos*ncells**3, unitcell=cubic_cell)
    structure.reindex([ncells,ncells,ncells,5])
    structure.get_atomic_numbers()[4::5] = BaK
    structure.update_atom_symbols()
    return structure


def to_bvo(structure = None): 
    shape = structure.atoms.index.to_numpy().max() + np.array(1)
    str = bvo.Crystal(shape[0], shape[-1],structure.get_cell())
    str.load(structure.get_atomic_numbers(), structure.xyz)
    return str
    
def from_bvo(str = None): 
    structure = Structure(numbers = str.get_atoms() , positions=str.get_positions(), 
    unitcell=str.get_cell()[0][0])
    structure.reindex(str.get_shape())
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

    E_BiO =  MorseEnergy(D = 10,  b = 0.318/lc(dop), desired = 1/2 ,         atom_type1 = 83, atom_type2 = 8)
    E_ko =   MorseEnergy(D = 2 ,  b = 0.398/lc(dop), desired = np.sqrt(2)/2, atom_type1 = 19, atom_type2 = 8)
    E_BaO =  MorseEnergy(D = 2,   b = 0.406/lc(dop), desired = np.sqrt(2)/2, atom_type1 = 56, atom_type2 = 8)


    # Yukawa interactions -> Effective charges taken to be proportional to average oxidation states 
    k = 40
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


def BvsOptimization(structure = None): 
    str = to_bvo(structure)
    doping = 1-structure.get_average_site(4)['Ba']['occ']
    Bi_bonds = structure.get_neighbors(0,1) + structure.get_neighbors(0,2) + structure.get_neighbors(0,3)
    Ba_bonds = structure.get_neighbors(4,1) + structure.get_neighbors(4,2) + structure.get_neighbors(4,3)
    O1_bonds = structure.get_neighbors(1,0) + structure.get_neighbors(1,4)
    O2_bonds = structure.get_neighbors(2,0) + structure.get_neighbors(2,4)
    O3_bonds = structure.get_neighbors(3,0) + structure.get_neighbors(3,4)

    OO_1 = structure.get_neighbors(1,2) + structure.get_neighbors(1,3) 
    OO_2 = structure.get_neighbors(2,1) + structure.get_neighbors(2,3) 
    OO_3 = structure.get_neighbors(3,1) + structure.get_neighbors(3,2) 


    str.set_bonds(0, Bi_bonds.values)
    str.set_bonds(4, Ba_bonds.values)
    str.set_bonds(1, O1_bonds.values)
    str.set_bonds(2, O2_bonds.values)
    str.set_bonds(3, O3_bonds.values)

    str.set_parameters(83,8,2.050,0.318)
    str.set_parameters(56,8,2.223,0.406)
    str.set_parameters(19,8,2.047,0.398)

    str.set_nominals_and_weights([83,56,19,8] ,[4+doping,2,1,2],[0,0,0,0], [1,1,1,1] , [1,1,1,1])
    sigma_Bi = 0.02
    sigma_Ba = 0.01
    sigma_O =  0.06
    

    mod_Bi  = bvo.Generator(0, [0,0,0] ,[sigma_Bi, sigma_Bi, sigma_Bi])
    mod_Ba  = bvo.Generator(4, [0.5,0.5,0.5], [sigma_Ba, sigma_Ba, sigma_Ba])
    mod_1  =  bvo.Generator(1, [0.5,0,0] ,    [0, sigma_O, sigma_O])
    mod_2  =  bvo.Generator(2, [0,0.5,0] ,    [sigma_O, 0, sigma_O])
    mod_3  =  bvo.Generator(3, [0,0,0.5] ,    [sigma_O, sigma_O, 0])


    
    D = 600
    R_BaO = np.sqrt(2)*lc(doping)/2 - 0.1*lc(doping)
    R_KO =  np.sqrt(2)*lc(0.45)/2
    R_BiO = 0
    K = 1000
    d = lc(doping)*np.sqrt(2)/2


    calc1 = bvo.Calculator(0, (Bi_bonds).values)
    calc1.add_interaction(83,8,"SoftVolumeExclusion", [D, R_BiO])

    calc2 = bvo.Calculator(4, (Ba_bonds).values)
    calc2.add_interaction(56,8,"SoftVolumeExclusion", [D, R_BaO])
    calc2.add_interaction(19,8,"SoftVolumeExclusion", [D, R_KO])
    
    
    calc3 = bvo.Calculator(1, (O1_bonds + OO_1).values)
    calc3.add_interaction(8,56,"SoftVolumeExclusion", [D, R_BaO])
    calc3.add_interaction(8,19,"SoftVolumeExclusion", [D, R_KO])
    calc3.add_interaction(8,83,"SoftVolumeExclusion", [D, R_BiO]) 
    calc3.add_interaction(8,8,"Spring", [K, d])

    calc4 = bvo.Calculator(2, (O2_bonds + OO_2).values)
    calc4.add_interaction(8,56,"SoftVolumeExclusion", [D, R_BaO])
    calc4.add_interaction(8,19,"SoftVolumeExclusion", [D, R_KO])
    calc4.add_interaction(8,83,"SoftVolumeExclusion", [D, R_BiO]) 
    calc4.add_interaction(8,8,"Spring", [K, d])


    calc5 = bvo.Calculator(3, (O3_bonds + OO_3).values)
    calc5.add_interaction(8,56,"SoftVolumeExclusion", [D, R_BaO])
    calc5.add_interaction(8,19,"SoftVolumeExclusion", [D, R_KO])
    calc5.add_interaction(8,83,"SoftVolumeExclusion", [D, R_BiO]) 
    calc5.add_interaction(8,8,"Spring", [K, d])


    sim = bvo.Simulation()
    sim.add_generator(mod_1)
    sim.add_generator(mod_2)
    sim.add_generator(mod_3)
    sim.add_calculator(calc3)
    sim.add_calculator(calc4)
    sim.add_calculator(calc5)


    sim.optimize_vector_valence = False
    sim.optimize_valence =  True

    av1 = lc(doping)*np.sqrt(3)/2
    av2 = lc(doping)*np.sqrt(2)/2

    cycles = 200
    for i in range(cycles):
        #T = 20*np.exp(-i/40) + 1
        T=.1
        print("Cycle: ", i) 
        sim.run_cycle(str, 0.005 ,1, T)

        print("Global instability index:", )
        print("Ba valence:", str.get_average_valence(56))
        print("K  valence: ", str.get_average_valence(19)) 
        print("Bi valence: ",str.get_average_valence(83))
    return from_bvo(str)
   
""" Running the simualtion """
n_samples = 10
samples = [CreateStructure(15, 0.35) for i in range(n_samples)]

temp = []
results = []
count = 0
with mp.Pool() as pool:
    # call the function for each item in parallel
    for result in pool.map(Yukawa_simulation, samples):
        temp.append(result)
        count = count + 1

with mp.Pool() as pool:
    # call the function for each item in parallel
    for result in pool.map(BvsOptimization, temp):
        results.append(result)
        count = count + 1

structure = CreateStructure(10, 0.1)

from custom.inout import write_structure


"""Writing structures to a formatted .txt file"""
for n in range(n_samples): 
    write_structure(results[n],"Samples/yukawa_model"+str(n)+".txt")


"""Calculating standard deviations of displacements and pair distances. """


BiBa = structure.get_neighbors(0,4)
Ba_bonds = structure.get_neighbors(4,1) + structure.get_neighbors(4,2)+  structure.get_neighbors(4,3)

OO_1 = structure.get_neighbors(1,2) + structure.get_neighbors(1,3) 
OO_2 = structure.get_neighbors(2,1) + structure.get_neighbors(2,3) 
OO_3 = structure.get_neighbors(3,1) + structure.get_neighbors(3,2) 
OO = OO_1 + OO_2 + OO_3

U_Bi = 0.0 
U_Ba = 0.0
U_K = 0.0
U_O = 0.0
BiK_distances =  np.array([])
BiBa_distances = np.array([])
BaO_distances  = np.array([])
KO_distances  =  np.array([])
OO_distances =   np.array([]) 


import custom.analysis as analysis
for sample in results: 
    U_Bi += analysis.std_dev(sample, site = 0, atom = 83)
    U_Ba += analysis.std_dev(sample, site = 4, atom = 56)
    U_K  += analysis.std_dev(sample, site = 4, atom = 19)
    U_O +=  analysis.std_dev(sample, site = 1, atom = 8)

    BiK_distances = np.concatenate([BiK_distances,  analysis.get_distances_bw_atoms(sample,BiBa, 83, 19)], axis = 0)
    BiBa_distances= np.concatenate([BiBa_distances, analysis.get_distances_bw_atoms(sample,BiBa, 83, 56)], axis = 0)
    BaO_distances = np.concatenate([BaO_distances,  analysis.get_distances_bw_atoms(sample,Ba_bonds, 19, 8)], axis = 0)
    KO_distances=   np.concatenate([KO_distances,   analysis.get_distances_bw_atoms(sample,Ba_bonds, 56, 8)], axis = 0)
    OO_distances=   np.concatenate([OO_distances,   analysis.get_distances_bw_atoms(sample,OO,8, 8)], axis = 0)

U_Bi = U_Bi/n_samples 
U_Ba = U_Ba/n_samples
U_K  = U_K/n_samples
U_O =  U_O/n_samples

print("Bi rms displacement: ", U_Bi )
print("Ba rms displacement: ", U_Ba )
print("K rms displacement: ",  U_K )
print("O_1 rms displacement: ",U_O)


print("Mean BiBa distance: ",  np.mean(BiBa_distances))
print("Mean BiK  distance: ",  np.mean(BiK_distances))
print("Mean BaO  distance: ",  np.mean(BaO_distances))
print("Mean KO  distance: ",   np.mean(KO_distances))
print("Mean OO  distance: ",   np.mean(OO_distances))



""" Writing scatty input files for diffuse scattering calculation. """
from custom.inout import Scatty as SC 
sc = SC() 
sc.basis = [(0.0,0.0,0.0),(0.5,0.0,0.0),(0.0,0.5,0.0),(0.0,0.0,0.5),(0.5,0.5,0.5)]
sc.centre = (0,0,0)
sc.x_axis = (10,0,0,500)
sc.y_axis = (0,10,0,500)
sc.z_axis = (0,0,0,1)
sc.cell_centering = "P"
sc.remove_bragg = True
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
Plotting pair distance distributions.
"""
plt.figure(figsize=[8,6])
plt.hist(BiBa_distances, bins = 90 , density = True, color = "green", fill = False, 
label="Bi-Ba pairs", histtype = "step", lw = 3)
plt.hist(BiK_distances , bins = 90 , density = True, color = "orange",fill = False, 
label="Bi-K pairs",  histtype = "step", lw = 3)
plt.axvline(np.sqrt(3)/2 , 0, 1, label='', ls = ':', color = "firebrick", lw = 2.5 )
plt.legend(edgecolor = "b", fancybox = True ,fontsize = 15)
plt.tick_params(axis = 'both', labelsize=15)
plt.xlabel("Distance (l.u.)", fontsize = 15)
plt.ylabel("Pair distribution", fontsize = 15)
plt.ylim((0,50))
plt.show() 


plt.figure(figsize=[8,6])
plt.hist(BaO_distances, bins = 90 , density = True, color = "green", fill = False, 
label="Ba-O pairs", histtype = "step", lw = 3)
plt.hist(KO_distances , bins = 90 , density = True, color = "orange",fill = False, 
label="K-O pairs",  histtype = "step", lw = 3)
plt.axvline(np.sqrt(2)/2 , 0, 1, label='', ls = ':', color = "firebrick", lw = 2.5 )
plt.legend(edgecolor = "b", fancybox = True ,fontsize = 15)
plt.tick_params(axis = 'both', labelsize=15)
plt.xlabel("Distance (l.u.)", fontsize = 15)
plt.ylabel("Pair distribution", fontsize = 15)
plt.ylim((0,50))
plt.show() 


plt.figure(figsize=[8,6])
plt.hist(OO_distances, bins = 90 , density = True, color = "orange", fill = False, 
label="O-O pairs", histtype="barstacked")
plt.axvline(np.sqrt(2)/2 , 0, 1, label='', ls = ':', color = "firebrick", lw = 2.5 )
plt.legend(edgecolor = "b", fancybox = True ,fontsize = 15)
plt.tick_params(axis = 'both', labelsize=15)
plt.xlabel("Distance (l.u.)", fontsize = 15)
plt.ylabel("Pair distribution", fontsize = 15)
plt.ylim((0,50))
plt.show() 

