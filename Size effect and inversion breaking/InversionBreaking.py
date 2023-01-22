#imports
import numpy as np
from matplotlib import pyplot as plt
from javelin.structure import Structure
from javelin.unitcell  import UnitCell
from javelin.fourier   import Fourier 
from custom.inout import write_structure 
import BvsOptimizer as bvo
import multiprocessing as mp


"""
Modeling the polar (inversion breaking) distortion of the Bi-O octahedra
by a simple model based onstrained MC minimization of valence 
mismatch in BKBO (bond valence sum theory). Model is constrained by 
simple rigidity and volume exclusion constraints.  

This code peforms the simulation on samples with different potassium doping 
levels and "measures" the doping dependence of the centrosymmetry parameter. 

Output: Plot containing the doping dependence of the centrosymmetry parameter. 

"""


def lc(x): 
    return 4.3548 - 0.1743*x

def CreateStructure(ncells = 10, doping = 0): 
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
        T=1
        print("Cycle: ", i) 
        sim.run_cycle(str, 0.005 ,1, T)
        print("Global instability index:", str.global_instability_index())
        print("O rms displcement: ", str.rms_displacement(1))
        print("Ba valence:", str.get_average_valence(56))
        print("K  valence: ", str.get_average_valence(19)) 
        print("Bi valence: ",str.get_average_valence(83))
        print("Ba-O average distance: ", (str.average_distance(56,8, Ba_bonds.values) - av2)/lc(doping)) 
        print("K-O: average distance: ", (str.average_distance(19,8, Ba_bonds.values) - av2)/lc(doping))
        print("O-O: average distance: ",str.average_distance(8,8, (OO_1 + OO_2 + OO_3).values)/lc(doping)) 
        cs = str.CS_parameter(Bi_bonds.values)
        print("CS_parameter:", cs)
 
    return from_bvo(str)


structure = CreateStructure(10, 0)
Bi_bonds = structure.get_neighbors(0,1) + structure.get_neighbors(0,2) + structure.get_neighbors(0,3)

  
""" Running the simualtion """
doping_samples = 10
measurement_samples = 10
centrosymmetry = []
x = np.linspace(0.2, 0.5, doping_samples )
for i in range(doping_samples): 
    samples = [CreateStructure(15, x[i]) for j in range(measurement_samples)]
    results = []
    count = 0
    with mp.Pool() as pool:
        for result in pool.map(BvsOptimization, samples):
            results.append(result)
            count = count + 1
    results = [to_bvo(r) for r in results]
    results = [r.CS_parameter(Bi_bonds.values) for r in results] 
    centrosymmetry.append(sum(results)/measurement_samples)

centrosymmetry = np.array(centrosymmetry)

plt.figure(figsize = [8,4])
plt.plot(x, centrosymmetry, color="orange", marker='o', linestyle='-', mew=2, lw = 2, mfc='w', markersize=8)
plt.xlabel("Doping", fontsize = 12)
plt.ylabel("Centrosymmetry parameter (arb. units)", fontsize = 12)
plt.show() 
