from javelin.structure import Structure
import numpy as np
import copy


def get_average_site_positions(structure=None, digits = 4): 
    """
    Returns the list of overage site positions. 
    """
    str_average = structure.get_average_structure()
    x,y,z = 0.0,0.0,0.0
    average_pos=[]
    count_atoms = 0
    for site in str_average.keys():

        x,y,z = 0.0,0.0,0.0
        count_atoms = 0

        for atom in str_average[site].keys(): 
            x += str_average[site][atom]['x']
            y += str_average[site][atom]['y']
            z += str_average[site][atom]['z']
            count_atoms += 1
            
        x = x/count_atoms
        y = y/count_atoms
        z = z/count_atoms
        
        average_pos.append(
        (round(x, digits), round(y,digits),round(z,digits)))
    return average_pos

def get_select_positions(structure = None, site = 0, atom = None, scaled = False): 
    """
    Returns an array of unit cell positions for a specified 
    atom and site. If atom is not given, returns positions 
    of all atoms at a given site. 
    """
    N_sites = len(structure.get_average_structure().keys()) 
    if scaled == True: 
        pos_by_site = structure.get_scaled_positions()[site::N_sites]
    else:
        pos_by_site = structure.xyz[site::N_sites]
    Zs = structure.get_atomic_numbers()[site::N_sites]
    if atom == None:  
        return pos_by_site
    else:
        mask = (Zs == atom) 
        pos_by_atoms = pos_by_site[mask]
        return pos_by_atoms 

def std_dev(structure = None, site = 0, atom = None): 
    """
    Returns an array of standard deviations of atomic positions
    in x,y,z directions for a specified crystal site and atom. 
    Standard deviations are in fractional coordinats. 
    """
    pos = get_select_positions(structure,site,atom)
    return np.std(pos, axis = 0)


def mean_position(structure = None, site = 0, atom = None): 
    """
    Returns the mean unit-cell atomic position of a
    given atom and site. 
    """
    pos = get_select_positions(structure,site,atom)
    return np.mean(pos, axis = 0)



def position_covariance_matrix(structure = None, site = 0, atom = None):
    """
    Returns the covariance matrix of positions for a given 
    atom and site. 
    """ 
    pos = get_select_positions(structure,site,atom)
    return np.cov(pos, rowvar = False)



def get_uncorrelated(structure = None): 
    """
    Returns a new structure with uncorrelated positions
    that have the same mean values and covariance matrices
    as the original structure. 
    """
    structure = copy.deepcopy(structure)
    str_average = structure.get_average_structure()
    N_sites = len(str_average.keys())
    for site in str_average.keys(): 
        for atom in structure.get_atom_Zs():
            pos = get_select_positions(structure, site = site, atom = atom)
            size = pos.shape[0]
            if size != 0:
                mean = np.mean(pos, axis = 0)
                cov =  np.cov(pos, rowvar = False)
                pos =  np.random.multivariate_normal(mean=mean,cov=cov,size=size)
                structure.xyz[site::N_sites] = pos
    return structure 


def get_distances_bw_atoms(structure, neighbors, atom_type1, atom_type2): 
    """ 
    Returns an array of distances between pairs of sites in a given structure
    specified by Neighborlist neighbors and occupied by given atom types.
    """
    shape = (len(structure.atoms.index.levels[0]),
            len(structure.atoms.index.levels[1]),
            len(structure.atoms.index.levels[2]),
            len(structure.atoms.index.levels[3]))
    x = structure.x.reshape(shape)
    y = structure.y.reshape(shape)
    z = structure.z.reshape(shape)
    Z = structure.get_atomic_numbers().reshape(shape)
    neighbors = np.asarray(neighbors)
    distances = []
    for i in range(shape[0]): 
        for j in range(shape[1]): 
            for k in range(shape[2]): 
                for n in neighbors: 
                    u = (i+n[2])%shape[0]
                    v = (j+n[3])%shape[1]
                    w = (k+n[4])%shape[2]
                    a1  = Z[i,j,k,n[0]]
                    a2  = Z[u,v,w,n[1]]
                    if((a1 == atom_type1 and a2 == atom_type2) or (a2 == atom_type1 and a1 == atom_type2)): 
                        dx = n[2] + x[u,v,w, n[1]] - x[i,j,k,n[0]]
                        dy = n[3] + y[u,v,w, n[1]] - y[i,j,k,n[0]]
                        dz = n[4] + z[u,v,w, n[1]] - z[i,j,k,n[0]]
                        distances.append(np.sqrt(dx**2 + dy**2 + dz**2))                 
    return np.array(distances)

def get_distances(structure, neighbors): 
    """ 
    Returns an array of distances between pairs of sites in a given structure 
    specified by Neighborlist neighbors.
    """
    shape = (len(structure.atoms.index.levels[0]),
            len(structure.atoms.index.levels[1]),
            len(structure.atoms.index.levels[2]),
            len(structure.atoms.index.levels[3]))
    x = structure.x.reshape(shape)
    y = structure.y.reshape(shape)
    z = structure.z.reshape(shape)
    neighbors = np.asarray(neighbors)
    distances = []
    for i in range(shape[0]): 
        for j in range(shape[1]): 
            for k in range(shape[2]): 
                for n in neighbors: 
                    u = (i+n[2])%shape[0]
                    v = (j+n[3])%shape[1]
                    w = (k+n[4])%shape[2]
                    dx = n[2] + x[u,v,w, n[1]] - x[i,j,k,n[0]]
                    dy = n[3] + y[u,v,w, n[1]] - y[i,j,k,n[0]]
                    dz = n[4] + z[u,v,w, n[1]] - z[i,j,k,n[0]]
                    distances.append(np.sqrt(dx**2 + dy**2 + dz**2))                 
    return np.array(distances)