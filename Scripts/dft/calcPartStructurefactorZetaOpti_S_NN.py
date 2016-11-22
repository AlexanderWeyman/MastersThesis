#!/usr/bin/python
"""This script computes the partial structure factors from dft approach after doi:10.1063/1.1853371 with RPM model"""

import numpy as np
import scipy
from scipy import math as smath
from scipy.optimize import curve_fit
import sys, os.path

import matplotlib.pyplot as plt

from progressBar import printProgress

#usage = "Usage: {} N delta_s eta l_Bjerrum q_min q_max q_step [useSympy=0] [zeta=1.0]".format(os.path.basename(sys.argv[0]))
usage = "Usage: {} sim_filename fit_q_min fit_q_max N delta_s eta l_Bjerrum".format(os.path.basename(sys.argv[0]))

try:
    sim_filename = sys.argv[1]
    fit_q_min = float(sys.argv[2])
    fit_q_max = float(sys.argv[3])
    N = int(sys.argv[4])
    delta_s = int(sys.argv[5])
    eta = float(sys.argv[6])
    l_b = float(sys.argv[7])
except:
    print "Given parameters are not valid. Exit."
    print usage
    exit()


zeta = 1.0
    

#q_min = 0.1
#q_max = 20
#q_step = 0.1


output_path = "output_dft/"
#create folder if necessary
if not os.path.isdir(output_path):
  os.makedirs(output_path)



#eta = 0.3 #volume fraction
#N = 30 #degree of polymerization
gamma0 = -(1.+2.*eta)**2/(1.-eta)**4
gamma1 = 6.*eta*(1.+0.5*eta)**2/(1.-eta)**4
gamma3 = -eta*(1.+2.*eta)**2/(2*(1.-eta)**4)

sigma = 1.0
#l_b = 7.0 #Bjerrum length

kappa_sigma = np.sqrt(48./(delta_s+2)*eta*l_b/sigma) #kappa * sigma, eq. 47
B = (kappa_sigma + 1. - np.sqrt(1.+2*kappa_sigma))/kappa_sigma if l_b else 0.0#see eq. 5

MPC = N #monomers per chain
#delta_s = 2 #aka sclen, number of side chain beads attached to the backbone
delta_b = 1 #aka side_bead_distance, if 1, every backbone monomer carries a sidechain

rho_tot = 6.*eta/(np.pi*sigma**3)
rho_p = 6.*eta/((delta_s+2)*N*np.pi*sigma**3)

#charge fraction of the polymer chain
f = 1./(delta_s + 1.)
#number of chain segments
N_seg = (delta_s + 1)*N

lambda_homogeneous = rho_p
S3_homogeneous = eta / (1. + f)

d_lny_d_S3 = 1./(S3_homogeneous-2.) - 3./(S3_homogeneous-1.)
d2_lny_d_S32 = 3./(S3_homogeneous-1.)**2 - 1./(S3_homogeneous-2.)**2


type_neutral_backbone = 1
type_charged_backbone = 2
type_neutral_sidebead = 3
type_charged_sidebead = 4
type_counterion = 5

def npFactorial(input_arr):
    ret = np.array([])
    
    for i in input_arr:
        ret = np.append(ret, smath.factorial(i))
    
    return ret

def plot_data(q, function, *args):
    #q = np.linspace(0.,30.,100)
    plot_vals = np.array([])
    for i in q:
        plot_vals = np.append(plot_vals, function(i,*args))
        #print i, plot_vals[-1]
    plt.figure()
    plt.plot(q, plot_vals)
    plt.show()

def primitive_Function(x, y, n=0):
    """Return values of the primitive function (dt. 'Stammfunktion'), see eq. 39"""
    #x, y must be scalar
    k = np.arange(n+1) #index variable
    
    ret = (-1.) * (smath.factorial(n)/npFactorial(n-k) * x**(n-k)/y**(k+2.) * np.cos(x*y+k*np.pi/2.)).sum() #sum over k from 0 to n
    return ret


def integral_function(y, n=0):
    """Calculation of I_n(y), eq. 39-43"""
    
    if y == 0.:
        return 1./(n+2.)
    
    if n == 1:
        return 1./y**3 * (np.sin(y) - y*np.cos(y))
    elif n == 2:
        return 1./y**4 * (-2+2*y*np.sin(y)+(2-y**2)*np.cos(y))
    elif n == 4:
        return 1./y**6 * (24 + (-24*y+4*y**3)*np.sin(y)+(-24+12*y**2-y**4)*np.cos(y))
    else:
        #y must be scalar
        return (primitive_Function(1., y, n) - primitive_Function(0., y, n))


def c_hat_py(q):
    """Calculation of the fourier transformed Percus-Yerick correlation function, eq. 38"""
    #q must be scalar
    
    ret = 4.*np.pi*sigma**3 * (gamma0*integral_function(q*sigma,n=1) + gamma1*integral_function(q*sigma, n=2) + gamma3*integral_function(q*sigma, n=4))
    return ret
#plot_data(np.linspace(0.,30.,100), c_hat_py)


def delta_c_hat(q):
    """Calculation of the coulomb part of the fourier transformed correlation function, eq. 44"""
    #q must be scalar
    
    ret = 4.*np.pi*l_b*sigma**2 * (2.*B*integral_function(q*sigma,n=1) - B**2*integral_function(q*sigma,n=2) + np.cos(q*sigma)/(q*sigma)**2.)
    return ret
#plot_data(np.linspace(1e-5,1,100), delta_c_hat)



### Begin setup PIL geometry ##################################################

def add_neighbor(part_id_a, part_id_b):
    """Register neighboring particles in neighbors dictionary, avoid double occurence of particle ids"""
    global neighbors
    
    if not part_id_a in neighbors[part_id_b]:
        neighbors[part_id_b] = np.append(neighbors[part_id_b], part_id_a)
    
    if not part_id_b in neighbors[part_id_a]:
        neighbors[part_id_a] = np.append(neighbors[part_id_a], part_id_b)


def add_particle(part_type, part_id=None):
    """Add a particle by defining the particle type in parttype array"""
    global parttype, n_part
    
    if part_id == None:
        part_id = n_part
    
    parttype[part_id] = part_type
    n_part+=1


def change_ptype(part_type, part_id):
    """Changes the particle type of the particle with the given particle id"""
    global parttype
    
    parttype[part_id] = part_type


#create lists/dictionaries that contain structural information of the polymer chains: valence[partID] (np array), neighbors[partID] (dict that returns np array - different dimensionality prevents using 2d np arrays), parttype[partID] (np array)
start_particle_id = 0 #defines the id of the first monomer of the backbone that is connected with a side bead
N_side_chains_per_polymer = int(smath.ceil(float(MPC - start_particle_id)/delta_b))
val_cM = 1.0 #valence of charged monomer beads
val_CI = -1.0 #valence of counterions

#init with zeros respectively empty np arrays
valence = np.zeros((delta_s+1)*N+1)
parttype = np.zeros((delta_s+1)*N+1,dtype=int)
neighbors = {}
for i in xrange((delta_s+1)*N+1):
    neighbors[i] = np.array([], dtype=int)

n_part = 0 #contains the number of particles 'added' to the system

first_monomer_id = n_part

# 'add' backbone of current polymer
# first particle
add_particle(type_neutral_backbone, first_monomer_id)

# second to last backbone particle
for backbone_id in xrange(first_monomer_id+1,first_monomer_id+MPC):
    add_particle(type_neutral_backbone, backbone_id)
    add_neighbor(backbone_id-1, backbone_id)

#counts the number of side chains (-1) beeing attachted to the backbone
side_chain_counter = 0

#current_particle_id runs over the backbone monomer id's where a side chain should be attached
for current_particle_id in xrange(first_monomer_id+start_particle_id, first_monomer_id+MPC, delta_b):
    current_particle_id_charged = current_particle_id
    
    #sub-loop for attaching side monomers
    for side_monomer_counter in xrange(0, delta_s):
        #print "attaching side monomer no",side_monomer_counter
        new_side_bead_id = MPC+side_monomer_counter*N_side_chains_per_polymer+side_chain_counter+first_monomer_id
        
        add_particle(type_neutral_sidebead, new_side_bead_id)
        add_neighbor(new_side_bead_id, current_particle_id_charged)
        
        current_particle_id_charged = new_side_bead_id
    
    #change particle type and charge last (side) bead
    if delta_s:
        change_ptype(type_charged_sidebead, current_particle_id_charged)
        valence[current_particle_id_charged] = val_cM
    else:
        change_ptype(type_charged_backbone, current_particle_id_charged)
        valence[current_particle_id_charged] = val_cM
    
    side_chain_counter += 1


#add counterion
valence[n_part] = val_CI
add_particle(type_counterion)

"""
for idx, val in enumerate(parttype):
    print "particle {0} of type {1} has following neighbors: {2}".format(idx, val, neighbors[idx])
"""
### End setup PIL geometry ####################################################


def is_neighboring(part_id_a, part_id_b):
    """Implementation of delta_alpha<->alpha' """
    if part_id_a in neighbors[part_id_b]:
        return 1.0
    else:
        return 0.0


def n_neighbors(part_id):
    """Calculates the number of connected particles, see eq. 16"""
    return len(neighbors[part_id])


def c_hat_assoc(q, part_id_a, part_id_b):
    """Calculation of the fourier transformed association correlation function, eq. 46"""
    #q must be scalar
    if parttype[part_id_a] == type_counterion or parttype[part_id_b] == type_counterion:
        ret = 0.0
    
    else:
        ret = np.pi/4. * sigma**3 * (n_neighbors(part_id_a)+n_neighbors(part_id_b)) * d_lny_d_S3 * integral_function(q*sigma,n=1) \
            + is_neighboring(part_id_a, part_id_b) * 1./lambda_homogeneous * np.sin(q*sigma)/(q*sigma) \
            + np.pi**2/4. * sigma**6 * lambda_homogeneous * (N_seg-1.) * d2_lny_d_S32 * integral_function(q*sigma,n=1)**2 \
            - 1./(2.*lambda_homogeneous) * (part_id_a == part_id_b) * n_neighbors(part_id_a) * (np.sin(q*sigma)/(q*sigma))**2
        
        #ret = np.pi/2 * sigma**3 * (n_neighbors(part_id_a)+n_neighbors(part_id_b))*(-12.*eta+40)/(9.*eta**2-36.*eta+32) * integral_function(q*sigma,n=1) \
            #+ is_neighboring(part_id_a, part_id_b) * (2*np.pi*N*sigma**3)/(3.*eta) * np.sin(q*sigma)/(q*sigma) \
            #+ 12*np.pi*sigma**3*(3-1./N) * (9.*eta**3-60.*eta**2+88.*eta)/(9.*eta**2-36.*eta+32.)**2 * integral_function(q*sigma,n=1)**2 \
            #- (part_id_a == part_id_b) * (np.pi*N*sigma**3)/(3.*eta) * n_neighbors(part_id_a) * (np.sin(q*sigma)/(q*sigma))**2
    
    return ret
#plot_data(np.linspace(1e-5,1,100), c_hat_assoc, 102, 100)


def c_hat(q, part_id_a, part_id_b):
    """Summation of all fourier transformed correlation functions"""
    ret = c_hat_py(q) - valence[part_id_a] * valence[part_id_b] * delta_c_hat(q) + zeta * c_hat_assoc(q, part_id_a, part_id_b)
    return ret


def mole_fraction(part_id):
    """Returns the mole fraction (dt. Stoffmengenanteil) for a given particle id, see eq. 26"""
    if parttype[part_id] == type_counterion:
        return 1./(delta_s+2)
    else:
        return 1./((delta_s+2)*N)


def c_matrix_entry(q, part_id_a, part_id_b):
    """Returns one entry of the curly C matrix, see eq. 28"""
    ret = rho_tot * np.sqrt(mole_fraction(part_id_a)*mole_fraction(part_id_b)) * c_hat(q, part_id_a, part_id_b)
    return ret


def calc_c_matrix(q):
    """Returns the entire curly C matrix for one q value, see eq. 28"""
    c_matrix = np.zeros(((delta_s+1)*N+1, (delta_s+1)*N+1))
    for part_id_a in xrange((delta_s+1)*N+1):
        for part_id_b in xrange((delta_s+1)*N+1):
            c_matrix[part_id_a][part_id_b] = c_matrix_entry(q, part_id_a, part_id_b)
    
    return c_matrix


def calc_s_matrix(q):
    """Calculates the curly S matrix, see eq. 29"""
    matr = np.identity((delta_s+1)*N+1) - calc_c_matrix(q)
    return np.linalg.inv(matr)




#get type ids
backbone_ids = np.where(np.in1d(parttype, [type_neutral_backbone, type_charged_backbone]))[0]
#terminal_ids = np.where(parttype == type_charged_sidebead)[0]
terminal_ids = np.where(np.in1d(parttype, [type_charged_sidebead, type_charged_backbone]))[0]
counterion_ids = np.where(parttype == type_counterion)[0]
particle_ids = range((delta_s+1)*N+1)


def get_S(q, zetaOpt):
    global zeta
    zeta = zetaOpt
    
    s_matrix = calc_s_matrix(q)
    
    S_BB = S_TT = S_BT = S_TC = S_CC = S_NN = 0.
    
    """
    for part_id_a in backbone_ids:
        for part_id_b in backbone_ids:
            S_BB += s_matrix[part_id_a][part_id_b]
        for part_id_b in terminal_ids:
            S_BT += s_matrix[part_id_a][part_id_b]
    
    
    for part_id_a in terminal_ids:
        for part_id_b in terminal_ids:
           S_TT += s_matrix[part_id_a][part_id_b] 
        for part_id_b in counterion_ids:
            S_TC += s_matrix[part_id_a][part_id_b]
    
    for part_id_a in counterion_ids:
        for part_id_b in counterion_ids:
            S_CC += s_matrix[part_id_a][part_id_b]
    """
    for part_id_a in particle_ids:
        for part_id_b in particle_ids:
            S_NN += np.sqrt(mole_fraction(part_id_a)*mole_fraction(part_id_b)) * s_matrix[part_id_a][part_id_b]
    """
    #normalization
    if len(backbone_ids):
        S_BB *= mole_fraction(backbone_ids[0])
    if len(backbone_ids) and len(terminal_ids):
        S_BT *= np.sqrt(mole_fraction(backbone_ids[0])*mole_fraction(terminal_ids[0]))
    if len(terminal_ids):
        S_TT *= mole_fraction(terminal_ids[0])
    if len(terminal_ids) and len(counterion_ids):
        S_TC *= np.sqrt(mole_fraction(terminal_ids[0])*mole_fraction(counterion_ids[0]))
    if len(counterion_ids):
        S_CC *= mole_fraction(counterion_ids[0])
    
    S_QQ = S_TT + S_CC - 2*S_TC
    """
    
    return S_NN


#loop over q values, calculate S matrix and therefrom derive partial structure factors of interest
#for current_q_idx, current_q in enumerate(q_vals):
def fit_S(q, zetaOpt):  
    
    if len(q) > 1:
        print "zeta={0}, {1} q values".format(zetaOpt, len(q))
        ret_arr = np.zeros(len(q))
        
        printProgress(0, len(q))
        for idx_q, curr_q in enumerate(q):
            ret_arr[idx_q] = get_S(curr_q, zetaOpt)
            printProgress(idx_q+1, len(q))
        
        return ret_arr
    
    else:
        return get_S(q, zetaOpt)
    


#print fit_S(0.5, 1.0)
#exit()
#read simulation results
sim_data = np.genfromtxt(sim_filename, skip_header=1)

#get xrange
def get_x_idx(x):
    for idx_x, x_val in enumerate(sim_data[:,0]):
        if x_val > x:
            return idx_x-1
    return idx_x

min_idx = get_x_idx(fit_q_min)
max_idx = get_x_idx(fit_q_max)
print "q[{0}]={1}".format(min_idx, sim_data[:,0][min_idx])
print "q[{0}]={1}".format(max_idx, sim_data[:,0][max_idx])

print curve_fit(fit_S, sim_data[:,0][:max_idx], sim_data[:,7][:max_idx], p0=0.85)



