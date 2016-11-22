#!/usr/bin/python
"""This script calculates the sum as a function of t in eq. 5 from doi:10.1021/jp501075g in order to get the (ionic) conductivity and the uncorrelated (ionic) conductivity. 
Tau is given in frame distances (for current simulation data this corresponds to 1 simulation time unit).
Before calculation this script checks if there is already a result for the current parameters in the output file."""

#sum_{i,j} z_i z_j <[\vec{R_i}(t+t0) - \vec{R_i}(t0)] \cdot [\vec{R_j}(t+t0) - \vec{R_j}(t0)]>_t0
#i,j: particles
#z_i: valency of particle i (e.g. 1, -1 or 0)

import numpy as np
import sys, os.path, pickle
import multiprocessing
from multiprocessing.managers import BaseManager
from progressBar import printProgress
import signalCatcher
import socket


print "Hostname: {0}".format(socket.gethostname())

limit_subiterations = True
max_subiterations = 50000

n_processes = multiprocessing.cpu_count()

output_path = "output_sigmasum/"
#create folder if necessary
if not os.path.isdir(output_path):
  os.makedirs(output_path)

usage = "Usage: {} filename tau_min tau_max tau_step (max_subiterations)".format(os.path.basename(sys.argv[0]))

#check if file exists  
try:
    filename = sys.argv[1]
    datafile_npart = open(filename, "r")
except:
    print "File '{}' not found.".format(filename)
    print usage
    exit()

try:
    tau_min = int(sys.argv[2])
    tau_max = int(sys.argv[3])
    tau_step = int(sys.argv[4])
except:
    print "No valid values for tau range found. Exit."
    print usage
    exit()

if limit_subiterations:
    try:
        sub_limit = int(sys.argv[5])
        max_subiterations = sub_limit
    except:
        print "No valid value for the subiteration limit found. Use {} instead.\n".format(max_subiterations)


# read in number of particles from header
n_part = 0
n_headerlines = 0
datafile = datafile_npart
for idx, line in enumerate(datafile):
    if idx == 0:
        continue
    
    current_line = line.rstrip()
    split_line = current_line.split()
    
    if split_line[0] == "atom":
        for part_ranges in split_line[1].split(","):
            end_range = int(part_ranges.split(":")[1])
            n_part = end_range if end_range > n_part else n_part
    elif line == "timestep ordered\n" or line =="t\n":
        n_headerlines = idx
        break

datafile_npart.close()
n_part += 1
print "Number of particles =", n_part
print "Number of header lines =", n_headerlines
print ""

ptypes = np.zeros(n_part, dtype=int)
list_ptypes = np.array([])
valencies = np.zeros(n_part, dtype=float)

tau_vals = np.arange(tau_min, tau_max, tau_step)
tau_vals = np.append(tau_vals, tau_max) #tau max should be the last datapoint in tau_vals
len_tau_vals = len(tau_vals)

# read in particle types from header
datafile_ptypes = open(filename, "r")
datafile = datafile_ptypes
for idx, line in enumerate(datafile):
    if idx >= n_headerlines:
        break
    
    current_line = line.rstrip()
    split_line = current_line.split()
    
    if split_line[0] == "atom":
        for part_ranges in split_line[1].split(","):
            start_range = int(part_ranges.split(":")[0])
            end_range = int(part_ranges.split(":")[1])
            ptype = int(split_line[7])
            ptypes[start_range:end_range+1] = ptype
            if not ptype in list_ptypes:
                list_ptypes = np.append(list_ptypes, ptype)
            valency = float(split_line[9])
            valencies[start_range:end_range+1] = valency
            
            print "Found type {0} with valency {1} for particle ids {2}-{3}.".format(ptype, valency, start_range, end_range)

#sort particle types
list_ptypes = np.sort(list_ptypes)

datafile_ptypes.close()

datafile_countframes = open(filename, "r")
def count_frames(datafile):
    counter_frames = 0
    
    for line in datafile:
        if line == "timestep ordered\n" or line =="t\n":
            counter_frames += 1    
    
    return counter_frames
print "\nCount total frames..."
n_timesteps = count_frames(datafile_countframes)
print "Found {} frames.\n".format(n_timesteps)
datafile_countframes.close()
"""
for idx, val in enumerate(ptypes):
    print "ptypes[{0}] = {1}".format(idx, val)
"""


#read in particle coordinates to multidimensional numpy array
# format: coords[timestep][particle_id] = [x, y, z]
datafile_numpy = open(filename, "r")

coords = np.zeros((n_timesteps, n_part, 3))
timestep_coords = np.zeros((n_part, 3))
part_id = 0
timestep = 0

print "Read in coordinates..."

for idx, line in enumerate(datafile_numpy):
    #skip headerlines plus first timestep line
    if idx < n_headerlines+1:
        continue
    
    #search for timestep string
    if line == "timestep ordered\n" or line =="t\n":
        if timestep%1000==0:
            #print "Read in timestep", timestep
            printProgress(timestep, n_timesteps)
        
        coords[timestep] = timestep_coords
        timestep_coords = np.zeros((n_part, 3))
        
        part_id = 0
        timestep += 1
        continue
    
    #with correct format this should not happen
    if line == "\n":
        continue
    
    current_line = line.rstrip()
    split_line = current_line.split()

    timestep_coords[part_id] = np.array(split_line, dtype=float)
    part_id += 1

#coordinates from last frame
coords[timestep] = timestep_coords
len_coords = len(coords)

printProgress(n_timesteps, n_timesteps)

datafile_numpy.close()

print ""
print "\nShape of coordinate matrix:", coords.shape
print "Memory usage of coordinate matrix:", round(coords.nbytes/1024.**2, 2), "mb\n"



#i of interest only for charged particles
i_charged = np.where(valencies)[0]
n_charged_part = len(i_charged)

#number of loop iterations
n_ij_iterations = n_charged_part/2. * (n_charged_part+1)


def ij_t0_average(tau, i, j):
    """Calculates average over different t0 times (which should correspond to the ensemble average) for particles i and j."""
    
    sigma_ij_addend = 0.
    sigma_ij_counter = 0
    
    for t0 in xrange(max_subiterations):
        #check first if coordinates for timestep=t0+tau do exist (avoid out of bounds)
        if t0+tau > len_coords-1:
            break
        
        r_i_last = coords[t0+tau][i] #returns numpy array [x,y,z]
        r_i_first = coords[t0][i]
        
        r_j_last = coords[t0+tau][j]
        r_j_first = coords[t0][j]
        
        #current_t0_addend = valencies[i]*valencies[j] * ((r_i_last - r_i_first)*(r_j_last - r_j_first)).sum()
        current_t0_addend = valencies[i]*valencies[j] * np.dot(r_i_last - r_i_first, r_j_last - r_j_first)
        
        #if i == j:
        #    #correlated part
        #    sigma_ii_counter += 1
        #    sigma_ii_addend = (sigma_ii_counter-1.0)/sigma_ii_counter * sigma_ii_addend + 1.0/sigma_ii_counter * current_t0_addend
        
        #uncorrelated (general) part
        sigma_ij_counter += 1
        sigma_ij_addend = (sigma_ij_counter-1.0)/sigma_ij_counter * sigma_ij_addend + 1.0/sigma_ij_counter * current_t0_addend
    
    return sigma_ij_addend




#file to save data
split_filename = os.path.splitext(filename)
if limit_subiterations:
    output_filename = "{0}{1}_SigmaSumMulticore_limit{2}_tau{3}-{4}-{5}.txt".format(output_path, os.path.basename(split_filename[0]), max_subiterations, tau_min, tau_max, tau_step)
else:
    max_subiterations = len_coords
    output_filename = "{0}{1}_SigmaSumMulticore_tau{2}-{3}-{4}.txt".format(output_path, os.path.basename(split_filename[0]), tau_min, tau_max, tau_step)


#checkpointing
calculated_taus = []
if os.path.isfile(output_filename):
    with open(output_filename, "r") as checkpoint_file:
        for line in checkpoint_file:
            try:
                old_tau = int(line.split()[0])
                calculated_taus.append(old_tau)
            except:
                pass
    output_file = open(output_filename, "a")

else:
    output_file = open(output_filename, "a")
    output_file.write("#tau\tsigma_sum\tsigma_uncorr_sum\n")
    output_file.flush()


class sharedProgress:
    """Shared object for managing the progress bar."""
    def __init__(self, maxIterations):
        self.iterationCounter = 0
        self.maxIterations = maxIterations
        #printProgress(0, maxIterations)
    
    def incrProgress(self, incr):
        self.iterationCounter += incr
        #if self.iterationCounter > self.maxIterations:
            #print "\nError: Iteration {0} > maxIterations {1}".format(self.iterationCounter, self.maxIterations)
        #else:
            #printProgress(self.iterationCounter, self.maxIterations)
        printProgress(self.iterationCounter, self.maxIterations)
    
    def resetProgress(self,printStatus=False):
        self.iterationCounter = 0
        if printStatus:
            printProgress(0, self.maxIterations)

BaseManager.register("sharedProgress", sharedProgress)
manager = BaseManager()
manager.start()
sharedProgressBar = manager.sharedProgress(n_ij_iterations)


def i_worker(tau, i_values):
    """worker function to calculate uncorrelated ii sum for given i values with fixed tau"""
    sigma_ii_sum = 0.
    for i in i_values:
        sigma_ii_sum += ij_t0_average(tau, i, i)
    
    return sigma_ii_sum


def ij_worker(tau, ij_dict):
    """universal worker function to sum up any desired ij combinations.
    ij_dict: {i1:[j11,j12,j13,...], i2:[j21,j22,...], ...}
    Addends are not doubled."""
    
    worker_sigma_ij_sum = 0.
    worker_iteration_counter = 0
    
    for i in ij_dict:
        for j in ij_dict[i]:
            worker_sigma_ij_sum += ij_t0_average(tau, i, j)
            worker_iteration_counter += 1
        
        sharedProgressBar.incrProgress(worker_iteration_counter)
        worker_iteration_counter = 0
    
    return worker_sigma_ij_sum


#calculate sigma sum as a function of tau from multidimensional coordinate numpy array
for current_tau_idx, current_tau in enumerate(tau_vals):
    if current_tau in calculated_taus:
        print "\nSigma sum for tau = {} already calculated. Skip.".format(current_tau)
        continue
    
    print "\nCalculate sigma sum for tau = {}...".format(current_tau)
    
    sigma_sum = 0. #sigma sum variable for the current tau value    
    sigma_uncorr_sum = 0. #sigma uncorrelated sum variable for the current tau value
    
    #distribute ij values equally to worker dictionaries
    sharedProgressBar.resetProgress(printStatus=True)
    
    n_addends_per_worker = int(np.ceil((n_ij_iterations-n_charged_part)/float(n_processes)))
    current_worker = 0
    distribution_counter = 0
    worker_ij_dicts = {} #worker_ij_dicts[worker_no] = {i1:[j11,j12,j13,...], i2:[j21,j22,...], ...}
    current_ij_dict = {} #{i1:[j11,j12,j13,...], i2:[j21,j22,...], ...}
    
    for i in i_charged:
        for j in i_charged:
            if j>=i:
                break
            
            if distribution_counter == n_addends_per_worker:
                #desired number of addends reached, now fill dictionary entry for current worker
                worker_ij_dicts[current_worker] = current_ij_dict
                current_ij_dict = {}
                current_worker += 1
                distribution_counter = 0
            
            if not i in current_ij_dict:
                current_ij_dict[i] = []
            
            current_ij_dict[i].append(j)
            distribution_counter += 1

    #last worker (probably with less addends)
    if n_processes-1 not in worker_ij_dicts:
        #this if condition should always be fulfilled
        worker_ij_dicts[n_processes-1] = current_ij_dict
    
    # DEBUG
    #print ""
    #total_entries = 0
    #for worker_idx in worker_ij_dicts:
        #n_entries = 0
        #for i_val in worker_ij_dicts[worker_idx]:
            #n_entries += len(worker_ij_dicts[worker_idx][i_val])
        #total_entries += n_entries
        #print "worker {0} has {1} i values and {2} total addends to compute.".format(worker_idx, len(worker_ij_dicts[worker_idx]), n_entries)
    #print "total entries:{0}, n_ij_iterations:{1}".format(total_entries, n_ij_iterations)
    #exit()
    
    pool = multiprocessing.Pool(processes=n_processes)
    res = {}
    
    for worker_idx in worker_ij_dicts:
        #distribute work (i.e. i,j dictionary) to the worker function
        res[worker_idx] = pool.apply_async(ij_worker, (current_tau, worker_ij_dicts[worker_idx]))
    
    #close pool and wait for the workers to finish
    pool.close()
    pool.join()
    
    ij_worker_sum = 0.
    for worker_idx in res:
            ij_worker_sum += res[worker_idx].get()
        
    sigma_sum += 2*ij_worker_sum #factor 2 because of entry_ij = entry_ji, avoid double computation
    
    
    #finally sum up all ii (correlated sum) in parallel
    if len(i_charged):
        chunk_entries = int(np.ceil(len(i_charged)/float(n_processes)))
        i_worker_distr = [i_charged[k:k+chunk_entries] for k in xrange(0, len(i_charged), chunk_entries)] #divide i of interest into n_processes parts
        
        pool = multiprocessing.Pool(processes=n_processes)
        res = {}
        
        for worker_idx, i_values in enumerate(i_worker_distr):
            #distribute work (i.e. i values for averaging) to the worker function
            res[worker_idx] = pool.apply_async(i_worker, (current_tau, i_values))
        
        #close pool and wait for the workers to finish
        pool.close()
        pool.join()
        
        sharedProgressBar.incrProgress(n_charged_part)
        
        #sum up all worker results (sum for all j<i for given i)
        for worker_idx in res:
            sigma_uncorr_sum += res[worker_idx].get()
        
        sigma_sum += sigma_uncorr_sum
    
    
    #append result for current tau to file
    output_file.write("{0}\t{1}\t{2}\n".format(current_tau, sigma_sum, sigma_uncorr_sum))
    output_file.flush()
    calculated_taus.append(current_tau)

output_file.close()
print "\nSuccess. Wrote all data to '{}'.".format(output_filename)
