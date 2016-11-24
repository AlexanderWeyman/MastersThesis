#!/usr/bin/python
"""Simplified version of calculateSigmaSumMulticore.py for the uncorrelated conductivity"""


import numpy as np
import sys, os.path, pickle

limit_subiterations = True
max_subiterations = 50000

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

for idx, line in enumerate(datafile_numpy):
    #skip headerlines plus first timestep line
    if idx < n_headerlines+1:
        continue
    
    #search for timestep string
    if line == "timestep ordered\n" or line =="t\n":
        if timestep%1000==0:
            print "Read in timestep", timestep
        
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

datafile_numpy.close()

print ""
print "\nShape of coordinate matrix:", coords.shape
print "Memory usage of coordinate matrix:", round(coords.nbytes/1024.**2, 2), "mb\n"



#file to save data
split_filename = os.path.splitext(filename)
if limit_subiterations:
    output_filename = "{0}{1}_SigmaSumII_limit{2}_tau{3}-{4}-{5}.txt".format(output_path, os.path.basename(split_filename[0]), max_subiterations, tau_min, tau_max, tau_step)
else:
    max_subiterations = len_coords
    output_filename = "{0}{1}_SigmaSumII_tau{2}-{3}-{4}.txt".format(output_path, os.path.basename(split_filename[0]), tau_min, tau_max, tau_step)
output_file = open(output_filename, "a")

output_file.write("#tau\tsigma_ii_sum\n")

#calculate sigma sum as a function of tau from multidimensional coordinate numpy array
for current_tau_idx, current_tau in enumerate(tau_vals):
    print "Calculate sigma sum for tau = {}...".format(current_tau)
      
    sigma_uncorr_sum = 0. #sigma uncorrelated sum variable for the current tau value
    
    #i,i sum
    for i in xrange(n_part):
        #continue if at least one of the two particles is neutral
        if not valencies[i]:
            continue
        
        sigma_ii_addend = 0.
        sigma_ii_counter = 0
        
        #calculate average over different t0 times (which should correspond to the ensemble average)
        for t0 in xrange(max_subiterations):
            #check first if coordinates for timestep=t0+current_tau do exist (avoid out of bounds)
            if t0+current_tau > len_coords-1:
                break
            
            r_i_last = coords[t0+current_tau][i] #returns numpy array [x,y,z]
            r_i_first = coords[t0][i]
            
            current_t0_addend = np.dot(r_i_last - r_i_first, r_i_last - r_i_first)
            

            #correlated part
            sigma_ii_counter += 1
            sigma_ii_addend = (sigma_ii_counter-1.0)/sigma_ii_counter * sigma_ii_addend + 1.0/sigma_ii_counter * current_t0_addend
            
            
        
        sigma_uncorr_sum += sigma_ii_addend
    
    #append result for current tau to file
    output_file.write("{0}\t{1}\n".format(current_tau, sigma_uncorr_sum))
    output_file.flush()

output_file.close()
print "Success. Wrote all data to '{}'.".format(output_filename)
