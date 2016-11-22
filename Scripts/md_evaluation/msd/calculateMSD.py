#!/usr/bin/python
"""This script calculates the mean square displacement from vtf files. Tau is given in frame distances (for current simulation data this corresponds to 1 simulation time unit)"""
import numpy as np
import sys, os.path, pickle

limit_subiterations = True
max_subiterations = 20000

output_path = "output_msd/"
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
msd = {}

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
                # also initialize msd dictionary that returns np array
                msd[ptype] = np.zeros(len_tau_vals)
            
            print "Found type {0} for particle ids {1}-{2}.".format(ptype, start_range, end_range)

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



#file to save msd data
split_filename = os.path.splitext(filename)
if limit_subiterations:
    output_filename = "{0}{1}_MSD_limit{2}_tau{3}-{4}-{5}.txt".format(output_path, os.path.basename(split_filename[0]), max_subiterations, tau_min, tau_max, tau_step)
else:
    output_filename = "{0}{1}_MSD_tau{2}-{3}-{4}.txt".format(output_path, os.path.basename(split_filename[0]), tau_min, tau_max, tau_step)
output_file = open(output_filename, "a")

output_file.write("#tau")
for t in list_ptypes:
    output_file.write("\tMSD_{}".format(int(t)))
output_file.write("\n")


#calculate msd from numpy array
for current_tau_idx, current_tau in enumerate(tau_vals):
    print "Calculate MSD for tau = {}...".format(current_tau)
    msd_tau = np.zeros(n_part) #contains the MSD(particle_id) for the current tau value
    msd_tau_counter = 0 # number of configurations averaged in msd_tau
    
    #average msd for current tau over different t0
    for t0 in range(max_subiterations):
        #check first if coordinates for timestep=t0+current_tau do exist (avoid out of bounds)
        if t0+current_tau > len_coords-1:
            break
        
        coords_start = coords[t0]
        coords_stop = coords[t0+current_tau]
        msd_tau_frame = (coords_stop[:,0] - coords_start[:,0])**2 + (coords_stop[:,1] - coords_start[:,1])**2 + (coords_stop[:,2] - coords_start[:,2])**2
        msd_tau_counter += 1
        # average msd_tau_frame into msd_tau
        msd_tau = (msd_tau_counter-1.0)/msd_tau_counter * msd_tau + 1.0/msd_tau_counter * msd_tau_frame
    
    tmp_msd_type_counter = {}
    for part_id in xrange(n_part):
        #initialize tmp counter
        if ptypes[part_id] not in tmp_msd_type_counter:
            tmp_msd_type_counter[ptypes[part_id]] = 0
        
        tmp_msd_type_counter[ptypes[part_id]] += 1
        msd[ptypes[part_id]][current_tau_idx] = (tmp_msd_type_counter[ptypes[part_id]] - 1.0)/tmp_msd_type_counter[ptypes[part_id]] * msd[ptypes[part_id]][current_tau_idx] + 1.0/tmp_msd_type_counter[ptypes[part_id]] * msd_tau[part_id]
    
    #write msd result for current tau to output file
    output_file.write(str(current_tau))
    for t in list_ptypes:
        output_file.write("\t{}".format(msd[t][current_tau_idx]))
    output_file.write("\n")
    output_file.flush()

output_file.close()

print "Success. Wrote all data to '{}'.".format(output_filename)
