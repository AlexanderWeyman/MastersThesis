#!/usr/bin/python
"""This script calculates the averaged velocities in a given direction for the different particle types."""

import numpy as np
import sys, os.path, pickle
from progressBar import printProgress

velocity_component = 0 #0,1,2 corresponds to x,y,z component
velocity_component_dict = {0:'x', 1:'y', 2:'z'}

output_path = "output_meanvelocities/"
#create folder if necessary
if not os.path.isdir(output_path):
  os.makedirs(output_path)

usage = "Usage: {0} filename (velocity_component={1})".format(os.path.basename(sys.argv[0]), velocity_component)

#check if file exists  
try:
    filename = sys.argv[1]
    datafile_npart = open(filename, "r")
except:
    print "File '{}' not found.".format(filename)
    print usage
    exit()

try:
    v_comp = int(sys.argv[2])
    if not v_comp in [0,1,2]:
        print "No valid velocity component found. Use {} instead\n".format(velocity_component)
    else:
        velocity_component = v_comp
except:
    print "No valid value for the velocity component found. Use {} instead.\n".format(velocity_component)


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

mean_velocities = {} #mean_velocities[part_id] = mean value

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
                mean_velocities[ptype] = 0.
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
        if timestep%100==0:
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

datafile_numpy.close()

print ""
print "\nShape of coordinate matrix:", coords.shape
print "Memory usage of coordinate matrix:", round(coords.nbytes/1024.**2, 2), "mb\n"



#file to save data
split_filename = os.path.splitext(filename)
output_filename = "{0}{1}_MeanVelocities.txt".format(output_path, os.path.basename(split_filename[0]))
output_file = open(output_filename, "a")

output_file.write("#")
for t in list_ptypes:
    output_file.write("v_{0}_mean_{1}\t".format(velocity_component_dict[velocity_component], int(t)))
output_file.write("\n")
output_file.flush()



#calculate mean velocities
for idx, particle_type in enumerate(list_ptypes):
    #indices of interest
    indices_current_ptype = np.where(ptypes == particle_type))[0]
    #average over coordinates of interest
    mean_velocities[particle_type] = np.mean(np.take(coords[:,:,velocity_component], indices_current_ptype, axis=1), dtype=np.float64))

    output_file.write("{}\t".format(mean_velocities))
    output_file.flush()
    printProgress(idx, len(list_ptypes))


output_file.close()
print "Success. Wrote all data to '{}'.".format(output_filename)
