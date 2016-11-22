#!/usr/bin/python
"""This script calculates the kinetic energy distribution of particle types along the z axis from vtf files."""

import numpy as np
import sys, os.path
from progressBar import printProgress
from readVtfToNumpy import readVtf

usage = "Usage: {0} filename_coordinates filename_velocities".format(os.path.basename(sys.argv[0]))

n_bins = 100
eta = 0.3
max_frames = -1 #for limitation define the maximum frame number (positive), -1 corresponds to no limitation

#check if files exist
try:
    filename_coordinates = sys.argv[1]
    datafile_coordinates = open(filename_coordinates, "r")
    datafile_coordinates.close()
except:
    print "File '{}' not found.".format(filename_coordinates)
    print usage
    exit()

try:
    filename_velocities = sys.argv[2]
    datafile_velocities = open(filename_velocities, "r")
    datafile_velocities.close()
except:
    print "File '{}' not found.".format(filename_velocities)
    print usage
    exit()

output_path = "output_velocitydistribution/"
#create folder if necessary
if not os.path.isdir(output_path):
  os.makedirs(output_path)


#coords[timestep][particle_id] = [x, y, z]
coords, list_ptypes, ptypes, valencies = readVtf(filename_coordinates)
vels, vel_list_ptypes, vel_ptypes, vel_valcencies = readVtf(filename_velocities) #dummy variables start with vel_, they should be equal to the ones from the coordinate file

if coords.shape != vels.shape:
    print "The given velocity file does not fit the coordinate file. Exit."
    exit()

n_frames = coords.shape[0]
n_part = coords.shape[1]
box_l = (n_part*np.pi/(6.*eta))**(1./3.)

#indices/ids for different particle types
indices_ptype = {}
kin_energy_distributions = {}
avg_counters = {} #avg_counter[type] = numpy array with length n_bins
for t in list_ptypes:
    indices_ptype[t] = np.where(ptypes == t)[0]
    kin_energy_distributions[t] = np.zeros(n_bins)
    avg_counters[t] = np.zeros(n_bins)


print "Calculate velocity distribution along z axis..."

printProgress(0, n_frames)

for current_frame in xrange(n_frames):
    if current_frame == max_frames:
        break
    
    for t in list_ptypes:
        #z coordinates of interest: np.take(coords[current_frame,:,2], indices_ptype[t])
        #velocities of interest (field is applied in x direction): np.take(vels[current_frame,:,0], indices_ptype[t])
        current_coords = np.take(coords[current_frame,:,2], indices_ptype[t])
        current_vels = np.take(vels[current_frame,:,:], indices_ptype[t], axis=0)
        
        for i in xrange(current_vels.shape[0]):
            #bin index is calculates from z coordinate
            #bin_idx = int(np.floor(current_coords[i]*n_bins/box_l))
            bin_idx = int(np.floor(current_coords[i]%box_l*n_bins/box_l))
            bin_idx = 0 if bin_idx == n_bins else bin_idx
            kin_energy_distributions[t][bin_idx] += np.dot(current_vels[i], current_vels[i])
            avg_counters[t][bin_idx] += 1.0
    
    printProgress(current_frame+1, n_frames)


#averaging velocities
print "Calculate mean velocities in bins"
for t in list_ptypes:
    for bin_idx, counter in enumerate(avg_counters[t]):
        if counter > 0.0:
            kin_energy_distributions[t][bin_idx] /= counter


#file to save data
split_filename = os.path.splitext(filename_velocities)
output_filename = "{0}{1}_Kinetic_Energy_Distribution_{2}bins.txt".format(output_path, os.path.basename(split_filename[0]), n_bins)
output_file = open(output_filename, "w")

bin_edges = np.linspace(0.,box_l,n_bins+1)

output_file.write("#bin_edges")
for t in list_ptypes:
    output_file.write("\tekin_x_type{0}".format(int(t)))
output_file.write("\n")

for i in xrange(n_bins):
    output_file.write("{}".format(bin_edges[i]))
    for t in list_ptypes:
        output_file.write("\t{}".format(kin_energy_distributions[t][i]))
    output_file.write("\n")

output_file.close()
print "Success. Wrote all data to '{}'.".format(output_filename)

