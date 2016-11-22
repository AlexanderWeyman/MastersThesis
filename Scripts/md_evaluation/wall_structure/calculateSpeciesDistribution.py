#!/usr/bin/python
"""This script calculates the distribution of particle types along the z axis from vtf files."""

import numpy as np
import sys, os.path
from progressBar import printProgress
from readVtfToNumpy import readVtf

usage = "Usage: {0} filename_coordinates".format(os.path.basename(sys.argv[0]))

n_bins = 10
eta = 0.3

#check if file exists  
try:
    filename = sys.argv[1]
    datafile_npart = open(filename, "r")
except:
    print "File '{}' not found.".format(filename)
    print usage
    exit()

output_path = "output_speciesdistribution/"
#create folder if necessary
if not os.path.isdir(output_path):
  os.makedirs(output_path)


#coords[timestep][particle_id] = [x, y, z]
coords, list_ptypes, ptypes, valencies = readVtf(filename)

n_frames = coords.shape[0]
n_part = coords.shape[1]
box_l = (n_part*np.pi/(6.*eta))**(1./3.)

#indices/ids for different particle types
indices_ptype = {}
z_distributions = {}
for t in list_ptypes:
    indices_ptype[t] = np.where(ptypes == t)[0]
    z_distributions[t] = np.zeros(n_bins)

print "Calculate distribution along z axis..."

printProgress(0, len(list_ptypes))


#for current_frame in xrange(len(coords)):
for idx, t in enumerate(list_ptypes):
    #z_distributions[t], bin_edges = np.histogram(np.take(coords[:,:,2], indices_ptype[t], axis=1), bins=n_bins, range=(0.,box_l), density=True)
    z_distributions[t], bin_edges = np.histogram(np.take(coords[:,:,2], indices_ptype[t], axis=1)%box_l, bins=n_bins, range=(0.,box_l), density=True)
    
    printProgress(idx+1, len(list_ptypes))
        
#file to save data
split_filename = os.path.splitext(filename)
output_filename = "{0}{1}_Species_Distribution_{2}bins.txt".format(output_path, os.path.basename(split_filename[0]), n_bins)
output_file = open(output_filename, "w")

output_file.write("#bin_edges")
for t in list_ptypes:
    output_file.write("\tdistr_type{0}".format(int(t)))
output_file.write("\n")

for i in xrange(len(bin_edges)-1):
    output_file.write("{}".format(bin_edges[i]))
    for t in list_ptypes:
        output_file.write("\t{}".format(z_distributions[t][i]))
    output_file.write("\n")

output_file.close()
print "Success. Wrote all data to '{}'.".format(output_filename)

