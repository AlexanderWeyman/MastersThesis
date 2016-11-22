#!/usr/bin/python
"""This script finds corrupt frames in a given vtf file."""
import numpy as np
import sys, os
from progressBar import printProgress

usage = "Usage: {} filename".format(os.path.basename(sys.argv[0]))

#check if file exists  
try:
    filename = sys.argv[1]
    datafile_npart = open(filename, "r")
except:
    print "File '{}' not found.".format(filename)
    print usage
    exit()


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


### count total number of frames ###
def count_frames(datafile):
    counter_frames = 0
    
    for line in datafile:
        #print "-begin-{0}-end-".format(line)
        #if line.rstrip() == "timestep ordered":
        if line == "timestep ordered\n" or line =="t\n":
            counter_frames += 1    
    
    return counter_frames

print "Count frames..."
with open(filename, "r") as datafile:
    n_frames = count_frames(datafile)
print "Found {} frames.\n".format(n_frames)


### check if coordinates are correctly formated (number of coordinates per frame and datatype) ###
print "Check validity of coordinate file..."

validity_buffer = ""

datafile_countframes = open(filename, "r")
counter_frames = 0
counter_particles = 0

#particle_id_deviation = 0 #deviations will be evaluated only for this given particle
particle_id_deviation = n_part - 1 #the last particle is a counterion
coord_part_last_frame = np.array([0.,0.,0.])
dx = np.zeros(n_frames)

for idx, line in enumerate(datafile_countframes):
    if idx < n_headerlines:
        continue
    
    #starts with first "timestep ordered" or "t" line
    
    if line == "timestep ordered\n" or line =="t\n":
        if not counter_frames % 100:
            printProgress(counter_frames, n_frames)
        if counter_particles != n_part and counter_frames:
            validity_buffer += "Frame {0} corrupt: found only {1}/{2} particles.\n".format(counter_frames, counter_particles, n_part)
        
        counter_particles = 0
        counter_frames += 1
    
    elif line == "\n":
        continue
    
    else:
        #should find only coordinates
        current_line = line.rstrip()
        split_line = current_line.split()
        try:
            current_coords = np.array(split_line, dtype=float)
            
            if counter_particles == particle_id_deviation:
                if counter_frames == 1:
                    coord_last_frame = np.copy(current_coords)
                else:
                    dx[counter_frames-1] = np.linalg.norm(current_coords-coord_last_frame)
                    coord_last_frame = np.copy(current_coords)
            
            if len(current_coords) != 3:
                raise Exception("Number of coordinates differs from 3.")
            counter_particles += 1
        except Exception as e:
            validity_buffer += "Frame {0} corrupt: coordinates of particle {1} faulty.\nLine: {2}\nException: {3}\n\n".format(counter_frames, counter_particles+1, current_line, e)

#check last frame
if counter_particles != n_part:
    validity_buffer += "Frame {0} corrupt: found only {1}/{2} particles.\n".format(counter_frames, counter_particles, n_part)

printProgress(counter_frames, n_frames)
datafile_countframes.close()

if validity_buffer == "":
    print "Coordinate file seems to be valid."
else:
    print validity_buffer


### check for gaps ###
print "\nChecking for gaps by observing particle {0}...".format(particle_id_deviation)
dx_mean = dx.mean()
dx_std = dx.std()

#find high deviations from mean value
std_criterion = 5.
dx_deviations = np.abs(dx-dx_mean)
idx_deviation = np.where(dx_deviations>std_criterion*dx_std)[0]

max_deviation = np.amax(dx_deviations)
idx_max_deviation = np.where(dx_deviations == max_deviation)[0]
print "Maximum deviation found: {0} sigma (between frames {1} and {2})\n".format(max_deviation/dx_std, idx_max_deviation, idx_max_deviation+1)

for idx in idx_deviation:
    print "Significant coordinate deviations: {0} sigma (between frames {1} and {2})".format(abs(dx[idx]-dx_mean)/dx_std, idx, idx+1)

if not len(idx_deviation):
    print "No gaps between frames (with significant deviations > {0} sigma) were found.".format(std_criterion)


