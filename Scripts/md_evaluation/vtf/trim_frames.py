#!/usr/bin/python
"""This script trims a vtf file for all frames after a given value. 
The index of the first frame is 1, first_frame and last_frame are included in the exported vtf file."""
import numpy as np
import sys, os

usage = "Usage: {} filename first_frame last_frame".format(os.path.basename(sys.argv[0]))

#filename = "R_GH_tester_N30_rho0.005_lB7.0_dist1_sclen3_polymer_sidechains_equi.vtf"
try:
    filename = sys.argv[1]
except:
    print usage
    exit()

if sys.argv[2] == "" or sys.argv[2] == "first":
    first_frame = 0
else:
    first_frame = int(sys.argv[2])

if sys.argv[3] == "" or sys.argv[3] == "last":
    last_frame = "last"
else:
    last_frame = int(sys.argv[3])

#check if file exists  
try:
    datafile = open(filename, "r")
except:
    print "File '{}' not found.".format(filename)
    exit()

def trim_vtf_file(datafile, last_frame, dst_filename):
    dst_file = open(dst_filename, "w")
    verbose_mode = False
    
    counter_frames = 0
    
    for line in datafile:        
        if line == "timestep ordered\n" or line =="t\n":
            counter_frames += 1
        #if counter_frames < last_frame and (counter_frames > first_frame or counter_frames == 0):
        #note that int < string is always true
        if counter_frames <= last_frame:
            if counter_frames >= first_frame or counter_frames == 0:
                dst_file.write(line)
                if verbose_mode and counter_frames > 0:
                    print counter_frames
        else:
            break
        
        
    
    dst_file.close()


pathname = os.path.dirname(filename)
if pathname != "":
    pathname += "/"
dest_filename = pathname + "frames_" + str(first_frame) + "-" + str(last_frame) + "_" + os.path.basename(filename)

trim_vtf_file(datafile, last_frame, dest_filename)

datafile.close()

print "Success. Exported frames {0}-{1} to file '{2}'".format(first_frame, last_frame, dest_filename)

#print "# of frames: {}".format(frame+1)
