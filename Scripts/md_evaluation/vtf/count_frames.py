#!/usr/bin/python
"""This script counts the number of frames for a given vtf file."""
import numpy as np
import sys, os


filename = sys.argv[1]

#check if file exists  
try:
    datafile = open(filename, "r")
except:
    print "File '{}' not found.".format(filename)
    exit()

def count_frames(datafile):
    counter_frames = 0
    
    for line in datafile:
        #if line.rstrip() == "timestep ordered":
        if line == "timestep ordered\n" or line =="t\n":
            counter_frames += 1    
    
    return counter_frames


print "{0} Frames were found in file '{1}'.".format(count_frames(datafile), os.path.basename(filename))
        
datafile.close()
