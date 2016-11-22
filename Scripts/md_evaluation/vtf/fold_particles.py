#!/usr/bin/python
"""This script folds the particle coordinates from vtf files back to the (original) simulation box."""
import numpy as np
import sys, os


#check if file exists  
try:
    filename = sys.argv[1]
    datafile = open(filename, "r")
except:
    print "File '{}' not found.".format(filename)
    exit()

try:
    boxsize = float(sys.argv[2]) #will be adjusted if contained in vtf file
except:
    boxsize=1e9

#faulty frame numbers are stored in np array
faulty_frames = np.array([])

#read lines in file and write to array
frame = -1

pathname = os.path.dirname(filename)
if pathname != "":
    pathname += "/"
dest_filename = pathname + "folded_" + os.path.basename(filename)
dest_file = open(dest_filename, "w")

for line in datafile:
    #current_line = line.rstrip()
    coords = True
    if line == "timestep ordered\n" or line =="t\n":
        frame += 1
        coords = False
        dest_file.write(line)
        #print "Folding Frame {}...".format(frame)
        if frame%1000 == 0:
            print frame
     
    if line[0] == "p":
        current_coords = (line.rstrip()).split()
        boxsize = float(current_coords[1])
        print "Set boxsize to {}.".format(boxsize)
        continue
    
    if frame == -1 or line == "\n":
        coords = False
        dest_file.write(line)
    
    if coords:
        current_coords = (line.rstrip()).split()
        #x_str, y_str, z_str = (line.rstrip()).split()
        if current_coords[0] == "pbc":
            boxsize = float(current_coords[1])
            print "Set boxsize to {}.".format(boxsize)
            continue
        
        try:
            x_folded = float(current_coords[0]) % boxsize
            y_folded = float(current_coords[1]) % boxsize
            z_folded = float(current_coords[2]) % boxsize
            
            dest_file.write("{0} {1} {2}\n".format(x_folded, y_folded, z_folded))
        except:
            faulty_frames = np.append(faulty_frames, frame)
            print "error in frame {}".format(frame)
            
datafile.close()
dest_file.close()

print "Faulty Frames:"
print faulty_frames
