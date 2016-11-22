#!/usr/bin/python
"""This script calculates the average for all analyze data (format of all rdf-, distribution-, formfactor data)"""
import numpy as np
import sys, os.path

filename = sys.argv[1]

#check if file exists  
try:
    datafile_check = open(filename, "r")
except:
    print "File '{}' not found.".format(filename)
    exit()

def write_coords(xdata, ydata, t_first, t_last, frames):
    y_datasets = len(ydata)
    
    #for two y datasets, ydata looks something like [[0.1,0.11,0.1], [0.09,0.1,0.1]] (len(xdata) would be 3)
    
    split_filename = os.path.splitext(filename)
    output_filename = split_filename[0] + "_average" + str(file_counter) + split_filename[1]
    output_file = open(output_filename, "w")
    
    output_file.write("# t_first={0}, t_last={1}, frames={2}\n".format(t_first, t_last, frames))
    
    #header
    output_file.write("# x")
    for i in xrange(y_datasets):
        output_file.write("\ty{}".format(i+1))
    output_file.write("\n")
    
    #data
    for i in xrange(len(xdata)):
        output_file.write("{}".format(xdata[i]))
        for j in xrange(y_datasets):
            output_file.write("\t{}".format(ydata[j][i]))
        output_file.write("\n")
    output_file.close()



#read first frame to check format
format_bins = 0
format_ydatasets = 0 #no of y datasets
found_xdata = False
found_yindices = np.array([])

for line in datafile_check:
    current_line = line.rstrip()
    split_line = current_line.split()

    if not found_xdata and split_line[0] != 'x':
        print "Cannot find x data in first line. Exit."
        exit()
    
    if found_xdata and split_line[0] == 'x':
        #case 1: next frame starts with x coordinates
        print "Initial format check of first frame finished."
        print "Found {0} datasets with {1} bins.".format(format_ydatasets, format_bins)
        
        break
    
    if split_line[0] == 'x':
        print "Read x data from first line."
        x_coords = np.array(split_line[2:], dtype=float)
        format_bins = len(x_coords)
        found_xdata = True
        
        continue
    
    #check y data
    current_y_dataset = int(split_line[0][1:])
    #print "found y dataset", split_line[0]
    
    if current_y_dataset not in found_yindices:
        #check correct order
        if current_y_dataset-1 != format_ydatasets:
            print "Wrong order of y data. Exit."
            exit()
        
        #correct order and new (unknown) y index ensured
        #add to known indices and increase y format counter
        found_yindices = np.append(found_yindices, current_y_dataset)
        format_ydatasets += 1
    
    else:
        #case 2: next frame starts with y coordinates
        print "Initial format check of first frame finished."
        print "Found {0} datasets with {1} bins.".format(format_ydatasets, format_bins)
        
        break

datafile_check.close()


#exit()


#open file again and start with averaging process
datafile = open(filename, "r")

#initialize data arrays and counter
x_coords_prev = np.zeros(format_bins)
y_data = np.zeros([format_ydatasets, format_bins])
y_datacounter = np.zeros(format_ydatasets)

x_stop_time = 0.
file_counter = 0
frame_counter = 0

for line in datafile:
    current_line = line.rstrip()
    
    split_line = current_line.split()
    
    if split_line[0] == 'x':
        #print "read in x coordinates ({} bins)".format(len(split_line[2:]))
        #x_start_time = float(split_line[1])
        x_coords = np.array(split_line[2:], dtype=float)
        
        #print x_coords
        
        if not np.array_equal(x_coords, x_coords_prev):
            #write old data to file
            print "Found new x values that differ from old ones. Write new file."
            if file_counter > 0:
                write_coords(x_coords_prev, y_data, x_start_time, x_stop_time, frame_counter)
            x_start_time = float(split_line[1])
            x_coords_prev = x_coords.copy()
            
            #reset y data and frame counter, increase file counter
            y_data = np.zeros([format_ydatasets, format_bins])
            frame_counter = 0
            y_datacounter = np.zeros(format_ydatasets)
            file_counter += 1
    
        continue

    #read y data
    current_y_dataset = int(split_line[0][1:])
    x_stop_time = float(split_line[1])
    current_y_coords = np.array(split_line[2:], dtype=float)
    
    if current_y_dataset == 1:
        frame_counter += 1 #the first frame starts with 1
        print "averaging dataset {0}.{1}".format(file_counter, frame_counter)
    
    y_datacounter[current_y_dataset-1] += 1 #separate counter for catching defect y data
    
    #print split_line[0]
    #print x_stop_time
    try:
        y_data[current_y_dataset-1] = (y_datacounter[current_y_dataset-1] - 1)/float(y_datacounter[current_y_dataset-1]) * y_data[current_y_dataset-1] + 1/float(y_datacounter[current_y_dataset-1]) * current_y_coords
    except:
        y_datacounter[current_y_dataset-1] -= 1
        print "Found defect y data. Adjust average and continue."
    
#write last frame
write_coords(x_coords, y_data, x_start_time, x_stop_time, frame_counter)

datafile.close()
