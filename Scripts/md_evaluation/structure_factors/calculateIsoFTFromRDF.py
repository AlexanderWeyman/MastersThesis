#!/usr/bin/python
"""This script calculates the isotropic fourier transform of pair correlation function h(r)=g(r)-1. It is meant to use files generated from averageAllAnalyzeData.py"""

#calculation of
# 4 pi/q int_0^infty r h(r) sin(q r) dr
import numpy as np
import sys, os.path


#check if file exists  
try:
    filename = sys.argv[1]
    avg_data = np.genfromtxt(filename)
except:
    print "File not found."
    exit()



y_datasets = len(avg_data[0])-1 #number of y colums

xdata = avg_data[:,0]
#ydata = avg_data[:,ydataset]

x_min = xdata[1]-xdata[0]
x_max = xdata[-1]

q_min = 2*np.pi/x_max
q_max = 2*np.pi/x_min
q_bins = 1000

q_max /=2. #some artificial factor to get rid of oscillations for large q values due to binning

q_vals = np.linspace(q_min, q_max, q_bins)


def integrate_iso_FT(q, y_vals, delta_x=x_min):
    return 4*np.pi/q * delta_x * (xdata*(y_vals-1.)*np.sin(q*xdata)).sum()
    

iso_FT_vals = np.zeros([y_datasets,q_bins])
for i in xrange(q_bins):
    for j in xrange(y_datasets):
        iso_FT_vals[j][i] = integrate_iso_FT(q_vals[i],avg_data[:,j+1])

"""
for i in xrange(len(q_vals)):
    print "{0}\t{1}".format(q_vals[i], iso_FT_vals[i])
"""

split_filename = os.path.splitext(filename)
output_filename = split_filename[0] + "_isoFT" + split_filename[1]
output_file = open(output_filename, "w")

output_file.write("#q")
for i in xrange(y_datasets):
    output_file.write("\ty{}".format(i+1))
output_file.write("\n")

#output_file.write("#q\tisoFT(h(r))(q)\n")
for i in xrange(q_bins):
    output_file.write("{}".format(q_vals[i]))
    for j in xrange(y_datasets):
        output_file.write("\t{}".format(iso_FT_vals[j][i]))
    output_file.write("\n")
    #output_file.write("{0}\t{1}\n".format(q_vals[i], iso_FT_vals[i]))
output_file.close()

print "Success. Wrote to '{}'".format(output_filename)

