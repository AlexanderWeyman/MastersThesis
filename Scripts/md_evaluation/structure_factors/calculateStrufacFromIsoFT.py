#!/usr/bin/python
"""This script calculates the partial structure factors from the h(q) values generated with calculateIsoFTFromRDF.py"""

#calculation of
# S_ij(q) = N_i/N (delta_i,j + N_j/V h_ij(q))
import numpy as np
import sys, os.path


#check if file exists  
try:
    filename = sys.argv[1]
    avg_data = np.genfromtxt(filename)
except:
    print "File not found."
    exit()

try:
    delta_s = int(sys.argv[2])
except:
    print "No valid value for the sidechain length found."
    exit()


N_P = 30 # number of polymer chains
N_M = 30 # degree of polymerization
box_l = np.array([14.645918834126679, 16.765391982241994, 18.45270143416798, 19.877570047286294]) #box_l[delta_s] for eta=0.3
V = box_l**3

N_B = N_P * N_M
N_T = N_P * N_M
N_C = N_P * N_M
N = (delta_s + 2) * N_P * N_M


q_vals = avg_data[:,0]

if delta_s > 0:
    h_BB = avg_data[:,1] # backbone - backbone
    h_TT = avg_data[:,2] # terminal - terminal
    h_BT = avg_data[:,3] # backbone - terminal
    h_TC = avg_data[:,5] # terminal - counterion
    h_CC = avg_data[:,6] # counterion - counterion
    h_NN = avg_data[:,8] # particle - particle
else:
    h_BB = avg_data[:,2]
    h_TT = avg_data[:,2]
    h_BT = avg_data[:,2]
    h_TC = avg_data[:,5]
    h_CC = avg_data[:,6]
    h_NN = avg_data[:,8]


S_BB = N_B/float(N) * (1 + N_B/float(V[delta_s]) * h_BB)
S_TT = N_T/float(N) * (1 + N_T/float(V[delta_s]) * h_TT)
S_BT = N_B/float(N) * N_T/float(V[delta_s]) * h_BT
S_TC = N_T/float(N) * N_C/float(V[delta_s]) * h_TC
S_CC = N_C/float(N) * (1 + N_C/float(V[delta_s]) * h_CC)
S_QQ = S_TT + S_CC - 2*S_TC # charge - charge structure factor
S_NN = 1 + N/float(V[delta_s]) * h_NN


split_filename = os.path.splitext(filename)
output_filename = split_filename[0] + "_Strufac" + split_filename[1]
output_file = open(output_filename, "w")

#header
output_file.write("#q\tS_BB\tS_TT\tS_BT\tS_TC\tS_CC\tS_QQ\tS_NN\n")

for i in xrange(len(q_vals)):    
    output_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(q_vals[i], S_BB[i], S_TT[i], S_BT[i], S_TC[i], S_CC[i], S_QQ[i], S_NN[i]))

output_file.close()

print "Success. Wrote to '{}'".format(output_filename)
