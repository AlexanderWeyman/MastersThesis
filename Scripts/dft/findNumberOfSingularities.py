#!/usr/bin/python
"""This script finds singularities by checking if the sign of the slope changes between previous, current and next datapoint. 
Besides the relation between absolute slopes has to be greater alpha."""

import sys, os.path
import numpy as np


defaultAlpha = 1.001
minSlopeSum = 10.

limit_range = True
qmin = 0.1
qmax = 1000.

usage = "Usage: {0} filename [alpha={1}]".format(os.path.basename(sys.argv[0]), defaultAlpha)

verboseOutput = False

#check if file exists  
try:
    filename = sys.argv[1]
    data = np.genfromtxt(filename, skip_header=1)
except:
    print "Error reading data. Exit."
    print usage
    exit()

try:
    alpha = float(sys.argv[2])
except:
    alpha = defaultAlpha


xColumn = 0
yColumn = 7

xData = data[:,xColumn]
yData = data[:,yColumn]
yDataMean = yData.mean()

deltaX = xData[1]-xData[0]
minSlopeSum *= deltaX

xSingularities = np.array([])


if len(xData) != len(yData):
    print "x and y data differ in the number of datapoints. Exit."
    exit()

#all(xData[i] <= xData[i+1] for i in xrange(len(xData)-1))

for idx in xrange(1, len(yData)-1):
    if xData[idx] < xData[idx-1]:
        print "Found unsorted x data. Exit."
        exit()
    
    if limit_range and (xData[idx] < qmin or xData[idx] > qmax):
        continue
        
    #print np.abs(yData[idx] - yData[idx-1])
    #if np.abs(yData[idx] - yData[idx-1]) > epsilon:
    #if np.abs(yData[idx] - yData[idx-1]) > epsilon and yData[idx]*yData[idx-1] < 0:
    #if (yData[idx+1]-yData[idx])*(yData[idx]-yData[idx-1]) and np.abs(yData[idx+1] - 2*yData[idx] + yData[idx-1]) > epsilon:
    
    if (yData[idx]-yData[idx-1])*(yData[idx+1]-yData[idx]) < 0 and np.abs((yData[idx+1]-yData[idx])/float(yData[idx]-yData[idx-1])) > alpha:
        appVal = xData[idx] if np.abs(yData[idx]-yDataMean) > np.abs(yData[idx+1]-yDataMean) else xData[idx+1]
        #if appVal not in xSingularities and yData[idx-1:idx+2].std() > 1e-6:
        #print np.abs(yData[idx+1]-yData[idx])+np.abs(yData[idx]-yData[idx-1]), appVal
        if appVal not in xSingularities and (np.abs(yData[idx+1]-yData[idx])+np.abs(yData[idx]-yData[idx-1])) > minSlopeSum:
            xSingularities = np.append(xSingularities, appVal)
            #print yData[idx-1:idx+2].std()
if verboseOutput:
    print "Found {0} singularities in interval [{1}:{2}]:\n".format(len(xSingularities), xData[0], xData[-1])
    for q in xSingularities:
        print q

print len(xSingularities)
