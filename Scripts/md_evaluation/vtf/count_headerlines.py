#!/usr/bin/python
"""This script counts the number of header lines of a vtf file"""
import sys

def getHeaderlines(vtfFilename):
    n_headerlines = 0
    datafile = open(vtfFilename, "r")
    for idx, line in enumerate(datafile):
        if idx == 0 and line != "timestep ordered\n" and line != "t\n":
            continue
        
        current_line = line.rstrip()
        split_line = current_line.split()
        
        if line == "timestep ordered\n" or line =="t\n":
            n_headerlines = idx
            break

    datafile.close()
    return n_headerlines


if __name__ == "__main__":
    print getHeaderlines(sys.argv[1])
