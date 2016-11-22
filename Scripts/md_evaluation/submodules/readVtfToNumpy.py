import numpy as np
from progressBar import printProgress

def readVtf(vtfFilename):
    """Read in vtf file to multidimensional numpy array with format coords[timestep][particle_id] = [x, y, z]
    Return values are
    coords[timestep][particle_id] = [x, y, z]
    list_ptypes: sorted list of all particle types
    ptypes[particle_id] = particle_type
    valencies[particle_id] = valency
    
    n_part and n_timesteps can be taken from shape of coords
    """
    
    # read in number of particles from header
    n_part = 0
    n_headerlines = 0
    datafile = open(vtfFilename, "r")
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

    datafile.close()
    n_part += 1
    print "Number of particles =", n_part
    print "Number of header lines =", n_headerlines
    print ""
    
    
    ptypes = np.zeros(n_part, dtype=int)
    list_ptypes = np.array([])
    valencies = np.zeros(n_part, dtype=float)


    # read in particle types from header
    datafile = open(vtfFilename, "r")
    for idx, line in enumerate(datafile):
        if idx >= n_headerlines:
            break
        
        current_line = line.rstrip()
        split_line = current_line.split()
        
        if split_line[0] == "atom":
            for part_ranges in split_line[1].split(","):
                start_range = int(part_ranges.split(":")[0])
                end_range = int(part_ranges.split(":")[1])
                ptype = int(split_line[7])
                ptypes[start_range:end_range+1] = ptype
                if not ptype in list_ptypes:
                    list_ptypes = np.append(list_ptypes, ptype)
                valency = float(split_line[9])
                valencies[start_range:end_range+1] = valency
                
                print "Found type {0} with valency {1} for particle ids {2}-{3}.".format(ptype, valency, start_range, end_range)

    #sort particle types
    list_ptypes = np.sort(list_ptypes)

    datafile.close()

    datafile_countframes = open(vtfFilename, "r")
    def count_frames(datafile):
        counter_frames = 0
        
        for line in datafile:
            if line == "timestep ordered\n" or line =="t\n":
                counter_frames += 1    
        
        return counter_frames
    print "\nCount total frames..."
    n_timesteps = count_frames(datafile_countframes)
    print "Found {} frames.\n".format(n_timesteps)
    datafile_countframes.close()
    """
    for idx, val in enumerate(ptypes):
        print "ptypes[{0}] = {1}".format(idx, val)
    """


    #read in particle coordinates to multidimensional numpy array
    # format: coords[timestep][particle_id] = [x, y, z]
    datafile_numpy = open(vtfFilename, "r")

    coords = np.zeros((n_timesteps, n_part, 3))
    timestep_coords = np.zeros((n_part, 3))
    part_id = 0
    timestep = 0
    
    print "Read in coordinates..."

    for idx, line in enumerate(datafile_numpy):
        #skip headerlines plus first timestep line
        if idx < n_headerlines+1:
            continue
        
        #search for timestep string
        if line == "timestep ordered\n" or line =="t\n":
            if timestep%100==0:
                #print "Read in timestep", timestep
                printProgress(timestep, n_timesteps)
            
            coords[timestep] = timestep_coords
            timestep_coords = np.zeros((n_part, 3))
            
            part_id = 0
            timestep += 1
            continue
        
        #with correct format this should not happen
        if line == "\n":
            continue
        
        current_line = line.rstrip()
        split_line = current_line.split()

        timestep_coords[part_id] = np.array(split_line, dtype=float)
        part_id += 1

    #coordinates from last frame
    coords[timestep] = timestep_coords
    len_coords = len(coords)

    printProgress(n_timesteps, n_timesteps)

    datafile_numpy.close()

    print ""
    print "\nShape of coordinate matrix:", coords.shape
    print "Memory usage of coordinate matrix:", round(coords.nbytes/1024.**2, 2), "mb\n"
    
    
    #return values
    return coords, list_ptypes, ptypes, valencies
