
# Homework 2
# G Besla

# Load Modules
import numpy as np
import astropy.units as u


def Read(filename):
    """ Function to read in our data file
    
    Input:  
        filename: str
            e.g. "MW_000.txt"
        
    Outputs: 
        time: astropy quantity
            Time of snapshot in Myr
        total: float
            Total number of particles 
        data: array of floats
            An array with the particle data, including position 
            vectors, velocity vectors and mass
            
    Example usage:  time, total, data = Read("filename")
    """
    
    
    # open the file 
    file = open(filename,'r')
    
    #read header info line by line (line will be a string)
    # read first two lines FIRST and store as variable
    
    # read and store time
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Myr

    # read and store total number of particles
    line2 = file.readline()
    label, value = line2.split()
    total = float(value)
    
    # close file
    file.close()

    # read the remainder of the file, 
    # "dtype=None" specifies data type. None is default float
    # default delimiter is line is split using white spaces
    # "skip_header=3"  skipping the first 3 lines 
    # the flag "names=True" creates arrays to store the date
    #       with the column headers given in line 4 like "m", "x"
    
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    
    # Note, another option is loadtxt, skipping the first 3 rows.  
    # data = np.loadtxt(filename,skiprows=3)
    # But this loses the information in the headers
    
    # this will return the time of the snapshot, 
    #total number of particles 
    #and an array that stores the remainder of the data. 
    return time, total, data


# Checking to see if the code works: (uncomment lines below)

#time, total, data =Read("MW_000.txt")

#print("Time of the snapshot", time)
#print("The total number of particles is", total)

#Mass of first particle, if you used loadgenfromtxt
#m1 = np.round(data['m'][0]*u.Msun*1e10)
#print("The mass of the first particle is", m1)

# if you wanted to use loadtxt, 
# first store the mass of all particles in a new array
# i.e. store the 2nd column
#Mass = data[:,1]
# Print the Mass of the first particle
#Mass[0]*u.Msun*1e10

# Type of first particle
#print("Data type of the first Particle is", data['type'][0])





