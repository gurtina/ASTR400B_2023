
# Homework 3 Solutions
# Computing the Mass of the Local Group
# G. Besla with Rixin Li  and Hayden Foote




# import relevant modules

# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# astropy provides unit system for astronomical calculations
import astropy.units as u
# import previous HW functions
from ReadFile import Read


def ComponentMass(filename, part_type):
    """ Function to read the data from a given snapshot and return the total mass
    of the specified particle type.
    
    INPUTS
    ------
    filename: 'str'
        Name of the snapshot file to read
    part_type: 'int: 1,2,3'
        Particle type that will be summed to return mass
        
        
    OUTPUTS
    ------
    mass: 'float'
        Total mass of teh specified particle type in 1e12 solar masses
    """
  
    # read teh particle data from the specified file
    time, total, data = Read(filename)
    
    # select particles with the same type and sum up the mass
    mass = np.sum(data[data['type'] == part_type]['m'])
    
    # round and return the result in the correct units (1e12)
    return np.round(mass*1e10/1e12, 3)



# --------------- MAIN --------------- #
if __name__ == '__main__' :


    # MW: Compute Mass for each component
    ######################################
    MW_halo = ComponentMass("MW_000.txt",1)
    MW_disk = ComponentMass("MW_000.txt",2)
    MW_bulge = ComponentMass("MW_000.txt",3)

    # Total MW Mass 
    MW_total = MW_halo + MW_disk + MW_bulge
    # MW Baryon Fraction
    MW_f_bar = (MW_disk + MW_bulge) / MW_total


    # M31: Compute Mass for each component
    ########################################
    M31_halo = ComponentMass("M31_000.txt",1)
    M31_disk = ComponentMass("M31_000.txt",2)
    M31_bulge = ComponentMass("M31_000.txt",3)

    # Total M31 Mass
    M31_total = M31_halo + M31_disk + M31_bulge
    # M31 Baryon Fraction 
    M31_f_bar = (M31_disk + M31_bulge) / M31_total


    # M33: Compute Mass for each component
    #####################################
    M33_halo = ComponentMass("M33_000.txt",1)
    M33_disk = ComponentMass("M33_000.txt",2)

    # Total M33 Mass
    M33_total = M33_halo + M33_disk
    # M33 Baryon Fraction 
    M33_f_bar = M33_disk / M33_total


    # Total mass for the Local Group
    ##################################
    LG_total = MW_total + M31_total + M33_total

    # Baryon fraction for the Local Group 
    LG_f_bar = (MW_disk + MW_bulge + M31_disk + M31_bulge + M33_disk) / LG_total


    # print out table
    print()
    print("Galaxy Name  | Halo Mass   |  Disk Mass   | Bulge Mass  | Total Mass  | f_bar ")
    print("             | [1e12 Msun] |  [1e12 Msun] | [1e12 Msun] | [1e12 Msun] |       ")
    print("-------------|-------------|--------------|-------------|-------------|-------")
    print(" Milky Way   | {:<8.3f}    | {:<8.3f}     | {:<8.3f}    | {:<8.3f}    | {:<8.3f}".format(MW_halo, MW_disk, MW_bulge, MW_total, MW_f_bar))
    print(" M31         | {:<8.3f}    | {:<8.3f}     | {:<8.3f}    | {:<8.3f}    | {:<8.3f}".format(M31_halo, M31_disk, M31_bulge, M31_total, M31_f_bar))
    print(" M33         | {:<8.3f}    | {:<8.3f}     | -           | {:<8.3f}    | {:<8.3f}".format(M33_halo, M33_disk, M33_total, M33_f_bar))
    print(" Local Group | -           | -            | -           | {:<8.3f}    | {:<8.3f}".format(LG_total, LG_f_bar))
    print()


# In[ ]:




