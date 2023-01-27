
# GB : ASTR 400B Solutions to PSet 2 
# ParticleInfo

# load modules
import numpy as np
import astropy.units as u
from ReadFile import Read




def ParticleInfo(PType, PNum, filename):
  
    """ Function to return properties of a particle of a given type
    
    Input: 
        PType: int
            particle type, e.g. Halo: 1, Disk: 2, Bulge: 3
        PNum: int 
            particle number, e.g. 100)
        filename: str
            e.g. "MW_000.txt")
        
    Output: 
        R3D: astropy quantity
            Magnitude of 3D Pos vector (kpc)
        V3D: astropy quantity
            Magnitude of 3D Velocity vector (km/s)
        Mass: astropy quantity
            Mass of the Particle (Msun)
    """


    # read in the file 
    time, total, data = Read(filename)

    
    #create an array to store indexes of particles of desired Ptype
    index = np.where(data['type'] == PType)

    # create new arrays with the m, x, y, z, 
    # vx, vy, vz of just the Ptype desired
    # Add units using Astropy
    # Recall Mass was stored in units of Msun/1e10
    mnew = data['m'][index]*1e10*u.Msun
    xnew = data['x'][index]*u.kpc
    ynew = data['y'][index]*u.kpc
    znew = data['z'][index]*u.kpc
    vxnew = data['vx'][index]*u.km/u.s
    vynew = data['vy'][index]*u.km/u.s
    vznew = data['vz'][index]*u.km/u.s
    
    # Compute the Magnitude of the 3D position
    # Value is rounded to 3 decimal places.
    R3D = np.round(np.sqrt(xnew[PNum-1]**2 + ynew[PNum-1]**2 + znew[PNum-1]**2),3)
    
    # Compute the magnitude of the 3D velocity
    # Value is rounded to 3 decimal places.
    V3D = np.round(np.sqrt(vxnew[PNum-1]**2 + vynew[PNum-1]**2 + vznew[PNum-1]**2),3)
    
    # Mass
    # Value is rounded to 3 decimal places
    Mass = np.round(mnew[PNum-1],3)
        
    return R3D, V3D, Mass
    
    



R3D, V3D, Mass = ParticleInfo(2,100,"MW_000.txt")

print("Magnitude of the distance of the 100th particle is",R3D)
print("Magnitude of the velocity of the 100th particle is",V3D)
print("Mass of the 100th particle is", Mass)
print("Distance of the 100th particle in lightyears is", np.round(R3D.to(u.lyr),3))




