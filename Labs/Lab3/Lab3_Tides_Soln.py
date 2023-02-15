
# In Class Lab 3
# G. Besla 

# import relevant modules 
import astropy.units as u
import numpy as np
from astropy import constants as const # import astropy constants


# The Large Magellanic Cloud is at a distance of 50 kpc from the Galactic Center. 
# It is observed to have a stellar disk that extends to a radius of at least 18.5 kpc.
# 
# ![LMC](./Lab3_Tidal.png)
# Deep photometric imaging reveals the faint stellar outskirts of the LMC. 
# Outskirts: DECam data Mackey+2016 MNRAS 459, 239. 
# Inner: shallower imaging from robotic telescopes Besla+2016 APJ 825.
# 
# In this lab we will determine
# the minimum mass required for the LMC so that it maintains the observed radius 
# in the face of the Milky Way's tidal field. 

# # Part A
# 
# We define the mass profile of the Milky Way using a Hernquist profile.
# 
# 
# $\rho(r) =  \frac{M_{halo}}{2\pi} \frac{a}{r(r+a)^3} \qquad M(r) =  \frac{M_{halo} r^2}{(a+r)^2}$ 
# 
# 

# ## #1
# 
# Create a function `hernquist_mass` that returns the dark matter halo mass at a given radius in units of solar mass.
# This function should take as input:  the distance from the Galactic center $r$, the scale radius $a$, and the halo mass $M_{halo}$.
# 
# 
# For the Hernquist scale radius for the Milky Way, use $a=60$ kpc. 
# 
# For $M_{halo}$ use your answer for the total mass of the simulated Milky Way you computed in Homework 3. 




def hernquist_mass(r,a=60*u.kpc, m_halo=1.96):
    """ Function that defines the Hernquist 1990 mass profile 
    Inputs:
        r: astropy quantity
            Galactocentric distance in kpc
        a: astropy quantity
            scale radius of the Hernquist profile in kpc
        m_halo: float
            total halo mass in units of 1e12 Msun 
        
    Ouputs:
        mass:  astropy quantity
            total mass within the input radius r in Msun
    """
    mass = m_halo*1e12*r**2/(a+r)**2*u.Msun # Hernquist mass  
    return mass




print(f"{hernquist_mass(10000*u.kpc):.2e}")




print(f"{hernquist_mass(260*u.kpc):.2e}")




print(f"{hernquist_mass(50*u.kpc):.2e}")


# ## #2
# 
# Compute the total mass of the Milky Way within 50 kpc, including it's baryonic component (Dark Matter + Bulge + Disk), in units of solar mass. Use your answers from Homework 3 for the Bulge and Disk Masses. 
# Store this as a variable called `mass_MW50`.
# 



# MW disk mass and bulge mass
# using answers from assignment 3 
mdisk =  7.5e10*u.Msun
mbulge = 1e10*u.Msun




# determine the total mass of the MW within 50 kpc
# Kochanek+1996 find 4.9e11
mass_MW50 = hernquist_mass(50*u.kpc) + mdisk + mbulge
print(f"{mass_MW50:.2e}")
mass_MW50


# # Part B
# 
# The Jacobi Radius for a satellite on a circular orbit about an extended host, where 
# the host is assumed to be well modeled as an isothermal sphere halo:
# 
# 
# $R_j = r  \bigg( \frac{M_{sat}}{2 M_{host}(<r)} \bigg)^{1/3}$
# 
# 
# The Isothermal Sphere approximation is not a bad one within 50 kpc.
# 
# Note also that the LMC is not on a circular orbit, but it is very close to its pericentric approach, where the velocity is all in the tangential component. So this isn't a terrible approximation either. 
# 
# ## #1
# Create a function called `jacobi_mass` that returns the total mass of a satellite galaxy in units of Msun, 
# such that it has a given size 
# 
# Do this by rearranging the Jacobi Radius equation to solve for the satellite mass. 
# 



# Rj = r*(Msat/2Mhost)**(1/3)
# Msat =  (Rj/r)**(3)*2*Mhost

def jacobi_mass(rj,r,m_host):
    """ Function that determines the minimum satellite mass needed
    to maintain a the size of a given satellite using the Jacobi Radius
    
    Inputs:
        rj : astropy quantity
            Jacobi Radius or the stellar radius of the satellite in kpc
        r : astropy quantity 
            Distance of the satellite from the host in kpc
        m_host: astropy quantity 
            Mass of the host galaxy in Msun within r in Msun
        
    Outputs:
        m_min: astropy quantity
            Minimum satellite mass in Msun
    """
    
    m_min = (rj/r)**3*2*m_host # satellite mass from Jacobi Radius 
    
    return m_min
    


# ## #2 
# 
# Determine the minimum total mass of the LMC needed to maintain its radius of 18.5 kpc in the face of the Milky Way's tidal 
# field at its current distance of 50 kpc. Store this as a variable called `LMC_jacobiM`.




sizeL = 18.5*u.kpc # Observed radius of/the LMC disk (Mackey+2016)

distL = 50.0*u.kpc # Galactocentric Distance to the LMC




# LMC minimum mass in maximal MW halo model (from Simulation)
LMC_jacobiM = jacobi_mass(sizeL,distL,mass_MW50)
print(f"{LMC_jacobiM:.2e}")
LMC_jacobiM


# ## #3
# 
# Recall that, ignoring centrifugal forces and assuming the host is a point mass, the tidal radius is given as :
# 
# $r_{tide} = r\left (\frac{m_{sat}}{4M_{host} } \right)^{1/3} $
# 
# Create a function to compute the total mass the must LMC possess to have a disk with radius 18.5 kpc.
# 
# The function should take as input, the current size of the satellite (kpc), this distnce ot the host(kpc) and the mass of the host (in Msun)
# 
# Use the function to determine the needed LMC mass and store it as a variable called `LMC_tidalM`. 




# rtide = r(msat/4/Mhost)**1/3
#msat = (rtide/r)**(3)*4*Mhost


def tidal_mass(r_tide, r, m_host):
    
    """ Function to compute the mass of the satellite needed such that
    the tidal radius is teh current size of the galaxy, ignoring
    centrifugal forces
    
    Inputs:
        r_tide: astropy quantity
            The tidal radius in kpc
        r : astropy quantity 
            The distance from the host galaxy in kpc
        mhost : astropy quantity 
            Mass of the host galaxy in Msun
        
    Output:
        sat_mass:  astropy quantity
            The mass of the satellite in Msun
    """
    sat_mass = 4*m_host*(r_tide/r)**3
    return sat_mass




LMC_tidalM = tidal_mass(sizeL,distL,mass_MW50) # Tidal mass of the LMC
print(f"{LMC_tidalM:.2e}")
LMC_tidalM


# ## #4
# 
# Compare `LMC_tidalM` to the calculation using the Jacobi Radius.
# How does the total mass of the LMC compare to its stellar mass (M$_\ast = 3 \times 10^9$ M$_\odot$)? 
# 



print(LMC_tidalM/LMC_jacobiM)
# Because of centrifugal forces the minimum mass is smaller
# using the Jacobi Radius. 




LMC_mstar = 3e9*u.Msun # LMC stellar mass




print(np.round(LMC_jacobiM/LMC_mstar))
# Mass ratio is ~ factor of 20. 


# # Part C: Consistency Check
#  NOT REQUIRED
# 
# The dynamical mass of the LMC at any radius can be determined by assuming a flat rotation curve.  "Dynamical mass" means mass needed to explain the rotation curve. 
# 
# $V_c^2 = \frac{G M}{r} = constant$
#  
#  The rotation curve of the LMC is observed to flatten at a value of 91.7 +/- 18.8 km/s  (van der Marel & Kallivayalil 2014 ApJ 781)
#  
#  
# Use `lambda` to create a function called `dyn_mass` that takes as input Vc (km/s) and distance to from the center of the galaxy (r) and returns the maximal dynamical mass in Msun. 
#  
#  $G$ = 4.4985e-6 kpc$^3$/Gyr$^2$/M$_\odot$
# 



const.G




# Gravitational Constant in the desired units

Grav = const.G.to(u.kpc**3/u.Gyr**2/u.Msun)




# Maximal LMC mass
# Assuming LMC has a flat rotation curve to 18.5 kpc
# Vc = 91.7 +/- 18.8 km/s  van der Marel & Kallivayalil 2014
vc_LMC =(91.7*u.km/u.s).to(u.kpc/u.Gyr)
print(vc_LMC)




# MLMC Vc^2 = GM/R = constant, rearrange for M:
# M = Vc**2*R/G
#def dyn_mass(vc, r):
#    return vc**2*r/Grav
dyn_mass = lambda vc, r: vc**2*r/Grav


#  
# ##  #1  
# Compute the dynamical mass enclosed by the LMC within the observed radius. Store it as a variable called `LMC_dynM`. 



LMC_dynM = dyn_mass(vc_LMC,sizeL) 
print(f"{LMC_maxM:.2e}")
LMC_dynM


# ## #2
# 
# Is `LMC_maxM` consistent with `LMC_jacobiM`, the minimum mass needed to explain the observed radius of the LMC given the tidal field of the MW? If not, how could the numbers be reconciled?



print(LMC_dynM/LMC_jacobiM)


# The minimum mass needed seems larger than the maximal mass possible.
# Either the LMC rotation curve needs to be higher (which it could within the errors)
# Or MW halo mass within 50 kpc is smaller, e.g. $3\times 10^{11}$ M$_\odot$.
# 
# Although note, the values are pretty close to being the same. 



# Try increasing the LMC circular speed by 1 sigma
v_new = vc_LMC + 18.8*u.kpc/u.Gyr
print(f"{dyn_mass(v_new,sizeL):.2e}")




# LMC minimum mass in lower mass MW halo model  
# (recall Hernquist model gives 4e11 Msun within 50 kpc)
min_MW = 3e11*u.Msun
print(f"{jacobi_mass(sizeL,distL,min_MW):.2e}")







