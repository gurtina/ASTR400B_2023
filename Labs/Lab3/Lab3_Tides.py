

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
# For the Hernquist scale radius for the Milky Way, use the default value of $a=60$ kpc. 
# 
# For $M_{halo}$ use your answer for the total mass of the simulated Milky Way you computed in Homework 3 as the default value (in units of 1e12). 






# ## #2
# 
# Compute the total mass of the Milky Way within 50 kpc, including its baryonic components (Dark Matter + Bulge + Disk). Use your answers from Homework 3 for the Bulge and Disk Masses. 
# Store this as a variable called `mass_MW50`.
# 






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






# ## #2 
# 
# Determine the minimum total mass of the LMC needed to maintain its radius of 18.5 kpc in the face of the Milky Way's tidal 
# field at its current distance of 50 kpc. Store this as a variable called `LMC_jacobiM`.






# ## #3
# 
# Recall that, ignoring centrifugal forces and assuming the host is a point mass, the tidal radius is given as :
# 
# $r_{tide} = r\left (\frac{m_{sat}}{4M_{host} } \right)^{1/3} $
# 
# Create a function to determine the total mass the must LMC possess to have a disk with radius 18.5 kpc.
# 
# The function should take as input, the current size of the satellite (kpc), this distance to the host (kpc) and the mass of the host (in Msun)
# 
# Use the function to determine the needed LMC mass and store it as a variable called `LMC_tidalM`. 






# ## #4
# 
# a) Compare `LMC_tidalM` to the calculation using the Jacobi Radius.
# 
# b)How does the total mass of the LMC compare to its stellar mass (M$_\ast = 3 \times 10^9$ M$_\odot$)? 
# 






# # Part C: Consistency Check
# 
# 
# The dynamical mass of the LMC at any radius can be determined by assuming a flat rotation curve.  "Dynamical mass" means mass needed to explain the rotation curve. 
# 
# $V_c^2 = \frac{G M}{r} = constant$
#  
#  The rotation curve of the LMC is observed to flatten at a value of 91.7 +/- 18.8 km/s  (van der Marel & Kallivayalil 2014 ApJ 781)
# 
#   
#    Create a function called `dyn_mass` that takes as input Vc (km/s) and distance to from the center of the galaxy (r) and returns the maximal dynamical mass in Msun. 
#  
# $G$ = 4.4985e-6 kpc$^3$/Gyr$^2$/M$_\odot$






#   
# ## #1  
# Compute the dynamical mass enclosed by the LMC within the observed radius. Store it as a variable called `LMC_dynM`. 






# ## #2
# 
# Is `LMC_dynM` consistent with `LMC_jacobiM`, the minimum mass needed to explain the observed radius of the LMC given the tidal field of the MW? If not, how can the numbers be reconciled?





