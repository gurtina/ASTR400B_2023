

# In Class Lab 6
# Surface Brightness Profiles




# Load Modules
import numpy as np
import astropy.units as u

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib
get_ipython().run_line_magic('matplotlib', 'inline')

# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile
from GalaxyMass import ComponentMass


# # Lab 6: Sersic Profiles
# 
# In this lab we will use Homework 5 solutions to compute the mass profile of the Milky Way's bulge. 
# We will turn the mass profile into a density profile and see if we can fit it reasonably well 
# with a sersic profile. 

# # Part A : 
# 
# Create a function called `sersicE` that returns the Sersic Profile in terms of the effective radius $R_e$ (i.e. the half light radius).
# 
# $I(r) = I_e exp^{-7.67 ( (r/R_e)^{1/n} - 1)}$
# 
# Where 
# 
# $ L = 7.2 I_e \pi R_e^2$
# 
# and  $R_e$ is the half light radius.  We will assume a mass to light ratio for the stellar bulge of 1, so this is also the half mass radius.
# 
# The function should take as input: the radius, $R_e$, $n$ and the total stellar mass of the system.
# 



  


# # Part B
# 
# a) Create an instance of the MassProfile Class for M31. Store it as a variable `M31`. 
# 





# b) Create an array of radii from 0.1 kpc to 30 kpc in increments of 0.1
# 




# c) Define a new array called `bulge_mass`, that uses the function `MassEnclosed` within MassProfile to compute the mass profile of the bulge.  Get rid of astropy units in `bulge_mass` by adding `.value` 
# 






# d) Compute the surface mass density profile for the simulated bulge and store it as an array called `bulge_I`. Assuming M/L ~ 1 this is also the surface brightness profile in Lsun/kpc^2






# # Part C
# 
# Compute $R_e$, the half mass radius, for the bulge






# # Part D
# 
# a) Plot the surface density profile of the simulated bulge
# 
# b) Plot the Sersic profile, assuming a de Vaucouleurs Profile.
# 
# c) If the profiles don't match, try changing either $R_e$ or $n$

# In[ ]:


# Plot the Bulge density profile vs 
# the Sersic profile
####################################


fig = plt.figure(figsize=(8,8))
ax = plt.subplot(111)


# plot the bulge mass density as a proxy for surface brighntess


# YOU ADD HERE: Sersic fit to the surface brightness Sersic fit
# Sersic 



#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size


# Add axis labels
plt.xlabel('Radius (kpc)', fontsize=22)
plt.ylabel('Log(I)  $L_\odot/kpc^2$', fontsize=22)



# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')


plt.savefig('Lab6.png')






