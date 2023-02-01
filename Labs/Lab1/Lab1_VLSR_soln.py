#!/usr/bin/env python
# coding: utf-8

# # In Class Lab 1
# 
# ## Must be uploaded to your 'Labs/Lab1' Github repository by 5 PM Jan 31st
# 
# ## Part A:  The Local Standard of Rest
# Proper motion of Sgr A* from Reid & Brunthaler 2004
# $\mu = 6.379$ mas/yr 
# 
# Peculiar motion of the sun, $v_\odot$ = 12.24 km/s  (Schonrich 2010)
# 
# 
# $v_{tan} = 4.74 \frac{\mu}{\rm mas/yr} \frac{R_o}{\rm kpc} = V_{LSR} + v_\odot$
# 
# 
# ### a)
# 
# Create a function called VLSR to compute the local standard of res (V$_{LSR}$).
# 
# The function should take as input: the solar radius (R$_o$), the proper motion (mu)
# and the peculiar motion of the sun in the $v_\odot$ direction.
# 
# Compute V$_{LSR}$ using three different values R$_o$: 
# 1. Water Maser Distance for the Sun :  R$_o$ = 8.34 kpc   (Reid 2014 ApJ 783) 
# 2. GRAVITY Collaboration Distance for the Sun:  R$_o$ = 8.178 kpc   (Abuter+2019 A&A 625)
# 3. Value for Distance to Sun listed in Sparke & Gallagher : R$_o$ = 7.9 kpc 
# 

# In[17]:


# Import Modules 
import numpy as np # import numpy
import astropy.units as u # import astropy units
from astropy import constants as const # import astropy constants


# In[2]:


def VLSR(Ro, mu=6.379, vsun=12.24*u.km/u.s):
    """This function will compute the velocity at the local standard of rest
                VLSR = 4.74*mu*Ro - vsun
                
        Inputs:  
            Ro:  'astropy quantity'
                The distance from teh Sun to the Galactic Center in kpc
            mu:  'float'
                The proper motion of Sag A* in mas/yr. 
                Default is from Reid & Brunthaler 2004
            vsun: 'astropy quantity'
                The peculiar motion of the Sun in the v direction (km/s)
                Default is from Schonrich + 2010
         
         Outputs:
             VLSR: 'astropy quantity'
                 The velocity of the local standard of rest (km/s)
        
    """
    
    return  4.74*mu*(Ro/u.kpc)*u.km/u.s - vsun


# In[3]:


# Define our distances to the Galactic Center from the Sun

RoReid = 8.34*u.kpc # Distance from Reid+2014 in kpc
RoGravity = 8.178*u.kpc # Distance from the Gravity Collab Abuter+2019 in kpc 
RoSG = 7.9*u.kpc # Distance from the textbook Sparke & Gallagher


# In[4]:


# Compute VLSR using Ro from Reid 2014
VLSR_Reid = VLSR(RoReid)
print(VLSR_Reid)


# In[5]:


# Compute VLSR using Ro from Gravity Collab
VLSR_Gravity = VLSR(RoGravity)
print(np.round(VLSR_Gravity))


# In[6]:


# Compute VLSR using Ro from Sparke & Gallagher
VLSR_SG = VLSR(RoSG)
print(np.round(VLSR_SG))


# ### b)
# 
# compute the orbital period of the sun in Gyr using R$_o$ from the GRAVITY Collaboration (assume circular orbit)
# 
# Note that 1 km/s $\sim$ 1kpc/Gyr

# In[7]:


def TorbSun(R, V):
    """  This Function will compute the orbital period of the Sun 
                T = 2 pi R / V
        
        Inputs:
            R: 'astropy quantity'
                Distance in kpc (Distance to the Galactic center)
            V: 'astropy quantity'
                Velocity in km/s (Velocity of the sun in v direction)
        
        Outputs:
            T:  'astropy quantity'
            Orbital period in Gyr
    """
    
    VkpcGyr = V.to(u.kpc/u.Gyr) # converting v from km/s to kpc/Gyr
    T = 2*np.pi*R/VkpcGyr # Orbital Period 
    
    return T
    
    
    


# In[8]:


# Velocity of the sun = VLSR + peculiar motion 
VsunPeculiar = 12.24*u.km/u.s
VSun = VLSR_Gravity  + VsunPeculiar


# In[9]:


# Compute the orbital period of the Sun! 
# Use Ro from Gravity Collaboration 
T_Grav = TorbSun(RoGravity, VSun)
print(T_Grav) 


# ### c)
# 
# Compute the number of rotations about the GC over the age of the universe (13.8 Gyr)

# In[10]:


# Age of Universe/ Orbital Period

Age  = 13.8*u.Gyr # Age of the Universe
print(Age/T_Grav)


# ## Part B:  Dark Matter Density Profiles
# 
# ### a)
# Try out Fitting Rotation Curves 
# [here](http://wittman.physics.ucdavis.edu/Animations/RotationCurve/GalacticRotation.html)
# 
# 
# ### b)
# 
# 
# In the Isothermal Sphere model, what is the mass enclosed within the solar radius (R$_o$) in units of M$_\odot$? 
# 
# Recall that for the Isothermal sphere :
# $\rho(r) = \frac{V_{LSR}^2}{4\pi G r^2}$
# 
# Where $G$ = 4.4985e-6 kpc$^3$/Gyr$^2$/M$_\odot$, r is in kpc and $V_{LSR}$ is in km/s
# 
# What about at 260 kpc (in units of  M$_\odot$) ? 

# In[18]:


# Gravitational Constant in the desired units

Grav = const.G.to(u.kpc**3/u.Gyr**2/u.Msun)


# In[19]:


# density profile  rho = VLSR^2 / (4*pi*G*R^2)
# Mass = Integrate rho dV 
#      = rho  4*pi*r**2 dr 
#      = VLSR**2 /G/(4*pi*r**2)  * (4*pi*r**2) dr 
#      =  VLSR**2 *r / G 

def MassIso(r, VLSR):
    """ This function will comput the dark matter mass enclosed within a given distance
        assuming an Isothermal Sphere Model for the dark matter.
        M = VLSR**2 / G * r
        
        Inputs:
        r : 'astropy quantity'
            Distance to the Galactic Center (kpc)
        VLSR:  'astropy quantity'
            Velocity of the Local Standard of Rest (km/s)
        
        Outputs:
        M : Mass enclosed within r in units of Msun 
    """

    VLSRkpcGyr = VLSR.to(u.kpc/u.Gyr) #  converting km/s to kpc/Gyr
    
    M = VLSRkpcGyr**2 / Grav * r # Mass for isothermal sphere
    
    return M 
    


# In[20]:


MIsoSolar = MassIso(RoGravity, VLSR_Gravity)
print(MIsoSolar)


# In[21]:


# print in scientific notation

print(f"{MIsoSolar:.2e}")


# In[22]:


# Compute mass within 260 kpc
MIso260 = MassIso(260*u.kpc, VLSR_Gravity)
print(f"{MIso260:.2e}")


# ## c) 
# 
# The Leo I satellite is one of the fastest moving satellite galaxies we know. 
# 
# 
# It is moving with 3D velocity of magnitude: Vtot = 196 km/s at a distance of 260 kpc (Sohn 2013 ApJ 768)
# 
# If we assume that Leo I is moving at the escape speed:
# 
# $v_{esc}^2 = 2|\Phi| = 2 \int G \frac{\rho(r)}{r}dV $ 
# 
# and assuming the Milky Way is well modeled by a Hernquist Sphere with a scale radius of $a$= 30 kpc, what is the minimum mass of the Milky Way (in units of M$_\odot$) ?  
# 
# How does this compare to estimates of the mass assuming the Isothermal Sphere model at 260 kpc (from your answer above)

# In[23]:


# Potential for a Hernquist Profile
# Phi = -G *M / (r+a)

#  Using the potential for a Hernquist Profile, the equation for the escape speed becomes:
#   vesc**2 =  2*G*M / (r+a)

# Rearrange the escape speed equation for M 
#   M = vesc**2 /2 /G *(r+a)

def MassFromVesc(vesc, r, a): 
    """ This function determines the total mass needed for a given escape speed
        assuming a Hernquist profile for the dark matter halo
            M = vesc**2*(r+a) /2/G
    
    Inputs: 
        vesc: 'astropy quantity'
            The scape speed in km/s (or the speed of the satellite)
        r: 'astropy quantity'
            The distance from the Galactic Center (kpc)
        a:  'astropy quantity'
            The Hernquist scale length (kpc)
            
    Outputs:
        M : 'astropy quantity'
            Total mass within r in Msun    
    """
    
    vescKpcGyr = vesc.to(u.kpc/u.Gyr) # Converting velocity units to kpc/Gyr
    
    M = vescKpcGyr**2/2/Grav*(r+a) # Required mass
    
    return M
    
    
    


# In[24]:


VLeoI = 196*u.km/u.s # Speed of Leo I from Sohn+2013 ApJ 768
a = 30*u.kpc # Scale Radius for the Hernquist Halo 
r = 260*u.kpc # Galactocentric distance of Leo I


# In[25]:


# Compute the mass needed to keep Leo I bound! 
MLeoI = MassFromVesc(VLeoI, r, a)
print(MLeoI)


# In[27]:


print(f"{MLeoI:.2e}")


# In[28]:


MIso260/MLeoI


# In[ ]:




