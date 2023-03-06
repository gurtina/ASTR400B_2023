
# In Class Lab 7 Solutions
# G. Besla 
# with code from R. Hoffman and E Patel


# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile


# # This Lab:
# 
# A) Using Contour plots to identify substructure within the stellar disk of M31.
# 
# B) Rotating the disk so that we are seeing it edge on
# 
# C) Create plots to examine the kinematics of the disk



# Code for plotting contours
#from https://gist.github.com/adrn/3993992

import scipy.optimize as so

def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

# Info about **kwargs, *args 
#https://book.pythontips.com/en/latest/args_and_kwargs.html

def density_contour(xdata, ydata, nbins_x, nbins_y, ax=None, **contour_kwargs):
    """ Create a density contour plot.
    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
        i.e. unknown number of keywords 
        
    Example Usage
    -------------
     density_contour(x pos, y pos, contour res, contour res, axis, colors for contours)
     e.g.:
     density_contour(xD, yD, 80, 80, ax=ax, colors=['red','orange', 'yellow', 'orange', 'yellow'])

    """

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))  
    
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T  # transpose of distribution fxn
    fmt = {}
    
    ### Adjust Here #### 
    
    # Contour Levels Definitions
    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68)) 
    #brentq is root finding method
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    
    # You might need to add a few levels
    # I added a few between 1 and 2 sigma to better highlight the spiral arm
    one_sigma1 = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.80))
    one_sigma2 = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.90))

    # Array of Contour levels. Adjust according to the above
    levels = [one_sigma, one_sigma1, one_sigma2, two_sigma, three_sigma][::-1]
    
    # contour level labels  Adjust accoding to the above.
    strs = ['0.68', '0.8','0.9','0.95', '0.99'][::-1]

    
    ###### 
    
    if ax == None:
        contour = plt.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        plt.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)

    else:
        contour = ax.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        ax.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)
    
    return contour


# Use the CenterOfMass code to compute the positions and velocities of all particles in M31's disk relative to its center of mass position and motion.



# Create a COM of object for M31 Disk Using Code from Assignment 4
COMD = CenterOfMass("M31_000.txt",2)




# Compute COM of M31 using disk particles
COMP = COMD.COM_P(0.1)
COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])




# Determine positions of disk particles relative to COM 
xD = COMD.x - COMP[0].value 
yD = COMD.y - COMP[1].value 
zD = COMD.z - COMP[2].value 

# total magnitude
rtot = np.sqrt(xD**2 + yD**2 + zD**2)

# Determine velocities of disk particles relatiev to COM motion
vxD = COMD.vx - COMV[0].value 
vyD = COMD.vy - COMV[1].value 
vzD = COMD.vz - COMV[2].value 

# total velocity 
vtot = np.sqrt(vxD**2 + vyD**2 + vzD**2)

# Vectors for r and v 
r = np.array([xD,yD,zD]).T # transposed for the Rotate Function later
v = np.array([vxD,vyD,vzD]).T


# # Part A:
# 
# Create plot of M31's disk density, using 2D Histograms 



# 1) Make plots 

# M31 Disk Density 
fig, ax= plt.subplots(figsize=(10, 10))

# plot the particle density for M31 using a 2D historgram
# can modify bin number to make the plot smoother
plt.hist2d(xD, yD, bins=150, norm=LogNorm(), cmap='magma')
plt.colorbar()

# Add axis labels
plt.xlabel('x (kpc)', fontsize=22)
plt.ylabel('y (kpc)', fontsize=22)

#set axis limits
plt.ylim(-40,40)
plt.xlim(-40,40)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.
density_contour(xD, yD, 80, 80, ax=ax, colors=['red','orange', 'yellow', 'orange', 'yellow'])


# Save to a file
plt.savefig('Lab7_M31Disk.png')


# In this projection, if the spiral arms are trailing it looks like M31 is rotating CW.
# But we don't know the inclination of the galaxy! This can change the interpretation.

# # Part B
# 
# Utilize the below code to rotate the M31 disk and plot it edge on and face on.
# 
# What is the sense of rotation of M31 ? 




def RotateFrame(posI,velI):
    """a function that will rotate the position and velocity vectors
    so that the disk angular momentum is aligned with z axis. 
    
    PARAMETERS
    ----------
        posI : `array of floats`
             3D array of positions (x,y,z)
        velI : `array of floats`
             3D array of velocities (vx,vy,vz)
             
    RETURNS
    -------
        pos: `array of floats`
            rotated 3D array of positions (x,y,z) such that disk is in the XY plane
        vel: `array of floats`
            rotated 3D array of velocities (vx,vy,vz) such that disk angular momentum vector
            is in the +z direction 
    """
    # compute the angular momentum
    L = np.sum(np.cross(posI,velI), axis=0)
    # normalize the vector
    L_norm = L/np.sqrt(np.sum(L**2))


    # Set up rotation matrix to map L_norm to z unit vector (disk in xy-plane)
    
    # z unit vector
    z_norm = np.array([0, 0, 1])
    
    # cross product between L and z
    vv = np.cross(L_norm, z_norm)
    s = np.sqrt(np.sum(vv**2))
    
    # dot product between L and z 
    c = np.dot(L_norm, z_norm)
    
    # rotation matrix
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
    R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2

    # Rotate coordinate system
    pos = np.dot(R, posI.T).T
    vel = np.dot(R, velI.T).T
    
    return pos, vel





# determine the rotated velocity vectors
rn, vn = RotateFrame(r,v)




# Rotated M31 Disk - EDGE ON

# M31 Disk Density 
fig, ax= plt.subplots(figsize=(10, 10))

# plot the particle density for M31 
# can modify bin number (bin =100 smoothest)
plt.hist2d(rn[:,0], rn[:,2], bins=150, norm=LogNorm(), cmap='magma')
plt.colorbar()

# Add axis labels
plt.xlabel('x (kpc)', fontsize=22)
plt.ylabel('y (kpc)', fontsize=22)

#set axis limits
plt.ylim(-40,40)
plt.xlim(-40,40)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.
#density_contour(rn[:,0], rn[:,2], 80, 80, ax=ax, 
#          colors=['red','orange', 'yellow', 'orange', 'yellow'])

plt.savefig('Lab7_EdgeOn_Density.png')




# Rotated M31 Disk - FACE ON

# M31 Disk Density 
fig, ax= plt.subplots(figsize=(10, 10))

# plot the particle density for M31 
# can modify bin number (bin =100 smoothest)
plt.hist2d(rn[:,0], rn[:,1], bins=150, norm=LogNorm(), cmap='magma')
plt.colorbar()

# Add axis labels
plt.xlabel('x (kpc)', fontsize=22)
plt.ylabel('y (kpc)', fontsize=22)

#set axis limits
plt.ylim(-40,40)
plt.xlim(-40,40)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.
density_contour(rn[:,0], rn[:,1], 80, 80, ax=ax, 
                colors=['red','white', 'white', 'white', 'white'])

plt.savefig('Lab7_FaceOn_Density.png')


# If the arms are trailing, then the disk should be rotating towards  NEGATIVE Y  along NEGATIVE X  (CCW)

# # Part C
# 
# a) Create a scatter plot of the edge on disk particles, weighted by velocity.
# 
# b) Create a phase diagram for the edge on disk (Position vs. Velocity) using a 2D Histogram.



# Plot velocity weighted EDGE ON DISK

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)


# plot position of disk particles color coded by velocity along the 3rd axis
plt.scatter(rn[:,0],rn[:,2],c=vn[:,1])

#colorbar
cbar = plt.colorbar()
cbar.set_label(r'Vy (km/s)', size=22)

# Add axis labels
plt.xlabel('x (kpc)', fontsize=22)
plt.ylabel('z (kpc)', fontsize=22)


#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

#set axis limits
plt.ylim(-40,40)
plt.xlim(-40,40)

# Save file
plt.savefig('Lab7_EdgeOn_Vel.png')


# Indeed, along NEGATIVE X the disk is rotating into the page, which is in the direction of NEGATIVE Y. So M31 is rotating Counter Clock Wise



# a)  Create a mass profile object for the Milky Way using solution to Homework 5
M31 = MassProfile("M31",0)




R = np.arange(0.01,40,0.1)
Vcirc = M31.CircularVelocityTotal(R)




# Make a phase diagram
# MW Disk Velocity Field edge on.

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# looking at MW edge on along x axis, vy is line of sight velocity

plt.hist2d(rn[:,0], vn[:,1], bins=150, norm=LogNorm())
plt.colorbar()

# Add the circular velocity
plt.plot(R, Vcirc, color="red")
plt.plot(-R, -Vcirc, color="red")

# Add axis labels
plt.xlabel('x (kpc)', fontsize=22)
plt.ylabel('Velocity Y Direction (km/s)', fontsize=22)

#set axis limits
#plt.ylim(0,200)
#plt.xlim(-30,30)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size


# Save file
plt.savefig('Lab7_PhaseDiagram.png')






