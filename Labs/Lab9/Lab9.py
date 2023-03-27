
# # In Class Lab 9
# 
# Tutorial to make some interesting plots with widgets and the simulaton data ! 
# 
# Note that the images in this lab will look a lot better if you use the high res files. 

# Graphical widgets -- helpful functions to make a "graphical user interface", or GUI.
# https://matplotlib.org/stable/api/widgets_api.html
# 
# These widgets need to be able to take input from the mouse and keyboard while the program is running. The most common way this is achieved is to have the code run in an infinite loop which is interrupted whenever input is provided. Some action is taken according to the input, and then the loop starts running again. This sort of algorithm is known as an *event loop* -- the code loops until a user event occurs.
# 
# `matplotlib` provides a number of simple widgets which automatically create an event loop for us. One can create a widget instance, and then tell the widget what function to run when something happens to the widget. Such a function is called a *callback* -- the event loop calls back to the function we give it in order to take some action before starting up again.
# 


import matplotlib.widgets as mwidgets  # get access to the widgets


# external modules
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import numpy as np

# my modules 
from ReadFile import Read
from CenterOfMass2 import CenterOfMass
from MassProfile import MassProfile

# I took the code from Lab 7 for making contours and made it into a separate script
# NOTE: it is more organized and easier to debug if you keep functions in separate scripts 
# and then call them when you want to e.g. make plots or do some analysis. 
from contour import density_contour


# # PART A : M-$\sigma$ Relation
# 
# - Compute the velocity dispersion of the bulge in the Simulated M31.
# 
# - create a function that defines the  M-$\sigma$ relation, from  Kormendy & Ho 2013 (eqn 7)
# 
# $\left (\frac{M_{BH}}{10^9 M_\odot} \right ) = \alpha \left( \frac{\sigma}{\sigma_0}\right)^\beta $  
# 
# $\alpha$ =  0.309, 
# $\beta$ = 4.38, 
# $\sigma_0$ = 200 km/s
# 
# Comparing to the Simulation:
# 
# - What mass black hole would you expect from the data? How does this compare to the mass of the SMBH in the real M31 (5-23$\times 10^7$ M$_\odot$) 



# Create a COM of object for MW Using Code from Assignment 4
# Here we are going to use bulge particles explicitly since we care about 
# the dynamics of the bulge




# Compute COM of M31 bulge




# Determine positions of bulge particles relative to COM 


# Determine the total magnitude of position vector from COM


# Determine the velocities of bulge particles relative to COM motion




# average velocity of bulge (within 5 kpc) 




# Standard Deviation. ( 1 sigma)




# Define Msigma 


# # PART B:  Widgets

# # PART 1




# Load in disk particles centered on the MW
# this is from the HighRes data files on nimoy so it might take a bit of time to load
COM = CenterOfMass("MW_000.txt",2)




# Compute COM of the MW at the new position using disk particles
COMP = COM.COM_P(0.1, 2)
COMV = COM.COM_V(COMP[0],COMP[1],COMP[2])
# Determine positions of disk particles relative to COM 
MW_Disk_x = COM.x - COMP[0].value 
MW_Disk_y = COM.y - COMP[1].value 

# Also store the disk velocity in the y direction
MW_Disk_vy = COM.vy - COMV[1].value




# Plot the disk of the MW with contours. 


# MW Disk Density 
fig, ax= plt.subplots(figsize=(10, 10))

## ADD HERE
# plot the particle density for MW using plt.hist2d 
# plt.hist2d(pos1,pos2, bins=, norm=LogNorm(), cmap= )
# cmap options: https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html  
#   e.g. 'magma', 'viridis'




#### ADD HERE 
# call density_contour to add contours
# density_contour(x pos, y pos, contour res, contour res, axis, colors=[colors,colors])





# Add axis labels
plt.xlabel('x (kpc)', fontsize=22)
plt.ylabel('y (kpc)', fontsize=22)

#set axis limits
plt.ylim(-30,30)
plt.xlim(-30,30)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size


# # Part 2  Zooming in on a plot with widgets
# 
# 
# We can catch characters typed on the keyboard -- *keypress events* -- by connecting a "key_press_event" to a callback function which takes an event as an argument.
# The event object contains a variety of data. The most useful being:
# 
#     event.key       # the key which was pressed
#     event.xdata     # the mouse x-position when the key was pressed
#     event.ydata     # the mouse y-position when the key was pressed
#     
# Another useful widget allows the user to select a rectangular region in some axes object, and then calls a callback function with the bounding coordinates (the extent) of the region selected. This is the RectangleSelector widget.
# 
# Note that click and release are not really that! Click contains the more-negative values and release the more positive values of both x and y coordinates.




def callbackRectangle1( click, release ):
    print( f"button {click.button} pressed" )
    print( f"button {release.button} released" )
    extent = [ click.xdata, release.xdata, click.ydata, release.ydata ]
    print( f"box extent is {extent}") 
    
    # ADD - to zoom in, reset the axes to the clicked box size


    # save the file as a .png
    # comment this out if your code is giving you problems





# add the  ability to reset the image
def onKeyPressed(event):

if event.key in ['R', 'r']:
    # ADD - to zoom back out reset the axes




# plot the particle density for the MW Disk and then zoom in on a region of the disk 

fig, ax = plt.subplots(figsize =(10 ,10))                             # create an axes

plt.hist2d(MW_Disk_x,MW_Disk_y, bins=300, norm=LogNorm(), cmap='magma')
plt.colorbar(label='Number  of  stars  per  bin')


# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.
density_contour(MW_Disk_x, MW_Disk_y, 80, 80, ax=ax,                 colors=['white'])
   

    
## NEW: Rectangle Selector.     
rs = mwidgets.RectangleSelector1( ax,                        # the axes to attach to
                           callbackRectangle,         # the callback function
                           button=[1, 3],             # allow us to use left or right mouse button
                                                      #button 1 is left mouse button
                           minspanx=5, minspany=5,    # don't accept a box of fewer than 5 pixels
                           spancoords='pixels' )      # units for above



#set axis limits
ax.set_ylim(-30,30)
ax.set_xlim(-30,30)


# Add axis labels
plt.xlabel('x (kpc)', fontsize=22)
plt.ylabel('y (kpc)', fontsize=22)

# ADDED THIS
# Press 'R' key to reset AND THEN
# to detect the 'R' key, click escape to reset the image
plt.connect("key_press_event", onKeyPressed)


# # Part C    Connecting Morphology to Kinematics
# 
# 
# Make a two panel plot.
# Left Panel:  Density 
# Right Panel: Phase Diagram 
# 
# Relect a section of the density plot and see where the particles are on the phase diagram



# ADD MassProfile Object.





# Add an array for radii up to 40 kpc


# Store Vcirc 





def callbackRectangle2( click, release ):
    print( f"button {click.button} pressed" )
    print( f"button {release.button} released" )
    extent = [ click.xdata, release.xdata, click.ydata, release.ydata ]
    print( f"box extent is {extent}") 
    
    # Add a yellow rectangle to where we selected a region rather than zooming in
    # xy need bottom left corner
    width = np.abs(release.xdata - click.xdata)
    height = np.abs(click.ydata-release.ydata )
    
    # Create a rectangle 
    Rect = plt.Rectangle( (click.xdata,click.ydata),  width, height, fill=False, color='yellow', linewidth=3)
    # xy, width, height, angle=0.0, **kwargs
    ax[0].add_patch(Rect) # add the rectangle to the axis object.
    
    
    # Plotting the corresponding points on the left panel
    # make sure pick rectangle from bottom left corner upwards 
    index = np.where( (MW_Disk_x > click.xdata) & (MW_Disk_x < release.xdata)& (MW_Disk_y > click.ydata )         & (MW_Disk_y < release.ydata))
    ax[1].scatter(MW_Disk_x[index],MW_Disk_vy[index])
    
    # save the file as a .png
    # comment this out if your code is giving you problems
    plt.savefig("Lab9_Position_Velocity.png")

    




# Now also want to see the corrsponding phase diagram for that region


fig,ax = plt.subplots(nrows=1, ncols=2, figsize=(30,10), constrained_layout=True)    
        #   ax[0] for Position
        #   ax[1] for Velocity 

                                                  
ax[0].hist2d(MW_Disk_x,MW_Disk_y, bins=200, norm=LogNorm(), cmap='magma')

density_contour(MW_Disk_x, MW_Disk_y, 80, 80, ax=ax[0],                 colors=['red','white','white', 'white', 'white', 'white', 'white', 'white'])
#set axis limits
ax[0].set_ylim(-30,30)
ax[0].set_xlim(-30,30)  

# Add axis labels
ax[0].set_xlabel('x (kpc)', fontsize=15)
ax[0].set_ylabel('y (kpc)', fontsize=15)


# Phase Diagram : X vs. VY 

ax[1].hist2d(MW_Disk_x,MW_Disk_vy, bins=500, norm=LogNorm(), cmap='magma')
ax[1].set_xlim(-30,30)

# Add axis labels
ax[1].set_xlabel('x (kpc)', fontsize=15)
ax[1].set_ylabel('Velocity Y Component (km/s)', fontsize=15)

# Add the circular velocity
ax[1].plot(R, -Vcirc, color="blue")
ax[1].plot(-R, Vcirc, color="blue")
                        
    
rs = mwidgets.RectangleSelector( ax[0],                        # the axes to attach to
                           callbackRectangle2,         # the callback function
                           button=[1, 3],             # allow us to use left or right mouse button
                           minspanx=5, minspany=5,    # don't accept a box of fewer than 5 pixels
                           spancoords='pixels' )      # units for above


# # Part D:  Flip it around : connect kinematics to morphology
# 
# Now Pick based on kinematics and find out where they are in the disk.
# This would be a good way to find e.g. high velocity particles. or particles that are really not obeying the normal kinematics of the disk at the current time.



# Copy over the Call back function and the onkeypressed function from Part C
# flip the axes ax[0] < --- > ax[1]









# Copy over the Density and phase diagram code
# flip the axes ax[0]<--> ax[1]







# # Part E : Connecting particles across snapshots
# 



# Load in a different Snapshot






# Compute COM of M31 using disk particles

# Determine positions of disk particles relative to COM 




# Copy over the Call back function from Part C
# Edit so that it overplots a scatter of particles from the new snapshot




# Copy over the plotting script from Part C
# Instead of the phase plot, have the second panel be the MW at a different snapshot

