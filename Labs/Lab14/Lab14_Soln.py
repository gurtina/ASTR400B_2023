
# # In Class Lab 14 Solution
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


# # PART A



# Load in disk particles centered on the MW
# this is from the HighRes data files on nimoy so it might take a bit of time to load
COM = CenterOfMass("MW_000.txt",2)
 # "MWa" was teh low res files




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

# plot the particle density for M31 
# can modify bin number (e.g. bin =100 for low res files)
plt.hist2d(MW_Disk_x,MW_Disk_y, bins=200, norm=LogNorm(), cmap='magma')
plt.colorbar(label='Number  of  stars  per  bin')


#### ADD HERE 
# call density_contour to add contours
# density_contour(x pos, y pos, contour res, contour res, axis, colors=[colors,colors])
density_contour(MW_Disk_x, MW_Disk_y, 80, 80, ax=ax,                colors=['white'])


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


# # Part B  Zooming in on a plot with widgets
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

# In[10]:



def callbackRectangle1( click, release ):
    print( f"button {click.button} pressed" )
    print( f"button {release.button} released" )
    extent = [ click.xdata, release.xdata, click.ydata, release.ydata ]
    print( f"box extent is {extent}") 
    
    # ADD - too zoom in reset the axes to the clicked box size
    ax.set_xlim( click.xdata, release.xdata ) # added
    ax.set_ylim( click.ydata, release.ydata ) # added

    # save the file as a .png
    # comment this out if your code is giving you problems
    plt.savefig("Lab14_Density_Zoom.png")


    # add the  ability to reset the image
def onKeyPressed(event):
    
    if event.key in ['R', 'r']:
        # ADD - to zoom back out reset the axes
        ax.set_xlim(-30,30 ) # added
        ax.set_ylim(-30,30) # added




# plot the particle density for the MW Disk and then zoom in on a region of the disk 

fig, ax = plt.subplots(figsize =(10 ,10))                             # create an axes

plt.hist2d(MW_Disk_x,MW_Disk_y, bins=300, norm=LogNorm(), cmap='magma')
plt.colorbar(label='Number  of  stars  per  bin')


# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.
density_contour(MW_Disk_x, MW_Disk_y, 80, 80, ax=ax,                 colors=['white'])
   

    
## NEW: Rectangle Selector.     
rs = mwidgets.RectangleSelector( ax,                        # the axes to attach to
                           callbackRectangle1,         # the callback function
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
# to detect the 'R' key, hit escape to reset the image
plt.connect("key_press_event", onKeyPressed)


# # Part C    Connecting Morphology to Kinematics
# 
# 
# Make a two panel plot.
# Left Panel:  Density 
# Right Panel: Phase Diagram 
# 
# Relect a section of the density plot and see where the particles are on the phase diagram

# In[13]:


# ADD MassProfile Object.

MWCirc = MassProfile("MW",0)




R = np.arange(0.001,30,0.5)
Vcirc = MWCirc.circularVelocityTotal(R)





def callbackRectangle2( click, release ):
    print( f"button {click.button} pressed" )
    print( f"button {release.button} released" )
    extent = [ click.xdata, release.xdata, click.ydata, release.ydata ]
    print( f"box extent is {extent}") 
    
    # No longer need to zoom  
    #ax[0].set_xlim( click.xdata, release.xdata ) # added
    #ax[0].set_ylim( click.ydata, release.ydata ) # added
    
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
    plt.savefig("Lab14_Position_Velocity.png")

    




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

#button 1 is left mouse button

# to detect the 'R' key,  press escape to reset the image
plt.connect("key_press_event", onKeyPressed)



# # Part D:  Flip it around : connect kinematics to morphology
# 
# Now Pick based on kinematics and find out where they are in the disk.
# This would be a good way to find e.g. high velocity particles. or particles that are really not obeying the normal kinematics of the disk at the current time.




def callbackRectangle3( click, release ):
    print( f"button {click.button} pressed" )
    print( f"button {release.button} released" )
    extent = [ click.xdata, release.xdata, click.ydata, release.ydata ]
    print( f"box extent is {extent}") 
    

    # xy need bottom left corner
    width = np.abs(release.xdata - click.xdata)
    height = np.abs(release.ydata - click.ydata )
    Rect = plt.Rectangle( (click.xdata,click.ydata),  width, height, fill=False, color='yellow', linewidth=3)
    # xy, width, height, angle=0.0, **kwargs
    ax[0].add_patch(Rect)
    
    # CHANGE HERE
    # make sure pick rectangle from bottom left corner upwards 
    # here switch y --> vy 
    index = np.where( (MW_Disk_x > click.xdata) & (MW_Disk_x < release.xdata)& (MW_Disk_vy > click.ydata )         & (MW_Disk_vy < release.ydata))
    # JUST NEED TO MODIFY THIS LINE
    ax[1].scatter(MW_Disk_x[index],MW_Disk_y[index])
    
    
    # save the file as a .png
    # CHANGE Filename
    plt.savefig("Lab14_Velocity_Position.png")

    




# flip the axes ax[0]<--> ax[1]


# plot the particle density for MW 
fig,ax = plt.subplots(nrows=1, ncols=2, figsize=(30,10), constrained_layout=True) 
        #   ax[0] for Velocity  # HERE
        #   ax[1] for Position    
                                                      
ax[1].hist2d(MW_Disk_x,MW_Disk_y, bins=200, norm=LogNorm(), cmap='magma')
# plt.colorbar()
# ax[0].scatter_density(MW_Disk.x, MW_Disk.z, norm=norm)
density_contour(MW_Disk_x, MW_Disk_y, 80, 80, ax=ax[1],                 colors=['red','white','white', 'white', 'white', 'white', 'white', 'white'])
#set axis limits
ax[1].set_ylim(-30,30)
ax[1].set_xlim(-30,30)  

# Add axis labels
ax[1].set_xlabel('x (kpc)', fontsize=15)
ax[1].set_ylabel('y (kpc)', fontsize=15)


# Phase Diagram : X vs. VY 

ax[0].hist2d(MW_Disk_x,MW_Disk_vy, bins=500, norm=LogNorm(), cmap='magma')
ax[0].set_xlim(-30,30)

# Add axis labels
ax[0].set_xlabel('x (kpc)', fontsize=15)
ax[0].set_ylabel('Velocity Y Component (km/s)', fontsize=15)

# Add the circular velocity
ax[0].plot(R, -Vcirc, color="blue")
ax[0].plot(-R, Vcirc, color="blue")
                        
    
rs = mwidgets.RectangleSelector( ax[0],                        # the axes to attach to
                           callbackRectangle3,         # the callback function
                           button=[1, 3],             # allow us to use left or right mouse button
                                                      #button 1 is left mouse button
                           minspanx=5, minspany=5,    # don't accept a box of fewer than 5 pixels
                           spancoords='pixels' )      # units for above



# to detect the 'R' key press to reset the image
#plt.connect("key_press_event", onKeyPressed)



# # Part E : Connecting particles across snapshots
# 



# Load in a different snapshot
#COM_2 = CenterOfMass("MW_450.txt",2)
COM_2 = CenterOfMass("MW_005.txt",2)





# Compute COM of M31 using disk particles
COMP_2 = COM_2.COM_P(0.1, 2)
# Determine positions of disk particles relative to COM 
MW_Disk_2_x = COM_2.x - COMP_2[0].value 
MW_Disk_2_y = COM_2.y - COMP_2[1].value 





def callbackRectangle4( click, release ):
    print( f"button {click.button} pressed" )
    print( f"button {release.button} released" )
    extent = [ click.xdata, release.xdata, click.ydata, release.ydata ]
    print( f"box extent is {extent}") 
    
    # ADD - too zoom in reset the axes
    #ax[0].set_xlim( click.xdata, release.xdata ) # added
    #ax[0].set_ylim( click.ydata, release.ydata ) # added
    
    # xy need bottom left corner
    width = np.abs(release.xdata - click.xdata)
    height = np.abs(click.ydata-release.ydata )
    Rect = plt.Rectangle( (click.xdata,click.ydata),  width, height, fill=False, color='yellow', linewidth=3)
    # xy, width, height, angle=0.0, **kwargs
    ax[0].add_patch(Rect)
    
    index = np.where( (MW_Disk_x > click.xdata) & (MW_Disk_x < release.xdata)& (MW_Disk_y > click.ydata )         & (MW_Disk_y < release.ydata))
    # ADJUST HERE 
    ax[1].scatter(MW_Disk_2_x[index],MW_Disk_2_y[index])

    # save the file as a .png
    # CHANGE Filename
    plt.savefig("Lab14_MW0_vs_MW5.png")

    
  



# What if we compared to another point in time


    
fig,ax = plt.subplots(nrows=1, ncols=2, figsize=(25,10))    
        #   ax[0] for Position
        #   ax[1] for Velocity 
  
                                                      
ax[0].hist2d(MW_Disk_x,MW_Disk_y, bins=200, norm=LogNorm(), cmap='magma')
# plt.colorbar()
density_contour(MW_Disk_x, MW_Disk_y, 80, 80, ax=ax[0],                 colors=['white'])
#set axis limits
ax[0].set_ylim(-30,30)
ax[0].set_xlim(-30,30)  

# Add axis labels
ax[0].set_xlabel('x (kpc)', fontsize=15)
ax[0].set_ylabel('y (kpc)', fontsize=15)


# SNAPSHOT 400  MODIFIED HERE 

ax[1].hist2d(MW_Disk_2_x, MW_Disk_2_y, bins=200, norm=LogNorm(), cmap='magma')
density_contour(MW_Disk_2_x, MW_Disk_2_y, 80, 80, ax=ax[1],                 colors=['white'])

# Set axis limits
ax[1].set_xlim(-30,30)
ax[1].set_ylim(-30,30)

#ax[1].set_xlim(-80,80)
#ax[1].set_ylim(-80,80)


# Add axis labels
ax[1].set_xlabel('x (kpc)', fontsize=15)
ax[1].set_ylabel('y (kpc)', fontsize=15)

                        
    
rs = mwidgets.RectangleSelector( ax[0],                        # the axes to attach to
                           callbackRectangle4,         # the callback function
                           button=[1, 3],             # allow us to use left or right mouse button
                           minspanx=5, minspany=5,    # don't accept a box of fewer than 5 pixels
                           spancoords='pixels' )      # units for above

#button 1 is left mouse button

# to detect the 'R' key press to reset the image
#plt.connect("key_press_event", onKeyPressed)






