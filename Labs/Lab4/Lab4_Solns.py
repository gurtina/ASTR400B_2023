

#In Class Lab 4 Solutions
# G Besla ASTR 400B

# Load Modules
import numpy as np
import astropy.units as u

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib
get_ipython().run_line_magic('matplotlib', 'inline')


# The Figure illustrates the color magnitude diagram (CMD) for the Carina Dwarf along with the interpreted 
# star formation history from isochrone fitting to the CMD.
# The image is from Tolstoy+2009 ARA&A 47 review paper about dwarf galaxies
# 
# NOTE: 
# The right hand plot x axis is Gyr SINCE the BIG BANG
# The data file naming convention lists the # of Gyr ago from today. 
# 
# ![Iso](./Lab4_Isochrones.png)
# 

# # This Lab:
# 
# Modify the template file of your choice to plot isochrones that correspond to the inferred star formation episodes (right panel of Figure 1) to recreate the CMD of Carina (left panel of Figure 1). 




# Some Notes about the Isochrone Data
# DATA From   http://stellar.dartmouth.edu/models/isolf_new.html
# files have been modified from download.  ( M/Mo --> M;   Log L/Lo --> L)
# removed #'s from all lines except column heading
# NOTE SETTINGS USED:  
# . Y = 0.245 default   [Fe/H] = -2.0  alpha/Fe = -0.2
# Colors:  UBVRIc + 2MASS + Kepler
# These could all be changed and it would generate a different isochrone





# Filename for data with Isochrone fit for 1 Gyr
# These files are located in the folder IsochroneData
filename1="./IsochroneData/Isochrone1.txt"





# READ IN DATA
# "dtype=None" means line is split using white spaces
# "skip_header=8"  skipping the first 8 lines 
# the flag "names=True" creates arrays to store the date
#       with the column headers given in line 8 

# Read in data for an isochrone corresponding to 1 Gyr
data1 = np.genfromtxt(filename1,dtype=None,names=True,skip_header=8)





## Solns: Step 1
# The major peak has ages of about 2-3 Gyr after the big bang --  11 to 10 Gyr ago from today
filename11= "./IsochroneData/Isochrone11.txt"
filename10= "./IsochroneData/Isochrone10.txt"

# Read in data for isochrones corresponding to 11 -10 Gyr ago from today
# This is the Dominant Peak
data11 = np.genfromtxt(filename11,dtype=None,names=True,skip_header=8)
data10 = np.genfromtxt(filename10,dtype=None,names=True,skip_header=8)





## Solns: Step 2
# The next peak is about 7-8 Gyr after the big bang - about 6-7 Gyr ago from today
filename6= "./IsochroneData/Isochrone6.txt"
filename7= "./IsochroneData/Isochrone7.txt"

# Read in data for isochrones corresponding to 6-7 Gyr ago from today
# This is the 2nd Dominant Peak
data6 = np.genfromtxt(filename6,dtype=None,names=True,skip_header=8)
data7 = np.genfromtxt(filename7,dtype=None,names=True,skip_header=8)





## Solns: Step 3

# The last peak is old, about 13.5 Gyr ago 
filename13_5= "./IsochroneData/Isochrone13_5.txt"

# Read in data for isochrones corresponding to 13.5 Gyr ago
# This is the last peak
data13_5 = np.genfromtxt(filename13_5,dtype=None,names=True,skip_header=8)





## Q2 Solns: Step 4 Try younger 

# The next peak 
#  Try 3 Gyr 
filename3= "./IsochroneData/Isochrone3.txt"
data3 = np.genfromtxt(filename3,dtype=None,names=True,skip_header=8)
# Try 2 Gyr
filename2= "./IsochroneData/Isochrone2.txt"
data2 = np.genfromtxt(filename2,dtype=None,names=True,skip_header=8)




# Plot Isochrones 
# For Carina

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot Isochrones

# Isochrone for 1 Gyr
#plt.plot(data1['B']-data1['R'], data1['R'], color='blue', linewidth=5, label='1 Gyr')
###EDIT Here, following the same format as the line above 


# Most Dominant Peak
# Isochrone for 11 Gyr
plt.plot(data11['B']-data11['R'], data11['R'], color='blue', linewidth=5, label='11 Gyr')
# Isochrone for 10 Gyr
plt.plot(data10['B']-data10['R'], data10['R'], color='blue',linestyle=":", linewidth=5, label='10 Gyr')

# 2nd Most Dominant Peak (there is a divide at 8 Gyr)
# Isochrone for 6 Gyr
plt.plot(data6['B']-data6['R'], data6['R'], color='magenta', linewidth=5, label='6 Gyr')
# Isochrone for 8 Gyr
plt.plot(data7['B']-data7['R'], data7['R'], color='magenta', linestyle=":", linewidth=5, label='7 Gyr')

# 3rd Peak
# Isochrone for 13.5 Gyr
plt.plot(data13_5['B']-data13_5['R'], data13_5['R'], color='red', linewidth=5, label='13.5 Gyr')

# Q2 : The Tolstoy Paper is a bit out of date. 
#  there might be even younger ages, between 2-8 Gyr ago. 
# Isochrone for 3 Gyr and 2 Gyr ago 
plt.plot(data3['B']-data3['R'], data3['R'], color='orange', linewidth=5, label='3 Gyr')
plt.plot(data2['B']-data2['R'], data2['R'], color='orange', linestyle=":", linewidth=5, label='2 Gyr')



# Add axis labels
plt.xlabel('B-R', fontsize=22)
plt.ylabel('M$_R$', fontsize=22)

#set axis limits
plt.xlim(-0.5,2)
plt.ylim(5,-2.5)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

#add figure text
plt.figtext(0.6, 0.15, 'Isochrone Carina', fontsize=22)


# Save to a file
ax.set_rasterized(True)
plt.savefig('IsochroneLabCarina.png') 


# # Q2
# 
# Could there be younger ages than suggested in the Tolstoy plot?
# Try adding younger isochrones to the above plot.
# 
# # Q3
# 
# What do you think might cause the bursts of star formation?
# 

# Star bursts could be induced by close passages between Carina and other galaxies, like the Milky Way or other satellites like the LMC.  
# The tidal field from those close encounters could funnel gas into the center of Carina, increasing the local density and therefore the star formation rate.





