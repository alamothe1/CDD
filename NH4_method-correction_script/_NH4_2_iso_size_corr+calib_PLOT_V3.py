# -*- coding: cp1252 -*-

### IMPORT MODULES ##############################################################################
# We here import commonly used modules
from __future__ import division     # true division
import numpy as np                  # array module
try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3                                 
from pylab import *
from time import *                  # module to measure time

from scipy import stats
import matplotlib.pyplot as plt

from matplotlib.font_manager import FontProperties
#################################################################################################




### Main inputs to the calculation
input_standards_filename = 'input_measured_iso_stds.txt'             # Name of the file with all the isotopic data of the standards
input_samples_filename   = 'input_measured_samples.txt'              # Name of the file with all the isotopic data of the samples
input_blanks_filename    = 'input_measured_blanks.txt'               # Name of the file with all the isotopic data of the blanks
# input_P_stds_filename    = 'input_measured_iso_stds_P_only.txt'      # Name of the file with all the isotopic data of the P standards only


### Import functions
import _ALL_functions_and_ref_material_data as f


### indicate which species is measured
Measured_species = 'NH\u2084\u207A'


### Import the size correction and calibration parameters
Abk, d15Nbk, d15Naz, slope, inter = f.read_file_size_corr_calib_params_d15NNH4('output_params_for_corr_NH4_Nitrogen.txt')
# and print the values
print ('Optimized blank size : ' + ' '*(5-len(str(round(Abk, 2))))    + str(round(Abk, 2))    + ' Vs')
print ('Optimized blank d15N : ' + ' '*(5-len(str(round(d15Nbk, 1)))) + str(round(d15Nbk, 1)) + ' per mill')
print ('Optimized azide d15N : ' + ' '*(5-len(str(round(d15Naz, 1)))) + str(round(d15Naz, 1)) + ' per mill')
print ('Slope                : ' + ' '*(5-len(str(round(slope , 3))))      + str(round(slope, 3)))
print ('Intercept            : ' + ' '*(5-len(str(round(inter , 3))))      + str(round(inter, 3)))
print ('---')


 
#################################################################################################
### Read inputs          ########################################################################
#################################################################################################
### Imports data for ALL ISOTOPIC STANDARDS
# Create object
IS_all = f.Isotopic_data()              # "IS" means isotopic standards
IS_all.read_file_without_header(input_standards_filename, 1)        # reads the file

# Convert data to arrays
# Only N for now
ID      = IS_all.ID
amount  = array(IS_all.amount)
Amea    = array(IS_all.m28)
d15Nmea = array(IS_all.d15Nmea)
d15Nstd = array(IS_all.d15Nstd)

# Calculate calibrated values, and size corrected + calibrated values and residues for the standards
d15Ncal  = slope*d15Nmea + inter
d15Nscc  = f.mixing_model_d15NNH4_REV(Amea, d15Ncal, Abk, d15Nbk, d15Naz)
d15Nres  = d15Nscc - d15Nstd
print ('d15Nres : - mean  = ' + str(round(mean(d15Nres), 3)))
print ('          - stdev = ' + str(round(std(d15Nres), 3)))

### Imports sample data
# Create object
samples = f.Isotopic_data()
samples.read_file_without_header(input_samples_filename, 1)        # reads the file

# Convert data to arrays
# Only N for now
samples.Amea    = array(samples.m28)
samples.d15Nmea = array(samples.d15Nmea)

### Imports blank data
# Create object
blanks = f.Isotopic_data()
blanks.read_file_without_header(input_blanks_filename, 1)        # reads the file

# Convert data to arrays
# Only N for now
blanks.Amea    = array(blanks.m28)
blanks.d15Nmea = array(blanks.d15Nmea)

##If P stds are shown :
##### Imports P stds data
### Create object
##Ponly = f.Isotopic_data()
##Ponly.read_file_without_header(input_P_stds_filename, 1)        # reads the file
##
### Convert data to arrays
### Only N for now
##Ponly.Amea    = array(Ponly.m28)
##Ponly.d15Nmea = array(Ponly.d15Nmea)
#################################################################################################


#################################################################################################
### Plot vs size measured             ###########################################################
#################################################################################################
Lstds  = []        # Create list to enumerate available stds for calibration (just returns the letters)
for i in range(len(ID)):
    # Append Lstds if ID is not there
    if ID[i] not in Lstds:
        Lstds.append(ID[i])



name0 = 'plot1_versus_Amea.png'
pad = 15

# fig = figure(figsize=(7, 10))        # first figure is the width, second the height
# subplots_adjust(hspace=0.35)
# subplots_adjust(wspace=0.4)

fig, axs = plt.subplots(3, sharex=True, figsize = (7,10))
# fig.suptitle('Data ICMS', fontsize = 20)
plt.subplots_adjust(hspace=.0) #(no blank between plots)


###### First PANEL
####################
ax1 = plt.subplot2grid((9,1), (0,0), rowspan=3)
### Plot data
# All NH4+ standards
for i in range(len(amount)):
    ax1.plot(Amea[i], amount[i], f.style[ID[i]], color=f.color[ID[i]], ms=f.ms[ID[i]])
# Measured blanks
for i in range(len(blanks.Amea)):
    ax1.plot(blanks.Amea[i], 0, f.style['blk'], color=f.color['blk'], ms=f.ms['blk'])
    
### Fake plots to add lines in legend
# Computed blank
ax1.plot([], [], f.style['comp_blk'], color=f.color['comp_blk'], ms=f.ms['comp_blk'], label='Blank, modeled')
# Standards
for element in Lstds:
    ax1.plot([], [], f.style[element], color=f.color[element], ms=f.ms[element], label=element + ', ' + f.species[element])
# Samples
ax1.plot([], [], f.style['spl'], color=f.color['spl'], ms=f.ms['spl'], label = 'Samples')
# Measured blank
ax1.plot([], [], f.style['blk'], color=f.color['blk'], ms=f.ms['blk'], label = 'Blank, measured')

### Ylabel
Measured_species = Measured_species[0:2] + Measured_species[2] + Measured_species[3]
ylabel(Measured_species + ' amount /nmol', fontsize=14, ha="center", va="center", labelpad=20)

### Legend
ax1.legend(numpoints=1, loc=2, prop=FontProperties(size='small'))

### Axis parameters
axis([0, int(max(max(Amea), max(samples.Amea))/10+1)*10, 0, int(max(amount)/10+1)*10])
tick_params(axis='both', pad=10)



###### Second PANEL
####################
ax2 = subplot2grid((9,1), (3,0), rowspan=4)
### Plot data
# All standards
for i in range(len(amount)):
    x = Amea[i]
    y = d15Nmea[i]
    ax2.plot(x, y, f.style[ID[i]], color=f.color[ID[i]], ms=f.ms[ID[i]], label=ID[i]+', ' + f.species[ID[i]])
# Position of computed blank
x = 2*Abk
mea_val = d15Nbk
y = f.calibration_d15NNH4_REV(f.mixing_model_d15NNH4(x, mea_val, Abk, d15Nbk, d15Naz), slope, inter)
ax2.plot(x, y, f.style['comp_blk'], color=f.color['comp_blk'], ms=f.ms['comp_blk'])
print(x,y)

# Measured blanks
x = blanks.Amea
y = blanks.d15Nmea
ax2.plot(x, y, f.style['blk'], color=f.color['blk'], ms=f.ms['blk'])
# All samples
ax2.plot(samples.Amea, samples.d15Nmea, f.style['spl'], color=f.color['spl'], ms=f.ms['spl'])
##If P stds are shown :
### All P stds
##ax2.plot(Ponly.Amea, Ponly.d15Nmea, f.style['P'], color=f.color['P'], ms=f.ms['P'])

### Plot model grey lines
win1 = 0.92
win2 = 0.97
LAmea0 = arange(0.1, int(max(max(Amea), max(samples.Amea))/10+1)*10, 0.1)
LAmea1 = arange(0.1, win1*int(max(max(Amea), max(samples.Amea))/10+1)*10, 0.1)
LAmea2 = arange(win2*int(max(max(Amea), max(samples.Amea))/10+1)*10, 1*int(max(max(Amea), max(samples.Amea))/10+1)*10, 0.1)
for val in arange(-200, 200, 10):
    # If a specific line (in the following list)
    if val in [-40, -20, 0, 20, 40, 60]:
        # Plot the two portions of the line
        x = LAmea1
        mea_val = val
        y = f.calibration_d15NNH4_REV(f.mixing_model_d15NNH4(x, mea_val, Abk, d15Nbk, d15Naz), slope, inter)
        ax2.plot(x, y, '-', color='grey', lw=0.5)
        x = LAmea2
        mea_val = val
        y = f.calibration_d15NNH4_REV(f.mixing_model_d15NNH4(x, mea_val, Abk, d15Nbk, d15Naz), slope, inter)
        ax2.plot(x, y, '-', color='grey', lw=0.5)
        # Write the ID of the line (the true value of d15NNH4+
        x1 = (win1 + win2)/2
        x2 = x1 * int(max(max(Amea), max(samples.Amea))/10+1)*10
        mea_val = val
        y = f.calibration_d15NNH4_REV(f.mixing_model_d15NNH4(x2, mea_val, Abk, d15Nbk, d15Naz), slope, inter)
        ax2.text(x2, y, val, ha="center", va="center", size=10, color='grey')
    else:
        # If not a specific line
        # Write it in one portion
        x = LAmea0
        mea_val = val
        y = f.calibration_d15NNH4_REV(f.mixing_model_d15NNH4(x, mea_val, Abk, d15Nbk, d15Naz), slope, inter)
        ax2.plot(x, y, '-', color='grey', lw=0.5)
# Plot model line for each available std
LAmea = arange(0.1, int(max(max(Amea), max(samples.Amea))/10+1)*10, 0.1)
for element in Lstds:
    x = LAmea
    mea_val = f.val.d15Nstd[element]
    y = f.calibration_d15NNH4_REV(f.mixing_model_d15NNH4(x, mea_val, Abk, d15Nbk, d15Naz), slope, inter)
    ax2.plot(x, y, '--', color=f.color[element], lw=1.5)

### Ylabel
ylabel(u'\u03B4\u00B9\u2075N(' + Measured_species + ') meas.', fontsize=14, ha="center", va="center", labelpad=20)

### Axis parameters
axis([0, int(max(max(Amea), max(samples.Amea))/10+1)*10, int(min(d15Nmea)/10-2)*10, int(max(d15Nmea)/10+2)*10])
##If P stds are shown :
##axis([0, int(max(max(Amea), max(samples.Amea))/10+1)*10, int(min(min(d15Nmea), min(Ponly.d15Nmea))/10-2)*10, int(max(d15Nmea)/10+2)*10])
tick_params(axis='both', pad=10)


###### Third PANEL
####################
ax3 = subplot2grid((9,1), (7,0), rowspan=2)
### Plot the residues
for i in range(len(amount)):
    ax3.plot(Amea[i], d15Nres[i], f.style[ID[i]], color=f.color[ID[i]], ms=f.ms[ID[i]])

### Plot mean + stdev values
# Mean
ax3.plot([-1000, 1000], [mean(d15Nres), mean(d15Nres)], 'k-', lw=1.5, label='mean')
# +/- std value
ax3.plot([-1000, 1000], [mean(d15Nres)+std(d15Nres), mean(d15Nres)+std(d15Nres)], 'k--', lw=1.5, label='mean 1 sigma')
ax3.plot([-1000, 1000], [mean(d15Nres)-std(d15Nres), mean(d15Nres)-std(d15Nres)], 'k--', lw=1.5, label='')

### Ylabel
ylabel(u'\u03B4\u00B9\u2075N(' + Measured_species + ') res.', fontsize=14, ha="center", va="center", labelpad=20)
axis([0, int(max(max(Amea), max(samples.Amea))/10+1)*10, (min(d15Nres)/0.5 - 1)*0.5, (max(d15Nres)/0.5 + 1)*0.5])

### Axis parameters
tick_params(axis='both', pad=10)

### Xlabel
xlabel('Amea / V.s (measured N2 peak size)', fontsize=14, ha="center", va="center", labelpad=20)


xticklabels = ax1.get_xticklabels() + ax2.get_xticklabels()
setp(xticklabels, visible=False)

plt.savefig(name0, bbox_inches='tight', dpi=600)            # save figure in a high quality PNG file
plt.close()
#################################################################################################




#################################################################################################
### Apply size correction + calibration to samples + export data in a file   ####################
#################################################################################################
### Calculate the size corr and calibrated samples data
samples.d15Nscc  = f.mixing_model_d15NNH4_REV(samples.Amea, f.calibration_d15NNH4(samples.d15Nmea, slope, inter), Abk, d15Nbk, d15Naz)


### Write output data
filename = 'output_SAMPLES_size_corrected+calibrated.txt'
out = open(filename, "w" )
out.write("'ID1 \t ID2 \t Area32 \t Area28 \t ratio28/32 \t d17Omea \t d18Omea \t D17Omea \t d15Nmea \t ******* \t d15Nscc \t d15Nerr")
out.write("\n")
for i in range(len(samples.d15Nscc)):
    out.write(samples.ID1[i])
    out.write("\t")
    out.write(samples.ID2[i])
    out.write("\t")
    out.write(str(samples.m32[i]).replace('.', ','))
    out.write("\t")
    out.write(str(samples.m28[i]).replace('.', ','))
    out.write("\t")
    out.write(str(samples.ratio[i]).replace('.', ','))
    out.write("\t")
    out.write(str(samples.d17Omea[i]).replace('.', ','))
    out.write("\t")
    out.write(str(samples.d18Omea[i]).replace('.', ','))
    out.write("\t")
    out.write(str(samples.D17Omea[i]).replace('.', ','))
    out.write("\t")
    out.write(str(samples.d15Nmea[i]).replace('.', ','))
    out.write("\t ******* \t")
    out.write(str(samples.d15Nscc[i]).replace('.', ','))
    out.write("\t")
    out.write(str(std(d15Nres)).replace('.', ','))
    out.write( "\n" )
out.close()
#################################################################################################
