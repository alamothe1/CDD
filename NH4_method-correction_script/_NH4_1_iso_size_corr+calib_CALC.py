
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

from matplotlib.font_manager import FontProperties
#################################################################################################



### Main inputs to the calculation
input_standards_filename = 'input_measured_iso_stds.txt'      # Name of the file with all the isotopic data of the standards
# input_samples_filename   = 'input_measured_samples.txt'       # Name of the file with all the isotopic data of the samples - inutile ici




import _ALL_functions_and_ref_material_data as f



#################################################################################################
### Read inputs          ########################################################################
#################################################################################################
### Imports data for ALL ISOTOPIC STANDARDS
# Create object
IS_all = f.Isotopic_data()              # "IS" means isotopic standards
IS_all.read_file_without_header(input_standards_filename, 1)        # reads the file

# Convert data to arrays
# Only N data here
Amea    = array(IS_all.m28)
d15Nmea = array(IS_all.d15Nmea)
d15Nstd = array(IS_all.d15Nstd)
#################################################################################################





#################################################################################################
### Determine BLANK and azide PROPERTIES   ######################################################
#################################################################################################
step = 0
condition = 0 # Condition to come out of the while loop

while condition == 0:
    step = step + 1
    print ('STEP ' + str(step))
    print ('--------------------------------------------------')
    ## Create list for results
    L = []

    if step == 1:       
        ### Find size and isotopic composition of the blank
        Abk_min,        Abk_max,        Abk_stp         =   0,      2,     0.1 #Vs
        d15Nbk_min,     d15Nbk_max,     d15Nbk_stp      =   -15,   -5,     0.5 #per mill
        d15Naz_min,     d15Naz_max,     d15Naz_stp      =   -15,   -5,     0.1 #per mill
    else:
        pass

    ## Write the values of each parameter
    print ('Search for :')
    print ('Abk in    [ ' + ' '*(5-len(str(round(Abk_min, 2)))) + str(round(Abk_min, 2))        + ', ' + ' '*(5-len(str(round(Abk_max, 2)))) + str(round(Abk_max, 2))       + '] @ ' + ' '*(5-len(str(round(Abk_stp, 2)))) + str(round(Abk_stp, 2)) + ' Vs')
    print ('d15Nbk in [ ' + ' '*(5-len(str(round(d15Nbk_min, 1)))) + str(round(d15Nbk_min, 1))  + ', ' + ' '*(5-len(str(round(d15Nbk_max, 1)))) + str(round(d15Nbk_max, 1)) + '] @ ' + ' '*(5-len(str(round(d15Nbk_stp, 1)))) + str(round(d15Nbk_stp, 1)) + ' per mill')
    print ('d15Naz in [ ' + ' '*(5-len(str(round(d15Naz_min, 1)))) + str(round(d15Naz_min, 1))  + ', ' + ' '*(5-len(str(round(d15Naz_max, 1)))) + str(round(d15Naz_max, 1)) + '] @ ' + ' '*(5-len(str(round(d15Naz_stp, 1)))) + str(round(d15Naz_stp, 1)) + ' per mill')

    ## Warn user if ranges are not correct
    if Abk_max - Abk_min < 2:
        print ('********************** Problem : Abk range is too small.')
        print ('********************** It must be at least 2 Vs wide. Please, adjust the values')
        print ('********************** Program stopped!!!')
        break
    elif d15Nbk_max - d15Nbk_min < 10:
        print ('********************** Problem : d15Nbk range is too small.')
        print ('********************** It must be at least 10 per mill wide. Please, adjust the values')
        print ('********************** Program stopped!!!')
        break
    elif d15Naz_max - d15Naz_min < 10:
        print ('********************** Problem : d15Naz range is too small.')
        print ('********************** It must be at least 10 per mill wide. Please, adjust the values')
        print ('********************** Program stopped!!!')
        break

    ## List all possible values for each parameter
    L_Abk    = arange(Abk_min, Abk_max, Abk_stp)
    L_d15Nbk = arange(d15Nbk_min, d15Nbk_max, d15Nbk_stp)
    L_d15Naz = arange(d15Naz_min, d15Naz_max, d15Naz_stp)


    print ('--------------------------------------------------')
    print ('Number of cases to compute = ' + str(len(L_Abk)*len(L_d15Nbk)*len(L_d15Naz)))
    print ('--------------------------------------------------')

    ctot  = 0
    pos   = 0
    comp  = 1e6


    for dz in L_d15Naz:
        for Ai in L_Abk:
            for dj in L_d15Nbk:
                # For each standard value, apply the mixing model to calculate the theoretical values from the accepted values according to the mixing model
                d15Nmix = f.mixing_model_d15NNH4(Amea, d15Nstd, Ai, dj, dz)
                # Calibrate the system by comparing these computed values with the measured ones to account for isotopic fractionation in the N2O line or in the MS
                slope, inter = polyfit(d15Nmea, d15Nmix, 1)

                ### Now reverse the previous approach to calculate the calibrated values obtained for each standard
                # Apply the calibration 
                d15Ncal = f.calibration_d15NNH4(d15Nmea, slope, inter)
                # Reverse the mixing model to obtain the size corrected and calibrated values from the measured ones
                d15Nscc = f.mixing_model_d15NNH4_REV(Amea, d15Ncal, Ai, dj, dz)   # "scc" stands for size corrected and calibrated
                errnorm_mean = mean(d15Nscc - d15Nstd)
                errnorm_std  = std(d15Nscc - d15Nstd)

                if  errnorm_std < comp:
                    comp = errnorm_std
                    pos  = ctot


                L.append(f.Sample())
                L[ctot].Amea         = Amea
                L[ctot].d15Nstd      = d15Nstd
                L[ctot].d15Nscc      = d15Nscc
                L[ctot].errnorm_mean = errnorm_mean
                L[ctot].errnorm_std  = errnorm_std
                L[ctot].slope        = slope
                L[ctot].inter        = inter
                L[ctot].Abk          = Ai
                L[ctot].d15Nbk       = dj
                L[ctot].d15Naz       = dz
                
                ctot = ctot + 1
    print(pos)
             
    #################################################################################################


    print ('Optimized blank size : ' + ' '*(5-len(str(round(L[pos].Abk, 2))))    + str(round(L[pos].Abk, 2))    + ' Vs')
    print ('Optimized blank d15N : ' + ' '*(5-len(str(round(L[pos].d15Nbk, 1)))) + str(round(L[pos].d15Nbk, 1)) + ' per mill')
    print ('Optimized azide d15N : ' + ' '*(5-len(str(round(L[pos].d15Naz, 1)))) + str(round(L[pos].d15Naz, 1)) + ' per mill')
    print ('---')
    print ('Slope                : ' + ' '*(5-len(str(round(L[pos].slope , 3))))      + str(round(L[pos].slope, 3)))
    print ('Intercept            : ' + ' '*(5-len(str(round(L[pos].inter , 3))))      + str(round(L[pos].inter, 3)))
    print ('---')
    print ('Error :')
    print ('- mean               : ' + ' '*(5-len(str(round(L[pos].errnorm_mean, 2))))      + str(round(L[pos].errnorm_mean, 2)) + ' per mill')
    print ('- stdev              : ' + ' '*(5-len(str(round(L[pos].errnorm_std , 2))))      + str(round(L[pos].errnorm_std, 2))  + ' per mill')

    print ('==================================================')

    ## Decide whether one more step is necessary
    range_shift = 0
    
    if L[pos].Abk == Abk_min:
        print ('********************** Problem : check the Abk_min value')
        print ('********************** Program stopped!!!')
        break
    elif L[pos].Abk >= Abk_max - 0.3:
        print ('********************** Abk range is shifted by 1 Vs')
        Abk_max = Abk_max + 1
        range_shift = 1
    elif L[pos].d15Nbk <= d15Nbk_min + 1:
        print ('********************** d15Nbk range is shifted by -5 per mill')
        d15Nbk_min = d15Nbk_min - 5
        d15Nbk_max = d15Nbk_max - 5
        range_shift = 1
    elif L[pos].d15Nbk >= d15Nbk_max - 1:
        print ('********************** d15Nbk range is shifted by +5 per mill')
        d15Nbk_max = d15Nbk_max + 5
        d15Nbk_min = d15Nbk_min + 5
        range_shift = 1
    elif L[pos].d15Naz <= d15Naz_min + 1:
        print ('********************** d15Naz range is shifted by -5 per mill')
        d15Naz_min = d15Naz_min - 5
        d15Naz_max = d15Naz_max - 5
        range_shift = 1
    elif L[pos].d15Naz >= d15Naz_max - 1:
        print ('********************** d15Naz range is shifted by +5 per mill')
        d15Naz_max = d15Naz_max + 5
        d15Naz_min = d15Naz_min + 5
        range_shift = 1

    if range_shift != 1:
        if Abk_stp > 0.05:
            print ('********************** Abk step range is set to 0.05 Vs')
            Abk_stp = 0.05
        elif d15Nbk_stp > 0.1:
            print ('********************** d15Nbk step range is set to 0.1 per mill')
            d15Nbk_stp = 0.1
        elif d15Naz_stp > 0.1:
            print ('********************** d15Naz step range is set to 0.1 per mill')
            d15Naz_stp = 0.1
        else:
            condition = 1
#################################################################################################

            


#################################################################################################
### Generate file with parameters for size correction and calibration     #######################
#################################################################################################
print ('')
print ('')
print ('')
print ('=======================================================')
print ('=======================================================')
print ('Parameters for size correction and calibration EXPORTED')
print ('=======================================================')
print ('=======================================================')
f.export_size_corr_calib_params_d15NNH4(L[pos].Abk, L[pos].d15Nbk, L[pos].d15Naz, L[pos].slope, L[pos].inter, 'output_params_for_corr_NH4_Nitrogen.txt')
#################################################################################################
