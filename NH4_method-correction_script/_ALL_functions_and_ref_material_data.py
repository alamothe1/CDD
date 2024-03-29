### IMPORT MODULES ##############################################################################
import numpy as np                  # array module
import re
#################################################################################################

# In this file, we attribute reference values and formats to each piece of reference material
# All stds : NO3-, NO2- and NH4+
# We also add the functions used to read the data and to calculate the (size corrected) raw values


species = {'A' : 'NO3-' ,
           'C' : 'NO3-' ,
           'E' : 'NO3-' ,
           'F' : 'NO3-' ,
           'Z' : 'NO3-' ,
           'J' : 'NH4+' ,
           'K' : 'NH4+' ,
           'L' : 'NH4+' ,
           'P' : 'NO2-' ,
           'Q' : 'NO2-' , 
           'R' : 'NO2-'     }

class Reference_material_data(object):
    "Define the reference material data class"

val = Reference_material_data()

val.d15Nstd = {'A' :    2.7 ,
               'C' :   0.45 ,
               'E' :   -1.8 ,
               'F' :   90.7 ,
               'Z' :  180.0 ,
               'J' :    0.4 ,
               'K' :  -30.4 ,
               'L' :   53.7 ,
               'P' :  -79.6 ,
               'Q' :    2.8 ,
               'R' :    3.7 }

val.d17Ostd = {'A' :  51.46 ,
               'C' :  18.33 ,
               'E' :  -14.8 ,
               'F' :  -0.72 ,
               'Z' :   13.4 ,
               'J' :  -9999 ,
               'K' :  -9999 ,
               'L' :  -9999 ,
               'P' :  -9999 ,
               'Q' :  -9999 ,
               'R' :  -9999 }

val.d18Ostd = {'A' :   57.5 ,
               'C' :  14.79 ,
               'E' :  -27.9 ,
               'F' :  -1.12 ,
               'Z' :   25.7 ,
               'J' :  -9999 ,
               'K' :  -9999 ,
               'L' :  -9999 ,
               'P' :    4.5 ,
               'Q' :   88.5 ,
               'R' :   11.4 }

val.D17Ostd = {'A' :  21.56 ,
               'C' :  10.64 ,
               'E' :   -0.3 ,
               'F' :  -0.14 ,
               'Z' :    0.0 ,
               'J' :  -9999 ,
               'K' :  -9999 ,
               'L' :  -9999 ,
               'P' :  -9999 ,
               'Q' :  -9999 ,
               'R' :  -9999 }

color   = {'A'   : 'blue',
           'C'   : 'green',
           'E'   : 'red',
           'F'   : 'grey',
           'Z'   : 'black',
           'J'   : 'green' ,
           'K'   : 'blue' ,
           'L'   : 'red'  ,
           'P'   : 'black',
           'Q'   : 'green',
           'R'   : 'blue' ,
           'blk' : 'grey',
           'comp_blk' : 'yellow',
           'spl' : 'black'  }

mec     = {'A'   : 'blue',
           'C'   : 'green',
           'E'   : 'red',
           'F'   : 'grey',
           'Z'   : 'black',
           'J'   : 'green' ,
           'K'   : 'blue' ,
           'L'   : 'red'  ,
           'P'   : 'black',
           'Q'   : 'green',
           'R'   : 'blue' ,
           'blk' : 'black',
           'comp_blk' : 'black',
           'spl' : 'black'  }

mfc   = {'A'   : 'blue',
           'C'   : 'green',
           'E'   : 'red',
           'F'   : 'grey',
           'Z'   : 'black',
           'J'   : 'green' ,
           'K'   : 'blue' ,
           'L'   : 'red'  ,
           'P'   : 'black',
           'Q'   : 'green',
           'R'   : 'blue' ,
           'blk' : 'white',
           'comp_blk' : 'yellow',
           'spl' : 'black'  }

style   = {'A'   : 'o'    ,
           'C'   : 'o'    ,
           'E'   : 'o'    ,
           'F'   : 'o'    ,
           'Z'   : 'o'    ,
           'J'   : 'o'    ,
           'K'   : 'o'    ,
           'L'   : 'o'    ,
           'P'   : '*'    ,
           'Q'   : 's'    ,
           'R'   : '^'    ,
           'blk' : 'o'    ,
           'comp_blk' : '*',
           'spl' : '+'      }

ms      = {'A'   :   5  ,
           'C'   :   5  ,
           'E'   :   5  ,
           'F'   :   5  ,
           'Z'   :   5  ,
           'J'   :   5  ,
           'K'   :   5  ,
           'L'   :   5  ,
           'P'   :   7  ,
           'Q'   :   5  ,
           'R'   :   5  ,
           'blk' :   7  ,
           'comp_blk' : 12,
           'spl' :   9       }




#################################################################################################
### Functions            ########################################################################
#################################################################################################
class Sample(object):
    "Define the sample class"

    
class Isotopic_data(object):
    "Define the isotopic data class"

    def __init__(self):
        self.ID      = 0
        self.amount  = 0
        self.m28     = 0
        self.m32     = 0
        self.ratio   = 0
        self.d17Omea = 0
        self.d18Omea = 0
        self.D17Omea = 0
        self.d15Nmea = 0
        self.d17Ostd = 0
        self.d18Ostd = 0
        self.D17Ostd = 0
        self.d15Nstd = 0

    def read_file_without_header(self, filename, nb_line_header):
        "Reads a file and removes a header"

        file_op = open(filename, 'r')
        lines = file_op.readlines()
        file_op.close()

        if 'samples' in filename:
            self.ID1 = []
            self.ID2 = []
        elif 'stds' in filename:
            self.ID       = []
            self.amount   = []

        self.m28      = []
        self.m32      = []
        self.ratio    = []
        self.d17Omea  = []
        self.d18Omea  = []
        self.D17Omea  = []
        self.d15Nmea  = []
        self.d17Ostd  = []
        self.d18Ostd  = []
        self.D17Ostd  = []
        self.d15Nstd  = []

        for i in range(nb_line_header, len(lines)):
            line_read = lines[i].replace(',', '.')
            box = re.split(r"\t", line_read)

            if 'samples' in filename:
                self.ID1.append(box[0])
                self.ID2.append(box[1])
            elif 'stds' in filename:
                self.ID.append(box[0])
                self.amount.append(float(box[1]))
                self.d17Ostd.append(val.d17Ostd[box[0]])
                self.d18Ostd.append(val.d18Ostd[box[0]])
                self.D17Ostd.append(val.D17Ostd[box[0]])
                self.d15Nstd.append(val.d15Nstd[box[0]])
            self.m32.append(float(box[2]))
            self.m28.append(float(box[3]))
            self.ratio.append(float(box[4]))
            self.d17Omea.append(float(box[5]))
            self.d18Omea.append(float(box[6]))
            self.D17Omea.append(float(box[7]))
            self.d15Nmea.append(float(box[8]))


def mixing_model_d15NNH4(Amea, d15Nmea, Abk, d15Nbk, d15Naz):
    "Apply the mixing model to generate the array of computed values from the std values"
    
    return ((Amea/2 - Abk)*d15Nmea + Abk*d15Nbk + (Amea/2)*d15Naz) / Amea

def mixing_model_d15NNH4_REV(Amea, d15Ncal, Abk, d15Nbk, d15Naz):
    "Reverse the mixing model to generate the array of size corrected and calibrated values"

    return (Amea*d15Ncal - Abk*d15Nbk - (Amea/2)*d15Naz) / (Amea/2 - Abk)

def calibration_d15NNH4(d15Nmea, slope, inter):
    "Apply the calibration to generate calibrated values from measured values"

    return slope*d15Nmea + inter

def calibration_d15NNH4_REV(d15Ncal, slope, inter):
    "Reverse the calibration to compute measured values from calibrated values"

    return (d15Ncal - inter)/slope





### Export / Read size corr + calib parameters for d15NNH4+

def export_size_corr_calib_params_d15NNH4(Abk, d15Nbk, d15Naz, slope, inter, filename):
        "Export the size corr and calibration parameters in a file"

        output_file = open(filename, "w" )
        output_file.write("Blank_size_(Vs)______ \t " + str(Abk).replace('.', ',') + "\n")
        output_file.write("Blank_d15N_(per_mill) \t " + str(d15Nbk).replace('.', ',') + "\n")
        output_file.write("Azide_d15N_(per_mill) \t " + str(d15Naz).replace('.', ',') + "\n")
        output_file.write("Slope________________ \t " + str(slope).replace('.', ',') + "\n")
        output_file.write("Intercept____________ \t " + str(inter).replace('.', ',') + "\n")
        output_file.close()

def read_file_size_corr_calib_params_d15NNH4(filename):
        "Reads the file with the size corr and calibration parameters"

        file_op = open(filename, 'r')
        lines = file_op.readlines()
        file_op.close()

        Abk    = float(lines[0].split()[1].replace(',', '.'))
        d15Nbk = float(lines[1].split()[1].replace(',', '.'))
        d15Naz = float(lines[2].split()[1].replace(',', '.'))
        slope  = float(lines[3].split()[1].replace(',', '.'))
        inter  = float(lines[4].split()[1].replace(',', '.'))
        return Abk, d15Nbk, d15Naz, slope, inter





