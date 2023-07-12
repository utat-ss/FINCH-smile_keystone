# # # CRITICALLY IMPORTANT! In order to acccess these variables, you have two options:
# # # 1. from config import * (Which will automatically let you access these as though they were defined in your current file)
# # # 2. import config (Which will let you access these as config.variable_name ONLY. In other words you may need to rewrite your code)
# # # The first option needs less revision to make work, but the second makes the job of your code linter easier. 
# # # While testing, use option 1.

# # package imports
import numpy as np
import matplotlib.pyplot as plt
import math

from scipy import stats
from IPython.display import clear_output

# Recursively import the whole lot of functions
import load_datacube_npy

# # indian pine array data
indian_pine_array_filepath = 'data\indian_pine_array.npy'
indian_pine_array = np.load(indian_pine_array_filepath)
indian_pine_wavelength_filepath = 'data\indian_pine_wavelength.txt'
indian_pine_wavelength = np.load(indian_pine_array_filepath)

# # Global Vars
wavelength_source, radianceData, g_data_dim, wavelength, wavelength_increment = load_datacube_npy.ldn(indian_pine_array_filepath, indian_pine_wavelength_filepath)

#wavelength_source, radianceData, g_data_dim = load_datacube('/pavia.npy', '/content/Science/fall_2021_onboarding/wavelength.txt')

#number of spectral pixels we're assuming our instrument has. 

#The g_data_dim[0] bands of reference spectra will be resampled to g_num_of_bands points
g_num_of_bands = 70 

g_num_shifts_1D = 5
g_shift_increment = .2 #nanometers
g_total_shifts = int(2*(g_num_shifts_1D/g_shift_increment)+1)
