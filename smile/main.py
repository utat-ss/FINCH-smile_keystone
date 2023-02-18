# # # UPDATE IMPERATIVES:
# 1. Conglomerate all config variables into the main function. DONE
# 2. Group one-function files into larger, single files. 
# 3. Handle imports of relevant packages in the scope of those files. 
# 4. Have each function take in global variables as parameters (So they don't rely
#    on the global scope, which is bad practice).

# main is operational; but there's a *lot* we could be doing to make this cleaner. Todo. - Andy


import numpy as np
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


# Function imports (DEPRECATED)
  # from create_ref_and_test_spectra import create_ref_and_test_spectra
  # from data_matrix_collapse import data_matrix_collapse
  # from determine_min_sa import determine_min_sa
  # from load_datacube_npy import ldn
  # import optical_sensor
  # from spectral_angle_calculation import spectral_angle_calculation
  # from spectral_response import test_spectral_response, make_random_SRFs
  # from smile_correction import smile_correction
  # from spline_interpolation import  spline_interpolation_1_pixel, spline_interpolation_all

#Function imports
from Quantificationfns import *
from Correctionfns import *

print("Imports completed, no issues.")

#Andy
if __name__ == '__main__':

  data = radianceData
  wavelength_input = wavelength
  wavelength_increment_input = wavelength_increment
  # Globals further defined here in case we want to pass in other data

  # # # Quantification
  column_averaged_spectra = data_matrix_collapse(data) # Step 1: Generate Column Averaged Spectra.
  demo_SRF = test_spectral_response # Step 2: Generate SRFs.
  ref_spectra, test_spectra = create_ref_and_test_spectra((0,1000), column_averaged_spectra, demo_SRF, wavelength, g_num_of_bands, g_num_shifts_1D, g_shift_increment) # Step 3, 4: Generate Reference and Test spectra. (0,100) is also a placeholder
  sa_deg = spectral_angle_calculation(test_spectra, ref_spectra, g_data_dim) # Step 5: Calculate spectral angle from test and reference spectra.
  min_spectral_angle = determine_min_sa(sa_deg, g_data_dim) # Step 6: Determine minimum spectral angle.

  print("Quantification complete, no issues.")
  
  # # # Correction
  # Step 7: Apply reverse of Quantified shifts to SRFs. DEPRECATED. STEP 7 IS NOW BUNDLED INTO STEP 9.
  # reverse_shifted_SRFS = reverse_shift(data, demo_SRF, min_spectral_angle) 
  spectra_wav, spectra_rad = spline_interpolation_all(data, wavelength_input, wavelength_increment_input) # Step 8: Spline interpolation of test spectra. INPUTS ALSO MAY NOT BE CORRECT
  # Step 9: Generate smile corrected spectra for each pixel.
  smile_correction(spectra_rad, min_spectral_angle, test_spectral_response)

  # The function runs through from beginning to end. Whether it's functional is yet to be seen.