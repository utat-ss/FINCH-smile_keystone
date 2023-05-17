# # # UPDATE IMPERATIVES:
# 1. Conglomerate all config variables into the main function. DONE
# 2. Group one-function files into larger, single files. 
# 3. Handle imports of relevant packages in the scope of those files. 
# 4. Have each function take in global variables as parameters (So they don't rely
#    on the global scope, which is bad practice).

# main is operational; but there's a *lot* we could be doing to make this cleaner. Todo. - Andy
# TODO: (tentative) Print the time it takes to run each step in a .txt file. - Shuhan
# TODO: (After project completion) Check each function's docstrings - Shuhan
# TODO: Feed cropped MODTRAN and Indian Pine data into quantification, and feed quantified smile distortion and FULL indian pine data into correction. - Shuhan


# Program starts here
import numpy as np
import timeit

import load_datacube_npy
import config

start = timeit.default_timer()

# # Load indian pine array data
indian_pine_array = np.load(config.indian_pine_array_filepath)
# indian_pine_wavelength = np.loadtxt(config.indian_pine_wavelength_filepath)
indian_pine_wavelength = np.linspace(400, 2500, 220)

# # Load MODTRAN data
from extract_from_MODTRAN import extract_from_MODTRAN

MODTRAN_x, MODTRAN_data= extract_from_MODTRAN(config.Reference_data_filepath)
np.savez_compressed(f'{config.data_folder_path}MODTRAN_data', MODTRAN_wl = MODTRAN_x, MODTRAN_data = MODTRAN_data)
print("MODTRAN data loaded, no issues.")

Reference_wl = MODTRAN_x
Reference_data = MODTRAN_data

# # Global Vars
# The rest of the global variables are defined in config.py
wavelength_source, radianceData, g_data_dim, wavelength, wavelength_increment = load_datacube_npy.ldn(config.indian_pine_array_filepath, config.indian_pine_wavelength_filepath)

g_width_of_band = len(wavelength)/config.g_num_of_bands

# Function imports (DEPRECATED)
    # from create_ref_and_test_spectra import create_ref_and_test_spectra
    # from data_matrix_collapse import data_matrix_collapse
    # from min_sa_shift import min_sa_shift
    # from load_datacube_npy import ldn
    # import optical_sensor
    # from spectral_angle_calculation import spectral_angle_calculation
    # from spectral_response import test_spectral_response, make_random_SRFs
    # from smile_correction import smile_correction
    # from spline_interpolation import  spline_interpolation_1_pixel,   spline_interpolation_all

#Function imports
from Quantificationfns import *
from Correctionfns import *

print("Imports completed, no issues.")
print(wavelength[1])

#Andy
if __name__ == '__main__':

    data = radianceData
    wavelength_input = wavelength
    wavelength_increment_input = wavelength_increment
    # Globals further defined here in case we want to pass in other data

    # # # Quantification
    # Step 1: Generate Column Averaged Spectra.
    column_averaged_spectra = data_matrix_collapse(data)
    print ("Step 1 Done, no issues.")
    np.savez_compressed(f'{config.data_folder_path}column_averaged_spectra', cas = column_averaged_spectra)

    # Step 2: Generate SRFs.
    demo_SRF = test_spectral_response
    print ("Step 2 Done, no issues.")

    # Step 3, 4: Generate Reference and Test spectra.
    # But first, select only the section corresponding to the feature range.
    if config.feature is not None:
        feature_start, feature_end = config.get_feature_index(wavelength_input, config.feature)
        column_averaged_spectra = column_averaged_spectra[:, feature_start:feature_end]
        wavelength_input = wavelength_input[feature_start:feature_end]
    
    ref_spectra, test_spectra = create_ref_and_test_spectra(column_averaged_spectra, demo_SRF, wavelength_input, ref_spectra=Reference_data)
    print ("Step 3, 4 Done, no issues.")
    np.savez_compressed(f'{config.data_folder_path}ref_and_test_spectra', ref = ref_spectra, test = test_spectra)

    # Step 5: Calculate spectral angle from test and reference spectra.
    sa_deg = spectral_angle_calculation(test_spectra[0], ref_spectra[0], g_data_dim)
    print ("Step 5 Done, no issues.")
    np.savez_compressed(f'{config.data_folder_path}sa_deg', sa_deg = sa_deg)

    # Step 6: Determine minimum spectral angle.
    min_spectral_angle = determine_min_sa_shift(sa_deg, g_data_dim)
    print ("Step 6 Done, no issues.")
    np.savez_compressed(f'{config.data_folder_path}min_spectral_angle', msa = min_spectral_angle)

    print(f"Quantification complete, no issues. Data saved to {config.data_folder_path}")

    # # # Correction
    # Step 7: Apply reverse of Quantified shifts to SRFs. DEPRECATED. STEP 7 IS NOW BUNDLED INTO STEP 9.
    # reverse_shifted_SRFS = reverse_shift(data, demo_SRF, min_spectral_angle)
    print ("There is no step 7, because it was deprecated.")

    # Step 8: Spline interpolation of test spectra. INPUTS ALSO MAY NOT BE CORRECT. TODO: Check if imputs are still not correct.
    spectra_wav, spectra_rad = spline_interpolation_all(data, wavelength, wavelength_increment_input, g_data_dim, wavelength)
    print("Step 8 Done, no issues.")
    np.savez_compressed(f'{config.data_folder_path}spline_interpolated', wav = spectra_wav, rad = spectra_rad)

    # Step 9: Generate smile corrected spectra for each pixel. (This turns out to be the most computationally expensive step.)
    corrected_datacube = smile_correction(spectra_rad, min_spectral_angle, test_spectral_response, config.g_num_of_bands, config.g_shift_increment, wavelength)
    np.savez_compressed(f'{config.data_folder_path}corrected_datacube', corrected_data = corrected_datacube)
    print(f"Correction complete, no issues. Data saved to {config.data_folder_path}.")

    # # # Clear all variables
    column_averaged_spectra = None
    demo_SRF = None
    [ref_spectra, test_spectra] = None, None
    sa_deg = None
    min_spectral_angle = None
    [spectra_wav, spectra_rad] = None, None
    corrected_datacube = None

    # The function runs through from beginning to end. Whether it's functional is yet to be seen.

end = timeit.default_timer()

print(f"Total time taken: {end - start} seconds.")