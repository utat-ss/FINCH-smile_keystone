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
from extract_from_MODTRAN import extract_from_MODTRAN
from config import *

#Function imports
from GenerateSmileFns import *
from Quantificationfns import *
from Correctionfns import *

# Log file
log_file_path = f'{data_folder_path}execution_times.txt'

# Function to log messages
def log_message(message):
    with open(log_file_path, 'a') as log_file:
        log_file.write(message + '\n')

start = timeit.default_timer()

# Start logging
with open(log_file_path, 'w') as log_file:
    log_file.write("Execution log:\n")

# Step-by-step timings
step_start = timeit.default_timer()

# # Load indian pine array data
full_indian_pine_array = indian_pine_array
# indian_pine_wavelength = np.loadtxt(indian_pine_wavelength_filepath)
indian_pine_wavelength = np.linspace(400, 2500, 220)


feature_start, feature_end = get_feature_index(indian_pine_wavelength, feature)
cropped_indian_pine_array = full_indian_pine_array[:, :, feature_start:feature_end]  # Cropped data
cropped_wavelength = indian_pine_wavelength[feature_start:feature_end]
step_end = timeit.default_timer()
log_message(f"Step: Load and crop Indian Pine data - {step_end - step_start:.4f} seconds")

# # Load MODTRAN data
step_start = timeit.default_timer()
MODTRAN_x, MODTRAN_data= extract_from_MODTRAN(Reference_data_filepath)

# Crop MODTRAN data for quantification
cropped_MODTRAN_x = MODTRAN_x[feature_start:feature_end]  # Cropped MODTRAN x
cropped_MODTRAN_data = MODTRAN_data[feature_start:feature_end]  # Cropped MODTRAN data
np.savez_compressed(f'{data_folder_path}cropped_MODTRAN_data', MODTRAN_wl=cropped_MODTRAN_x, MODTRAN_data=cropped_MODTRAN_data)
print("MODTRAN data loaded, no issues.")
step_end = timeit.default_timer()
log_message(f"Step: Load and crop MODTRAN data - {step_end - step_start:.4f} seconds")

Reference_wl = MODTRAN_x
Reference_data = MODTRAN_data

# # Global Vars
# The rest of the global variables are defined in py
# wavelength_source, radianceData, g_data_dim, wavelength, wavelength_increment = load_datacube_npy.ldn(indian_pine_array_filepath, indian_pine_wavelength_filepath)

# NOTE: this doesn't seem to be used anywhere.
g_width_of_band = len(wavelength)/g_num_of_bands

if feature is not None:
    feature_start, feature_end = get_feature_index(wavelength, feature)
    Reference_wl = Reference_wl[feature_start:feature_end]
    Reference_data = Reference_data[feature_start:feature_end]
    wavelength_input = wavelength[feature_start:feature_end]

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


print("Imports completed, no issues.")
print(wavelength[1])

#Andy
if __name__ == '__main__':

    data = radianceData
    wavelength_input = wavelength
    wavelength_increment_input = wavelength_increment
    # Globals further defined here in case we want to pass in other data

    # # # Quantification
    # Step 0: Generate artificil SMILE shift in the radianceData
    step_start = timeit.default_timer()
    cropped_data = generate_smile_shift(cropped_indian_pine_array)  # Quantify using cropped data
    print ("Step 0 Done, no issues.")
    step_end = timeit.default_timer()
    log_message(f"Step 0: Generate artificial SMILE shift - {step_end - step_start:.4f} seconds")

    # Step 1: Generate Column Averaged Spectra.
    step_start = timeit.default_timer()
    column_averaged_spectra = data_matrix_collapse(data)
    print ("Step 1 Done, no issues.")
    np.savez_compressed(f'{data_folder_path}column_averaged_spectra', cas = column_averaged_spectra)
    step_end = timeit.default_timer()
    log_message(f"Step 1: Generate Column Averaged Spectra - {step_end - step_start:.4f} seconds")

    # Step 2: Generate SRFs.
    step_start = timeit.default_timer()
    demo_SRF = test_spectral_response
    print ("Step 2 Done, no issues.")
    step_end = timeit.default_timer()
    log_message(f"Step 2: Generate SRFs - {step_end - step_start:.4f} seconds")

    # Step 3, 4: Generate Reference and Test spectra.
    step_start = timeit.default_timer()
    if feature is not None:
        column_averaged_spectra = column_averaged_spectra[:, feature_start:feature_end]
    ref_spectra, test_spectra = create_ref_and_test_spectra(column_averaged_spectra, demo_SRF, cropped_wavelength, ref_spectra=cropped_MODTRAN_data)
    np.savez_compressed(f'{data_folder_path}ref_and_test_spectra', ref = ref_spectra, test = test_spectra)
    print ("Step 3, 4 Done, no issues.")
    step_end = timeit.default_timer()
    log_message(f"Step 3, 4: Generate Reference and Test spectra - {step_end - step_start:.4f} seconds")

    # Step 5: Calculate spectral angle from test and reference spectra.
    step_start = timeit.default_timer()
    sa_deg = spectral_angle_calculation(test_spectra["spectra"], ref_spectra["spectra"], g_data_dim)
    print ("Step 5 Done, no issues.")
    np.savez_compressed(f'{data_folder_path}sa_deg', sa_deg = sa_deg)
    step_end = timeit.default_timer()
    log_message(f"Step 5: Calculate spectral angle from test and reference spectra - {step_end - step_start:.4f} seconds")

    # Step 6: Determine minimum spectral angle.
    step_start = timeit.default_timer()
    min_spectral_angle = determine_min_sa_shift(sa_deg, g_data_dim)
    print(min_spectral_angle, "at length = ", len(min_spectral_angle))
    print ("Step 6 Done, no issues.")
    np.savez_compressed(f'{data_folder_path}min_spectral_angle', msa = min_spectral_angle)
    print(f"Quantification complete, no issues. Data saved to {data_folder_path}")
    step_end = timeit.default_timer()
    log_message(f"Step 6: Determine minimum spectral angle - {step_end - step_start:.4f} seconds")

    # # # Correction
    # Step 7: Apply reverse of Quantified shifts to SRFs. DEPRECATED. STEP 7 IS NOW BUNDLED INTO STEP 9.
    # reverse_shifted_SRFS = reverse_shift(data, demo_SRF, min_spectral_angle)
    print ("There is no step 7, because it was deprecated.")

    # Step 8: Spline interpolation of test spectra. INPUTS ALSO MAY NOT BE CORRECT. TODO: Check if imputs are still not correct.
    step_start = timeit.default_timer()
    spectra_wav, spectra_rad = spline_interpolation_all(full_indian_pine_array, indian_pine_wavelength, wavelength_increment, g_data_dim, indian_pine_wavelength)
    print("Step 8 Done, no issues.")
    np.savez_compressed(f'{data_folder_path}spline_interpolated', wav = spectra_wav, rad = spectra_rad)
    step_end = timeit.default_timer()
    log_message(f"Step 8: Spline interpolation of test spectra - {step_end - step_start:.4f} seconds")

    # Step 9: Generate smile corrected spectra for each pixel. (This turns out to be the most computationally expensive step.)
    step_start = timeit.default_timer()
    interpolated_wl = np.transpose(spectra_wav, ((1, 2, 0)))
    corrected_datacube = smile_correction(spectra_rad, min_spectral_angle, test_spectral_response, np.transpose(cropped_wavelength, ((1, 2, 0))))
    np.savez_compressed(f'{data_folder_path}corrected_datacube', corrected_data = corrected_datacube)
    print(f"Correction complete, no issues. Data saved to {data_folder_path}.")
    step_end = timeit.default_timer()
    log_message(f"Step 9: Generate smile corrected spectra for each pixel - {step_end - step_start:.4f} seconds")

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

print(f"Total time taken: {end - start:.4f} seconds")
print(f"Execution timings logged in {log_file_path}")
