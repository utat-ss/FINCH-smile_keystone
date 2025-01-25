"""
Main Script for SMILE Quantification and Correction

This script performs SMILE quantification and correction on hyperspectral data. It involves data preprocessing, 
quantification of spectral misregistration, and application of corrective measures to generate smile-corrected 
spectra. The process is logged with step-by-step execution times for performance tracking.

Key Steps:
1. **Load and Crop Indian Pine Data**:
   - Loads the hyperspectral datacube and wavelength array.
   - Crops the data based on the specified wavelength feature range.

2. **Load and Crop MODTRAN Data**:
   - Extracts MODTRAN data for the specified wavelength range for use as a reference.

3. **Quantification of SMILE Effects**:
   - Generates artificial SMILE shifts in the data.
   - Computes column-averaged spectra for quantification.
   - Generates Spectral Response Functions (SRFs).
   - Produces reference and test spectra by comparing MODTRAN data with adjusted column-averaged spectra.
   - Calculates spectral angles between test and reference spectra to measure SMILE.

4. **Correction of SMILE Effects**:
   - Uses the minimum spectral angles to adjust the hyperspectral data.
   - Interpolates the data to correct wavelength shifts.
   - Produces smile-corrected spectra for all pixels.

5. **Logging and Outputs**:
   - Logs execution times for each step to an output file for performance analysis.
   - Saves intermediate and final results to compressed `.npz` files.

Dependencies:
- **Modules**:
  - `GenerateSmileFns`: Functions for generating artificial SMILE shifts.
  - `Quantificationfns`: Functions for calculating spectral angles and generating SRFs.
  - `Correctionfns`: Functions for applying SMILE correction to the datacube.
  - `extract_from_MODTRAN`: Extracts MODTRAN reference data for comparisons.
  - Additional utility modules for spectral response testing, spline interpolation, and spectral analysis.
- **Inputs**:
  - Hyperspectral data cube (Indian Pine dataset).
  - MODTRAN reference data.
  - Configurations from a `config` file.

Outputs:
- Cropped MODTRAN data.
- Column-averaged spectra.
- Reference and test spectra.
- Spectral angle results.
- Smile-corrected datacube.

Notes:
- The script heavily relies on external modules and global configurations.
- Several steps involve saving intermediate data to allow for step-by-step debugging.
- Certain inputs are configurable via the `config` file.

"""

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

if feature is not None:
    feature_start, feature_end = get_feature_index(wavelength, feature)
    Reference_wl = Reference_wl[feature_start:feature_end]
    Reference_data = Reference_data[feature_start:feature_end]
    wavelength_input = wavelength[feature_start:feature_end]

# Function imports
from create_ref_and_test_spectra import create_ref_and_test_spectra
from data_matrix_collapse import data_matrix_collapse
from min_sa_shift import min_sa_shift
from spectral_angle_calculation import spectral_angle_calculation
from spectral_response import test_spectral_response
from smile_correction import smile_correction
from spline_interpolation import  spline_interpolation_all


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

    # Step 7: Spline interpolation of test spectra.
    step_start = timeit.default_timer()
    spectra_wav, spectra_rad = spline_interpolation_all(full_indian_pine_array, indian_pine_wavelength, wavelength_increment)
    print("Step 7 Done, no issues.")
    np.savez_compressed(f'{data_folder_path}spline_interpolated', wav = spectra_wav, rad = spectra_rad)
    step_end = timeit.default_timer()
    log_message(f"Step 7: Spline interpolation of test spectra - {step_end - step_start:.4f} seconds")

    # Step 8: Generate smile corrected spectra for each pixel. (This turns out to be the most computationally expensive step.)
    step_start = timeit.default_timer()
    interpolated_wl = np.transpose(spectra_wav, ((1, 2, 0)))
    corrected_datacube = smile_correction(spectra_rad, min_spectral_angle, test_spectral_response, np.transpose(cropped_wavelength, ((1, 2, 0))))
    np.savez_compressed(f'{data_folder_path}corrected_datacube', corrected_data = corrected_datacube)
    print(f"Correction complete, no issues. Data saved to {data_folder_path}.")
    step_end = timeit.default_timer()
    log_message(f"Step 8: Generate smile corrected spectra for each pixel - {step_end - step_start:.4f} seconds")

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

# Log total time taken
print(f"Total time taken: {end - start:.4f} seconds")
print(f"Execution timings logged in {log_file_path}")
