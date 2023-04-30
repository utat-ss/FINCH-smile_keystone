"""This is a temporary file that contains code snippets that are used for testing purposes. 

DO NOT IMPORT THIS FILE INTO ANY OTHER FILES.
"""
from main import *

def get_feature_index(wavelength_source, feature_range):
    """Finds the index of the feature range in the wavelength_source array.
    
    Args:
        wavelength_source (array): array of wavelengths
        feature_range (tuple): tuple of the start and end of the feature range in nm
    
    Returns:
        feature_index (tuple): tuple of the start and end index of the feature range
    """
    feature_start, feature_end = feature_range
    
    start_index = np.where(abs(wavelength_source - feature_start) == min(abs(wavelength_source - feature_start)))[0][0]

    end_index = np.where(abs(wavelength_source - feature_end) == min(abs(wavelength_source - feature_end)))[0][0]

    return(start_index, end_index)
    
data = radianceData
wavelength_input = wavelength
wavelength_increment_input = wavelength_increment

start_index, end_index = get_feature_index(wavelength_input, (1300, 1500))

# Globals further defined here in case we want to pass in other data
# # # Quantification
# Step 1: Generate Column Averaged Spectra.
column_averaged_spectra = data_matrix_collapse(data)

demo_SRF = test_spectral_response

# ref_spectra, test_spectra = create_ref_and_test_spectra(column_averaged_spectra, demo_SRF, wavelength_input, g_num_of_bands, g_num_shifts_1D, g_shift_increment, ref_spectra=Reference_data)

print(column_averaged_spectra.shape)

feature_crop = column_averaged_spectra[:, start_index:end_index]

print(feature_crop.shape)