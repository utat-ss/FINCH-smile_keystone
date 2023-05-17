"""This is the config file for defining the constants and global variables. Not to be confused with config_legacy.py, which is the old config file."""
import numpy as np

# File paths
data_folder_path = 'data/TempData/'
indian_pine_array_filepath = 'data\indian_pine_array.npy'
indian_pine_wavelength_filepath = 'data\indian_pine_wavelength.txt'

Reference_data_filepath = "data/MODTRANdata.json"

# Global Vars
g_num_of_bands = 70
g_num_shifts_1D = 5
g_shift_increment = .2 #nanometers
g_total_shifts = int(2*(g_num_shifts_1D/g_shift_increment)+1)

# Feature range
feature = (1300, 1500) #nm
# feature = None

# Make the wavelength array, rather than using the provided .txt file.
with open(indian_pine_wavelength_filepath, 'r') as f:
    wavelength_temp = [float(i) for i in f.read().splitlines()]

wavelength = np.linspace(min(wavelength_temp), max(wavelength_temp), 200) # NOTE: bandaid solution, 200 should be np.shape(indian_pine_array)[1], but that can't be obtained without loading the array first.

# Some constants that must be calculated
def get_feature_index(wavelength_source, feature_range):
    """Finds the index of the feature range in the wavelength_source array.
    
    Args:
        wavelength_source (array): array of wavelengths
        feature_range (tuple): tuple of the start and end of the feature range in nm
    
    Returns:
        feature_index (tuple): tuple of the start and end index of the feature range
    """
    feature_start, feature_end = feature_range
    
    # start_index = np.where(abs(wavelength_source - feature_start) == min(abs(wavelength_source - feature_start)))[0][0]
    # # end_index = np.where(abs(wavelength_source - feature_end) == min(abs(wavelength_source - feature_end)))[0][0]
    # diff_start, diff_end = abs(wavelength_source - feature_start), abs(wavelength_source - feature_end)
    diff_start, diff_end = [abs(i - feature_start) for i in wavelength_source], [abs(i - feature_end) for i in wavelength_source]
    start_index = diff_start.index(min(diff_start))
    end_index = diff_end.index(min(diff_end))

    return(start_index, end_index)

# Resampled wavelength to create an array indicating the spectral boudnaries of each band
band_bounds = np.round(np.linspace(min(wavelength), max(wavelength), g_num_of_bands+1), 2)

if feature is not None:
    # wavelength indices covering the feature range
    feature_index_wavelength = get_feature_index(wavelength, feature)
    # band numbers covering the feature range
    feature_index_band = get_feature_index(band_bounds, feature)
    # number of bands in the feature
    feature_num_of_bands = feature_index_band[1] - feature_index_band[0] + 1

    # If there is specfied feature range, overrite it over default parameters
    num_of_bands = feature_num_of_bands
    wavelength_input = wavelength[feature_index_wavelength[0]:feature_index_wavelength[1]+1]
else:
    num_of_bands = g_num_of_bands
    wavelength_input = wavelength

# For each number in config.resampled_wavelength, find the index of the closest value in config.wavelength_input
band_index = [np.where(abs(wavelength_input - i) == min(abs(wavelength_input - i)))[0][0] for i in band_bounds]