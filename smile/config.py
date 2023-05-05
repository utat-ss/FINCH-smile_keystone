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
# feature = (1300, 1500) #nm
feature = None

# Some constants that must be calculated
with open(indian_pine_wavelength_filepath, 'r') as f:
    wavelength = [float(i) for i in f.read().splitlines()]

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

resampled_wavelength = np.linspace(min(wavelength), max(wavelength), g_num_of_bands+1)

if feature is not None:
    # wavelength indices covering the feature range
    feature_index_wavelength = get_feature_index(wavelength, feature)
    # band numbers covering the feature range
    feature_index_band = get_feature_index(resampled_wavelength, feature)
    # number of bands in the feature
    feature_num_of_bands = feature_index_band[1] - feature_index_band[0] + 1

    # If there is specfied feature range, overrite it over default parameters
    num_of_bands = feature_num_of_bands
    wavelength_input = resampled_wavelength[feature_index_band[0]:feature_index_band[1]+1]
else:
    num_of_bands = g_num_of_bands
    wavelength_input = wavelength


# NOTE: From this point on, I just need to pass feature bands into optical_sensor. Execute optical_sensor with the feature only.
# NOTE: Trust the result.