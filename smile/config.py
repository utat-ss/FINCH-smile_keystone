"""This is the config file for defining the constants and global variables. Not to be confused with config_legacy.py, which is the old config file."""

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