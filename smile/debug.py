"""This is a temporary file that contains code snippets that are used for testing purposes. 

DO NOT IMPORT THIS FILE INTO ANY OTHER FILES.
"""

import config
import test
import numpy as np
from matplotlib import pyplot as plt
import main
import load_datacube_npy

wavelength_source, radianceData, g_data_dim, wavelength, wavelength_increment = load_datacube_npy.ldn(config.indian_pine_array_filepath, config.indian_pine_wavelength_filepath)

reference_and_test =np.load(f'{config.data_folder_path}ref_and_test_spectra.npz', allow_pickle=True)

# # Unpack data
reference_, test_ = reference_and_test['ref'], reference_and_test['test']


sa_deg = main.spectral_angle_calculation(test_.item()["spectra"], reference_.item()['spectra'], g_data_dim)


