"""This is a temporary file that contains code snippets that are used for testing purposes. 

DO NOT IMPORT THIS FILE INTO ANY OTHER FILES.
"""

import config
import test
import numpy as np
from matplotlib import pyplot as plt
import main

min_spectral_angle = np.load(config.data_folder_path + "min_spectral_angle.npz")["msa"]

print(np.shape(min_spectral_angle))