"""This is a temporary file that contains code snippets that are used for testing purposes. 

DO NOT IMPORT THIS FILE INTO ANY OTHER FILES.
"""

import config
import test
import numpy as np
from matplotlib import pyplot as plt

wavelength_index = [np.where(abs(config.wavelength_input - i) == min(abs(config.wavelength_input - i)))[0][0] for i in config.resampled_wavelength]

output = [config.wavelength_input[i] for i in wavelength_index]

print(output)