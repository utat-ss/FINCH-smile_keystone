# Author: Shuhan
"""One script to display the plots for a given coordinate in the datacube."""
# TODO: fix plot naming
# TODO: Figure out wby spectral angle is still wrong.

import config
from plot_config import *

# Step 0: MODTRAN
plot_MODTRAN_data()

# Step 1: Generate the Column Averaged Spectra
plot_column_average_spectra(row_to_plot=0)

# Step 2: Generate the SRFs (there is no need to plot this)
# Step 3, 4: Generate Reference and Test Spectra
plot_resampled_ref_and_test(to_be_plotted=10)

# Step 5: Calculate Spectral Angle
plot_spectral_angle(50) # TODO: fix the x axis
# NOTE: Investigate why the spectral angle is not a smile curve
# Potential causes:
# 1. the reference spectra is not correct (fix this by using MODTRAN)
# 2. mathematical error with the spectral angle calculation

# Step 6: Determine Minimum Spectral Angle
plot_min_sa_shift()

# Final step: 
plot_corrected_datacube(0)