# Author: Shuhan
# One script to display the plots for a given coordinate in the datacube.

from plot_config import *

# Step 1: Generate the Column Averaged Spectra
plot_column_average_spectra()

# Step 2: Generate the SRFs (there is no need to plot this)

# Step 3, 4: Generate Reference and Test Spectra
plot_resampled_ref_and_test(wavelength, to_be_plotted='reference')
plot_resampled_ref_and_test(wavelength, to_be_plotted=10)

# Step 5: Calculate Spectral Angle
plot_spectral_angle()

# Step 6: Determine Minimum Spectral Angle
plot_min_sa()