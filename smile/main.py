from config import *

# main is operational; but there's a *lot* we could be doing to make this cleaner. Todo. - Andy

# Function imports
from create_ref_and_test_spectra import create_ref_and_test_spectra
from data_matrix_collapse import data_matrix_collapse
from determine_min_sa import determine_min_sa
from load_datacube_npy import ldn
import optical_sensor
from spectral_angle_calculation import spectral_angle_calculation
from spectral_response import test_spectral_response, make_random_SRFs
from smile_correction import smile_correction
from spline_interpolation import  spline_interpolation_1_pixel, spline_interpolation_all

print("Imports completed, no issues.")

#Andy
if __name__ == '__main__':
  # skeleton setup, currently nonfunctional. Continuing to work on inputs.

  data = radianceData
  wavelength_input = wavelength
  wavelength_increment_input = wavelength_increment
  # Globals further defined here in case we want to pass in other data

  # # # Quantification
  column_averaged_spectra = data_matrix_collapse(data) # Step 1: Generate Column Averaged Spectra.
  demo_SRF = test_spectral_response # Step 2: Generate SRFs.
  ref_spectra, test_spectra = create_ref_and_test_spectra((0,1000), column_averaged_spectra, demo_SRF) # Step 3, 4: Generate Reference and Test spectra. (0,100) is also a placeholder
  sa_deg = spectral_angle_calculation(test_spectra, ref_spectra) # Step 5: Calculate spectral angle from test and reference spectra.
  min_spectral_angle = determine_min_sa(sa_deg) # Step 6: Determine minimum spectral angle.

  print("Quantification complete, no issues.")

  # # # Correction
  # Step 7: Apply reverse of Quantified shifts to SRFs. DEPRECATED. STEP 7 IS NOW BUNDLED INTO STEP 9.
  # reverse_shifted_SRFS = reverse_shift(data, demo_SRF, min_spectral_angle) 
  spectra_wav, spectra_rad = spline_interpolation_all(data, wavelength_input, wavelength_increment_input) # Step 8: Spline interpolation of test spectra. INPUTS ALSO MAY NOT BE CORRECT
  # Step 9: Generate smile corrected spectra for each pixel.
  smile_correction(spectra_rad, min_spectral_angle, test_spectral_response)

  # The function runs through from beginning to end. Whether it's functional is yet to be seen.