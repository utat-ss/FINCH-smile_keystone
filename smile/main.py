from config import *

# NEED FUNCTION IMPORTS. CURRENTLY NOT WORKING.

#Andy
if __name__ == '__main__':
  # This block is currently intended to be a skeleton setup, calling all relevant functions in order. Inputs will likely need to be changed as circumstances do.

  data = radianceData
  wavelength_input = wavelength
  wavelength_increment_input = wavelength_increment
  # Globals further defined here in case we want to pass in other data

  # # # Quantification
  # Quantification seems to work fine, the code runs properly.
  column_averaged_spectra = data_matrix_collapse(data) # Step 1: Generate Column Averaged Spectra.
  demo_SRF = test_spectral_response([100]) # Step 2: Generate SRFs. 100 IS A PLACEHOLDER VALUE
  ref_spectra, test_spectra = create_ref_and_test_spectra((0,1000)) # Step 3, 4: Generate Reference and Test spectra. (0,100) is also a placeholder
  sa_deg = spectral_angle_calculation(test_spectra, ref_spectra) # Step 5: Calculate spectral angle from test and reference spectra.
  min_spectral_angle = determine_min_sa(sa_deg) # Step 6: Determine minimum spectral angle.

  # # # Correction
  # Double check inputs; should be functional. Confirm inputs to be correct first and foremost.
  # Step 7: Apply reverse of Quantified shifts to SRFs. DEPRECATED. STEP 7 IS NOW BUNDLED INTO STEP 9.
  # reverse_shifted_SRFS = reverse_shift(data, demo_SRF, min_spectral_angle) 
  spectra_wav, spectra_rad = spline_interpolation_all(data, wavelength_input, wavelength_increment_input) # Step 8: Spline interpolation of test spectra. INPUTS ALSO MAY NOT BE CORRECT
  # Step 9: Generate smile corrected spectra for each pixel.
  smile_correction(spectra_rad, min_spectral_angle, test_spectral_response)