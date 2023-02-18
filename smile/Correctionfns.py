  # # # Correction
  # Step 7: Apply reverse of Quantified shifts to SRFs. DEPRECATED. STEP 7 IS NOW BUNDLED INTO STEP 9.
  # reverse_shifted_SRFS = reverse_shift(data, demo_SRF, min_spectral_angle) 
  spectra_wav, spectra_rad = spline_interpolation_all(data, wavelength_input, wavelength_increment_input) # Step 8: Spline interpolation of test spectra. INPUTS ALSO MAY NOT BE CORRECT
  # Step 9: Generate smile corrected spectra for each pixel.
  smile_correction(spectra_rad, min_spectral_angle, test_spectral_response)