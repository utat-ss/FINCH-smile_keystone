# # # Correction
# Step 7: Apply reverse of Quantified shifts to SRFs. DEPRECATED. STEP 7 IS NOW BUNDLED INTO STEP 9.
# reverse_shifted_SRFS = reverse_shift(data, demo_SRF, min_spectral_angle)
# spectra_wav, spectra_rad = spline_interpolation_all(data, wavelength_input, wavelength_increment_input) # Step 8: Spline interpolation of test spectra. INPUTS ALSO MAY NOT BE CORRECT
# Step 9: Generate smile corrected spectra for each pixel.
# smile_correction(spectra_rad, min_spectral_angle, test_spectral_response)

# Import the necessary modules

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from IPython.display import clear_output

from optical_sensor import run_resampling_spectra

# Edit the following snippet
# Author: REDIET
# Step 8 of Smile. Function that performs spline interpolation on one pixel's spectra.
# The original code also generates plots, but the plot-creating lines of code have been inhibited for now.
# Call function spline_interpolation_all to utilize this file's codes.


def spline_interpolation_1_pixel(test_spectra_rad, test_spectra_wav, interp_step, wavelength, make_plot=False):

    """ Generates spline interpolation of test spectral data
      
        Args:
            test_spectra_rad: radiance values of test data
            test_spectra_wav: corresponding wavelength values of each radiance value
            interp_step: new interval between wavelengths of interpolated spectra

        Variables used:
            tck: coefficients of spline curve (spectral respone function)

        Returns:
            new_test_spectra_wav: arrays containing the interpolated spectra wavelenth values
        new_test_spectra_rad: arrays containing the interpolated spectra radiance values

    """
    #Getting spline coefficients
    tck = interpolate.splrep(test_spectra_wav, test_spectra_rad)

    #Applying interpolation
    new_test_spectra_wav = np.arange(min(wavelength), max(wavelength)+interp_step, interp_step)
    new_test_spectra_rad = interpolate.splev(new_test_spectra_wav, tck)

    #optional plotting
    if make_plot:
        plt.figure(figsize=(20,5))
        plt.plot(test_spectra_wav, test_spectra_rad, marker='o', linestyle='None')
        plt.plot(new_test_spectra_wav, new_test_spectra_rad)
        plt.title('Spline interpolation')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Radiance')
        plt.show()

    return new_test_spectra_wav, new_test_spectra_rad


# Author: Julia
# The following function execute Radiet's spline interpolation function
def spline_interpolation_all(test_spectra_rad_all,test_spectra_wav, interp_step, g_data_dim, wavelength):
    """ Generates spline interpolation of test spectral data for all pixels 
      
        Args:
            test_spectra_rad_all: radiance values of test data, data cube that contain all radiance values for all pixels. Shape = (spectral dimension, num of rows, num of cols)
            test_spectra_wav: corresponding wavelength values of each radiance value
            interp_step: new interval between wavelengths of interpolated spectra

        Returns:
            new_test_spectra_wav_all: arrays containing the interpolated spectra wavelenth values for all pixels, shape = (??? TBD, should be dependent on the step size, num of rows, num of cols)
            new_test_spectra_rad_all: arrays containing the interpolated spectra radiance values for all pixels, shape = (??? TBD should be dependent on the step size, num of rows, num of cols)

    """
    # getting the length of the array for interpolated wavelength and radiance values. This step is needed because we need to figure out the shape of the array before we do anything,
    # the length of the interpolated values array is dependent on the step size and i can't figure it out lol
    dim_1 = np.shape(spline_interpolation_1_pixel(test_spectra_rad_all[:,1, 1], test_spectra_wav, interp_step, wavelength)[0])[0]
    # Creating zeros arrays to store the values of spiline interpolation for all pixels.
    new_test_spectra_wav_all = np.zeros((dim_1, g_data_dim[1], g_data_dim[2]))
    new_test_spectra_rad_all = np.zeros((dim_1, g_data_dim[1], g_data_dim[2]))

    # looping through each pixel by going through each row and column, and implement the spline_interpolation_1_pixel to all pixels.
    for num_of_row in range(g_data_dim[1]):
        for num_of_col in range(g_data_dim[2]):
            new_test_spectra_wav_all[:, num_of_row, num_of_col] = spline_interpolation_1_pixel(test_spectra_rad_all[:,num_of_row, num_of_col], test_spectra_wav, interp_step, wavelength)[0]
            new_test_spectra_rad_all[:, num_of_row, num_of_col] = spline_interpolation_1_pixel(test_spectra_rad_all[:,num_of_row, num_of_col], test_spectra_wav, interp_step, wavelength)[1]

    return new_test_spectra_wav_all, new_test_spectra_rad_all

# Author: Shuhan and Andy
# The final step of Smile. This file applies the reverse of the calculated shift to apply smile correction

def smile_correction(interpolated_data, calculated_shift, srf, g_num_of_bands, g_shift_increment, wavelength):
    """Use run_resampling_spectra with the reverse of the calculated shift to apply smile correction.
    First, take out a 2D array of data that has the same row

    Args:
        interpolated_data: spline-interpolated version of the original spectral data cube. 
        calculated_shift: A 1D array detailing the shift constants to each column. Since we are 
            operating under the assumption that each column has a uniform shift throughout itself, 
            this can remain a 1D array. Consider expanding this in the future.
        srf: spectral response function. 

    Outputs: 
        corrected_data_cube: smile-corrected version of the original data cube, with dimensions of 
            (spectra, row, column)
    """
    # calculated_shift spans in the row direction.
    data_shape = np.shape(interpolated_data)

    if len(calculated_shift) != data_shape[2]:
        raise ValueError("len(calcualted_shift) and np.shape(interpolated_data)[2] don't match!")
    if len(np.shape(srf)) > 1:
        print ("WARNING: Accomodation for multiple distinct spectra response functions has not been implemented yet. Only the first row of srfs are accepted")
        srf = srf[0]

    corrected_spectra_collection = []
    for i in range(data_shape[2]):
        # Take out a slice of the interpolated data, such that all the data came from the same column. data_slice starts with shape (rows, spectral)
        clear_output(wait=True)
        print (f"Progress: {i}/{data_shape[2]}")
        data_slice = np.transpose(interpolated_data[:, :, i])

        output_temp = []

        for count, row in enumerate(data_slice):
            # For each row, identify a shift constant and isolate a spectra from a colum
            shift_constant = calculated_shift[count]
            isolated_spectra = row

            corrected_spectra, _, _ = run_resampling_spectra(isolated_spectra, srf, shift_constant, g_num_of_bands, g_shift_increment, wavelength, show_progress=False)
            output_temp.append(corrected_spectra[0])

        corrected_spectra_collection.append(output_temp)

    corrected_data_cube = np.transpose(corrected_spectra_collection, (2, 1, 0))

    return corrected_data_cube