# Author: REDIET
# Step 8 of Smile. Function that performs spline interpolation on one pixel's spectra.
# The original code also generates plots, but the plot-creating lines of code have been inhibited for now.
# Call function spline_interpolation_all to utilize this file's codes. 
from config import *

def spline_interpolation_1_pixel(test_spectra_rad,test_spectra_wav, interp_step):

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
  from scipy import interpolate
  #Getting spline coefficients
  tck = interpolate.splrep(test_spectra_wav, test_spectra_rad)

  #Applying interpolation
  new_test_spectra_wav = np.arange(min(wavelength), max(wavelength)+interp_step, interp_step)
  new_test_spectra_rad = interpolate.splev(new_test_spectra_wav, tck)

  #optional plotting

  # plt.figure(figsize=(20,5))
  # plt.plot(test_spectra_wav, test_spectra_rad, marker='o', linestyle='None')
  # plt.plot(new_test_spectra_wav, new_test_spectra_rad)
  #plt.show()

  return new_test_spectra_wav, new_test_spectra_rad

#testing on one pixel
#x = spline_interpolation_1_pixel(radianceData[:,1, 1], wavelength, 1)[0]
#y = spline_interpolation_1_pixel(radianceData[:,1, 1], wavelength, 1)[1]
#plt.figure(figsize=(20,5))
#plt.plot(wavelength, radianceData[:,1, 1], marker='o', linestyle='None', label='test data points')
#plt.plot(x,y, label='interpolated test data')
#plt.ylabel('radiance'); plt.xlabel=('wavelength')
#plt.legend()
#plt.show()


# Author: Julia
# The following function execute Radiet's spline interpolation function
def spline_interpolation_all(test_spectra_rad_all,test_spectra_wav, interp_step):
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
    dim_1 = np.shape(spline_interpolation_1_pixel(test_spectra_rad_all[:,1, 1], test_spectra_wav, interp_step)[0])[0]
    # Creating zeros arrays to store the values of spiline interpolation for all pixels. 
    new_test_spectra_wav_all = np.zeros((dim_1,g_data_dim[1], g_data_dim[2]))
    new_test_spectra_rad_all = np.zeros((dim_1,g_data_dim[1], g_data_dim[2]))

    # looping through each pixel by going through each row and column, and implement the spline_interpolation_1_pixel to all pixels. 
    for num_of_row in range(g_data_dim[1]):
      for num_of_col in range(g_data_dim[2]):
        new_test_spectra_wav_all[:, num_of_row, num_of_col] = spline_interpolation_1_pixel(test_spectra_rad_all[:,num_of_row, num_of_col], test_spectra_wav, interp_step)[0]
        new_test_spectra_rad_all[:, num_of_row, num_of_col] = spline_interpolation_1_pixel(test_spectra_rad_all[:,num_of_row, num_of_col], test_spectra_wav, interp_step)[1]

    return new_test_spectra_wav_all, new_test_spectra_rad_all
