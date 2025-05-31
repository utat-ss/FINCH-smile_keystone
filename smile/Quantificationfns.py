# # # Quantification

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

from optical_sensor import *
import config

# Author: Andy
# Step 1 of Smile. Takes a 3-dimensional array of size a*b*c, and returns a 2-dimensional array of size b*c
#   with values equal to the averages of the columns of dimension a. In other words, replaces dimension a with
#   the mean of its values, collapsing the 3-D matrix into a 2-D one.
def data_matrix_collapse(image_file):
    """
    Takes a 3-dimensional array of shape (bands, rows, cols) and returns a 2-dimensional array
    of shape (bands, cols) by averaging over the rows dimensions. In other words, for each
    spectral band and column, it computes the mean across all rows, producing a spectrum for
    every column that is averaged vertically.

    Args:
        3-D array of shape (bands, rows, cols)
    Returns:
        2-D  array of shape (bands, cols) (or (cols, bands) if transposed)
    """
    # Get dimensions from the input image 
    AdimSpectral = len(image_file)    # Number of bands
    AdimLength = len(image_file[0])   # Number of rows
    AdimWidth = len(image_file[0][0]) # Number of columns

    # Read the input image and create a copy aimg with the same dimensions
    Aimg = image_file

    # Create an empty array to hold column-averaged spectra
    # Shape: (bands, cols)
    k = np.zeros((AdimSpectral, AdimWidth))

    # For each spectral band (slice), extract the 2D spatial image
    for band in range(AdimSpectral):
        b = np.reshape(Aimg[band, :, :], (AdimLength, AdimWidth))
        
        # Average over rows to get column-averaged values for this band
        collapsed_b = np.zeros(AdimWidth)
        for col_index in range(AdimWidth):
            collapsed_b[col_index] = np.mean(b[:, col_index])

        k[band, :] = collapsed_b

    # Tranpose to obtain shape (cols, bands)
    k_transpose = np.transpose(k)

    return k_transpose

# Author: Shuhan
# This fn contains mathematical functions that are currently placeholders for the sensor's spectral response functions. 
def test_spectral_response(x, mu = None, delta = None):
    """
    A test spectral response function. Basically a Gaussian function centered 
        at mu * (max(x) - min(x)) with std = sigma

    Args: 
        x: a 1D array

    Output: 
        A Gaussian function corresponding to the inxpxut. A 1D array
    """
    if mu is None:
        mu = x[int(0.5 * len(x))]

    if delta is None:
        delta = (x.max() - x.min()) / 10

    sigma = delta / (2 * np.sqrt(2 * np.log(2)))

    # Add a buffer on both sides of x in case the gaussian is at the edge
    buffer_width = 3 * sigma
    left_buffer = np.linspace(x.min() - buffer_width, x.min(), 50, endpoint=False)
    right_buffer = np.linspace(x.max(), x.max() + buffer_width, 50, endpoint=False)
    x_extended = np.concatenate([left_buffer, x, right_buffer])

    gaussian = stats.norm.pdf(x_extended, mu, sigma)
    normed_gaussian = gaussian / gaussian.max()
    
    indices = np.searchsorted(x_extended, x)
    normed_gaussian = normed_gaussian[indices]

    return normed_gaussian

# Author: Shuhan
# Step 3 of Smile. This fn creates reference and test spectra from the provided data cube
def create_ref_and_test_spectra(data_for_resampling:list, test_spectral_response, wavelength, ref_spectra = None):
    """
    Returns reference and test spectra. 

    Args:
        data_for_resampling: the data that is being resampled and processed for SMILE correction;
        wavelength: a 1D array of wavelengths;
        test_spectral_response: a mathematical function simulating each pixel's spectral response;
        g_num_shifts_1D: global variable, number of shifts in 1D;
        g_shift_increment: global variable, shift increment;
        reference_spectra: spot for an external reference data column. Without user input, the code automatically sources a data column from the data cube.

    Outputs:
        sampled_reference: cropped resampled reference spectra, which was originally generated using the first colum. 
            Its shape is (# of shifts, # of bands)
        sampled_test: cropped resampled test spectra. Its shape is (# of columns, # of shifts, # of bands)
    """

    shift_bound = config.g_num_shifts_1D * config.g_shift_increment
    shift_range = (-shift_bound, shift_bound)

    no_reference = False

    if ref_spectra is None:
        no_reference = True
        ref_spectra = data_for_resampling[0]

    if config.feature is tuple:
        feature_start, feature_end = config.feature
        ref_spectra = ref_spectra[:, feature_start:feature_end]
        data_for_resampling = data_for_resampling[:, feature_start:feature_end]

    sampled_ref_spec, ref_pos, ref_srf = run_resampling_spectra(ref_spectra,
                                                                test_spectral_response,
                                                                shift_range, wavelength,
                                                                data_is_feature=True)
    
    sampled_test_spec, test_pos, test_srf = run_resampling_spectra(data_for_resampling,
                                                                   test_spectral_response,
                                                                   0, wavelength,
                                                                   data_is_feature=True)
    
    # TODO: normalzie one spectra to the other
    
    
    sampled_reference = {"spectra": sampled_ref_spec, "pos": ref_pos, "srf": ref_srf}
    sampled_test = {"spectra": sampled_test_spec, "pos": test_pos, "srf": test_srf}
    if no_reference:
        print('Reference spectra not found, used first column by default.')

    return sampled_reference, sampled_test

# Author: Julia
# Step 5 of Smile. This fn computes the spectral angle from the test and reference spectra

def spectral_angle_calculation(test_spectra, ref_spectra, g_data_dim, plot_col=-1):
    """ 
    Calculates the spectral angle
    
    Arguments:
    test_spectra: shape = (number of columns, number of cropped bands). This array represents the collection of the dot products 
    of response functions and hyperspectral image data. 

    ref_spectra: shape = (number of shifts, number of cropped bands). This will be used to be compared with test spectra to calculate the 
    spectra angle. 
    
    Returns:
        sa_deg: (number of columns, number of shifts) np-array containing the spectral angle
    
    """
    # calculating the SA angle in degrees
    sa_deg = np.zeros((g_data_dim[2],np.shape(ref_spectra)[0]))
    for num_of_col in range(g_data_dim[2]):
        for num_of_shift, _ in enumerate (ref_spectra):
            dot1 = test_spectra[num_of_col][:]
            dot2 = ref_spectra[num_of_shift][:]

            norm1 = np.linalg.norm(test_spectra[num_of_col][:])
            norm2 = np.linalg.norm(ref_spectra[num_of_shift][:])
 
            dot = np.dot(dot1, dot2) / (norm1 * norm2)
            sa_deg[num_of_col][num_of_shift] = np.degrees(np.arccos(dot))

    if plot_col == -1:
        pass
    else:
        plt.plot(sa_deg[plot_col], marker='+',linestyle='None')
        plt.show()
        print("This is column ", plot_col)
        print(sa_deg[plot_col])

    return sa_deg

# Author: Rediet
# Step 6 of Smile. This fn calculates the minimum spectral angle.
def determine_min_sa_shift(sa_deg, g_data_dim):
    """Calculates the minimum spectral angle
        Variables used: 
        min_each_row[g_data_dim[2]] holds the minimum value of each row in sa_deg (comparing between shifts for each column of the data)
        min_col_num[g_data_dim[2]] holds the colmumn number of the minimum values of each row in sa_deg

    Args: 
        sa_deg: A 2D matrix (g_data_dim[2], g_total_shifts)containing SA in degrees 

    Returns: 
        min_col_num: A 1D matrix of size g_data_dim[2] containing the shifts of the best matched spectrum in sa_deg
    """
    # Create a 1d array of size g_data_dim[2] conta`ining the smallest value in each row of sa_deg (the spectral angle in degrees)
    # Create another 1D array of size g_data_dim[2] containing the column number of the smallest value in each row
    shift_bound = config.g_num_shifts_1D * config.g_shift_increment
    shift_range = np.linspace(-shift_bound, shift_bound, config.g_num_shifts_1D * 2 + 1)

    min_sa_value = np.zeros(g_data_dim[2])
    best_wavelength_shifts = np.zeros(g_data_dim[2])

    for i in range(g_data_dim[2]):
        col = sa_deg[i]
        # finding the first occurance of the index of min_each_row[i]
        # min_col_num[i] = np.where(col == min_row)[0][0]
        min_index = np.argmin(col)
        min_sa_value[i] = col[min_index]
        best_wavelength_shifts[i] = shift_range[min_index]

    return best_wavelength_shifts, min_sa_value
