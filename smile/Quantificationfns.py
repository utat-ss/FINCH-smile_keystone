# # # Quantification

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

from optical_sensor import *

# Author: Andy
# Step 1 of Smile. Takes a 3-dimensional array of size a*b*c, and returns a 2-dimensional array of size b*c
#   with values equal to the averages of the columns of dimension a. In other words, replaces dimension a with
#   the mean of its values, collapsing the 3-D matrix into a 2-D one.
def data_matrix_collapse(image_file):
    ''' 
    Takes a 3-dimensional array of size a*b*c, and returns a 2-dimensional array of size b*c
    with values equal to the averages of the columns of dimension a. In other words, replaces
    dimension a with the mean of its values, collapsing the 3-D matrix into a 2-D one. 

    Args: 
      3-D array of size a*b*c
    Returns: 
      2-D array of size a*c (Or size c*a if transposed)
    '''
    AdimSpectral = len(image_file) # spectral
    AdimLength = len(image_file[0]) # row
    AdimWidth = len(image_file[0][0]) # column

    # These are not expected values, but I used them to verify the output dimensions.

    # Read the input image and create a copy aimg with the same dimensions.
    Aimg = image_file

    # Create a new empty matrix k, which has one dimension less than aimg (Instead of being spectral x row x            col, it is spectral x row only)
    k = np.zeros((AdimSpectral, AdimWidth))
    

    # For each "slice"  of aimg, set b to be equal to a 2D version of that slice (isolating that slice).
    for slice in range(AdimSpectral):
        b = np.reshape(Aimg[slice, :, :], (AdimLength, AdimWidth))
        # In b, take the average of all spectral values in the spectral dimension
        collapsed_b = np.zeros(AdimWidth)
        for smallerslice in range(AdimWidth):
            collapsed_b[smallerslice] = np.mean(b[smallerslice, :])
        k[slice, :] = collapsed_b

    # Then, set the corresponding cell row in k to be equal to b.
    # Repeat the previous three steps for each of the initial slices.
    # In essence, k takes the averages of the each wavelength per row and "collapses" them into size slices of column averaged spectra.
    # Finally, tranpose k to get L, size 5x10.
    k_transpose = np.transpose(k) # k_transpose is unused but may be useful
    return k_transpose

# Author: Shuhan
# This fn contains mathematical functions that are currently placeholders for the sensor's spectral response functions. 
def test_spectral_response(x):
    """
    A test spectral response function. Basically a Gaussian function centered 
        at mu * (max(x) - min(x)) with std = sigma

    Args: 
        x: a 1D array

    Output: 
        A Gaussian function corresponding to the inxpxut. A 1D array
    """
    sigma = 0.25 * len(x)
    mu = 0.5 * len(x)
    gaussian = stats.norm.pdf(x, mu, sigma)
    normed_gaussian = gaussian / max(gaussian)
    
    return normed_gaussian

def make_random_SRFs(x, g_num_of_bands):
    """
    Generates a list of unique and random Gaussian SRFs for each pixel. The parameters of each
        generated Gaussian function is a random number between 0 and 1. These random numbers 
        follow normal(Gaussian) distribution
    Args:
        x: a 1D array

    Output: 
        A list of Gaussian functions. Interact with them like this: 
            func_list = make_random_SRFs
            
            test_SRF_1 = func_list[1](x's)
            test_SRF_2 = func_list[2](x's)
            ...
    """
    output = []
    for i in range(g_num_of_bands):
        sigma_temp = np.random.normal(loc = 0.5, scale = 0.1) * len(x)
        mu_temp = np.random.normal(loc = 0.5, scale = 0.1) * len(x)

        def temp (x):

            gaussian = stats.norm.pdf(x, mu_temp, sigma_temp)
            normed_gaussian = gaussian / max(gaussian)

            return normed_gaussian
          
        globals()[f'random_SRF_{i}'] = temp
    
        output.append(globals()[f'random_SRF_{i}'])

    return output
    # plt.plot(test_spectral_response(np.arange(0, 100)))
    # plt.plot(test_spectral_response(np.arange(0, 100) + 10))

# Author: Shuhan
# Step 3 of Smile. This fn creates reference and test spectra from the provided data cube
def create_ref_and_test_spectra(data_for_resampling:list, test_spectral_response, wavelength, g_num_of_bands, g_num_shifts_1D, g_shift_increment, reference_spectra = None):
    """
    Returns reference and test spectra. 

    Args:
        data_for_resampling: the data that is being resampled and processed for SMILE correction;
        wavelength: a 1D array of wavelengths;
        test_spectral_response: a mathematical function simulating each pixel's spectral response;
        g_num_of_bands: global variable, number of bands;
        g_num_shifts_1D: global variable, number of shifts in 1D;
        g_shift_increment: global variable, shift increment;
        reference_spectra: spot for an external reference data column. Without user input, the code automatically sources a data column from the data cube.

    Outputs:
        sampled_reference: cropped resampled reference spectra, which was originally generated using the first colum. 
            Its shape is (# of shifts, # of bands)
        sampled_test: cropped resampled test spectra. Its shape is (# of columns, # of shifts, # of bands)
    """
    shift_bound = g_num_shifts_1D * g_shift_increment
    shift_range = (-shift_bound, shift_bound)

    if type(reference_spectra) is not np.array:
        no_reference = True
        ref_spectra = data_for_resampling[80]

    sampled_reference = run_resampling_spectra(ref_spectra, test_spectral_response, shift_range, g_num_of_bands, g_shift_increment, wavelength)

    sampled_test = run_resampling_spectra(data_for_resampling, test_spectral_response, 0, g_num_of_bands, g_shift_increment, wavelength)

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
        for num_of_shift in range(np.shape(ref_spectra)[0]):
            sa_deg[num_of_col][num_of_shift] = np.degrees(np.arccos(np.dot(test_spectra[num_of_col][:],ref_spectra[num_of_shift][:])/(np.linalg.norm(test_spectra[num_of_col][:])*np.linalg.norm(ref_spectra[num_of_shift][:]))))
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
def determine_min_sa(sa_deg, g_data_dim):
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
    min_each_row = np.zeros(g_data_dim[2])
    min_col_num = np.zeros(g_data_dim[2])

    for i in range(0, g_data_dim[2]):
        col = sa_deg[i]
        #seaching each row for the minimum
        min_each_row[i] = min(col)
        #finding the first occurance of the index of min_each_row[i]
        min_col_num[i] = np.where(col == min_each_row[i])[0][0]

        return min_col_num.astype(int)
