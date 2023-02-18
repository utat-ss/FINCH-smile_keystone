# # # Quantification
  column_averaged_spectra = data_matrix_collapse(data) # Step 1: Generate Column Averaged Spectra.
  demo_SRF = test_spectral_response # Step 2: Generate SRFs.
  ref_spectra, test_spectra = create_ref_and_test_spectra((0,1000), column_averaged_spectra, demo_SRF) # Step 3, 4: Generate Reference and Test spectra. (0,100) is also a placeholder
  sa_deg = spectral_angle_calculation(test_spectra, ref_spectra) # Step 5: Calculate spectral angle from test and reference spectra.
  min_spectral_angle = determine_min_sa(sa_deg) # Step 6: Determine minimum spectral angle.




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


