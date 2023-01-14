# Author: Julia 
# Step 5 of Smile. This file computes the spectral angle from the test and reference spectra

def spectral_angle_calculation(test_spectra, ref_spectra, plot_col = -1):
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
          sa_deg[num_of_col][num_of_shift] = np.degrees(np.arccos(np.dot(test_spectra[num_of_col][:], ref_spectra[num_of_shift][:])/(np.linalg.norm(test_spectra[num_of_col][:])*np.linalg.norm(ref_spectra[num_of_shift][:]))))
    if plot_col == -1:
      pass
    else:
      plt.plot(sa_deg[plot_col], marker='+',linestyle='None');plt.show()
      print("This is column ", plot_col)
      print(sa_deg[plot_col])
      
    return sa_deg
