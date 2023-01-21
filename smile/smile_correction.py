# Author: Shuhan and Andy
# The final step of Smile. This file applies the reverse of the calculated shift to apply smile correction. 

def smile_correction(interpolated_data, calculated_shift, srf):
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
        """Take out a slice of the interpolated data, such that all the data came from the same column.
        data_slice starts with shape (rows, spectral)
        """
        clear_output(wait=True)
        print (f"Progress: {i}/{data_shape[2]}")
        data_slice = np.transpose(interpolated_data[:, :, i]) 

        output_temp = []
        
        for j in range(len(data_slice)):
            """For each row, identify a shift constant and isolate a spectra from a column"""
            shift_constant = calculated_shift[j]
            isolated_spectra = data_slice[j]

            corrected_spectra, _, _ = run_resampling_spectra(isolated_spectra, srf, shift_constant, show_progress=False)
            output_temp.append(corrected_spectra[0])
        
        corrected_spectra_collection.append(output_temp)

    corrected_data_cube = np.transpose(corrected_spectra_collection, (2, 1, 0))

    return corrected_data_cube
