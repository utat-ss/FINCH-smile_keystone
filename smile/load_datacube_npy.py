import numpy as np
def ldn(datacube_filepath, wavelength_filepath):
    """
    Load and process a datacube and wavelength data to convert them to radiance units.

    Parameters:
    datacube_filepath (str): File path to the datacube .npy file, containing hyperspectral data in a 3D array (rows x columns x spectral bands).
    wavelength_filepath (str): File path to the wavelength .txt file, containing the spectral wavelength data.

    Returns:
    tuple: A tuple containing:
        - wavelength_source (np.array): The original wavelength array from the wavelength file, excluding the first 20 entries.
        - radiance_data (np.array): Transformed radiance data from the datacube, reshaped to spectral bands x rows x columns and scaled for radiance conversion.
        - data_dim (tuple): Dimensions of the data array after transformation, in the form (spectral bands, rows, columns).
        - wavelength (np.array): Resampled wavelength array, using a linear spacing based on the original wavelength file's range.
        - wavelength_increment (float): The increment value between consecutive entries in the resampled wavelength array.
    """
    data = np.load(datacube_filepath)
    wavelength_source = np.loadtxt(wavelength_filepath)[20:]
    data = np.transpose(data, (2,0,1))
    radiance_data = (data-1000)/500
    data_dim = data.shape #this is the dimension of the dataset we are using as reference data. 
        #Since we do not have real test data, we will generate 'test spectra' by resampling the reference spectra
        # g_data_dim[0] = spectral, g_data_dim[1] = row, g_data_dim[2] = col
  
    # Hot patch to fix the original wavelength.txt file's indexing problem
    wavelength = np.linspace(min(wavelength_source), max(wavelength_source), len(wavelength_source))
    wavelength_increment = wavelength[1] - wavelength[0]
  
    return(wavelength_source, radiance_data, data_dim, wavelength, wavelength_increment)
