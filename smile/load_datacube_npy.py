import numpy as np
from config import *
def ldn(datacube_filepath, wavelength_filepath):
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
