import numpy as np
def ldn(datacube_filepath, wavelength_filepath):
  data = np.load(datacube_filepath)
  if data.ndim != 3: # Ensure data is 3D before transposing
    raise ValueError("Datacube must be a 3D numpy array, but got {data.shape}")
    
  wavelength_source = np.loadtxt(wavelength_filepath)[20:]
  if len(wavelength_source) <= 20:
    raise ValueError("Wavelength file must have at least 20 entries, but got {len(wavelength_source)}") 
    
  wavelength_source = wavelength_source[20:]
  if np.any(np.diff(wavelength_source) < 0):
    wavelength_source = np.sort(wavelength_source)
  data = np.transpose(data, (2,0,1))
  
  radiance_data = (data-1000)/500
  
  data_dim = data.shape #this is the dimension of the dataset we are using as reference data. 
      #Since we do not have real test data, we will generate 'test spectra' by resampling the reference spectra
      # g_data_dim[0] = spectral, g_data_dim[1] = row, g_data_dim[2] = col
  # Hot patch to fix the original wavelength.txt file's indexing problem
  wavelength = np.linspace(min(wavelength_source), max(wavelength_source), len(wavelength_source))
  wavelength_increment = wavelength[1] - wavelength[0]
  
  return(wavelength_source, radiance_data, data_dim, wavelength, wavelength_increment)
