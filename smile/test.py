"""This is a rewrite of optical_sensor"""

import numpy as np
import config

class optical_sensor:
    """This class object is used to resample spectra. It is the core of Smile's spectra resampling algorithm, which is used in many places."""

    def __init__(self, data_input, sensor_number:int, spectral_response_function1, shift_constant):
        """Note: domain of each sensor's spectral response function (SRF) is 
            currently unknown, here I'm assuming that it's based on the physical 
            boundaries limited by each sensor's size (i.e. no blooming or 
            smearing)

        Args: 
            sensor_number: from smaller wavelengths to larger wavelengths, from 
                0 to N.
            spectral_response_function: a mathematical function. The spectral 
                response function of this individual sensor

        Output: 
            self.output: resampled spectra [1D array with length = # of sensors]
            self.position: x position of each sensor's center[1D array with 
                length = # of sensors]
            self.spectral_response: the spectral response curve [1D array with 
                length = # of sensors]
        """
        self.name = f'sensor_{sensor_number}'
        self.spectral_response_function = spectral_response_function1

        # The left and right bounds define tha range of of sensor in terms of the actual wavelength in nm
        self.leftbound = config.resampled_wavelength[sensor_number]
        self.rightbound = config.resampled_wavelength[sensor_number + 1]
        

        self.shift_constant = shift_constant

