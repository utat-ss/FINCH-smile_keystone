# Author: Shuhan
# This file defines the optical_sensor class object, which is initially defined in step 3 of Smile.

import numpy as np

class optical_sensor:

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
        #data_input = data_for_resampling[col_num]

        # from shorter wavelengths to longer wavelengths, assign each sensor 
        # with a name, just to tell them apart
        self.name = f'sensor_{sensor_number}'
        # each sensor may have a different SRF, or they may not. 
        # SRFs are mathematical expressions describing how sensors respond to 
        # photons landing on different parts of the sensor (refer to how
        # hyperspecral cameras work)
        self.spectral_response_function = spectral_response_function1 #changed to 1 to avoid dupes

        # x_left_bound and x_right_bound are the given CCD's spectral range in
        # terms of the indices of the raw data array.
        # x_left_bound and x_right_bound are the single sensor's spectral range.
        # For convenience's sake, they are in terms of indices of the raw data 
        # array (0 ~ N for len(data) == N)
        sensor_width = len(data_input) / g_num_of_bands

        self.x_left_bound = int(sensor_number * sensor_width)
        self.x_right_bound = int((1 + sensor_number) * sensor_width)
        
        # an arbitrary x axis for the SRF, origin is seated on the left bound of
        # the sensor in question, ends at the sensor's right bound
        x_axis = np.arange(0, self.x_right_bound - self.x_left_bound)
        x_axis_demo = np.linspace(0, self.x_right_bound - self.x_left_bound, 100)

        # The spectral response curve
        self.shift_constant = shift_constant
        self.spectral_response = self.spectral_response_function(x_axis + self.shift_constant)
        self.spectral_response_demo = self.spectral_response_function(x_axis_demo + self.shift_constant)

        # Now that we have the spectral response curve, we simply need to do the
        # dot product between the spectral response curve and the corresponding 
        # segment of the actual data
        # The dot product = the sensor's total output, since in Python, the 
        # spectral response curve and the corresponding segment of data are 1xN 
        # matrices of the same length.
        self.raw_intensity = data_input[self.x_left_bound: self.x_right_bound]
        self.output = np.dot(self.raw_intensity, self.spectral_response)
        self.position = 0.5 * (self.x_left_bound + self.x_right_bound)
