# Author: Shuhan
# This file defines the optical_sensor class object, which is initially defined in step 3 of Smile.
# Following optical_sensor class object is the run_resampling_spectra function. This function executes the functionalities provided by optical_sensor, and is the only function that does so.

# ===IMPORTANT:===
# It is necessary to import run_resampling_spectra while importing this file.

import numpy as np
from matplotlib import pyplot as plt
from IPython.display import clear_output
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
        g_num_of_bands = 70
        #data_input = data_for_resampling[col_num]

        # from shorter wavelengths to longer wavelengths, assign each sensor 
        # with a name, just to tell them apart
        self.name = f'sensor_{sensor_number}'
        # each sensor may have a different SRF, or they may not. 
        # SRFs are mathematical expressions describing how sensors respond to 
        # photons landing on different parts of the sensor (refer to how hyperspecral cameras work)
        self.spectral_response_function = spectral_response_function1 #changed to 1 to avoid dupes

        # x_left_bound and x_right_bound are the given CCD's spectral range in terms of the indices of the raw data array.
        # x_left_bound and x_right_bound are the single sensor's spectral range.
        # For convenience's sake, they are in terms of indices of the raw data array (0 ~ N for len(data) == N)
        sensor_width = len(data_input) / g_num_of_bands

        self.x_left_bound = int(sensor_number * sensor_width)
        self.x_right_bound = int((1 + sensor_number) * sensor_width)
        
        # An arbitrary x axis for the SRF, origin is seated on the left bound of the sensor in question, ends at the sensor's right bound
        x_axis = np.arange(0, self.x_right_bound - self.x_left_bound)
        x_axis_demo = np.linspace(0, self.x_right_bound - self.x_left_bound, 100)

        # The spectral response curve
        self.shift_constant = shift_constant
        self.spectral_response = self.spectral_response_function(x_axis + self.shift_constant)
        self.spectral_response_demo = self.spectral_response_function(x_axis_demo + self.shift_constant)

        # Now that we have the spectral response curve, we simply need to do the dot product between the spectral response curve and the corresponding 
        # segment of the actual data
        # The dot product = the sensor's total output, since in Python, the  spectral response curve and the corresponding segment of data are 1xN matrices of the same length.
        self.raw_intensity = data_input[self.x_left_bound: self.x_right_bound]
        self.output = np.dot(self.raw_intensity, self.spectral_response)
        self.position = 0.5 * (self.x_left_bound + self.x_right_bound)

def run_resampling_spectra(data_input, srf_input:list, shift_range:tuple or int, g_num_of_bands, g_num_of_shifts_1D, wavelength, show_plots = False, show_progress = True):
    """
    One function that runs it all. If the SRF for each sensor is unique, compile
        them into a list in a low -> high wavelength order; if all sensors have 
        the same SRF, input that SRF. 
    
    Shift is 0 by default, enter a different value if it is not 0. 

    Args: 
        data_input = the 2D collapsed data array for column-by-column resampling. 
            Essentially output of Andy's code. There can also be only one column
        srf_input = a list of spectral resposne functions for all the sensors.
            If SRFs for all sensors are the same, then just enter one function.
        shift_range = a tuple indicating the minimum and maximum shifts. For 
            example, shift_range = (-5, 5) will produce shifts -5, -4, ..., 4, 5. 
            Alternatively, if you want a single shift, just put a single value
        show_plots = if true, shows plot. False by default. 
        show_progress = if true, shows progress. True by default

    Outputs: 
        output: resampled spectra, a 1D array with the length = the # of sensors
        sensor_pos: a 1D array indicating the centers of each sensor
        spectral_response_demo: a 1D array that is the combination of spectral
            response curves of each sensor
    """
    sampled_spectra_columns = []
    sensor_pos_columns = []
    srf_columns = []

    wavelength_increment = wavelength[1] - wavelength[0]
    input_shape = np.shape(data_input)

    # If data_input is 1D, num_of_columns = 1;
    # if data_input is 2D, num_of_columns = input_shape[0]
    if len(input_shape) == 1:
        num_of_columns = 1
    else:
        num_of_columns = input_shape[0]

    if isinstance(shift_range, tuple):
        min_shift, max_shift = shift_range/wavelength_increment
        shift_range = np.linspace(min_shift, max_shift, g_num_of_shifts_1D)

    else:
        shift_range = [shift_range]

    for column in range(num_of_columns):
        # Run resampling for each column
        if show_progress:
            clear_output(wait=True)
            print(f'Working through column {column}/{num_of_columns}')

        if num_of_columns > 1:
            # If there are more than 1 columns, pick out the column
            data_temp = data_input[column]
        else:
            # If there is only one column, the data is the column
            data_temp = data_input

        sampled_spectra_shift = []
        sensor_pos_shift = []
        srf_band_shift = []

        for shift in shift_range:
            # optical_sensor.shift_constant = shift
            sampled_spectra = []
            sensor_pos = []
            srf_band_temp = []

            for bands in range(g_num_of_bands):
                if srf_input is list:
                    srf = srf_input[bands]
                else:
                    srf = srf_input

                single_sensor = optical_sensor(data_temp, bands, srf, shift)
                single_sensor.shift_constant = shift

                sampled_spectra.append(single_sensor.output)
                sensor_pos.append(single_sensor.position)

                x_for_demo = np.arange(0, 100)
                srf_band_temp.append(single_sensor.spectral_response_function(x_for_demo))
                srf_band = np.concatenate(srf_band_temp)

            sampled_spectra_shift.append(sampled_spectra)
            sensor_pos_shift.append(sensor_pos)
            srf_band_shift.append(srf_band)

        sampled_spectra_shift = np.array(sampled_spectra_shift)
        sensor_pos_shift = np.array(sensor_pos_shift)
        srf_band_shift = np.array(srf_band_shift)

        sampled_spectra_columns.append(sampled_spectra_shift)
        sensor_pos_columns.append(sensor_pos_shift)
        srf_columns.append(srf_band_shift)

    sampled_spectra_columns = np.array(sampled_spectra_columns)
    sensor_pos_columns = np.array(sensor_pos_columns)
    srf_columns = np.array(srf_columns)

    # # Deprecated code for plotting
    # if show_plots:
    #     if num_of_columns > 1:
    #         print ("YouError: too many columns! This is on you.")

    #     else:
    #         fig, plot_for_show = plt.subplots(1, 1, figsize=(15, 7))

    #         plot_for_show.plot(wavelength, data_input, label = 'Input data')
    #         plot_for_show.set_xlabel('Wavelength [nm]')
    #         plot_for_show.set_ylabel('Radiance')

    #         for i in range(len(sampled_spectra_shift)): 
    #             plot_for_show.scatter(stretch_horizontal(sensor_pos_shift, wavelength), sampled_spectra_shift, label = f'Resampled band {i}')
    #             plot_for_show.plot(np.linspace(min(wavelength), max(wavelength), np.shape(srf_band_shift)[1]), srf_band_shift[i])

    #         plot_for_show.legend()
    #         fig.tight_layout()

    if num_of_columns == 1:
        return sampled_spectra_shift, sensor_pos, srf_band_shift

    return sampled_spectra_columns, sensor_pos, srf_columns
