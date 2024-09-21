"""This is a rewrite of optical_sensor"""

import numpy as np
import config

class optical_sensor:
    """This class object is used to resample spectra. It is the core of Smile's spectra resampling algorithm, which is used in many places."""

    def __init__(self, data_input, band_number:int, spectral_response_function1, shift_constant):
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
        self.name = f'sensor_{band_number}'
        self.spectral_response_function = spectral_response_function1

        # Compute the band center, band width, and band wavelength
        self.left_index = config.band_index[band_number]
        self.right_index = config.band_index[band_number+1]
        self.band_wavelength = config.wavelength_input[self.left_index:self.right_index]
        band_width = self.band_wavelength[-1] - self.band_wavelength[0]
        self.band_center = (band_width)/2

        # Compute the resampled spectra
        self.intensity = data_input[self.left_index:self.right_index]

        # Compute the statistical weights to assign via the spectral response function
        self.xpos = self.band_wavelength - self.band_center
        self.shift_constant = shift_constant
        self.spectral_response = self.spectral_response_function(self.xpos + self.shift_constant)

        # print(np.shape(self.intensity), np.shape(self.band_wavelength))
        # print(self.intensity, self.band_wavelength)

        try:
            self.output = np.dot(self.intensity, self.spectral_response)
        except ValueError:
            print(self.intensity, self.spectral_response, self.xpos, self.band_wavelength)

        # Generate a concatenated spectral response curve
        self.sr_demo = self.spectral_response_function(np.linspace(min(self.xpos),
                                                                   max(self.xpos),
                                                                   100) + self.shift_constant)

def run_resampling_spectra(data_input, srf_input:list, shift_range:tuple or int, wavelength,
                           show_progress = True):
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
        # min_shift, max_shift = shift_range/wavelength_increment
        min_shift = shift_range[0] / wavelength_increment
        max_shift = shift_range[1] / wavelength_increment
        shift_range = np.linspace(min_shift, max_shift, config.g_num_shifts_1D*2+1)
        # shift_range = np.linspace(min_shift, max_shift, config.g_num_shifts_1D)

    else:
        shift_range = [shift_range]

    for column in range(num_of_columns):
        # Run resampling for each column
        if show_progress:
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

            for bands in range(config.num_of_bands):
                if srf_input is list:
                    srf = srf_input[bands]
                else:
                    srf = srf_input

                try:
                    single_sensor = optical_sensor(data_temp, bands, srf, shift)
                except IndexError:
                    break
                single_sensor.shift_constant = shift

                sampled_spectra.append(single_sensor.output)
                sensor_pos.append(single_sensor.band_center)

                # srf_band_temp.append(single_sensor.spectral_response_function(x_for_demo))
                srf_band_temp.append(single_sensor.sr_demo)
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

    if num_of_columns == 1:
        return sampled_spectra_shift, sensor_pos, srf_band_shift

    return sampled_spectra_columns, sensor_pos, srf_columns