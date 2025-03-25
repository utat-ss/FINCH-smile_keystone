# Author: Shivesh Prakash, Patuan
# This file contains functions to generate artificial SMILE shift in the radianceData
# This file is only for internal testing purposes, it is not to be added to the algorithm


from scipy.interpolate import CubicSpline
import numpy as np


def generate_shift_val(column_length, row_length, maxShift=None):
    """This function calculates the shift corresponding to each column.
    It acts as an helper function for generate_smile_shift().

    Args:
        column_length (int) : The length of the column
        maxShift (int, optional): The maximum shift value. Defaults to 0.05 x row_length.

    Returns:
        numpy array[float] : List of shift for each column
    """
    if(maxShift is None):
        maxShift = 0.05 * row_length
    column_container = np.array([index_column for index_column in range(column_length)]) # Make list [1, 2, ..., column_length]
    a = - (4 * maxShift / (column_length**2))
    b = column_length
    return a * (column_container * (column_container - b))


def generate_smile_shift(data, wavelength, generateShiftFunction, maxShift = None):
    """This function induces artificial SMILE shift in the radianceData.
    It uses generateShiftFunction in the paramater as the function to generate the shift value for each column.
    It then uses CubicSpline interpolation to generate the shifted data.

    Args:
        data (numpy array): The radianceData with dimension row x column x spectral
        wavelength (numpy array) : The wavelength (spectral band) for the data. Used for interpolation
        generateShiftFunction (function): A function to generate the shift for each column
        maxShift: The maximum of shift on the peak (Center of data for quadratic function)
    Returns:
        numpy array(float32): The shifted radianceData
    """
    if (data.dtype == np.uint16):
        data = data.astype(np.float32) # To avoid overflow while doing interpolation
    

    num_cols = np.shape(data)[1]
    num_rows = np.shape(data)[0]
    shiftValue = generateShiftFunction(num_cols, num_rows, maxShift) # How the shift is generated (Simplest one is quadratic function)
    shiftedData = np.zeros_like(data)

    for col in range(num_cols):
        for row in range(num_rows):
            cs = CubicSpline(wavelength, data[row][col])
            shiftedData[row][col] = cs(wavelength + shiftValue[col])

    return shiftedData