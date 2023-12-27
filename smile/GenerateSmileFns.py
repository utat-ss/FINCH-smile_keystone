# Author: Shivesh Prakash
# This file contains functions to generate artificial SMILE shift in the radianceData
# This file is only for internal testing purposes, it is not to be added to the algorithm


from scipy.interpolate import CubicSpline
import numpy as np


def generate_shift_val(x, maxShift=2):
    """This function calculates the shift corresponding to each column.
    It acts as an helper function for generate_smile_shift().

    Args:
        x (int): The column number
        maxShift (int, optional): The maximum shift value. Defaults to 2.

    Returns:
        float: The shift value
    """
    a = -(maxShift / 672.4)
    b = 144 * maxShift / 6724
    c = 1540 * (maxShift / 6724)
    return (a * (x**2)) + (b * x) + c


def generate_smile_shift(data):
    """This function induces artificial SMILE shift in the radianceData.
    It uses the generate_shift_val() function to calculate the shift value for each column.
    It then uses CubicSpline interpolation to generate the shifted data.

    Args:
        data (numpy array): The radianceData

    Returns:
        numpy array: The shifted radianceData
    """
    num_cols = np.shape(data)[1]
    num_rows = np.shape(data)[2]
    x = [t for t in range(num_cols)]
    for col in range(num_cols):
        shift_subpixel = generate_shift_val(col)
        for row in range(num_rows):
            cs = CubicSpline(x, data[:][row][col])
            data[:][row][col] = [cs(x + shift_subpixel) for x in x]
    return data
