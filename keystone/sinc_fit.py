# Author: Shivesh Prakash
# This file returns the array of sinc function values

import numpy as np

def sinc_fit(x: np.ndarray, amplitude: float, frequency: float, phase_shift: float) -> np.ndarray:
    """
    Calculate the value of a sinc function for a given array of x values.

    Args:
        x (np.ndarray): Array of input values.
        amplitude (float): Amplitude of the sinc function.
        frequency (float): Frequency of the sinc function.
        phase_shift (float): Phase shift of the sinc function.

    Returns:
        np.ndarray: Array of sinc function values.
    """
    sinc_values = amplitude * np.sinc((frequency * x) - phase_shift)
    return sinc_values
