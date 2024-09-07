# Author: Shivesh Prakash
# This file applies the keystone distortion artifically to an image

import numpy as np
from scipy.interpolate import CubicSpline

def apply_keystone(img: np.ndarray, p: float) -> np.ndarray:
    """
    Apply keystone correction to an image.

    Args:
        img (np.ndarray): The input image as a 2D NumPy array.
        p (float): The intended subpixel shift.

    Returns:
        np.ndarray: The keystone-corrected image as a 2D NumPy array.
    """
    size = len(img)
    k = 1 + (p * 2 / size)
    x = np.arange(size)

    if k <= 1:
        # Calculate the keystone correction values for interpolation
        lst = np.linspace(size * (1 - k) / 2, size * (1 + k) / 2, size)
        # Perform cubic spline interpolation on all rows
        cs = CubicSpline(x, img, axis=1)
        new_img = cs(lst)
    else:
        # Calculate the number of pixels to pad on each side
        to_pad = int((k - 1) * size / 2)
        # Pad the image with edge values to handle keystone correction
        padded_img = np.pad(img, ((0, 0), (to_pad, to_pad)), mode='edge')
        y = np.arange(len(padded_img[0]))
        # Calculate the keystone correction values for interpolation
        lst = np.linspace(to_pad + (size * (1 - k) / 2), to_pad + (size * (1 + k) / 2), size)
        # Perform cubic spline interpolation on all rows
        cs = CubicSpline(y, padded_img, axis=1)
        new_img = cs(lst)

    return new_img