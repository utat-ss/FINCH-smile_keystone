# Author: Shivesh Prakash
# This file returns the peak of sinc function

import numpy as np
from scipy.optimize import curve_fit

def compute_sinc_fit_peak(img, p):
    """
    Calculate the maximum amplitude of a fitted sinc function.

    Args:
        img (np.ndarray): The input image as a 2D NumPy array.
        p (float): The intended subpixel shift.

    Returns:
        float: The maximum amplitude of the fitted sinc function.
    """
    x3 = keystone(img, p)
    r = fourier_t(img, x3)
    l = []
    
    n = len(img)
    for i in range(n):
        y = [r[j, i] for j in range(n)]
        l.append(np.mean(y))
    
    r3 = np.array(l)
    x = np.array([x for x in range(- (n // 2), n - (n // 2))])
    popt, pcov = curve_fit(sincfit, x, r3)
    x2 = np.linspace(-10, 10, 500)
    y2 = np.array([sincfit(x, popt[0], popt[1], popt[2]) for x in x2])
    
    return max(y2)