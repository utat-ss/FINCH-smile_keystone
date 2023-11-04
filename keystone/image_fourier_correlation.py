# Author: Shivesh Prakash, Juttada Satya Sathwik
# This file returns the result of cross-phase correlation

import numpy as np

def image_fourier_correlation(x1: np.ndarray, x3: np.ndarray) -> np.ndarray:
    """
    Perform Fourier Transform based cross-phase correlation of two images.
    
    Args:
        x1 (np.ndarray): The first input image.
        x3 (np.ndarray): The second input image.
        
    Returns:
        np.ndarray: The result of the cross-phase correlation in the spatial domain.
    """
    # Compute the Fourier Transform of the input images
    F_x1 = np.fft.fft2(x1)
    F_x3 = np.fft.fft2(x3)
    
    # Compute the complex conjugate of the Fourier Transform of x3
    Gconj = np.conj(F_x3)
    
    # Calculate the Cross Phase Correlation
    R = (F_x1 * Gconj) / np.abs(F_x1 * Gconj)
    
    # Shift the result back to the real space
    r = np.fft.fftshift(np.fft.ifft2(R))
    
    return r
