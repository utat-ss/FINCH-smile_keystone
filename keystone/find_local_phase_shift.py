#Author: Shivesh Prakash

import numpy as np

def find_local_phase_shift(x1: np.ndarray, x3: np.ndarray) -> int:
    """
    Calculates the local phase shift between two images along the x-axis using cross-correlation.

    Args:
    - x1: np.ndarray - Input image 1
    - x3: np.ndarray - Input image 2

    Returns:
    - int: The calculated x-shift representing the local phase shift between the images
    """

    # Compute cross-correlation using Fourier Transform
    r = fourier_t(x1, x3)

    # Find the index of the maximum correlation value
    max_row, max_col = np.unravel_index(np.argmax(r), r.shape)

    # Calculate the x-shift from the correlation peak
    # Shift relative to the center of the cross-correlation result
    x_shift = max_col - (r.shape[1] // 2)

    return int(x_shift)
