# Author: Shivesh Prakash
# This file contains the main function to correct a keystone distorted datacube

import numpy as np


def main(
    radianceData,
    adjustment=0.0001,
    threshold=0.005,
    max_iteration=10,
    first_guess_threshold=0.007,
):
    """
    Corrects the keystone distortion in a datacube.

    Args:
        radianceData (np.ndarray): The input datacube with keystone distortion.
        adjustment (float): Adjustment value for bisection.
        threshold (float): Phase shift threshold for termination.
        max_iteration (int): Maximum number of iterations allowed.
        first_guess_threshold (float): Threshold for the initial guess phase shift.

    Returns:
        np.ndarray: The corrected datacube.
    """
    num_bands = radianceData.shape[0]

    # Iterate over each pair of adjacent bands
    for band_index in range(num_bands - 1):
        # Step 1: Get the first guess for the keystone shift
        first_guess_value = first_guess(radianceData, band_index, first_guess_threshold)

        # Step 2: Refine the guess to get the final keystone shift
        final_shift = final_guess(
            radianceData,
            band_index,
            first_guess_value,
            adjustment,
            threshold,
            max_iteration,
        )

        # Step 3: Apply the final keystone correction to the band x+1
        radianceData[band_index + 1] = apply_keystone(
            radianceData[band_index + 1], -final_shift
        )

    return radianceData
