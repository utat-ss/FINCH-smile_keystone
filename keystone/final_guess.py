# Author: Shivesh Prakash
# This file contains code to generate the final guess.

from apply_keystone import apply_keystone
from compute_sinc_fit_peak import compute_sinc_fit_peak
from coarse_bisection_recursive import coarse_bisection_recursive


def final_guess(
    radianceData,
    band_index,
    first_guess,
    adjustment=0.0001,
    threshold=0.005,
    max_iteration=500,
):
    """
    Estimates the final keystone shift value using coarse bisection method.

    Args:
        radianceData: The radiance data containing spectral images.
        band_index: The index of the band to process.
        first_guess: The initial guess for the keystone shift.
        adjustment: Adjustment value for bisection.
        threshold: Phase shift threshold for termination.
        max_iteration: Maximum number of iterations allowed.

    Returns:
        float: The final keystone shift value.
    """
    # Get images for band x and band x+1
    imgx = radianceData[band_index]
    imgx1 = radianceData[band_index + 1]

    # Compute initial phase shift for the first guess
    initial_adjusted_img = apply_keystone(imgx1, -first_guess)
    initial_shift = compute_sinc_fit_peak(imgx, initial_adjusted_img)

    # Refine the guess using recursive coarse bisection
    final_guess = coarse_bisection_recursive(
        imgx,
        imgx1,
        -first_guess,
        initial_shift,
        threshold,
        0,
        max_iteration,
        adjustment,
    )

    return final_guess
