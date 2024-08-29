# Author: Shivesh Prakash
# This file contains code to generate the first guess for the keystone subpixel shift using a particular band of the hyperspectral cube.

import numpy as np
import cv2
from apply_keystone import apply_keystone
from compute_sinc_fit_peak import compute_sinc_fit_peak
from run_polynomial_regression import run_polynomial_regression
from invert_polynomial_regression import invert_polynomial_regression
from detect_sift_features import detect_sift_features

def first_guess(radianceData, band_index, threshold=0.006):
    """
    Generates an initial estimate of the keystone subpixel shift between two adjacent bands in a hyperspectral image cube.

    Parameters:
    radianceData (numpy.ndarray): The hyperspectral data cube, where each band is a 2D array of spectral values.
    band_index (int): The index of the current band in the hyperspectral cube.
    threshold (float): The maximum allowable error in phase shift estimation. If the estimated shift is within this
                       threshold, the function returns the shift. Otherwise, further refinement is performed using SIFT.

    Returns:
    float: The estimated phase shift between the two bands, refined using SIFT if necessary.
    """

    # Get the current band (imgx) and the next band (imgx1) from the hyperspectral cube.
    imgx = radianceData[band_index]
    imgx1 = radianceData[band_index + 1]

    # Step 1: Initial Guess Using Phase Shift
    # Degrade the next band (imgx1) with different keystone functions to simulate potential shifts.
    degraded_images = [apply_keystone(imgx1, p) for p in np.linspace(-0.05, 0.05, 10)]

    # Compute the phase shift between the current band and each degraded image using a sinc fit.
    phase_shifts = [
        compute_sinc_fit_peak(imgx, degraded_img) for degraded_img in degraded_images
    ]

    # Perform polynomial regression to relate the simulated keystone functions to the observed phase shifts.
    regression_model = run_polynomial_regression(
        np.linspace(-0.05, 0.05, 10), phase_shifts
    )

    # Invert the regression model to estimate the keystone shift between the original bands.
    estimated_pixel_shift = invert_polynomial_regression(
        regression_model, compute_sinc_fit_peak(imgx, imgx1)
    )

    # Apply the estimated keystone shift to correct the next band.
    fixed_x1 = apply_keystone(imgx1, -estimated_pixel_shift)

    # Compute the phase shift between the original band and the corrected band.
    estimated_phase_shift = compute_sinc_fit_peak(imgx, fixed_x1)

    # Step 2: Check if the estimated phase shift is within the acceptable threshold
    if abs(estimated_phase_shift) <= threshold:
        return estimated_phase_shift

    # Step 3: SIFT Method
    # If the initial estimate is not within the threshold, use SIFT to refine the estimate.
    points1, points2 = detect_sift_features(imgx, imgx1)

    # Perform polynomial regression on the matched SIFT feature points to predict the keystone shift.
    sift_regression_model = run_polynomial_regression(points1[:, 0], points2[:, 0])

    # Use the SIFT-based regression model to predict the pixel shift.
    predicted_pixel_shift = sift_regression_model(0)

    # Apply the predicted keystone shift to correct the next band.
    sift_fixed_x1 = apply_keystone(imgx1, -predicted_pixel_shift)

    # Compute the phase shift between the original band and the SIFT-corrected band.
    sift_estimated_phase_shift = compute_sinc_fit_peak(imgx, sift_fixed_x1)

    return sift_estimated_phase_shift
