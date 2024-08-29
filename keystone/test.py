# Author: Shivesh Prakash
# This file contains code to test our Keystone Algorithm

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error, mean_squared_error
from apply_keystone import apply_keystone
from compute_sinc_fit_peak import compute_sinc_fit_peak
from first_guess import first_guess
from final_guess import final_guess


def apply_keystone_to_datacube(radianceData, start_shift=-0.01, end_shift=0.01):
    """
    Applies a linearly varying keystone distortion to each band of a corrected datacube.

    Parameters:
    radianceData (np.ndarray): The hyperspectral datacube with shape (num_bands, height, width).
    start_shift (float): The starting value for the keystone distortion.
    end_shift (float): The ending value for the keystone distortion.

    Returns:
    np.ndarray: The datacube with keystone distortions applied to each band.
    """
    num_bands = radianceData.shape[0]
    distorted_datacube = np.empty_like(radianceData)

    # Generate keystone shift values that vary linearly across the bands
    shifts = np.linspace(start_shift, end_shift, num_bands)

    # Apply the corresponding keystone shift to each band
    for i in range(num_bands):
        distorted_datacube[i] = apply_keystone(radianceData[i], shifts[i])

    return distorted_datacube


def test(
    radianceData,
    realData,
    adjustment=0.0001,
    threshold=0.005,
    max_iteration=10,
    first_guess_threshold=0.007,
):
    """
    Tests the keystone correction algorithm by comparing the calculated keystone shifts to the known real shifts.

    Parameters:
    radianceData (np.ndarray): The distorted hyperspectral datacube to be corrected.
    realData (np.ndarray): The reference datacube with known correct shifts for comparison.
    adjustment (float): The adjustment step size used in the refinement of the shift estimate.
    threshold (float): The error threshold for determining convergence during the refinement process.
    max_iteration (int): The maximum number of iterations allowed during the refinement process.
    first_guess_threshold (float): The error threshold for accepting the initial guess without further refinement.

    Returns:
    None
    """
    calculated_shifts = [0]  # List to store the shifts calculated by the algorithm
    real_shifts = [0]  # List to store the known real shifts

    num_bands = radianceData.shape[0]

    # Iterate over each pair of adjacent bands in the datacube
    for band_index in range(num_bands - 1):
        # Get the current band and the reference band from the real data
        imgx = radianceData[band_index]
        imgx_real = realData[band_index]

        # Calculate the real shift between the current and next band using sinc fitting
        real_shifts.append(-compute_sinc_fit_peak(imgx_real, imgx))

        # Step 1: Obtain an initial estimate of the keystone shift using the first guess method
        first_guess_value = first_guess(radianceData, band_index, first_guess_threshold)

        # Step 2: Refine the initial guess to obtain the final keystone shift
        final_shift = final_guess(
            radianceData,
            band_index,
            first_guess_value,
            adjustment,
            threshold,
            max_iteration,
        )
        calculated_shifts.append(final_shift)

        # Step 3: Apply the final keystone correction to the next band in the datacube
        radianceData[band_index + 1] = apply_keystone(
            radianceData[band_index + 1], final_shift
        )

    # Calculate error metrics to evaluate the performance of the keystone correction algorithm
    mae = mean_absolute_error(real_shifts, calculated_shifts)
    rmse = np.sqrt(mean_squared_error(real_shifts, calculated_shifts))

    # Print error metrics
    print(f"Mean Absolute Error (MAE): {mae:.6f}")
    print(f"Root Mean Square Error (RMSE): {rmse:.6f}")

    # Plot real vs calculated shifts for visual comparison
    plt.figure(figsize=(10, 6))
    plt.plot(real_shifts, label="Real Shifts", marker="o")
    plt.plot(calculated_shifts, label="Calculated Shifts", marker="x")
    plt.xlabel("Band Index")
    plt.ylabel("Shift Value")
    plt.title("Real vs Calculated Shifts")
    plt.legend()
    plt.grid(True)
    plt.show()
