"""This module contains functions for plotting the results of the smile algorithm."""
from matplotlib import pyplot as plt
import numpy as np

from load_datacube_npy import ldn
import main

DataFolder = 'data/TempData/'
PlotFolder = "data/SavedPlots/"
indian_pine_array_filepath = 'data\indian_pine_array.npy'
indian_pine_wavelength_filepath = 'data\indian_pine_wavelength.txt'

indian_pine_array = np.load(indian_pine_array_filepath)
indian_pine_wavelength = np.load(indian_pine_array_filepath)

# # Global Variables
_, _, _, wavelength, _ = ldn(indian_pine_array_filepath, indian_pine_wavelength_filepath)

# column_average_spectra
def plot_column_average_spectra(save=True, row_to_plot = None):
    """Plots column average spectra. By default, it plots the whole picture of the 2D collapse result. However, it can also plot a single row of the 2D collapse result if plot_row is specified.
    Args: 
        save (bool): If True, saves the plot to the current directory.
        row_to_plot (int): If specified, plots a single row of the 2D collapse result.
    Returns:
        None
    """
    column_average_spectra = np.load(f'{DataFolder}column_averaged_spectra.npy')

    fig, ax = plt.subplots(1, 1, figsize=(9, 7))
    row_text = ''

    if row_to_plot is None:
        ax.imshow(column_average_spectra, aspect='auto')
        ax.set_ylabel('Row Number')
    else:
        ax.plot(wavelength, column_average_spectra[row_to_plot])
        ax.set_ylabel('Radiance')
        row_text = f"_row{row_to_plot}"

    ax.set_xlabel('Wavelength (nm)')

    if save:
        fig.savefig(f"{PlotFolder}ColumnAverageSpectra{row_text}.png")

# create_ref_and_test_spectra
def plot_ref_and_test_spectra(save=True, crop_range=None):
    """
    """
    crop_start, crop_end = crop_range
    band_length = int((max(wavelength) - min(wavelength)) / main.g_num_of_bands)
    crop_band_start = crop_start * band_length
    crop_band_end = (crop_end + 1) * band_length

    # Retrieve data
    ref, test = np.load(f'{DataFolder}ref_and_test_spectra.npy', mmap_mode='r')

    # Unpack data
    reference_spectra, sensors_position_reference, srf_columns_reference = ref
    test_spectra, sensors_position_test, srf_columns_test = test

    # Crop data
    cropped_reference = [i[crop_start:crop_end] for i in reference_spectra]
    cropped_test = [i[crop_start:crop_end] for i in test_spectra]

    # Plot data
    fig, ref_and_test_plot = plt.subplots(1, 1, figsize=(9, 7))

