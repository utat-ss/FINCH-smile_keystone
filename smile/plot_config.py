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

# # QoL Functions
def stretch_horizontal(to_be_stretched, target):
    """
    A linear transformation that stretches one array of any size, so it starts
        and ends at the same points as the target array.

    Args: 
        to_be_stretched: the array that needs to be stretched
        target: the array which length to_be_stretched needs to match

    Outputs:
        A 1D array, stretched input array
    """
    return np.array(to_be_stretched) * np.mean(np.diff(target)) + min(target)

# # Plotting Functions
# column_average_spectra
def plot_column_average_spectra(save=True, row_to_plot = None):
    """Plots column average spectra. By default, it plots the whole picture of the 2D collapse result. However, it can also plot a single row of the 2D collapse result if plot_row is specified.
    Args: 
        save (bool): If True, saves the plot to the current directory.
        row_to_plot (int): If specified, plots a single row of the 2D collapse result.
    Returns:
        None
    """
    column_average_spectra = np.load(f'{DataFolder}column_averaged_spectra.npz', mmap_mode='r')['cas']

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
def plot_resampled_ref_spectra(wl = wavelength, save=True, crop_range=None, to_be_plotted=None):
    """Plots resampled reference spectra along with the shifts.
    Args:
        save (bool): If True, saves the plot to the current directory.
        crop_range (tuple): If specified, crops the plot to the specified range.
        to_be_plotted (str) or (int) or None: If input is 'refrence', plots the reference spectra. If input is any integer within the range of the number of columns, plots the test spectra of the specified column. If input is None, plots nothing.
    Returns:
        None
    """
    shift_extent = main.g_num_shifts_1D * main.g_shift_increment
    shift_range = np.arange(-shift_extent, shift_extent, main.g_shift_increment)

    if crop_range is not None:
        crop_start, crop_end = crop_range
        band_length = int((max(wl) - min(wl)) / main.g_num_of_bands)
        crop_wl_start = crop_start * band_length
        crop_wl_end = (crop_end + 1) * band_length

    # # Retrieve data
    input_data = np.load(f'{DataFolder}column_averaged_spectra.npz',mmap_mode='r')['cas']
    reference_and_test = np.load(f'{DataFolder}ref_and_test_spectra.npz', allow_pickle=True)

    # # Unpack data
    # Reference_spectra = (num of shifts, num of bands)
    reference_, test_ = reference_and_test['ref'], reference_and_test['test']
    reference_spectra, sensors_position_reference, srf_columns_reference = reference_
    test_spectra, sensors_position_test, srf_columns_test = test_

    # # Plot data
    fig, reference_plot = plt.subplots(1, 1, figsize=(25, 7))

    reference_plot.set_xlabel('Wavelength (nm)', fontsize = 15)
    reference_plot.set_ylabel('Radiance', fontsize = 15)

    # Plot the original data
    if main.MODTRAN_data is None:
        if crop_range is not None:
            wl = wl[crop_wl_start:crop_wl_end]
            input_data = input_data[crop_wl_start:crop_wl_end]
        reference_plot.plot(wl, input_data[0], label='First Column as Reference Spectra')
        ref_ = "\n Reference: first comlumn"
    else:
        ref_data = main.MODTRAN_data
        if crop_range is not None:
            wl = wl[crop_wl_start:crop_wl_end]
            ref_data = input_data[crop_wl_start:crop_wl_end]
        reference_plot.plot(wl, ref_data, label='RefSerence Spectra from MODTRAN')

    fig_title = 'empty'

    # Plot the results
    if to_be_plotted == 'reference':
        sensors_positions = stretch_horizontal(sensors_position_reference[0], wl)
        fig_title = f"Reference Spectra with {main.g_num_shifts_1D} Shifts"
        for i, ref_shift in enumerate(reference_spectra):
            if crop_range is not None:
                sensors_positions = sensors_positions[crop_start:crop_end]
                ref_shift = ref_shift[crop_start:crop_end]
            reference_plot.scatter(sensors_positions, ref_shift, label=f'shift {shift_range[i]}')
            reference_plot.plot(np.linspace(min(wl), max(wl), main.g_num_of_bands*100), srf_columns_reference[i])

    elif to_be_plotted < len(test_spectra):
        # print(np.shape(sensors_position_test[to_be_plotted][0]), np.shape(test_spectra[to_be_plotted][0]))
        sensors_positions = stretch_horizontal(sensors_position_test[to_be_plotted][0], wl)
        reference_plot.scatter(sensors_positions, test_spectra[to_be_plotted][0], label=f'Test Spectra of Column {to_be_plotted}')
        reference_plot.plot(np.linspace(min(wl), max(wl), main.g_num_of_bands*100), srf_columns_test[to_be_plotted][0], label=f'SRF of Column {to_be_plotted}') # TODO: Crop this
        fig_title = f"Test Spectra of Column {to_be_plotted}"

    elif to_be_plotted >= len(test_spectra):
        raise UserWarning(f"Error: to_be_plotted is {to_be_plotted}, which is greater than the number of columns in the test spectra ({len(test_spectra)}).")
    
    elif to_be_plotted is None:
        raise UserWarning("Error: to_be_plotted is None. Please specify what to be plotted.")
    
    reference_plot.legend()
    reference_plot.set_title(fig_title + ref_, fontsize = 15)
    if save:
        fig.savefig(f"{PlotFolder}RefAndTestSpectra.png")

# spectral angle calculation

def plot_spectral_angle(column_to_plot = None, save=True):
    """Plots the spectral angle between the reference and test spectra.
    Args:
        column_to_plot (int): The column to plot the spectral angle of. If None, plots the all columns in the same figure.
        save (bool): If True, saves the plot to the current directory.
    Returns:
        None
    """
    # Retrieve data
    spectral_angle = np.load(f'{DataFolder}sa_deg.npz', allow_pickle=True)['sa_deg']

    # Plot data
    fig, spectral_angle_plot = plt.subplots(1, 1, figsize=(9, 7))

    if column_to_plot is None:
        for i, sa in enumerate(spectral_angle):
            spectral_angle_plot.plot(sa, label=f'Column {i}')
            columns = "all columns"
    elif column_to_plot < len(spectral_angle):
        spectral_angle_plot.plot(spectral_angle[column_to_plot], label=f'Column {column_to_plot}')
        columns = f"Column {column_to_plot}"
    elif column_to_plot >= len(spectral_angle):
        raise UserWarning(f"Error: column_to_plot is {column_to_plot}, which is greater than the number of columns in the test spectra ({len(spectral_angle)}).")

    spectral_angle_plot.set_xlabel('Shifts', fontsize = 15)
    spectral_angle_plot.set_ylabel('Spectral Angle (deg)', fontsize = 15)
    spectral_angle_plot.set_title(f"Spectral Angle of {columns}")
    spectral_angle_plot.grid()

    if save:
        fig.savefig(f"{PlotFolder}SpectralAngle.png")


# plot_column_average_spectra()
plot_resampled_ref_spectra(save=True, crop_range=(0, 20), to_be_plotted='reference')
# plot_spectral_angle()
