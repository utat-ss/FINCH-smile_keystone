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

# # Plotting configs
image_size = (9, 9)
spectrum_size = (25, 7)

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
# 0. Plot MODTRAN data

# 1. column_average_spectra
def plot_column_average_spectra(save=True, row_to_plot = None):
    """Plots column average spectra. By default, it plots the whole picture of the 2D collapse result. However, it can also plot a single row of the 2D collapse result if plot_row is specified.
    Args: 
        save (bool): If True, saves the plot to the current directory.
        row_to_plot (int): If specified, plots a single row of the 2D collapse result.
    Returns:
        None
    """
    column_average_spectra = np.load(f'{DataFolder}column_averaged_spectra.npz', mmap_mode='r')['cas']

    fig, ax = plt.subplots(1, 1, figsize=spectrum_size)
    row_text = ''

    if row_to_plot is None:
        ax.imshow(column_average_spectra, aspect='auto')
        ax.set_ylabel('Row Number')
    else:
        ax.plot(wavelength, column_average_spectra[row_to_plot])
        ax.set_ylabel('Radiance')
        row_text = f"_row{row_to_plot}"

    ax.set_xlabel('Wavelength (nm)')

    return_msg = ''
    if save:
        filename = f"{PlotFolder}ColumnAverageSpectra{row_text}.png"
        fig.savefig(filename)
        return_msg += f"Saved image to {filename}"
    print(f'column_average_spectra plotting done. {return_msg}')

# 2. create_ref_and_test_spectra
def plot_resampled_ref_and_test(wl = wavelength, to_be_plotted=None, crop_range=None, save=True):
    """Plots resampled reference spectra along with the shifts.
    Args:
        wl (1D array): The wavelength array.
        crop_range (tuple): If specified, crops the plot to the specified range.
        to_be_plotted (str) or (int) or None: If input is 'refrence', plots the reference spectra. If input is any integer within the range of the number of columns, plots the test spectra of the specified column. If input is None, plots nothing.
        save (bool): If True, saves the plot to PlotFolder.
    Returns:
        None§
    """
    shift_extent = main.g_num_shifts_1D * main.g_shift_increment
    shift_range = np.arange(-shift_extent, shift_extent, main.g_shift_increment)

    if crop_range is not None:
        crop_start, crop_end = crop_range
        band_length = round(len(wl) / main.g_num_of_bands)
        crop_wl_start = (crop_start * band_length)
        crop_wl_end = ((crop_end) * band_length)

    # # Retrieve data
    input_data = np.load(f'{DataFolder}column_averaged_spectra.npz',mmap_mode='r')['cas']
    reference_and_test = np.load(f'{DataFolder}ref_and_test_spectra.npz', allow_pickle=True)

    # # Unpack data
    reference_, test_ = reference_and_test['ref'], reference_and_test['test']
    reference_spectra, sensors_position_reference, srf_columns_reference = reference_
    test_spectra, sensors_position_test, srf_columns_test = test_

    # # Plot data
    fig, reference_plot = plt.subplots(1, 1, figsize=spectrum_size)

    reference_plot.set_xlabel('Wavelength (nm)', fontsize = 15)
    reference_plot.set_ylabel('Radiance', fontsize = 15)
    reference_plot.grid()

    # Plot the original data
    if main.MODTRAN_data is None:
        spectra_name = 'first column'

        ref_data = input_data[0]
    
    else:
        spectra_name = 'MODTRAN'

        ref_data = main.MODTRAN_data

    if crop_range is None:
        reference_plot.plot(wl, ref_data, label=spectra_name)
    else:
        reference_plot.plot(wl[crop_wl_start:crop_wl_end], ref_data[crop_wl_start:crop_wl_end], label=f"Reference Spectra ({spectra_name})")

    fig_title = 'empty'

    # Plot the results
    if to_be_plotted == 'reference':
        sensors_positions = stretch_horizontal(sensors_position_reference, wl)
        fig_title = f"ReferenceSpectra{2 * main.g_num_shifts_1D}Shifts"

        for i, ref_shift in enumerate(reference_spectra):
            # Determine the plot's range
            SRF_x = np.linspace(min(wl), max(wl), main.g_num_of_bands*100)

            if crop_range is None:
                reference_plot.scatter(sensors_positions, ref_shift, label=f'shift {round(shift_range[i], 2)}')

                # Plot the SRF
                reference_plot.plot(SRF_x, srf_columns_reference[i], label = f"SRF with shift {round(shift_range[i], 2)}")

            else:
                reference_plot.scatter(sensors_positions[crop_start:crop_end], ref_shift[crop_start:crop_end], label=f'shift {round(shift_range[i], 2)}')

                # Plot the SRF
                reference_plot.plot(SRF_x[crop_start * 100:crop_end * 100], srf_columns_reference[i][crop_start * 100:crop_end * 100], label = f"SRF with shift {round(shift_range[i], 2)}")

    elif to_be_plotted < len(test_spectra): # TODO: Edit this whole thing. Make sure it's up to date as the other sections.
        sensors_positions = stretch_horizontal(sensors_position_test, wl)
        fig_title = f"TestSpectraColumn {to_be_plotted}"

        reference_plot.scatter(sensors_positions, test_spectra[to_be_plotted], label=f'Resampled test spectra of Column {to_be_plotted}')
        SRF_x = np.linspace(min(wl), max(wl), main.g_num_of_bands*100)
        reference_plot.plot(SRF_x, srf_columns_test[to_be_plotted][0], label=f'SRF of Column {to_be_plotted}')
        

    elif to_be_plotted >= len(test_spectra):
        raise UserWarning(f"Error: to_be_plotted is {to_be_plotted}, which is greater than the number of columns in the test spectra ({len(test_spectra)}).")
    
    elif to_be_plotted is None:
        raise UserWarning("Error: to_be_plotted is None. Please specify what to be plotted.")
    
    reference_plot.legend(bbox_to_anchor=(1, 1), loc='upper left')
    reference_plot.set_title(fig_title + "\n Reference: " + spectra_name, fontsize = 15)
    fig.tight_layout()

    return_msg = ''
    if save:
        filename = f"{PlotFolder}{fig_title}.png"
        fig.savefig(filename)
        return_msg = f"Saved image to {PlotFolder}{fig_title}.png"

    print(f"resampled_ref_and_test plotting done. {return_msg}")


# 3. spectral angle calculation
def plot_spectral_angle(column_to_plot = None, save=True):
    """Plots the spectral angle between the reference and test spectra.
    Args:
        column_to_plot (int): The column to plot the spectral angle of. If None, plots the all columns in the same figure.
        save (bool): If True, saves the plot to PlotFolder.
    Returns:
        None
    """
    # Retrieve data
    spectral_angle = np.load(f'{DataFolder}sa_deg.npz', allow_pickle=True)['sa_deg']

    # Plot data
    fig, spectral_angle_plot = plt.subplots(1, 1, figsize=image_size)

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

    return_msg = ''
    if save:
        filename = f"{PlotFolder}SpectralAngle.png"
        fig.savefig(filename)
        return_msg = f"Saved image to {filename}"

    print(f"spectral_angle plotting done. {return_msg}")

# 4. determine minimum spectral angle
def plot_min_sa(save=True):
    """Plots the minimum spectral angle for each column.
    Args:
        save (bool): If True, saves the plot to PlotFolder.
    """
    # Retrieve data
    min_sa = np.load(f'{DataFolder}min_spectral_angle.npz', allow_pickle=True)['msa']

    # Plot data
    fig, min_sa_plot = plt.subplots(1, 1, figsize=image_size)

    min_sa_plot.plot(min_sa, label='Minimum Spectral Angle')
    min_sa_plot.set_xlabel('Columns', fontsize = 15)
    min_sa_plot.set_ylabel('Spectral Angle (deg)', fontsize = 15)
    min_sa_plot.set_title("Minimum Spectral Angle of all Columns")

    min_sa_plot.grid()
    min_sa_plot.legend()
    fig.tight_layout()

    return_msg = ''
    if save:
        filename = f"{PlotFolder}MinimunSpectralAngles.png"
        fig.savefig(filename)
        return_msg = f"Saved image to {filename}"

    print(f"min_sa plotting done. {return_msg}")

# corrected datacube
def plot_corrected_datacube(slice_number=None, spatial_coordinate=None, save=True):
    """Plots the corrected datacube at a specific wavelength or spatial coordinate.
    Args:
        slice_number (int): The wavelength index to plot. If None, plots the datacube at a specific spatial coordinate.
        spatial_coordinate (tuple): The spatial coordinate to plot. If None, plots the datacube at a specific wavelength. The coordinate is in the form (row, column).
        save (bool): If True, saves the plot to PlotFolder.
    Returns:
        None
    """
    # Retrieve data
    corrected_datacube = np.load(f'{DataFolder}corrected_datacube.npz', allow_pickle=True)['corrected_data']

    # Plot data
    if slice_number is not None and spatial_coordinate is None:
        # Plot the datacube at a specific wavelength
        fig, corrected_datacube_plot = plt.subplots(1, 1, figsize=image_size)

        corrected_datacube_plot.imshow(corrected_datacube[0])
        corrected_datacube_plot.set_xlabel('x [pixel]', fontsize = 15)
        corrected_datacube_plot.set_ylabel('y [pixel]', fontsize = 15)
        corrected_datacube_plot.set_title(f"Target's satellite image at {round(wavelength[slice_number], 2)}nm", fontsize=15)

    elif slice_number is None and spatial_coordinate is not None:
        # Plot the corrected datacube at a specific spatial coordinate
        corrected_wavelength = np.linspace(min(wavelength), max(wavelength), main.g_num_of_bands)

        fig, corrected_datacube_plot = plt.subplots(1, 1, figsize=spectrum_size)

        corrected_datacube_plot.plot(corrected_wavelength, corrected_datacube[:, spatial_coordinate[0], spatial_coordinate[1]])
        corrected_datacube_plot.set_xlabel('Wavelength [nm]', fontsize = 15)
        corrected_datacube_plot.set_ylabel('Reflectance', fontsize = 15)
        corrected_datacube_plot.set_title(f"Corrected reflectance at ({spatial_coordinate[0]}px, {spatial_coordinate[1]}px)", fontsize=15)

    # Error handling
    elif slice_number is None and spatial_coordinate is None:
        raise UserWarning("Error: slice_number and spatial_coordinate are both None. Please specify what to be plotted.")

    elif slice_number is not None and spatial_coordinate is not None:
        raise UserWarning("Error: can not plot both a slice and a spatial coordinate. Please specify what to be plotted.")

    fig.tight_layout()

    if save:
        fig.savefig(f"{PlotFolder}CorrectedDatacube.png")

    print("corrected_datacube plotting done")
