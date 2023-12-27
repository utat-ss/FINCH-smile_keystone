"""This module contains functions for plotting the results of the smile algorithm."""
from matplotlib import pyplot as plt
import numpy as np

from load_datacube_npy import ldn
import main
import config

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
    return np.array(to_be_stretched) * np.mean(np.diff(target))/np.mean(np.diff(target)) + min(target)


# # Plotting Functions
# 0. Plot MODTRAN data
def plot_MODTRAN_data(crop_range = None, save=True):
    """Plots MODTRAN data. By default, it plots the whole picture of the MODTRAN data. However, it can also plot a cropped version of the MODTRAN data if crop_range is specified.

    Args:
        crop_range (tuple): If specified, plots a cropped version of the MODTRAN data. The tuple should be in the form (min_wavelength, max_wavelength).
        save (bool): If True, saves the plot to the current directory.
    
    Outputs:
        A plot of the MODTRAN data.
    """
    #TODO: gneerate an overlay plot of the calibration feature wavelength range for MODTRAN vs the data. The goal is to check alignment. 
    MODTRAN_data = np.load(f'{DataFolder}MODTRAN_data.npz', mmap_mode='r')['MODTRAN_data']
    MODTRAN_wavelength = np.load(f'{DataFolder}MODTRAN_data.npz', mmap_mode='r')['MODTRAN_wl']

    fig, MODTRAN_plot = plt.subplots(1, 1, figsize=image_size)
    MODTRAN_plot.set_xlabel('Wavelength (nm)', fontsize = 15)
    MODTRAN_plot.set_ylabel('Radiance', fontsize = 15)
    MODTRAN_plot.grid()

    if crop_range is None:
        MODTRAN_plot.plot(MODTRAN_wavelength, MODTRAN_data, label='MODTRAN Data')
        MODTRAN_plot.set_title("Full MODTRAN Data", fontsize = 15)
    else:
        crop_start_index = np.argmin(np.abs(MODTRAN_wavelength - crop_range[0]))
        crop_end_index = np.argmin(np.abs(MODTRAN_wavelength - crop_range[1]))

        MODTRAN_plot.plot(MODTRAN_wavelength[crop_start_index:crop_end_index], MODTRAN_data[crop_start_index:crop_end_index], label='MODTRAN Data')

        MODTRAN_plot.set_title(f"MODTRAN Data from {crop_range[0]} to {crop_range[1]} nm", fontsize = 15)

    fig.tight_layout()
    return_msg = ''
    if save:
        filename = f"{PlotFolder}MODTRAN_data.png"
        fig.savefig(filename)
        return_msg += f"Saved image to {filename}"

    print(f'MODTRAN data plotting done. {return_msg}')


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
    ax.set_title(f"Column Average Spectra row {row_text}")
    ax.grid()

    fig.tight_layout()

    return_msg = ''
    if save:
        filename = f"{PlotFolder}ColumnAverageSpectra{row_text}.png"
        fig.savefig(filename)
        return_msg += f"Saved image to {filename}"
    print(f'column_average_spectra plotting done. {return_msg}')

# 2. create_ref_and_test_spectra
def plot_resampled_ref_and_test(to_be_plotted=None,
                                crop_range=None,
                                save=True):
    """Plots resampled reference spectra along with the shifts.
    Args:
        wl (1D array): The wavelength array.
        crop_range (tuple): If specified, crops the plot to the specified range.
        to_be_plotted (str) or (int) or None: If input is 'refrence', plots the reference spectra.
        If input is any integer within the range of the number of columns, plots the test spectra of the specified column.
        If input is None, plots nothing.
        save (bool): If True, saves the plot to PlotFolder.
    Returns:
        None
    """
    # # Retrieve data
    input_data = np.load(f'{DataFolder}column_averaged_spectra.npz',mmap_mode='r')['cas']
    reference_and_test = np.load(f'{DataFolder}ref_and_test_spectra.npz', allow_pickle=True)

    # # Unpack data
    reference_, test_ = reference_and_test['ref'], reference_and_test['test']
    reference_ = reference_.tolist()
    test_ = test_.tolist()

    ref_spectra = reference_['spectra']
    srf_columns_ref = reference_['srf']

    test_spectra = test_['spectra']
    srf_columns_test = test_['srf']

    # # Catch the error before the rest of the code runs
    if to_be_plotted is None or int(to_be_plotted) >= len(test_spectra):
        raise ValueError("Input {to_be_plotted} is not valid."
                         f"Please input an integer between 0 and {len(test_spectra)-1}.")

    if main.Reference_data is None:
        spectra_name = 'first column'

    else:
        spectra_name = 'MODTRAN' # Artifact from using MODTRAN data.

    # # Set up a plot
    fig, ref_plot = plt.subplots(1, 1, figsize=spectrum_size)

    ref_plot.set_xlabel('Wavelength (nm)', fontsize = 15)
    ref_plot.set_ylabel('Radiance', fontsize = 15)
    ref_plot.grid()

    if crop_range is not None:
        ref_plot.set_xlim(crop_range)


    # # Plot the un-resampled reference spectra from main.py
    sensors_position_ref = np.linspace(min(main.Reference_wl),
                                        max(main.Reference_wl),
                                        len(ref_spectra[config.g_num_shifts_1D]))
    sensors_position_test = np.linspace(min(main.Reference_wl),
                                        max(main.Reference_wl),
                                        len(test_spectra[to_be_plotted][0]))
    
    ref_plot.plot(main.Reference_wl,
                  main.Reference_data,
                  label="Un-resampled Reference")
    
    # # Plot the resampled reference spectra
    ref_plot.scatter(sensors_position_ref,
                  ref_spectra[config.g_num_shifts_1D],
                  label="Resampled Reference")

    # # Plot the resampled test spectra
    ref_plot.scatter(sensors_position_test,
                  test_spectra[to_be_plotted][0],
                  label="Resampled Test")

    # # Plot the SRFs
    srf_x = np.linspace(min(main.Reference_wl),
                        max(main.Reference_wl),
                        np.shape(srf_columns_test)[2])
    ref_plot.plot(srf_x, srf_columns_test[to_be_plotted][0], label="SRF Reference")


    fig_title = "Overlay of Reference and Test Spectra" + '\n'
    fig_title += f"Test Spectra: {to_be_plotted}th column" + '\n' 
    fig_title += f"Reference Spectra: {spectra_name}"
    ref_plot.legend(bbox_to_anchor=(1, 1), loc='upper left')
    ref_plot.set_title(fig_title,
                       fontsize = 15)
    fig.tight_layout()

    return_msg = ''
    if save:
        img_title = f"TestSpectraColumn{to_be_plotted}"
        filename = f"{PlotFolder}{img_title}.png"
        fig.savefig(filename)
        return_msg = f"Saved image to {PlotFolder}{img_title}.png"

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
def plot_min_sa_shift(save=True):
    """Plots the minimum spectral angle for each column.
    Args:
        save (bool): If True, saves the plot to PlotFolder.
    """
    # Retrieve data
    min_sa = np.load(f'{DataFolder}min_spectral_angle.npz', allow_pickle=True)['msa']

    # Plot data
    fig, min_sa_plot = plt.subplots(1, 1, figsize=image_size)

    # TODO: fix the hard code
    min_sa_plot.plot(min_sa, label='Minimum Spectral Angle')
    min_sa_plot.set_xlabel('Columns', fontsize = 15)
    min_sa_plot.set_ylabel('Shifts that minimize the spectral angle', fontsize = 15)
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

        datacube_shape = np.shape(corrected_datacube)
        img_center = (int(datacube_shape[1]/2), int(datacube_shape[2]/2))

        corrected_datacube_plot.imshow(corrected_datacube[slice_number][0:-1][0:-10], vmin=0, vmax=corrected_datacube[0][img_center[0]][img_center[1]])
        corrected_datacube_plot.set_xlabel('x [pixel]', fontsize = 15)
        corrected_datacube_plot.set_ylabel('y [pixel]', fontsize = 15)
        corrected_datacube_plot.set_title(f"Target's satellite image at {round(wavelength[slice_number], 2)}nm", fontsize=15)

        location = f"{round(wavelength[slice_number], 2)}nm"

    elif slice_number is None and spatial_coordinate is not None:
        # Plot the corrected datacube at a specific spatial coordinate
        corrected_wavelength = np.linspace(min(wavelength), max(wavelength), config.g_num_of_bands)

        fig, corrected_datacube_plot = plt.subplots(1, 1, figsize=spectrum_size)

        corrected_datacube_plot.plot(corrected_wavelength, corrected_datacube[:, spatial_coordinate[0], spatial_coordinate[1]])
        corrected_datacube_plot.set_xlabel('Wavelength [nm]', fontsize = 15)
        corrected_datacube_plot.set_ylabel('Reflectance', fontsize = 15)
        corrected_datacube_plot.set_title(f"Corrected reflectance at ({spatial_coordinate[0]}px, {spatial_coordinate[1]}px)", fontsize=15)
        fig.colorbar()

        location = f"({spatial_coordinate[0]}px, {spatial_coordinate[1]}px)"

    # Error handling
    elif slice_number is None and spatial_coordinate is None:
        raise UserWarning("Error: slice_number and spatial_coordinate are both None. Please specify what to be plotted.")

    elif slice_number is not None and spatial_coordinate is not None:
        raise UserWarning("Error: can not plot both a slice and a spatial coordinate. Please specify what to be plotted.")

    fig.tight_layout()

    return_msg = ''
    if save:
        filename = f"{PlotFolder}CorrectedDatacube{location}.png"
        fig.savefig(filename)
        return_msg = f"Saved image to {filename}"

    print(f"corrected_datacube plotting done. {return_msg}")
