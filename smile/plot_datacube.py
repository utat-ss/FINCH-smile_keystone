from config import *

from matplotlib import pyplot as plt
import numpy as np

DataFolder = 'data/TempData/'

original_datacube = np.copy(radianceData)
smiled_datacube = np.load(f'{DataFolder}smiled_data.npz')['smiled_data']
corrected_datacube = np.load(f'{DataFolder}corrected_datacube.npz')['corrected_data']

num_cols = original_datacube.shape[2]

def plot_pixel_spectrum(row, col):
    plt.subplot(2, 1, 1)
    plt.plot(wavelength, original_datacube[:, row, col], 'r', label='Original Radiance')
    plt.plot(wavelength, smiled_datacube[:, row, col], 'b', label='Smiled Radiance')
    plt.title(f'Radiance Data (Original + Smiled) at (row, col) = {row, col}')
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Radiance')
    plt.legend(loc='upper right')   

    plt.subplot(2, 1, 2)
    plt.plot(wavelength, original_datacube[:, row, col], 'r', label='Original Radiance')
    plt.plot(wavelength, corrected_datacube[:, row, col], 'g', label='Corrected Radiance')
    plt.title(f'Radiance Data (Original + Corrected) at (row, col) = {row, col}')
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Radiance')
    plt.legend(loc='upper right')

    plt.subplots_adjust(hspace=0.75)
    plt.show()

def plot_row(row, step):
    for col in range(0, num_cols, step):
        plot_pixel_spectrum(row, col)

if __name__ == '__main__':
    row = 50
    col = 80
    step = 10

    # plot_pixel_spectrum(row, col)
    plot_row(row, step)