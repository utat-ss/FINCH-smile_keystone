# Author: Shuhan
# Step 3 of Smile. This file creates reference and test spectra from the provided data cube

def create_ref_and_test_spectra(crop_range:tuple, show_plots = False):
    """
    Returns reference and test spectra. 

    Args:
        crop_range: just enter the band numbers as (int, int). For example, enter crop_range = (0, 5) for showing bands 0~5 
        show_plot: boolean, set it to True if you need to show plots

    Outputs:
        cropped_sampled_reference: cropped resampled reference spectra, which was originally generated using the first colum. 
            Its shape is (# of shifts, # of bands)
        cropped_sampled_test: cropped resampled test spectra. Its shape is (# of columns, # of shifts, # of bands)
    """
    band_length = int((max(wavelength) - min(wavelength)) / g_num_of_bands)
    crop_band_start = crop_range[0] * band_length
    crop_band_end = (crop_range[1] + 1) * band_length
    shift_bound = g_num_shifts_1D * g_shift_increment 
    shift_range = (-shift_bound, shift_bound)

    sampled_reference_spectra, sensors_position_reference, srf_columns_reference = run_resampling_spectra(data_for_resampling[0], test_spectral_response, shift_range)

    sampled_test_spectra, sensors_position_test, srf_columns_test = run_resampling_spectra(data_for_resampling, test_spectral_response, 0) 
    
    cropped_sampled_reference = [i[crop_range[0]:crop_range[1]] for i in  sampled_reference_spectra]
    cropped_sampled_test = [[j[crop_range[0]:crop_range[1]] for j in i] for i in sampled_test_spectra]
            

    if show_plots:
        fig, plot_for_show = plt.subplots(1, 1, figsize = (21, 7))
        plot_for_show.plot(wavelength, data_for_resampling[0])
        plot_for_show.set_xlabel('Wavelength [nm]')
        plot_for_show.set_xlim(crop_band_start + min(wavelength), crop_band_end + min(wavelength))

        i = 0
        for shift in range(len(sampled_reference_spectra)):
            plot_for_show.scatter(stretch_horizontal(sensors_position_reference[i], wavelength), sampled_reference_spectra[i], label = f"shift = {shift}")
            plot_for_show.plot(np.linspace(min(wavelength), max(wavelength), len(srf_columns_reference[i])), srf_columns_reference[i], label = f"shift = {shift}")
            i += 1
        
        plot_for_show.legend(loc = 'center left', ncol = 1, bbox_to_anchor = (1, 0.5))
        fig.savefig(f'Demo')

    return cropped_sampled_reference, cropped_sampled_test

