# Author: Shuhan
# This file contains mathematical functions that are currently placeholders for the sensor's spectral response functions. 

def test_spectral_response(input):
    """
    A test spectral response function. Basically a Gaussian function centered 
        at mu * (max(input) - min(input)) with std = sigma

    Args: 
        input: a 1D array

    Output: 
        A Gaussian function corresponding to the input. A 1D array
    """
    sigma = 0.25 * len(input)
    mu = 0.5 * len(input)
    gaussian = stats.norm.pdf(input, mu, sigma)
    normed_gaussian = gaussian / max(gaussian)
    
    return normed_gaussian

def make_random_SRFs(input):
    """
    Generates a list of unique and random Gaussian SRFs for each pixel. The parameters of each
        generated Gaussian function is a random number between 0 and 1. These random numbers 
        follow normal(Gaussian) distribution
    Args:
        input: a 1D array

    Output: 
        A list of Gaussian functions. Interact with them like this: 
            func_list = make_random_SRFs
            
            test_SRF_1 = func_list[1](inputs)
            test_SRF_2 = func_list[2](inputs)
            ...
    """
    output = []
    for i in range(g_num_of_bands):
        sigma_temp = np.random.normal(loc = 0.5, scale = 0.1) * len(input)
        mu_temp = np.random.normal(loc = 0.5, scale = 0.1) * len(input)

        def temp (input):

            gaussian = stats.norm.pdf(input, mu_temp, sigma_temp)
            normed_gaussian = gaussian / max(gaussian)

            return normed_gaussian
          
        globals()[f'random_SRF_{i}'] = temp
    
        output.append(globals()[f'random_SRF_{i}'])

    return output

# plt.plot(test_spectral_response(np.arange(0, 100)))
# plt.plot(test_spectral_response(np.arange(0, 100) + 10))
