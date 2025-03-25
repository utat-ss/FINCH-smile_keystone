import numpy as np
import matplotlib.pyplot as plt

# Import the functions here (Look for these functions in the smile folder)
from GenerateSmileFns import *
from Quantificationfns import *
from Correctionfns import *


def cheese(maxed:int, wavelength: int, height: int):
    '''
    Crazy unoptimized. Look at cook_a_line and implement it one day so I don't copy everything over
    '''
    data_shape = data.shape

    print(data_shape)

    to_graph_array = []
    numbers = []
    
    for pixel_number in range(data_shape[1]):
        first = pixel_information([height, pixel_number])
        pixel_informations = pixel_wavelength_information([height, pixel_number, wavelength])
        
        to_graph_array.append(pixel_informations)
        numbers.append(pixel_number)
        
    plt.title("Band Graph")
    plt.xlabel("Location")
    plt.ylabel("Brilliance")
    plt.plot(numbers, to_graph_array, zorder=1)

    print(numbers)
    for number in range(len(numbers)): 
        numbers[number] = numbers[number]+(-maxed)/((data_shape[1]/2) ** 2) * number * (number - data_shape[1])
    print(numbers)
    plt.plot(numbers, to_graph_array, zorder=2)
    plt.show()

    g_1 = []
    for _ in range(145):
        g_1.append(_)

    g_2 = numbers

    print(g_1)
    for number in range(len(g_2)): 
        g_2[number] = (-maxed)/((data_shape[1]/2) ** 2) * number * (number - data_shape[1])
    print(g_2)

    plt.plot(g_1, g_2)
    plt.show()

def cheese2(amplitude: int, wavelength: int, height: int, mean: float = None, std_dev: float = None):
    '''
    An improved version of the cheese function that uses a normal distribution.
    
    Parameters:
    -----------
    amplitude: int
        The maximum amplitude of the normal distribution curve. This part is the maximum shift, and it appears at the mean 
    wavelength: int
        The wavelength band you are selecting
    height: int
        The height position in the image
    mean: float, optional
        The mean of the normal distribution. If None, defaults to the center of the data
    std_dev: float, optional
        The standard deviation of the normal distribution
    '''
    import numpy as np
    
    data_shape = data.shape
    
    # Extract the original data
    to_graph_array = []
    numbers = []
    
    for pixel_number in range(data_shape[1]):
        pixel_informations = pixel_wavelength_information([height, pixel_number, wavelength])
        to_graph_array.append(pixel_informations)
        numbers.append(pixel_number)
    
    # Set default values for mean and std_dev if not provided
    if mean is None:
        mean = data_shape[1] / 2  # Center of the data
    
    if std_dev is None:
        std_dev = data_shape[1] / 6  # 1/6 of the data width
    
    # Plot the original data
    plt.figure(figsize=(10, 6))
    plt.title("Band Graph with Normal Distribution")
    plt.xlabel("Location")
    plt.ylabel("Brilliance")
    plt.plot(numbers, to_graph_array, label="Original Data", zorder=1)
    
    # Calculate normal distribution values
    normal_dist = []
    for x in numbers:
        # Normal distribution formula: f(x) = amplitude * exp(-(x-mean)²/(2*std_dev²))
        normal_value = amplitude * np.exp(-((x - mean) ** 2) / (2 * std_dev ** 2))
        normal_dist.append(normal_value)
    
    # Create modified data with normal distribution
    modified_numbers = []
    for i, x in enumerate(numbers):
        modified_numbers.append(x + normal_dist[i])
    
    # Plot the modified data
    plt.plot(modified_numbers, to_graph_array, label="Modified Data", zorder=2)
    plt.legend()
    plt.show()
    
    # Plot just the normal distribution curve
    plt.figure(figsize=(10, 6))
    plt.title("Normal Distribution Curve")
    plt.xlabel("Location")
    plt.ylabel("Value")
    plt.plot(numbers, normal_dist)
    plt.grid(True)
    plt.show()
    
    return modified_numbers, normal_dist