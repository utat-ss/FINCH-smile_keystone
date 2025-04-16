import numpy as np
import matplotlib.pyplot as plt 
import plotly.express as px
import plotly.graph_objects
import pandas as pd
import os

# a seperate file named data_from which stores location of whatever dataset you want analyzed
from data_from import indian_pine, indian_pine_wavelength_separation, spectral_monuments_min, spectral_monuments_max

metadata = indian_pine_wavelength_separation

# The data that is being analyzed (.npy file), and the file that the data will be written into 
data = np.load(indian_pine)

# ---------- Aux Functions (These functions probs won't be called, but they're used for the other functions) ---------- #
file_path = indian_pine

def normalize_pixels(array_to_normalize: np.array, maximum_pixel: int, new_max:int):
    '''
    I have to normalize the pixels because its in watts per square meter per steradian (whut) 
    
    Function takes two parameters: 
    array_to_normalize is a numpy array. It is the array with values in watts per square meter per steradian (this function is more or less attatched to the create_band_sheet function)
    maximum_pixel is an integer that is the maximum value in watts per square meter per steradian 
    '''
    
    # Gets the size of the data that will be created so we can iterate all the data
    data_shape = data.shape

    # Creates the array that we will be putting data into once the data is normalized. Creates an array that is the same size as the data but with all 0s
    normalized_array = np.zeros((data.shape[0], data.shape[1]))

    # normalizes the data 
    for number1 in range(0, data_shape[0]):
        for number2 in range(0, data_shape[1]):
            to_put = (array_to_normalize[number1][number2] / maximum_pixel) * new_max

            normalized_array[number1][number2] = to_put

            print(normalized_array[number1][number2])

    return normalized_array

# Specific band for a certain pixel
def pixel_wavelength_information(result=[int, int, int]):
    '''
    Takes a pixel and prints out a particular wavelength for it in the specified band. Takes 3 integers in a list, first two inputs (result[0] and result[1]) 
    correspond to the spot they take up on the image, and the third input (result[2]) corresponds to the wavelength band that you want to use
    '''
    
    # result[0] is the x component and result[1] is the y component and result [2] is the specific wavelength you want to see
    return (data[result[0]][result[1]][result[2]])

# Gets all the wavelength data for a particular pixel
def pixel_information(result=[int, int]):
    '''
    Takes a pixel and prints out all the wavelengths for that particular pixel. result[0] corresponds to x location and result[1] corresponds to y location on the image
    '''
    
    # result[0] is the x component and result[1] is the y component
    return (data[result[0]][result[1]])

def create_a_panda(wavelength: int): 
    '''
    Creates up a panda
    '''

    # Sets the thing to infinity so none of that "..." truncating funny business happens
    np.set_printoptions(threshold=np.inf)
    
    # Gets the shape of the data, only ever tried this on the indian pines dataset but hopefully it'll work for anything 
    data_shape = data.shape

    # Makes an array with data_shape[0] * data_shape[1] pixels to put numbers into 
    my_array = np.zeros((data_shape[0], data_shape[1]))
    
    # This gets the highest watt per square meter per steradian for the normalization function
    maximum_pixel_brilliance = 0
    
    # iterates through the elements in the array and puts it into my_array, units are in watts per square unit per steradian
    for number1 in range(0, data_shape[0]):
        for number2 in range(0, data_shape[1]):
                      
            # Specifies the pixel you want the information for, and the band you want it from
            to_get = [number1, number2, wavelength]
            
            # Passes above information into the pixel_wavelength_information, and it gets the (brilliance?) of a pixel
            to_put = pixel_wavelength_information(to_get)
            
            # To find the maximum brightness pixel, the max is compared
            maximum_pixel_brilliance = max(maximum_pixel_brilliance, to_put)
            
            # The numpy array at [number1][number2] becomes the value that was found in the variable to_put
            my_array[number1][number2] = to_put

    df = pd.DataFrame(my_array)

    print(df)
    
    df.to_csv(f"{wavelength}_panda.csv", index=False)

# ----- Real Functions (this is the stuff that'll get used I think) ----- #

def cook_a_line(wavelength:int, height: int):
    '''
    Gets you a graph for a line at a certain height at a certain wavelength. 
    
    wavelength: Represents the wavelength band you are selecting
    height: I think its from the top down???? so maybe it should be depth??? idk ill figure it out in the future
    '''
    
    # Gets the size of the data that will be created so we can iterate all the data
    data_shape = data.shape
    
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
    plt.plot(numbers, to_graph_array)
    plt.show()

def create_band_sheet(wavelength: int, selection: int, selection_number:int = 255):
    # WORKS PROPERLY!!
    '''
    Creates a big fat 145 * 145 numpy array of all the pixels at a certain band/wavelength, then creates an image using Pillow 
    
    Parameter 1: wavelength integer is the index of the wavelength you want to call (might change this over so you can put in the actual 
    wavelength, but that is a something for another time :D)
    
    Parameter 2: 1 indicates it is normalized to 255, 2 indicates use brilliance values
    '''
    
    # Sets the thing to infinity so none of that "..." truncating funny business happens
    np.set_printoptions(threshold=np.inf)
    
    # Gets the shape of the data, only ever tried this on the indian pines dataset but hopefully it'll work for anything 
    data_shape = data.shape
    
    # this and the "thing += 1" in the nested for loop below counts the number of pixels, should match up to data_shape[0] * data_shape[1], if it doesn't ya done goofed up 
    #thing = 0
    
    # Makes an array with data_shape[0] * data_shape[1] pixels to put numbers into 
    my_array = np.zeros((data_shape[0], data_shape[1]))
    
    # This gets the highest watt per square meter per steradian for the normalization function
    maximum_pixel_brilliance = 0
    
    # iterates through the elements in the array and puts it into my_array, units are in watts per square unit per steradian
    for number1 in range(0, data_shape[0]):
        for number2 in range(0, data_shape[1]):
                      
            # Specifies the pixel you want the information for, and the band you want it from
            to_get = [number1, number2, wavelength]
            
            # Passes above information into the pixel_wavelength_information, and it gets the (brilliance?) of a pixel
            to_put = pixel_wavelength_information(to_get)
            
            # To find the maximum brightness pixel, the max is compared
            maximum_pixel_brilliance = max(maximum_pixel_brilliance, to_put)
            
            # The numpy array at [number1][number2] becomes the value that was found in the variable to_put
            my_array[number1][number2] = to_put

    # Normalizes the array using the normalize_pixels function
    if selection == 1:
        normalized_array = normalize_pixels(my_array, maximum_pixel_brilliance, selection_number)
 
        plt.imshow(normalized_array)

        # img.show()
        plt.show()
    else: 
        plt.imshow(my_array)
        plt.show()

def interactive_3d_pixel_line_display(pixel_x:int, pixel_y:int):
    '''
    
    Gets all the wavelengths for a pixel and graphs it out
    
    pixel_x = x location of pixel 
    pixel_y = y location of pixel 
    
    Same as the pixel_graph function, except made in plotly, which gives more interactivitiy abilities and I can also graph the vertical lines
    
    Uses something called plotly to run, must be run in a interactive window (or jupyter)
    
    Import ipykernel and pip install --upgrade nbformat and hopefully it works
    
    
    '''

    # Gets the brilliance values for all the wavelengths of the pixel location provided
    thing = pixel_information([pixel_x, pixel_y])

    # Shape of the data 
    data_shape = thing.shape

    # Two lists that are appended to, creating the stuff that will be used to graph later
    thing_1 = []
    thing_2 = []

    # Goes over every wavelength and adds its brillaince value to thing_2, the wavelength is added to thing_1
    for thinint in range(data_shape[0]):
        thing_1.append(thinint * metadata[1] + metadata[0])
        thing_2.append(thing[thinint])

    # Makes a panda table with the lists created above to be made into a list later 
    df = pd.DataFrame({"Wavelength": thing_1, "Brilliance": thing_2})

    # Graph created 
    fig = px.line(df, x="Wavelength", y="Brilliance")

    # Iterates over all the minimum wavelengths and puts them onto the graph 
    for minimum_wavelength in spectral_monuments_min:

        # Used to make sure the list is longer than 1, it is is, a rectangle is added
        try:

            # DO NOT DELETE THIS, this triggers an error if its only a single element (variable isn't used that butat's how its supposed to be) 
            items = len(minimum_wavelength)
       
            # Rectangle added to the figure
            fig.add_vrect(x0=minimum_wavelength[0], x1=minimum_wavelength[1], opacity=0.1, line_width=0, fillcolor="red")
            
        # If the length of minimum wavelength gives an error, it is not a list, so it only makes a line 
        except TypeError:
            
            # Line added to the figure 
            fig.add_vline(x=minimum_wavelength, line_width=0.5, line_dash="dash", line_color="red", opacity=0.75)
            
    for maximum_wavelength in spectral_monuments_max:

        # Used to make sure the list is longer than 1, it is is, a rectangle is added
        try:
            
            # DO NOT DELETE THIS, this triggers an error if its only a single element (variable isn't used that butat's how its supposed to be) 
            items = len(maximum_wavelength)
            
            # Rectangle added to the figure
            fig.add_vrect(x0=maximum_wavelength[0], x1=maximum_wavelength[1], opacity=0.1, line_width=0, fillcolor="blue")
            
        # If the length of minimum wavelength gives an error, it is not a list, so it only makes a line 
        except TypeError:
            
            # Line added to the figure 
            fig.add_vline(x=maximum_wavelength, line_width=0.75, line_dash="dash", line_color="blue", opacity=0.75)

    # Figure is displayed 
    fig.show()

def interactive_3d_graph_display(wavelength: int):
    '''
    What does this function do o-O 

    Hmm, I wonder why it only takes in wavelength and nothing else, what magic does it perform :O 
    '''

    supposed_path = f"{wavelength}_panda.csv"

    if os.path.exists(supposed_path):
        pass
    else: 
        create_a_panda(wavelength)

    z_data = pd.read_csv(supposed_path)
    z = z_data.values 

    wavelength = 10 * wavelength + 400

    sh_0, sh_1 = z.shape
    x, y = np.linspace(0, 1, sh_0), np.linspace(0, 1, sh_1)
    fig = plotly.graph_objects.Figure(data=[plotly.graph_objects.Surface(z=z, x=x, y=y)])

    fig.update_traces(contours_z=dict(show=True, usecolormap=True,
                                  highlightcolor="limegreen", project_z=True))

    fig.update_layout(title=dict(text=f'Wavelength {wavelength}'), autosize=False,
                    width=800, height=800,
                    margin=dict(l=65, r=50, b=65, t=90))
    fig.show()

def pixel_graph(pixel_x:int, pixel_y:int): 
    # WORKS PROPERLY!!
    '''
    Gets all the wavelengths for a pixel and graphs it out
    
    pixel_x = x location of pixel 
    pixel_y = y location of pixel 
    '''

    # Gets the brilliance values for all the wavelengths of the pixel location provided
    thing = pixel_information([pixel_x, pixel_y])

    # Gets the shape of the thing
    data_shape = thing.shape
    
    # Lists to be graphed 
    thing_1 = []
    thing_2 = []

    # Graphs the wavelength in thing_1 and brilliance in thing_2
    for thinint in range(data_shape[0]):
        thing_1.append(thinint * metadata[1] + metadata[0])
        thing_2.append(thing[thinint])

    # Makes the stuff into a numpy to be graphed
    thing_1 = np.array(thing_1)
    thing_2 = np.array(thing_2)
    
    # Graphs the stuff 
    plt.title("Pixel Graph")
    plt.xlabel("Wavelength Band")
    plt.ylabel("Brilliance")
    plt.plot(thing_1, thing_2)
    plt.show()   

def cook_a_line(wavelength:int, height: int):
    '''
    Gets you a graph for a line at a certain height at a certain wavelength. 
    
    wavelength: Represents the wavelength band you are selecting
    height: I think its from the top down???? so maybe it should be depth??? idk ill figure it out in the future
    '''
    
    # Gets the size of the data that will be created so we can iterate all the data
    data_shape = data.shape
    
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
    plt.plot(numbers, to_graph_array)
    plt.show()

# ---------- START OF CHEESE FUNCTIONS ---------- #

def cheese1(maxed: int, wavelength: int, height: int):
    '''
    Parabolic Shift. Applies a parabolic smile shift to spectral data.
    
    Args:
        maxed (int): Maximum shift amount at the edges
        wavelength (int): The wavelength band to visualize
        height (int): The row position in the image to analyze
    '''
    data_shape = data.shape
    
    # Extract the original data
    my_data = []
    numbers = []
    
    for pixel_number in range(data_shape[1]):
        pixel_informations = pixel_wavelength_information([height, pixel_number, wavelength])
        my_data.append(pixel_informations)
        numbers.append(pixel_number)
    
    # Calculate the shift amount for each position
    shift_data = []
    for number in numbers:
        # Parabolic function that creates a smile/frown effect
        shift = (-maxed)/((data_shape[1]/2) ** 2) * number * (number - data_shape[1])
        shift_data.append(shift)
    
    # Apply the shift to create modified positions
    modified_numbers = []
    for i, number in enumerate(numbers):
        modified_numbers.append(number + shift_data[i])
    
    # Graph the results
    grapher(numbers, my_data, modified_numbers, shift_data, 1)


def cheese2(amplitude: int, wavelength: int, height: int, mean: float = None, std_dev: float = None, shift_direction: int = 1):
    '''
    An improved version of the cheese function that uses a normal distribution.
    
    Parameters:
    -----------
    amplitude: int
        The maximum height of the normal distribution curve. This part is the maximum shift, and it appears at the mean (This is the height of the tallest part of the normal distibution (the part that is shifted the most is shifted by this amount))
    wavelength: int
        The wavelength band you are selecting. It is based on this, you get the row of wavelengths, that is then shifted over 
    height: int
        The height position in the image. After you have selected the wavelegnth of the image, you also get to choose how high up on the image you want to gerneate smile shift for. (70 would be the 70th pixel from the bottom, either that, or the70th pixel from the top)
    mean: float, optional
        The mean of the normal distribution. If None, defaults to the center of the data. Where is the largest distibution located?
    std_dev: float, optional
        The standard deviation of the normal distribution. How quickly does the normol distibution drop off
    shift_direction: int, optional 
        Defines whether to implement smile shift forward or backwards, 1 implements it forward, -1 implements it backwards, I'll come up with a purpose for 0 one day
    '''
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
    
    # Calculate normal distribution values
    normal_dist = []
    for x in numbers:
        # Normal distribution formula: f(x) = amplitude * exp(-(x-mean)²/(2*std_dev²))
        normal_value = amplitude * np.exp(-((x - mean) ** 2) / (2 * std_dev ** 2))
        normal_dist.append(normal_value)
    
    # Create modified data with normal distribution
    modified_numbers = []
    for i, x in enumerate(numbers):
        if shift_direction == 1:
            modified_numbers.append(x - normal_dist[i])
        elif shift_direction == -1:
            modified_numbers.append(x + normal_dist[i])
    
    grapher(numbers, to_graph_array, modified_numbers, normal_dist, 2)

def cheese3(amplitude: int, wavelength: int, height: int, mean: float = None, std_dev_left: float = None, std_dev_right: float = None, shift_direction: int = 1):
    '''
    An improved version of the cheese function that uses a normal distribution.
    
    Parameters:
    -----------
    data: np.ndarray
        The input image or data array from which we extract pixel information.
    amplitude: int
        The maximum height of the normal distribution curve. Determines the maximum shift.
    wavelength: int
        The wavelength band to be processed.
    height: int
        The vertical position in the image where the shift is applied.
    mean: float, optional
        The mean of the normal distribution. Defaults to the center of the image width.
    std_dev_left: float, optional
        Standard deviation for the left side of the distribution.
    std_dev_right: float, optional
        Standard deviation for the right side of the distribution.
    shift_direction: int, optional 
        Direction of the shift: 1 for forward, -1 for backward.
    '''

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
    
    if std_dev_right is None:
        std_dev_right = data_shape[1] / 6  # Default: 1/6 of the data width

    if std_dev_left is None:
        std_dev_left = data_shape[1] / 6  # Default: 1/6 of the data width

    # Calculate normal distribution values
    normal_dist = []
    for x in numbers:
        if x < mean:
            normal_value = amplitude * np.exp(-((x - mean) ** 2) / (2 * std_dev_left ** 2))
        else: 
            normal_value = amplitude * np.exp(-((x - mean) ** 2) / (2 * std_dev_right ** 2))
        normal_dist.append(normal_value)

    # Create modified data with normal distribution
    modified_numbers = []
    for i, x in enumerate(numbers):
        if shift_direction == 1:
            modified_numbers.append(x - normal_dist[i])
        elif shift_direction == -1:
            modified_numbers.append(x + normal_dist[i])
        else:
            modified_numbers.append(x) 

    grapher(numbers, to_graph_array, modified_numbers, normal_dist, 3)


def grapher(numbers: list, my_data, modified_numbers, shift_data, version):
    """
    This function was created to reduce the amount of code that was used in graphing the cheese functions. By putting it into one function, I won't have to figure out how to graph things again for every version of smile shift generator I create. 
    
    This code only kind of works for SMILE Version 1, and I don't think I'll have to use that one, so I'm not gonna spend a bunch of time trying to figure that part out!
    Args:
        numbers (list): This is the original numbers that are passed into the function. By default, this is just counting up by 1
        my_data (list): Takes in a list of the spectral data, this is the same for both the modified and unmodified data, as smile shift does not occur in in the y axis (no reduction or increase in brilliance)
        modified_numbers (list): Same as numbers, but after smile shift has been applied
        shift_data (list): This is the graph of the amount of shift applied
        version (int): The version of smile shift that was used, allows you to pinpoint the function that generated the shift
    """
    '''
    There are way too many graphing functions for the cheese stuff, and I have to figure out which function it is, which is quite a pain. This function will be used as the one that graphs all the cheese information. 

    Takes in two lists, the original data, and the modified data, and graphs them. This makes it so I don't have to rewrite the plotting function every time (wastes a lot of code)
    
    The curve that is used to modify the data will be plotted under everything else. (might make a separate function for this, for only the affected, but we'll figure that out in the future)
    '''
    # Plot just the normal distribution curve
    plt.figure(figsize=(10, 6))
    plt.title(f"Spectra Graph: SMILE Version {version}")
    plt.xlabel("Location")
    plt.ylabel("Brilliance")
    plt.plot(modified_numbers, my_data, zorder=3, color="red", label="Modified Data")
    plt.plot(numbers, my_data, zorder=2, color="blue", label="Original Data")
    plt.legend(loc="upper left")
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.title(f"Change Curve: SMILE Version {version}")
    plt.xlabel("Location")
    plt.ylabel("Location Shift")
    plt.plot(numbers, shift_data)
    plt.grid(True)
    plt.show()

# ---------- END OF CHEESE FUNCTIONS ---------- #

create_band_sheet(3, 2)
cook_a_line(5, 70)
pixel_graph(30, 30)
interactive_3d_pixel_line_display(100, 100)

interactive_3d_graph_display(3)

# Figure out if the shifted amoutn is always forward or always backwards. Nobody has an idea! We ball!
cheese2(amplitude=2, wavelength=100, height=70, mean=100, std_dev=60, shift_direction=1)
cheese3(amplitude=3, wavelength=100, height=70, mean=100, std_dev_left=60, std_dev_right=30, shift_direction=1)