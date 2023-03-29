# Author: Shivesh Prakash
# This file extracts an array of data from MODTRAN corresponding to a specific feature

import numpy as np
import json as js


def extract_from_MODTRAN(
    json_file: str,
    dataset_no: int = 0,
    start_wavelength: int = 400,
    end_wavelength: int = 1700,
    deeper_by: int = 0,
    wider_by: int = 0,
) -> np.array:
    """This function returns an array of data points corresponding to the specific feature.

    Args:
        json_file = the location of the json file
        dataset_no = the number of the dataset in the json file
        start_wavelength = the start wavelength around which data is to be extracted
        end_wavelength = the end wavelength around which data is to be extracted
        deeper_by = the amount by which the y-axis is to be lowered in order to fit in with the satellite data
        wider_by = the amount by which the x-axis is to be widened in order to fit in with the satellite data

    Outputs:
        np.array consisting of x and y-axis data points, wavelength in nm on x-axis and transmittance on the y-axis
    """
    f = open(json_file, "r")
    data = js.load(f)
    json = data["datasets"][dataset_no]
    key = list(json.keys())[0]
    ind = 0
    for i in range(len(json[key]["roots"]["references"])):
        if "data" in json[key]["roots"]["references"][i]["attributes"]:
            ind = i

    xpoint = np.array(json[key]["roots"]["references"][ind]["attributes"]["data"]["x"])
    ypoints = np.array(json[key]["roots"]["references"][ind]["attributes"]["data"]["y"])
    xpoints = np.array([d * 1000 for d in xpoint])

    n1, n2 = start_wavelength, end_wavelength
    xp_temp = np.array([d for d in xpoints if n1 < d < n2])
    pos = [list(xpoints).index(d) for d in list(xpoints) if n1 < d < n2]
    optimised_y = []
    optimised_x = []
    for i in pos:
        if ypoints[i] < 0.91:
            optimised_y.append(ypoints[i] - deeper_by)
        else:
            optimised_y.append(ypoints[i])
    x_of_min = xp_temp[np.argmin(ypoints[pos])]
    for point in xp_temp:
        if point < x_of_min:
            optimised_x.append(point - wider_by)
        elif point > x_of_min:
            optimised_x.append(point + wider_by)
        else:
            optimised_x.append(point)
    yp = np.array(optimised_y)
    xp = np.array(optimised_x)
    datapoints = []
    for i in range(len(xp)):
        datapoints.append((xp[i], yp[i]))

    return np.array(datapoints)


extract_from_MODTRAN("..data/MODTRANdata.json")
