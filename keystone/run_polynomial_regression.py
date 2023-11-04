# Author: Shivesh Prakash
# This file returns a polynomial regression model

import numpy as np

def run_polynomial_regression(list1, list2, degree = 1):
    """
    Perform polynomial regression to predict list2 from list1.

    Args:
        list1 (list): List of values (independent variable).
        list2 (list): List of values (dependent variable).

    Returns:
        np.poly1d: A cubic polynomial regression model.
    """
    if len(list1) != len(list2):
        raise ValueError("Both lists must be of the same length.")
    
    coefficients = np.polyfit(list1, list2, degree)
    regression_model = np.poly1d(coefficients)
    
    return regression_model