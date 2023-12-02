# Author: Shivesh Prakash
# This file inverts a polynomial regression model

import numpy as np

def invert_polynomial_regression(regression_model, value):
    """
    Use the polynomial regression model to predict list1 value from list2.
    Usage: invert_polynomial_regression(model, value)

    Args:
        regression_model (np.poly1d): Cubic polynomial regression model obtained from `run_cubic_polynomial_regression`.
        value: A value from list2 for which you want to find the corresponding list1 value.

    Returns:
        float: The predicted value.
    """
    if not isinstance(regression_model, np.poly1d):
        raise ValueError("Invalid regression model. Please provide a valid cubic polynomial regression model.")
    
    # Find the roots (solutions) of the polynomial to invert the regression
    roots = np.roots(regression_model - value)
    
    # Select the real root (ignoring complex roots)
    real_root = np.real(roots[np.isreal(roots)])[0]
    
    return real_root
