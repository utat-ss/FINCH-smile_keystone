# Author: Evan Song
# This file contains a helper function to implement the coarse bisection technique

from apply_keystone import apply_keystone
from compute_sinc_fit_peak import compute_sinc_fit_peak

def coarse_bisection(guess, imagex, adjustment, imagex1):
    """
    Executes coarse bisection method to find the guess with the lowest phase shift.
    
    Args:
    - guess: Initial polynomial guess for coarse bisection
    - imagex: Image at band x
    - adjustment: User-defined amount to adjust the guess up and down to obtain two more guesses for coarse bisection
    - imagex1: Image at band x+1

    Returns:
    - Tuple: (Guess with the lowest phase shift, Associated phase shift)
    """
  
    # Adjust guesses
    up_adjusted_guess = guess + adjustment
    down_adjusted_guess = guess - adjustment
    
    # Compute phase shifts for adjusted guesses
    adjusted_1 = apply_keystone(imagex1, up_adjusted_guess)
    adjusted_2 = apply_keystone(imagex1, down_adjusted_guess)
    adjusted_3 = apply_keystone(imagex1, guess)
    
    shift_1 = compute_sinc_fit_peak(imagex, adjusted_1)
    shift_2 = compute_sinc_fit_peak(imagex, adjusted_2)
    shift_3 = compute_sinc_fit_peak(imagex, adjusted_3)
    
    # Get the guess with the lowest phase shift
    shift_list = [shift_1, shift_2, shift_3]
    minimum = min(shift_list)
    
    # Return the guess with the lowest phase shift
    if minimum == shift_1:
        return (up_adjusted_guess, shift_1)
    elif minimum == shift_2:
        return (down_adjusted_guess, shift_2)
    else:
        return (guess, shift_3)
