# Author: Shivesh Prakash
# This file contains a recursive function to implement the coarse bisection technique

from coarse_bisection import coarse_bisection

def coarse_bisection_recursive(imgx, imgx1, guess, shift, threshold, iteration, max_iteration):
    """
    Applies a recursive coarse bisection method to refine a guess until a phase shift threshold is reached 
    or maximum iterations are reached.

    Args:
    - imgx: Image at band x
    - imgx1: Image at band x+1
    - guess: Initial guess for coarse bisection
    - shift: Current phase shift value associated with the guess
    - threshold: Threshold value for the phase shift, determines termination condition
    - iteration: Current iteration count
    - max_iteration: Maximum number of iterations allowed

    Returns:
    - Updated guess or the last calculated guess if the threshold or maximum iterations are reached
    """

    # Termination condition: phase shift below threshold or reached maximum iterations
    if shift <= threshold or iteration >= max_iteration:
        return guess
    else:
        # Refine the guess using coarse bisection method
        guess, shift = coarse_bisection(guess, imgx, adjustment, imgx1)
        
        # Recursive call with updated guess and shift values
        return coarse_bisection_recursive(imgx, imgx1, guess, shift, threshold, iteration + 1, max_iteration)
