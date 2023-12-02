
import numpy as np
import matplotlib.pyplot as plt
import math

from scipy import stats
from IPython.display import clear_output
def coarse_bisection(guess, imagex, adjustment, imagex1):
  """
  Args:
  guess: our first initial polynomial guess for coarse bisection
  imagex: image at band x
  adjustment: amount we are adjusting the guess up and down by, to get 2 more guesses for coarse bisection (user defined)
  imagex1: image at band x+1

  Returns:
  the guess with the lowest phase shift, in the form: tuple(guess, associated phase shift)
  """
  
  
  up_adjusted_guess = guess + adjustment
  down_adjusted_guess = guess - adjustment
  adjusted_1 = apply_keystone(imagex1, up_adjusted_guess)
  adjusted_2 = apply_keystone(imagex1, down_adjusted_guess)
  adjusted_3 = apply_keystone(imagex1, guess)
  shift_1 = compute_sinc_fit_peak(imagex, adjusted_1)
  shift_2 = compute_sinc_fit_peak(imagex, adjusted_2)
  shift_3 = compute_sinc_fit_peak(imagex, adjusted_3)
  shift_list = [shift_1, shift_2, shift_3]
  minimum = min(shift_list)
  if minimum == shift_1:
    return (up_adjusted_guess, shift_1)
  elif minimum == shift_2:
    return (down_adjusted_guess, shift_2)
  else:
    return (guess, shift_3)
