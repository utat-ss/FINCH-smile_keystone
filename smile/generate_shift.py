# Author: Shivesh Prakash
# This file contains a function to induce artificial SMILE distortions
# Work in Progress, do not run yet

import numpy as np
import matplotlib.pyplot as plt

def gen_shift(x, y, lam, a, b, c):
    subpixel_length = (x[-1] - x[0]) / (len(x) - 1)
    for i in range(len(x)):
        X = x[i]
        shift = int((a * X**2 + b * X + c) / subpixel_length)
        for j in range(len(y)):
            if j + shift in y:
                lam[i][j] = lam[i][j + shift]
                lam_shift[i][j] = shift
    shifts = [lam_shift[i][0] for i in range(len(lam_shift))]
    plt.plot(x, shifts)
    plt.show()


gen_shift(x, y, lam, -1 / 14, 19 / 14, 0)
