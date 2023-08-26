"""This is a temporary file that contains code snippets that are used for testing purposes. 

DO NOT IMPORT THIS FILE INTO ANY OTHER FILES.
"""
from matplotlib import pyplot as plt
import numpy as np

test = np.zeros((200, 145, 145))

for i in range(np.shape(test)[2]):
    data_slice = np.transpose(test[:, :, i])
    print(np.shape(data_slice))

    for j, row in enumerate(data_slice):
        print(np.shape(row))
        break
    break

