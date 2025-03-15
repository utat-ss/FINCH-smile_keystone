import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from torch.utils.data import Dataset, DataLoader
import time
from datetime import datetime

# ----------------------- DEFINING PARAMETERS ----------------------- #

# Load the data
data = pd.read_csv('simpler_data.csv')

# Defining the device to be used, this will use the GPU if available
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# print(f"Using device: {device}")

# ----------------------- DATASET INIITIATOR ----------------------- #

# Initiating it with the Dataset class allows for seperation between normal models and dataset, makes it easier to do stuff 
class Spectral_Dataset(Dataset):
    def __init__(self, csv_data_file):
        """
        Initialize the dataset from a CSV file containing spectral data.
        
        Args:
            csv_file (str): Path to the CSV file containing spectral data
        """

        self.data = pd.read_csv(csv_data_file) 
        
        # Extract spectral bands (columns starting from the 6th column, after use)
        spectral_columns = self.data.columns[5:]
        self.spectra = self.data[spectral_columns].values.astype(np.float32)
        



# ----------------------- DEFINING GENERATOR ----------------------- #
class Generator(nn.Module):
    def __init__(self, input_dim, output_dim):
        super(Generator, self).__init__()
        self.main = nn.Sequential(
            nn.Linear(input_dim, 128),

            # the inplace=True reduces the amount of memory used (very slightly) by destroying the original input used. Might speed it up a tiny bit. It should be applicable in this use case because later layers never refer back to a previous layer
            # No idea which activation function to use, so LeakyReLu to test it out 
            nn.LeakyReLU(0.2, inplace=True),
            nn.Linear(128, 256),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Linear(256, 512),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Linear(512, 1024),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Linear(1024, output_dim),

            # Output layer 
            nn.Tanh()
        )

    def forward(self, x):
        return self.main(x)