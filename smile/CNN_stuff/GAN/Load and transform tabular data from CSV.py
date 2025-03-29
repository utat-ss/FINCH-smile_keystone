import os 
import pandas as pd
import torch
from torch.utils.data import TensorDataset, DataLoader
import torchvision.transforms as transforms 
from torchvision.utils import make_grid
from torchvision import datasets
import torch.nn as nn 
import torch.nn.functional as F
import torch 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.pyplot import figure 
from tqdm import tqdm
plt.ion() 
from IPython.display import clear_output


# Hyperparameters
batch_size = 64
lr = 0.00005
n_epochs = 50
clip_value = 1
n_critic = 5
z_dim = 256
input_size = 3  # gv, npv and soil

# Device configuration
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Data path
simpler_data_location = r"C:\Users\liche\OneDrive\Desktop\PycharmProjects\Django Todo\FINCH-smile_keystone\smile\CNN_stuff\GAN\simpler_data.csv"

# Load the data from CSV
def load_data():
    # Read the CSV file
    print(f"Loading data from {simpler_data_location}")
    df = pd.read_csv(simpler_data_location)

    df = df.drop(columns=['Spectra'])
    
    # Print some information about the data
    print(f"Data shape: {df.shape}")
    print(f"Column names: {df.columns.tolist()}")
    
    # Convert DataFrame to numpy array
    data_array = df.values.astype(np.float32)
    
    # Normalize the data to [-1, 1] range
    # This is equivalent to what transforms.Normalize([0.5], [0.5]) would do
    data_normalized = (data_array - 0.5) / 0.5
    
    # Convert numpy array to PyTorch tensor
    data_tensor = torch.tensor(data_normalized)
    
    # Create a TensorDataset
    dataset = TensorDataset(data_tensor)
    
    # Create a DataLoader
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=True,
        drop_last=True  # Drop the last batch if it's smaller than batch_size
    )
    
    print(f"Created dataloader with {len(dataloader)} batches")
    return dataloader, df.shape[1]  # Return dataloader and number of features

# Load data and get the number of features
dataloader, num_features = load_data()
print(f"Number of features in the dataset: {num_features}")


class Generator(nn.Module):
    def __init__(self):
        super(Generator, self).__init__()

        self.model = nn.Sequential(          
            # Input: random noise vector of size z_dim
            nn.Linear(in_features=z_dim, out_features=200), 
            nn.ReLU(inplace=True),
            nn.Dropout(0.05),
            
            nn.Linear(in_features=200, out_features=400),
            nn.ReLU(inplace=True),
            nn.Dropout(0.05),

            nn.Linear(in_features=400, out_features=400),  # Fixed input size
            nn.ReLU(inplace=True),
            nn.Dropout(0.05),

            nn.Linear(in_features=400, out_features=200),  # Fixed input size
            nn.ReLU(inplace=True),
            nn.Dropout(0.05),
            
            # Output layer: should match the number of features in our dataset
            nn.Linear(in_features=200, out_features=num_features),
            # Tanh to output values in the range [-1, 1], matching our normalized data
            nn.Tanh()
        )

    def forward(self, z):
        return self.model(z)
    
class Critic(nn.Module):
    def __init__(self):
        super(Critic, self).__init__()

        self.model = nn.Sequential(
            # Input: data with num_features dimensions
            nn.Linear(in_features=num_features, out_features=200),
            nn.ReLU(inplace=True),
            nn.Dropout(0.05),

            nn.Linear(in_features=200, out_features=400),
            nn.ReLU(inplace=True),
            nn.Dropout(0.05),

            nn.Linear(in_features=400, out_features=200),
            nn.ReLU(inplace=True),
            nn.Dropout(0.05),

            nn.Linear(in_features=200, out_features=100),
            nn.ReLU(inplace=True),
            nn.Dropout(0.05),

            nn.Linear(in_features=100, out_features=1),
            # No activation at the end for WGAN
        )

    def forward(self, x):
        return self.model(x)

# Initialize models
generator = Generator().to(device)
critic = Critic().to(device)

# Optimizers
optimizer_G = torch.optim.RMSprop(generator.parameters(), lr=lr)
optimizer_C = torch.optim.RMSprop(critic.parameters(), lr=lr)

def train_wgan(): 
    # Store generated samples
    gen_samples = []

    # Losses for the critic and the generator
    critic_losses = []
    generator_losses = []

    for epoch in range(n_epochs):
        # Track losses for this epoch
        epoch_critic_loss = 0
        epoch_gen_loss = 0
        
        # Progress bar for batches
        progress_bar = tqdm(enumerate(dataloader), total=len(dataloader), desc=f"Epoch {epoch+1}/{n_epochs}")
        
        for i, (real_data,) in progress_bar:
            # Move data to device
            real_data = real_data.to(device)
            batch_size = real_data.size(0)
            
            # ---------------------
            # Train Critic
            # ---------------------
            for _ in range(n_critic):
                optimizer_C.zero_grad()
                
                # Generate random noise
                z = torch.randn(batch_size, z_dim, device=device)
                
                # Generate fake data
                fake_data = generator(z)
                
                # Compute critic loss
                real_validity = critic(real_data)
                fake_validity = critic(fake_data.detach())
                
                # Wasserstein loss
                critic_loss = -torch.mean(real_validity) + torch.mean(fake_validity)
                
                # Backward pass
                critic_loss.backward()
                
                # Clip critic weights
                for p in critic.parameters():
                    p.data.clamp_(-clip_value, clip_value)
                
                # Update critic
                optimizer_C.step()
                
                epoch_critic_loss += critic_loss.item()
            
            # ---------------------
            # Train Generator
            # ---------------------
            optimizer_G.zero_grad()
            
            # Generate random noise
            z = torch.randn(batch_size, z_dim, device=device)
            
            # Generate fake data
            fake_data = generator(z)
            
            # Compute generator loss
            fake_validity = critic(fake_data)
            
            # Generator wants critic to think its outputs are real
            generator_loss = -torch.mean(fake_validity)
            
            # Backward pass
            generator_loss.backward()
            
            # Update generator
            optimizer_G.step()
            
            epoch_gen_loss += generator_loss.item()
            
            # Update progress bar
            progress_bar.set_postfix({
                'C Loss': critic_loss.item(),
                'G Loss': generator_loss.item()
            })
        
        # Average losses for this epoch
        avg_critic_loss = epoch_critic_loss / len(dataloader)
        avg_gen_loss = epoch_gen_loss / len(dataloader)
        
        critic_losses.append(avg_critic_loss)
        generator_losses.append(avg_gen_loss)
        
        print(f"Epoch {epoch+1}/{n_epochs} - Critic Loss: {avg_critic_loss:.4f}, Generator Loss: {avg_gen_loss:.4f}")
        
        # Generate and save a sample at the end of each epoch
        with torch.no_grad():
            z = torch.randn(1, z_dim, device=device)
            sample = generator(z).cpu().numpy()
            gen_samples.append(sample)
    
    # Plot losses
    plt.figure(figsize=(10, 5))
    plt.plot(critic_losses, label='Critic Loss')
    plt.plot(generator_losses, label='Generator Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    plt.title('Training Losses')
    plt.savefig('wgan_losses.png')
    plt.show()
    
    return gen_samples, critic_losses, generator_losses

# Run the training
if __name__ == "__main__":
    print(f"Starting WGAN training on device: {device}")
    gen_samples, critic_losses, generator_losses = train_wgan()
    
    # Save the trained models
    torch.save(generator.state_dict(), 'generator.pth')
    torch.save(critic.state_dict(), 'critic.pth')
    
    print("Training complete. Models saved.")
