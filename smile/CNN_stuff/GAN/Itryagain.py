import os 

# This is necessary because it repackages data into a format that the neural network can understand 
import torchvision.transforms as transforms 

# Creates a grid of images, used to display a bunch of images 
from torchvision.utils import make_grid

# Load in a lot of data, useful for parallel processing (parallel data loading, I think big improvements only in GPUs) 
from torch.utils.data import DataLoader

# pretty sure I don't need this cas the data is local, but useful for prototyping if I need to 
from torchvision import datasets

# Important building block for neural networks in here (super inheritance) 
import torch.nn as nn 

# provides a bunch of activation functions, not sure why this isn't included in torch.nn but it is what it is
import torch.nn.functional as F

# import the torch library
import torch 

import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.pyplot import figure 

from tqdm import tqdm

# Makes matplotlib plots interactive (did not know you could do that before this), good for real time visualization (does not block code execution) 
plt.ion() 
from IPython.display import clear_output


# Some of the hyperparameters (probably not optimal cas I'm just seeing what other people use)
batch_size = 64 # How many samples used per time
lr = 0.00005 # Learning rate
n_epochs = 5000 # Number of epochs, obviously, increase it if you're creating the final model
clip_value = 1 # This is the 1-Lipschitz thing, and in this case, it means the gradient cannot be greator than 1 (I'm pretty sure its 1-Lipschitz, but some people use 0.01???)
n_critic = 5  # Number of critic iterations per generator iteration. For WGANs, you need the critic to update more frequently than the generator (I'm guessing this is due to the much higher complexity of the critic compared to the genrator)
z_dim = 256  # Dimension of the latent, you can this of this as having multiple bags of marbles, the more bags, the more combinations of marbles you can make (this controls the compexity of the data that your generator can capture)

input_size = 3 # This is the number of abundancies that we can put in (gv, npv and soil) to predict the spectra. Might change it to use spectra to preduct abundancies later on, but this model will be simpler (I hope) 

# This is for when we use google colab to do this thing, I don't have a nvidia GPU though, so cpu for now
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Loading in the data 
simpler_data_location = r"C:\Users\liche\OneDrive\Desktop\PycharmProjects\FINCH-smile_keystone\smile\CNN_stuff\GAN\simpler_data.csv" # This is just the simpler_data file with the abundancies and spectra, I removed the "use" column because it serves no use :D

# Data Preprocessing Pipeline
data_transform = transforms.Compose([
    transforms.ToTensor(),

    # Normalization might help reduce the difficulty that the generator has in the training process, also ensures the real images and the generator are in the same format, I'm pretty sure its important
    transforms.Normalize([0.5], [0.5])   # Normalize to [-1, 1], all the data is from 0 to 1, so this should fit pretty nicely (the numbers resulting from this normalization will be anything from -1 to 1)
    # Same output range as the tanh function, so this is usable for the generator's last layer
])



# Technically, this should be on another python file, but I'm keeping it here for simplicity and to see if it even works
# IMPORTANT: The final output needs to be 213 columns, because there are 213 columns in total, in this GAN, I am lumping all the data together, so the generator comes up with the abundance and spectra in tandem, hopefully reinforces more connection between the two
# This 213 columns will then be passed into the critic (called the discriminator in traditional GANs) to see if it can tell the difference between the real and fake data, and this is used to train the generator
class Generator(nn.Module): # This is the actual generator model, there isn't much to this code, and its basically just a standard neural network, and its the training loop that does the cool stuff
    # No idea what this args kwargs stuff is, but it autocompleted so in it goes (its the inputs, never seen it represented it like that though)
    def __init__(self, *args, **kwargs):
        super(Generator, self).__init__(*args, **kwargs)

        # Trying out random crap in the neural network, will require a lot of tweaking so it better fits the data
        self.model = nn.Sequential(          
            # WGANs and GANs in general take random noise and try to generate something that makes sense from it, the random inpput is the z stuff  
            # Input is Z, going into the neural network 
            nn.Linear(in_features=z_dim, out_features=200), 

            # The inplace = True is something that reduces the amount of memory used (very slightly) by destroying the original input used. Might speed it up a tiny bit. It should be applicable in this use case because later layers never refer back to a previous layer
            # No idea which activation function to use, so LeakyReLu to test it out
            nn.ReLU(inplace=True),

            # This controls how many neurons are no longer considered, a lot of tweaking required for this value, but this is typically used to prevent overfitting (higher value, less overfitting)
            nn.Dropout(0.05),
            
            # 200 -> 400 -> 250 -> 212
            # Completely random in_features and out_features, change these depending on what works best, this is just a hopefully this works thing
            nn.Linear(in_features=200, out_features=400),
            nn.ReLU(inplace=True),
            nn.Dropout(0.05),

            nn.Linear(in_features=400, out_features=400),
            nn.ReLU(inplace=True),
            nn.Dropout(0.05),

            nn.Linear(in_features=400, out_features=250),
            nn.ReLU(inplace=True),
            nn.Dropout(0.05),
            
            
            nn.Linear(in_features=250, out_features=212),
            nn.ReLU(inplace=True),
            nn.Tanh(),  # Added tanh activation to pass all values through the tanh function (if you remember from earlier, we normalized it so this would match, I think it should make it easier now, but this could also be a really garbage approach, we'll figure it out in the future!)

            # This is the output layer, and it should be the same size as the real data, because the generator is trying to generate data that looks like the real data
        )

    # This is the forward function, and it is what is called when you pass data through the neural network
    def forward(self, z):
        return self.model(z)
    
# This is the critic, and it is the thing that tries to tell the difference between the real and fake data
class Critic(nn.Module):
    def __init__(self, *args, **kwargs):
        super(Critic, self).__init__(*args, **kwargs)

        # This is the neural network, and it is the thing that tries to tell the difference between the real and fake data
        self.model = nn.Sequential(
            # Input: 212 important columns, this is the same as the generator's output columns, and it is necessary to have the same input size as the generator
            nn.Linear(in_features=212, out_features=200),
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
            # No sigmoid at the end because WGANs don't use BCE loss, they use Earth Mover's Distance, which is not a value between 0 or 1, it is a random number who's difference with another number is what's important 
        )

    # This is the forward function, and it is what is called when you pass data through the neural network
    def forward(self, x):
        return self.model(x)
    
    # Three abundances are what the discriminator predicts 

# This part will be important when we have a GPU (or TPU I think?), but for now, we're just using the CPU, so this code doesn't do much (except send it to the CPU, which I guess is doing a lot...)
generator = Generator().to(device)
critic = Critic().to(device)


def gradient_penalty(critic, real, fake):
    """
    Calculate the gradient penalty for WGAN-GP.
    Absolutely no idea what this does, but its necessary for wgan-gp, so we ball

    Args:
        critic: The critic model
        real: Batch of real data
        fake: Batch of generated data
        
    Returns:
        The gradient penalty term
    """
    # Get batch size and device
    batch_size, _ = real.shape
    
    # Create random interpolation factors for each sample in the batch
    alpha = torch.rand(batch_size, 1).to(device)
    
    # Create interpolated samples between real and fake data
    interpolated = alpha * real + (1 - alpha) * fake
    interpolated = interpolated.requires_grad_(True)
    
    # Get critic scores for the interpolated samples
    critic_interpolated_scores = critic(interpolated)
    
    # Calculate gradients of critic scores with respect to the interpolated samples
    gradients = torch.autograd.grad(
        outputs=critic_interpolated_scores,
        inputs=interpolated,
        grad_outputs=torch.ones_like(critic_interpolated_scores).to(device),
        create_graph=True,
        retain_graph=True,
    )[0]
    
    # Calculate the gradient penalty
    # The L2 norm of the gradients should be close to 1
    gradients = gradients.view(batch_size, -1)
    gradient_norm = gradients.norm(2, dim=1)
    gradient_penalty = ((gradient_norm - 1) ** 2).mean()
    
    return gradient_penalty


# This is where the real magic happens!! This is the training function and implements the Wasserstein loss function, which is the thing that makes WGANs different from GANs (btw, I have no idea if I'm implementing it correctly) 
def train_wgan():
    # Optimizers for the critic and generator
    critic_optimizer = torch.optim.Adam(critic.parameters(), lr=lr, betas=(0.0, 0.9))
    generator_optimizer = torch.optim.Adam(generator.parameters(), lr=lr, betas=(0.0, 0.9))
    
    # Load the data using pandas to handle mixed data types
    import pandas as pd
    df = pd.read_csv(simpler_data_location)
    
    # Skip the first column and any non-numeric columns
    numeric_data = df.select_dtypes(include=['float64', 'int64']).iloc[:, 1:].values
    
    # Convert to tensor
    data_tensor = torch.tensor(numeric_data, dtype=torch.float32).to(device)
    
    # Print data shape for debugging
    print(f"Data shape after processing: {data_tensor.shape}")
    
    # Create DataLoader
    dataset = torch.utils.data.TensorDataset(data_tensor)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
    
    # Lambda for gradient penalty (standard value is 10)
    lambda_gp = 10
    
    # Losses for the critic and the generator are stored here for plotting at the end 
    critic_losses = []
    generator_losses = []
    
    # For tracking progress
    total_batches = len(dataloader)
    
    # How many training cycles will we run on the thing, n_epochs in the hyperparameters at the top
    for epoch in range(n_epochs):
        # Track losses for this epoch
        epoch_critic_losses = []
        epoch_generator_losses = []
        
        # Progress bar for this epoch
        progress_bar = tqdm(enumerate(dataloader), total=total_batches, desc=f"Epoch {epoch+1}/{n_epochs}")
        
        for i, (real_data,) in progress_bar:
            # Get batch size (might be smaller for the last batch)
            current_batch_size = real_data.size(0)
            
            # ---------------------
            # Train Critic
            # ---------------------
            for _ in range(n_critic):
                # Zero gradients
                critic_optimizer.zero_grad()
                
                # Generate noise
                z = torch.randn(current_batch_size, z_dim).to(device)
                
                # Generate fake data
                fake_data = generator(z)
                
                # Critic scores
                real_scores = critic(real_data)
                fake_scores = critic(fake_data.detach())  # Detach to avoid training generator
                
                # Compute WGAN loss with gradient penalty
                # WGAN loss: maximize E[critic(real)] - E[critic(fake)]
                # Or equivalently, minimize E[critic(fake)] - E[critic(real)]
                gp = gradient_penalty(critic, real_data, fake_data.detach())
                critic_loss = fake_scores.mean() - real_scores.mean() + lambda_gp * gp
                
                # Backpropagation
                critic_loss.backward()
                critic_optimizer.step()
                
                # Store loss
                epoch_critic_losses.append(critic_loss.item())
            
            # ---------------------
            # Train Generator
            # ---------------------
            # Zero gradients
            generator_optimizer.zero_grad()
            
            # Generate noise
            z = torch.randn(current_batch_size, z_dim).to(device)
            
            # Generate fake data
            fake_data = generator(z)
            
            # Critic scores for fake data
            fake_scores = critic(fake_data)
            
            # Generator loss: minimize -E[critic(fake)]
            # This is equivalent to maximizing E[critic(fake)]
            generator_loss = -fake_scores.mean()
            
            # Backpropagation
            generator_loss.backward()
            generator_optimizer.step()
            
            # Store loss
            epoch_generator_losses.append(generator_loss.item())
            
            # Update progress bar
            progress_bar.set_postfix({
                'C_Loss': f"{critic_loss.item():.4f}",
                'G_Loss': f"{generator_loss.item():.4f}"
            })
        
        # Average losses for this epoch
        avg_critic_loss = sum(epoch_critic_losses) / len(epoch_critic_losses)
        avg_generator_loss = sum(epoch_generator_losses) / len(epoch_generator_losses)
        
        # Store average losses
        critic_losses.append(avg_critic_loss)
        generator_losses.append(avg_generator_loss)
        
        # Print epoch results
        print(f"Epoch [{epoch+1}/{n_epochs}] | Critic Loss: {avg_critic_loss:.4f} | Generator Loss: {avg_generator_loss:.4f}")
        
        # Visualize some generated data every few epochs
        if (epoch + 1) % 100 == 0 or epoch == 0:
            with torch.no_grad():
                # Generate some samples
                z = torch.randn(5, z_dim).to(device)
                generated_samples = generator(z).cpu().numpy()
                
                # Plot the first 3 columns (abundances) and a few spectral bands
                fig, axs = plt.subplots(2, 1, figsize=(12, 10))
                
                # Plot abundances
                axs[0].bar(range(3), generated_samples[0, :3], color='blue', alpha=0.7)
                axs[0].set_title('Generated Abundances (First Sample)')
                axs[0].set_xticks(range(3))
                axs[0].set_xticklabels(['GV', 'NPV', 'Soil'])
                axs[0].set_ylim(-1, 1)  # Since we're using tanh
                
                # Plot spectral bands (just a subset)
                spectral_indices = range(3, min(50, generated_samples.shape[1]))
                for i in range(min(5, generated_samples.shape[0])):
                    axs[1].plot(spectral_indices, generated_samples[i, spectral_indices], 
                                label=f'Sample {i+1}', alpha=0.7)
                
                axs[1].set_title('Generated Spectral Bands (First 50 bands)')
                axs[1].set_xlabel('Band Index')
                axs[1].set_ylabel('Value')
                axs[1].legend()
                
                plt.tight_layout()
                plt.savefig(f'generated_samples_epoch_{epoch+1}.png')
                plt.close()
    
    # Plot the losses
    plt.figure(figsize=(10, 6))
    plt.plot(critic_losses, color='blue', label="Critic Loss", zorder=2)
    plt.plot(generator_losses, color='red', label="Generator Loss", zorder=1)
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('WGAN Training Losses')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('wgan_training_losses.png')
    plt.show()
    
    # Save the models
    torch.save(generator.state_dict(), 'wgan_generator.pth')
    torch.save(critic.state_dict(), 'wgan_critic.pth')
    
    return generator, critic, critic_losses, generator_losses

if __name__ == "__main__":
    print("Starting WGAN training...")
    generator, critic, critic_losses, generator_losses = train_wgan()
    print("Training completed successfully!")