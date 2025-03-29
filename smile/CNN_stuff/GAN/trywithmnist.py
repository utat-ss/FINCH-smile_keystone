import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from torchvision import datasets, transforms
import matplotlib.pyplot as plt
import numpy as np

# Set random seed for reproducibility
torch.manual_seed(42)

# Hyperparameters
batch_size = 64
lr = 0.00005
n_epochs = 50
clip_value = 0.01
n_critic = 5  # Number of critic iterations per generator iteration
z_dim = 100  # Dimension of the latent space
image_size = 28  # MNIST image size
image_channels = 1  # Grayscale images for MNIST

# Device configuration
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# Data Loading
transform = transforms.Compose([
    transforms.Resize(image_size),
    transforms.ToTensor(),
    transforms.Normalize([0.5], [0.5])  # Normalize to [-1, 1]
])

# MNIST dataset
mnist_dataset = datasets.MNIST(root='./data', train=True, transform=transform, download=True)
dataloader = DataLoader(mnist_dataset, batch_size=batch_size, shuffle=True, num_workers=2)

# Generator network
class Generator(nn.Module):
    def __init__(self):
        super(Generator, self).__init__()
        
        self.model = nn.Sequential(
            # Input is Z, going into a convolution
            nn.Linear(z_dim, 128 * 7 * 7),
            nn.ReLU(True),
            nn.Unflatten(1, (128, 7, 7)),
            
            nn.ConvTranspose2d(128, 64, 4, 2, 1, bias=False),  # Output: 64 x 14 x 14
            nn.BatchNorm2d(64),
            nn.ReLU(True),
            
            nn.ConvTranspose2d(64, image_channels, 4, 2, 1, bias=False),  # Output: 1 x 28 x 28
            nn.Tanh()  # Output range: [-1, 1]
        )

    def forward(self, z):
        return self.model(z)

# Critic network
class Critic(nn.Module):
    def __init__(self):
        super(Critic, self).__init__()
        
        self.model = nn.Sequential(
            # Input: 1 x 28 x 28
            nn.Conv2d(image_channels, 64, 4, 2, 1, bias=False),  # Output: 64 x 14 x 14
            nn.LeakyReLU(0.2, inplace=True),
            
            nn.Conv2d(64, 128, 4, 2, 1, bias=False),  # Output: 128 x 7 x 7
            nn.LayerNorm([128, 7, 7]),
            nn.LeakyReLU(0.2, inplace=True),
            
            nn.Flatten(),
            nn.Linear(128 * 7 * 7, 1)
            # No sigmoid at the end because WGAN doesn't use BCE loss
        )

    def forward(self, img):
        return self.model(img)

# Initialize the networks
generator = Generator().to(device)
critic = Critic().to(device)

# Optimizers
optimizer_G = optim.RMSprop(generator.parameters(), lr=lr)
optimizer_C = optim.RMSprop(critic.parameters(), lr=lr)

# Training function
def train_wgan():
    # Lists to store generated images for visualization
    gen_imgs = []
    
    # Lists to store losses for plotting
    critic_losses = []
    generator_losses = []
    
    for epoch in range(n_epochs):
        for i, (real_imgs, _) in enumerate(dataloader):
            # Configure input
            real_imgs = real_imgs.to(device)
            
            # ---------------------
            #  Train Critic
            # ---------------------
            optimizer_C.zero_grad()
            
            # Sample noise as generator input
            z = torch.randn(real_imgs.size(0), z_dim).to(device)
            
            # Generate a batch of images
            fake_imgs = generator(z)
            
            # Compute outputs
            real_validity = critic(real_imgs)
            fake_validity = critic(fake_imgs.detach())
            
            # Compute loss (Wasserstein loss)
            critic_loss = -torch.mean(real_validity) + torch.mean(fake_validity)
            
            # Backward pass
            critic_loss.backward()
            optimizer_C.step()
            
            # Clip weights of critic
            for p in critic.parameters():
                p.data.clamp_(-clip_value, clip_value)
            
            # Train generator every n_critic iterations
            if i % n_critic == 0:
                # ---------------------
                #  Train Generator
                # ---------------------
                optimizer_G.zero_grad()
                
                # Generate a batch of images
                fake_imgs = generator(z)
                
                # Compute output
                fake_validity = critic(fake_imgs)
                
                # Compute loss
                generator_loss = -torch.mean(fake_validity)
                
                # Backward pass
                generator_loss.backward()
                optimizer_G.step()
                
                # Store losses for plotting
                critic_losses.append(critic_loss.item())
                generator_losses.append(generator_loss.item())
                
                print(
                    f"[Epoch {epoch}/{n_epochs}] [Batch {i}/{len(dataloader)}] "
                    f"[D loss: {critic_loss.item():.4f}] [G loss: {generator_loss.item():.4f}]"
                )
        
        # Save a sample of generated images
        if epoch % 5 == 0:
            with torch.no_grad():
                z = torch.randn(16, z_dim).to(device)
                sample_imgs = generator(z).detach().cpu()
                gen_imgs.append(sample_imgs)
    
    return gen_imgs, critic_losses, generator_losses

# Visualization function
def visualize_results(gen_imgs, critic_losses, generator_losses):
    # Plot losses
    plt.figure(figsize=(10, 5))
    plt.plot(critic_losses, label='Critic loss')
    plt.plot(generator_losses, label='Generator loss')
    plt.xlabel('Iterations')
    plt.ylabel('Loss')
    plt.legend()
    plt.title('WGAN Training Losses')
    plt.savefig('wgan_losses.png')
    plt.show()
    
    # Plot generated images
    for i, imgs in enumerate(gen_imgs):
        epoch = i * 5
        plt.figure(figsize=(10, 6))
        for j in range(min(16, len(imgs))):
            plt.subplot(4, 4, j+1)
            plt.imshow(imgs[j][0].numpy(), cmap='gray_r')
            plt.axis('off')
        plt.tight_layout()
        plt.savefig(f'wgan_epoch_{epoch}.png')
        plt.show()

# Train the model
if __name__ == "__main__":
    gen_imgs, critic_losses, generator_losses = train_wgan()
    visualize_results(gen_imgs, critic_losses, generator_losses)