import numpy as np
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def perform_pca(data, n_components):
    """
    Perform Principal Component Analysis (PCA) on the given dataset.
    
    Parameters:
        data (numpy.ndarray): The input dataset with shape (n_samples, n_features).
        n_components (int): Number of principal components to retain.
    
    Returns:
        transformed_data (numpy.ndarray): Data projected onto principal components.
        explained_variance (numpy.ndarray): Variance explained by each principal component.
        components (numpy.ndarray): Principal component vectors.
    """
    # Standardize the dataset
    scaler = StandardScaler()
    standardized_data = scaler.fit_transform(data)
    
    # Perform PCA
    pca = PCA(n_components=n_components)
    transformed_data = pca.fit_transform(standardized_data)
    
    return transformed_data, pca.explained_variance_ratio_, pca.components_


def plot_pca(transformed_data, mode="2d", marker_size=2):
    """
    Create an interactive PCA plot using Plotly.
    
    Parameters:
        transformed_data (numpy.ndarray): Data projected onto principal components.
        mode (str): "2d" for a 2D scatter plot, "3d" for a 3D scatter plot.
        marker_size (int): Size of the markers in the scatter plot.
    """
    if mode == "3d" and transformed_data.shape[1] >= 3:
        fig = px.scatter_3d(
            x=transformed_data[:, 0], y=transformed_data[:, 1], z=transformed_data[:, 2],
            title="3D PCA Projection", labels={"x": "PC1", "y": "PC2", "z": "PC3"},
            # size_max=1  # This sets the maximum marker size
        )
        # Alternatively, you can use:
        fig.update_traces(marker=dict(size=marker_size))
    else:
        fig = px.scatter(
            x=transformed_data[:, 0], y=transformed_data[:, 1],
            title="2D PCA Projection", labels={"x": "PC1", "y": "PC2"},
            size_max=marker_size  # This sets the maximum marker size
        )
        # Alternatively, you can use:
        # fig.update_traces(marker=dict(size=marker_size))
    
    fig.show()

# Whats the data that is being given to the PCAing Function?

data = "SOMETHING GOES HERE"

transformed_data, explained_variance, components = perform_pca(data, 3)

print("Transformed Data:\n", transformed_data[:5])  # Print first 5 samples
print("Explained Variance Ratio:\n", explained_variance)
print("Principal Components:\n", components)

# Plot PCA results (user can choose "2d" or "3d")
plot_pca(transformed_data, mode="3d", marker_size=2)

