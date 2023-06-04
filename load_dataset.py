"""
Loads MATLAB datasets.
"""

from scipy.io import loadmat


def load_dataset(path):
    """
    Loads the 2D spatial coordinates, gene expressions,
    gene symbols, and cluster IDs from MATLAB .mat files.
    """
    name = path.split("/")[-1].split(".")[0]
    # print(f"Loading {name} dataset from {path}...")
    return loadmat(path)[name]
