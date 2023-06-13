"""
Convert a MATLAB .mat file to a CSV file.
"""

from scipy.io import loadmat
import numpy as np


def convert_mat_to_csv(mat_file, csv_file):
    """
    Convert a MATLAB .mat file to a CSV file.

    This function reads the MATLAB .mat file specified by 'mat_file' and converts
    its contents to a CSV file specified by 'csv_file'. The function saves each
    variable from the .mat file as a separate CSV file. The variable names are used
    as the file names for the CSV files.

    Parameters:
    mat_file (str): Path to the input .mat file.
    csv_file (str): Path to the output CSV file.

    Returns:
    None. The function saves the converted CSV files.

    Example:
    >>> convert_mat_to_csv('data.mat', 'data.csv')
    Converted data.mat to data.csv
    """

    data = loadmat(mat_file)

    for i in data:
        if "__" not in i and "readme" not in i:
            np.savetxt((csv_file), data[i], delimiter=",")

    print(f"Converted {mat_file} to {csv_file}")

    return
