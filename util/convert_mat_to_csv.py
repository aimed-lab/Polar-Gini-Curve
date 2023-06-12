from scipy.io import loadmat
import numpy as np


def convert_mat_to_csv(mat_file, csv_file):
    data = loadmat(mat_file)

    for i in data:
        if "__" not in i and "readme" not in i:
            np.savetxt((csv_file), data[i], delimiter=",")

    print(f"Converted {mat_file} to {csv_file}")

    return
