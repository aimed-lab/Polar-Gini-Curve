import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.io import loadmat
from joblib import Parallel, delayed
import matplotlib.pyplot as plt

# Load the data
idx = loadmat("data/clusterID.mat")["clusterID"]
express_mat = loadmat("data/Expression.mat")["Expression"]
gene_list = loadmat("data/geneList.mat")["geneList"]
coordinate = loadmat("data/coordinate.mat")["coordinate"]

express_mat = express_mat.T
max_idx = np.max(idx)
percen_exp = np.zeros((len(gene_list), max_idx))
matrix_result = np.empty((len(gene_list), max_idx))
matrix_result[:] = np.NaN


def compute_gini(pop, val, makeplot=False):
    """
    Compute the Gini coefficient and plot the Lorenz curve.

    Parameters:
    pop (array-like): Population vector.
    val (array-like): Value vector.
    makeplot (bool, optional): Whether to plot the Lorenz curve. Default is False.

    Returns:
    gini (float): Gini coefficient.
    lorenz (ndarray): Lorenz curve.
    curve_points (ndarray): Coordinate points of the Lorenz curve.
    """

    assert len(pop) == len(val), "compute_gini expects two equally long vectors."

    pop = np.append(0, pop)  # pre-append a zero
    val = np.append(0, val)  # pre-append a zero

    isok = ~np.isnan(pop) & ~np.isnan(val)  # filter out NaNs
    if np.sum(isok) < 2:
        print("Warning: Not enough data")
        return np.nan, np.nan, np.nan

    pop = pop[isok]
    val = val[isok]

    assert np.all(pop >= 0) and np.all(
        val >= 0
    ), "compute_gini expects nonnegative vectors."

    # process input
    weighted = val * pop
    sorted_indices = np.argsort(val)
    pop = pop[sorted_indices]
    weighted = weighted[sorted_indices]
    pop = np.cumsum(pop)
    weighted = np.cumsum(weighted)
    relpop = pop / pop[-1]
    relz = weighted / weighted[-1]

    # Gini coefficient
    gini = 1 - np.sum((relz[:-1] + relz[1:]) * np.diff(relpop))

    # Lorentz curve
    lorenz = np.column_stack([relpop, relz])
    curve_points = np.column_stack([pop, weighted])

    if makeplot:  # ... plot it?
        plt.fill_between(relpop, relz, color=[0.5, 0.5, 1.0])  # the Lorentz curve
        plt.plot([0, 1], [0, 1], "--k")  # 45 degree line
        plt.axis(
            "tight"
        )  # ranges of abscissa and ordinate are by definition exactly [0,1]
        plt.axis("equal")  # both axes should be equally long
        plt.grid()
        plt.title(f"Gini coefficient = {gini}")
        plt.xlabel("Share of population")
        plt.ylabel("Share of value")
        plt.show()

    return gini, lorenz, curve_points


def make_2d_gini(x, cluster_id, cluster_name=None):
    """
    Compute the Gini coefficient for 2D data points in different clusters.

    Parameters:
    x (ndarray): 2D data points.
    cluster_id (ndarray): Cluster ID for each data point.
    cluster_name (list, optional): Names of the clusters. Default is None.

    Returns:
    angle_list (ndarray): List of angles.
    all_gini (list): List of Gini coefficients for each cluster.
    """

    assert x.shape[1] == 2, "We need 2D data."

    num_cluster = len(np.unique(cluster_id))

    # automatically set the cluster name if needed
    if cluster_name is None:
        cluster_name = ["cluster " + str(i) for i in range(1, num_cluster + 1)]

    # get the angle list with resolution of 1000
    resolution = 1000
    angle_list = np.pi * np.linspace(0, 360, resolution) / 180

    all_gini = [None] * num_cluster

    # for each cluster, get the list of gini corresponding to each angle
    for cluster in range(1, num_cluster + 1):
        coordinate = x[cluster_id == cluster]
        gini = np.zeros(len(angle_list))

        for i, angle in enumerate(angle_list):
            value = np.dot(coordinate, [np.cos(angle), np.sin(angle)])
            value -= np.min(value)
            gini[i] = compute_gini(np.ones(len(value)), value)[0]

        all_gini[cluster - 1] = gini

    return angle_list, all_gini


def compute_matrix(i, gene):
    results = []
    for j in range(max_idx):
        cluster_index = np.where(idx == j + 1)[0]
        cluster_coor = coordinate[cluster_index, :]
        cluster_lbl = np.ones(len(cluster_index), dtype=int)

        gene_index = np.where((idx == j + 1) & (express_mat[:, i] > 0))[0]
        if len(gene_index) > 0:
            gene_coor = coordinate[gene_index, :]
            gene_lbl = np.full(len(gene_index), 2, dtype=int)

            angle_list, all_gini = make_2d_gini(
                np.concatenate((gene_coor, cluster_coor)),
                np.concatenate((gene_lbl, cluster_lbl)),
                [gene_list[i], "cluster " + str(j + 1)],
            )

            y = np.square(np.abs(all_gini[0] - all_gini[1]))
            results.append((i, j, np.sqrt(np.mean(y))))
    return results


# Use joblib to parallelize the computation
results = Parallel(n_jobs=-1)(
    delayed(compute_matrix)(i, gene) for i, gene in enumerate(gene_list)
)
for result in results:
    for i, j, value in result:
        matrix_result[i, j] = value


# Compute expressions and percentages
for i in range(len(gene_list)):
    for j in range(
        1, max_idx + 1
    ):  # Indices start from 1 in the cluster id array (idx)
        cluster_index = np.where(idx == j)[0]  # Get the indices for the current cluster
        percen_exp[i, j - 1] = len(
            np.where(express_mat[i, cluster_index] > 0)[0]
        ) / len(cluster_index)

# Initialize variables for normalization
normalized_result = np.ones(matrix_result.shape)
p_val = np.ones((len(gene_list), max_idx))

# Compute normalization and p-values
for j in range(max_idx):
    all_cluster_RSMD = matrix_result[:, j]
    index = np.where(percen_exp[:, j] > 0.1)[0]
    p_value_able = all_cluster_RSMD[index]
    max_score = np.max(p_value_able)

    normalized_result[index, j] = matrix_result[index, j] / max_score
    p_value_able = p_value_able / max_score

    mu = np.mean(p_value_able)
    sigma = np.std(p_value_able)
    p = norm.cdf(p_value_able, mu, sigma)
    p_val[index, j] = p

flattened_gene_list = [item[0][0] for item in gene_list]

# Create a dictionary to store the data
data_dict = {
    "Gene": np.repeat(flattened_gene_list, max_idx),
    "Cluster": np.tile(np.arange(1, max_idx + 1), len(gene_list)),
    "percentage_cell_expressing": percen_exp.flatten(order="F"),
    "raw_RSMD": matrix_result.flatten(order="F"),
    "normalized_RSMD": normalized_result.flatten(order="F"),
    "P_value": p_val.flatten(order="F"),
}

# Convert the dictionary to a DataFrame
final_output = pd.DataFrame(data_dict)

# Sort the DataFrame by Gene and Cluster
final_output = final_output.sort_values(by=["Gene", "Cluster"])

# Save the DataFrame to a CSV file
final_output.to_csv("RSMD_clusters.csv", index=False)
