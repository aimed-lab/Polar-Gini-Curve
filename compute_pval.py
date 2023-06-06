import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy.stats import norm
from matplotlib import path
from draw_pgc import make2DGini

# TODO: Debug this - I have no idea what I'm doing
def monte_carlo_simulation(cluster_id, coordinate, gene_lbl, RMSD, num_exp=200):
    """
    This function performs a Monte Carlo simulation to generate a new cell-point
    clusters with the same shape and number of points to a specified cluster.
    The function calculates the RMSD for each simulation and generates a histogram of the RMSDs.
    It also computes a p-value comparing the given RMSD to the simulated RMSDs.

    Parameters:
    cluster_id (numpy array): The ID for each cluster.
    coordinate (numpy array): The 2D spatial coordinates.
    gene_lbl (numpy array): The labels for cells expressing the marker gene.
    RMSD (float): The Root Mean Square Deviation of the actual data.
    num_exp (int, optional): The number of simulations to run. Default is 200.

    Returns:
    allRanRMSD (numpy array): The RMSDs for each simulation.
    p_val (float): The computed p-value.
    """

    point_idx = np.where(cluster_id == 1)
    cluster_pt = coordinate[point_idx]

    if len(cluster_pt.shape) != 2 or cluster_pt.shape[1] != 2:
        raise ValueError(
            f"Unexpected shape for cluster_pt: {cluster_pt.shape}. Expected (N, 2)."
        )

    hull = ConvexHull(cluster_pt)
    bound_pt = cluster_pt[hull.vertices]

    x_range = [np.min(bound_pt[:, 0]), np.max(bound_pt[:, 0])]
    y_range = [np.min(bound_pt[:, 1]), np.max(bound_pt[:, 1])]

    all_ran_RMSD = np.ones(num_exp)

    for expr in range(num_exp):
        back_pt = []
        while len(back_pt) < len(point_idx):
            rand_x = x_range[0] + (x_range[1] - x_range[0]) * np.random.rand()
            rand_y = y_range[0] + (y_range[1] - y_range[0]) * np.random.rand()

            if path.Path(bound_pt).contains_points([(rand_x, rand_y)]):
                back_pt.append([rand_x, rand_y])

        back_pt = np.array(back_pt)
        back_lbl = np.ones(len(point_idx))

        fore_pt = back_pt[: len(gene_lbl)]
        fore_lbl = 2 * np.ones(len(gene_lbl))

        _, all_gini = make2DGini(
            np.vstack((fore_pt, back_pt)),
            np.hstack((fore_lbl, back_lbl)),
            ["fore", "back"],
        )

        all_ran_RMSD[expr] = np.sqrt(np.mean(np.abs(all_gini[0] - all_gini[1]) ** 2))

    plt.hist(all_ran_RMSD, 100)
    plt.show()

    p_val = 1 - norm.cdf(RMSD, np.mean(all_ran_RMSD), np.std(all_ran_RMSD))

    return all_ran_RMSD, p_val
