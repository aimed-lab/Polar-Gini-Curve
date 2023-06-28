"""
Calculate the Gini coefficient and plot the Lorenz curve.
"""

import numpy as np
import matplotlib.pyplot as plt


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


def compute_rmsd(
    coordinate, cluster_id, expression, gene_list, marker="Actc1", cluster=1
):
    """
    This function computes the RMSD (Root Mean Square Deviation) between the Polar Gini Curves (PGC) of
    all cells in a cluster and cells in the same cluster expressing a given marker.

    Parameters:
    coordinate (numpy array): The 2D spatial coordinates.
    cluster_id (numpy array): The ID for each cluster.
    expression (numpy array): The gene expression data.
    gene_list (numpy array): The list of gene symbols.
    marker (str, optional): The marker gene to consider. Default is 'Actc1'.

    Returns:
    RMSD (float): Between the PGCs of all cells in the cluster and markered cells in the cluster.
    """

    cluster_index = np.where(cluster_id == cluster)[0]  # get all cells in cluster
    cluster_coor = coordinate[
        cluster_index, :
    ]  # get the spatial coordinate of cells in cluster
    cluster_lbl = np.ones(len(cluster_index))  # label these cells as '1'

    marker_index = np.where(gene_list == marker)[0]
    marker_expression = expression[:, marker_index]
    gene_index = np.where((cluster_id == cluster) & (marker_expression > 0))[
        0
    ]  # find cells expressing marker in cluster

    gene_coor = coordinate[
        gene_index, :
    ]  # get the spatial coordinate of cells expressing marker in cluster 1
    gene_lbl = 2 * np.ones(len(gene_index))  # label these cells as '2'

    # Draw the PGCs for these two sets of cells
    _, all_gini = make_2d_gini(
        np.vstack((gene_coor, cluster_coor)),
        np.hstack((gene_lbl, cluster_lbl)),
        [f"{marker} cluster {cluster} cells", f"All cluster {cluster} cells"],
    )

    try:
        rmsd = np.sqrt(
            np.mean(np.abs(all_gini[0] - all_gini[1]) ** 2)
        )  # the RMSD between two PGCs
    except:
        rmsd = np.nan

    # Round to 2 precision points
    return round(rmsd, 2)
