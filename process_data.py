"""
Identify clusters using t-SNE and DBSCAN.
"""

import numpy as np
from sklearn.manifold import TSNE
from sklearn.metrics import pairwise_distances
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt


def identify_clusters(expression):
    """
    This function takes in gene expression data, computes the spatial coordinates using t-SNE,
    calculates pairwise distances, and then uses DBSCAN for clustering. It also plots a
    k-distance graph to help determine the optimal DBSCAN parameters.

    Parameters:
    expression (numpy array): The gene expression data.

    Returns:
    clusterID (numpy array): The ID for each cluster computed using DBSCAN.
    coordinate (numpy array): The 2D spatial coordinates computed using t-SNE.

    Note:
    This function also plots a k-distance graph to visualize the 50th nearest distances of points
    and help determine the suitable parameters for the DBSCAN algorithm.
    """

    # compute the spatial coordinate
    coordinate = TSNE().fit_transform(expression)

    # compute pairwise distances
    p_dist = pairwise_distances(coordinate, metric="euclidean")

    # find the suitable parameter for dbscan algorithm, with the minimum of neighborhood is 50
    kd_sorted = np.sort(p_dist, axis=1)[:, :50]

    # plot the parameter graph. It looks like the elbow point 3 is a good parameter
    plt.plot(np.sort(kd_sorted[-1, :]))
    plt.title("k-distance graph")
    plt.xlabel("Points sorted with 50th nearest distances")
    plt.ylabel("50th nearest distances")
    plt.show()

    # cluster ID
    cluster_id = DBSCAN(eps=3, min_samples=50).fit_predict(coordinate)

    return cluster_id, coordinate
