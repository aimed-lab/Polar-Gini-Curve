"""
Micellaneous functions for plotting.
"""

from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import numpy as np

from draw_pgc import compute_gini


def draw_tsne(
    marker_gene,
    coordinate,
    cluster_id,
    target_cluster_id,
    expression_data,
    gene_list,
    random_state=0,
):
    """
    Plot t-SNE visualization for a marker gene in a specific cluster.

    Parameters:
        marker_gene (str): The marker gene to visualize.
        coordinate (ndarray): The 2D spatial coordinates.
        cluster_id (ndarray): The cluster ID for each cell.
        target_cluster_id (int): The ID of the target cluster.
        expression_data (ndarray): The gene expression data.
        gene_list (ndarray): The list of gene symbols.
        random_state (int, optional): Random state for t-SNE. Default is 0.

    Returns:
        None. The plot is displayed and saved as a PNG file.
    """

    # Create a t-SNE model
    model = TSNE(n_components=2, random_state=random_state)

    # Flatten cluster_id
    cluster_id = np.ravel(cluster_id)

    # Calculate t-SNE coordinates for all cells
    tsne_coordinates = model.fit_transform(coordinate)

    # Set up the plot
    _, ax = plt.subplots()

    # Plot target cluster in grey
    target_cluster = tsne_coordinates[cluster_id == target_cluster_id]
    ax.scatter(
        target_cluster[:, 0], target_cluster[:, 1], color="grey", s=1, alpha=0.25
    )

    # Find cells in the target cluster that express the marker gene
    marker_gene_index = np.where(gene_list == marker_gene)[
        0
    ]  # find index of marker gene
    marker_gene_expression = (
        expression_data[:, marker_gene_index].flatten() > 0
    )  # find cells that express the marker gene
    marker_gene_in_target_cluster = np.logical_and(
        cluster_id == target_cluster_id, marker_gene_expression
    )

    # Plot these cells in red
    marker_gene_coordinates = tsne_coordinates[marker_gene_in_target_cluster.ravel(), :]
    ax.scatter(
        marker_gene_coordinates[:, 0],
        marker_gene_coordinates[:, 1],
        color="red",
        s=1,
    )

    # Labels and title
    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2")
    ax.set_title(f"t-SNE for marker gene {marker_gene} in cluster {target_cluster_id}")

    # Save figure in png format
    plt.savefig(f"tsne_{marker_gene}_c-{target_cluster_id}.png")
    plt.show()


def plot_gini(
    coordinate, cluster_id, expression, gene_list, marker="Actc1", target_cluster=1
):
    """
    Plot Gini coefficients for each cluster.

    Parameters:
    coordinate (numpy array): The 2D spatial coordinates.
    cluster_id (numpy array): The ID for each cluster.
    expression (numpy array): The gene expression data.
    gene_list (numpy array): The list of gene symbols.
    marker (str, optional): The marker gene to consider. Default is 'Actc1'.
    target_cluster (int, optional): The cluster number to plot. Default is 1.

    Returns:
    None. The plot is displayed.
    """

    cluster_id = cluster_id.flatten()  # Flatten the cluster_id array

    cluster_index = np.where(cluster_id == target_cluster)[
        0
    ]  # get all cells in target_cluster
    cluster_coor = coordinate[
        cluster_index, :
    ]  # get the spatial coordinate of cells in target_cluster
    cluster_lbl = np.ones(len(cluster_index))  # label these cells as '1'

    marker_index = np.where(gene_list == marker)[0]
    marker_expression = expression[:, marker_index]
    gene_index = np.where((cluster_id == target_cluster) & (marker_expression > 0))[
        0
    ]  # find cells expressing marker in target_cluster

    gene_coor = coordinate[
        gene_index, :
    ]  # get the spatial coordinate of cells expressing marker in target_cluster
    gene_lbl = 2 * np.ones(len(gene_index))  # label these cells as '2'

    x = np.vstack((gene_coor, cluster_coor))
    cluster_id = np.hstack((gene_lbl, cluster_lbl))
    cluster_name = [
        f"{marker} cluster {target_cluster} cells",
        f"All cluster {target_cluster} cells",
    ]

    assert x.shape[1] == 2, "We need 2D data."

    num_cluster = len(np.unique(cluster_id))

    # automatically set the cluster name if needed
    if cluster_name is None:
        cluster_name = ["cluster " + str(i) for i in range(1, num_cluster + 1)]

    # get the angle list with resolution of 1000
    resolution = 1000
    angle_list = np.pi * np.linspace(0, 360, resolution) / 180

    # for each cluster, get the list of gini corresponding to each angle
    for cluster in range(1, num_cluster + 1):
        coordinate = x[cluster_id == cluster]
        gini = np.zeros(len(angle_list))

        for i, angle in enumerate(angle_list):
            value = np.dot(coordinate, [np.cos(angle), np.sin(angle)])
            value -= np.min(value)
            gini[i] = compute_gini(np.ones(len(value)), value)[0]

        # polar plot
        plt.polar(angle_list, gini)
        plt.title("Gini Coefficient")

    plt.legend(cluster_name)
    plt.savefig(f"gini_{marker}_c-{target_cluster}.png")
    plt.show()
