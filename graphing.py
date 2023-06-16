from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import numpy as np


def draw_tsne(
    marker_gene,
    coordinate,
    clusterID,
    target_cluster_id,
    expression_data,
    gene_list,
    random_state=0,
):
    # Create a t-SNE model
    model = TSNE(n_components=2, random_state=random_state)

    # Flatten clusterID
    clusterID = np.ravel(clusterID)

    # Calculate t-SNE coordinates for all cells
    tsne_coordinates = model.fit_transform(coordinate)

    # Set up the plot
    _, ax = plt.subplots()

    # Plot target cluster in grey
    target_cluster = tsne_coordinates[clusterID == target_cluster_id]
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
        clusterID == target_cluster_id, marker_gene_expression
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
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_title(f"t-SNE for marker gene {marker_gene} in cluster {target_cluster_id}")

    # Save figure in png format
    plt.savefig(f"tsne_{marker_gene}_c-{target_cluster_id}.png")
    plt.show()
