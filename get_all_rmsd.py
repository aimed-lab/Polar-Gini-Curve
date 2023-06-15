"""
Get all RMSD values for each gene in each cluster.
"""

import csv
import numpy as np
from draw_pgc import compute_rmsd


def get_all_rmsd(coordinate, cluster_id, expression, gene_list):
    """
    Compute and save the RMSD values for each gene in each cluster.

    Parameters:
    - coordinate (numpy.ndarray): The 2D spatial coordinates.
    - cluster_id (numpy.ndarray): The ID for each cluster.
    - expression (numpy.ndarray): The gene expression data.
    - gene_list (numpy.ndarray): The list of gene symbols.
    - output_file (str): The file path to save the results.

    Returns:
    None. The function saves the RMSD values to the output_file.
    """

    # Load the dataset
    num_cluster = int(np.max(cluster_id))

    with open("data/all_rmsd.csv", "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        header = ["Gene"]

        for i in range(1, num_cluster + 1):
            header.append("RMSD for cluster " + str(i))

        writer.writerow(header)

        for i, gene in enumerate(gene_list):
            matrix_result = np.empty(num_cluster)
            matrix_result[:] = np.nan

            for j, cluster in enumerate(range(1, num_cluster + 1)):
                rmsd = compute_rmsd(
                    coordinate,
                    cluster_id,
                    expression,
                    gene_list,
                    marker=gene,
                    cluster=cluster,
                )
                matrix_result[j] = rmsd

            writer.writerow([gene[0][0], *matrix_result])
