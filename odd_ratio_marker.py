"""
Calculates mean expression, percentage expression, p-value of 
fisher's exact test, and odd ratios of a gene in a cell cluster.
"""

import csv
import numpy as np
from scipy.stats import fisher_exact


def odd_ratio_marker(cluster_id, expression, gene_list):
    """
    Function to compute the mean expression, percentage expression, p-value of Fisher's exact test,
    and odd ratios for a set of genes in a cell cluster.

    Args:
        cluster_id (numpy.array): An array containing the cluster IDs.
        expression (numpy.array): A 2D array representing gene expression values.
                                  Columns represent cells and rows represent genes.
        gene_list (numpy.array): An array of gene symbols.

    Returns:
        None. The function writes the results into a
        .csv file named 'Stat.csv' in the 'data' directory.

    """

    expression = expression.T

    # Open file to write the results
    with open("data/odd_ratio_marker.csv", "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        # Write header
        header = ["Gene", "Cluster", "Mean expression", "Percentage expression", "p-value", "Odd ratio"]
        writer.writerow(header)

        # Loop through each gene
        for i in range(len(gene_list)):
            for j in range(1, np.max(cluster_id) + 1):
                cluster_index = np.nonzero(cluster_id == j)[0]
                mean_exp = np.mean(expression[i, cluster_index])
                percent_exp = len(np.nonzero(expression[i, cluster_index] > 0)[0]) / len(cluster_index)

                count1 = len(np.nonzero(expression[i, np.setdiff1d(np.arange(len(cluster_id)), cluster_index)] == 0)[0])
                count2 = len(np.nonzero(expression[i, cluster_index] == 0)[0])
                count3 = len(np.nonzero(expression[i, np.setdiff1d(np.arange(len(cluster_id)), cluster_index)] > 0)[0])
                count4 = len(np.nonzero(expression[i, cluster_index] > 0)[0])
                count_num = np.array([[count1, count2], [count3, count4]]) + 1

                p_value, _ = fisher_exact(count_num, "greater")
                odd_ratio = (count_num[0, 0] * count_num[1, 1]) / (count_num[0, 1] * count_num[1, 0])

                # Write results for current gene-cluster combination
                writer.writerow([gene_list[i][0][0], j, mean_exp, percent_exp, p_value, odd_ratio])
