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

    # mean expression and percentage expression
    mean_exp = np.zeros((len(gene_list), np.max(cluster_id)))
    percent_exp = np.zeros(
        (len(gene_list), np.max(cluster_id))
    )  # percentage of cell expressing the specific gene_list in the cluster
    p_fisher = np.ones(
        (len(gene_list), np.max(cluster_id))
    )  # p-value of fisher's exact test
    odd_ratio = np.ones((len(gene_list), np.max(cluster_id)))  # odd ratio

    for i in range(len(gene_list)):
        for j in range(np.max(cluster_id)):
            cluster_index = np.nonzero(cluster_id == j + 1)[0]
            mean_exp[i, j] = np.mean(expression[i, cluster_index])
            percent_exp[i, j] = len(
                np.nonzero(expression[i, cluster_index] > 0)[0]
            ) / len(cluster_index)

            count1 = len(
                np.nonzero(
                    expression[
                        i, np.setdiff1d(np.arange(len(cluster_id)), cluster_index)
                    ]
                    == 0
                )[0]
            )
            count2 = len(np.nonzero(expression[i, cluster_index] == 0)[0])
            count3 = len(
                np.nonzero(
                    expression[
                        i, np.setdiff1d(np.arange(len(cluster_id)), cluster_index)
                    ]
                    > 0
                )[0]
            )
            count4 = len(np.nonzero(expression[i, cluster_index] > 0)[0])
            count_num = np.array([[count1, count2], [count3, count4]]) + 1

            _, p_fisher[i, j] = fisher_exact(count_num, "greater")
            odd_ratio[i, j] = (count_num[0, 0] * count_num[1, 1]) / (
                count_num[0, 1] * count_num[1, 0]
            )

    # save result into file Stat.csv into data directory
    with open(
        "data/odd_ratio_marker.csv", "w", newline="", encoding="utf-8"
    ) as csvfile:
        writer = csv.writer(csvfile)
        header = [
            "Gene",
            "Mean expression",
            "Percentage expression",
            "p-value",
            "Odd ratio",
        ]
        writer.writerow(header)
        for i, gene in enumerate(gene_list):
            writer.writerow(
                [gene[0][0], *mean_exp[i], *percent_exp[i], *p_fisher[i], *odd_ratio[i]]
            )
