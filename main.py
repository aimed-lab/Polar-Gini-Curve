"""
Main driver code for the project.
"""

# import pandas as pd
# import numpy as np
from load_dataset import load_dataset

# # from process_data import identify_clusters, plot_clusters_and_expression
# from draw_pgc import make_2d_gini
# # from get_all_rmsd import get_all_rmsd
# # from odd_ratio_marker import odd_ratio_marker
# from compute_pval import compute_pval
# from graphing import draw_tsne, plot_gini

coordinate = load_dataset("data/coordinate.mat")  # 2D spatial coordinate
Expression = load_dataset("data/Expression.mat")  # gene expression
geneList = load_dataset("data/geneList.mat")  # gene symbol
clusterID = load_dataset("data/clusterID.mat")  # cluster ID

# draw_tsne(
#     "Actc1", coordinate, clusterID, 2, Expression, geneList, random_state=0
# )

# plot_gini(coordinate, clusterID, Expression, geneList)

# numCluster = np.max(clusterID)  # get the number of clusters
# print(numCluster)

# identify_clusters(Expression)
# plot_clusters_and_expression(coordinate, clusterID, Expression, geneList)

# Test case with Actc1 marker
# RSMD = compute_rmsd(coordinate, clusterID, Expression, geneList)

# print(RSMD)

# odd_ratio_marker(clusterID, Expression, geneList)
# get_all_rmsd(coordinate, clusterID, Expression, geneList)
# odd_ratio_marker = pd.read_csv("data/odd_ratio_marker.csv")
# all_rmsd = pd.read_csv("data/all_rmsd.csv")

# compute_pval(clusterID, geneList, odd_ratio_marker, all_rmsd)
