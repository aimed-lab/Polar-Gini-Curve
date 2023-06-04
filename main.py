"""
Main driver code for the project.
"""

from load_dataset import load_dataset
from process_data import identify_clusters, plot_clusters_and_expression
from draw_pgc import compute_rmsd

coordinate = load_dataset("data/coordinate.mat")  # 2D spatial coordinate
Expression = load_dataset("data/Expression.mat")  # gene expression
geneList = load_dataset("data/geneList.mat")  # gene symbol
clusterID = load_dataset("data/clusterID.mat")  # cluster ID

# numCluster = np.max(clusterID)  # get the number of clusters
# print(numCluster)

# print(identify_clusters(Expression))
# plot_clusters_and_expression(coordinate, clusterID, Expression, geneList)

# Test case with Actc1 marker
print(compute_rmsd(coordinate, clusterID, Expression, geneList))
