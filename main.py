import numpy as np

from load_dataset import load_dataset

coordinate = load_dataset("data/coordinate.mat")  # 2D spatial coordinate
Expression = load_dataset("data/Expression.mat")  # gene expression
geneList = load_dataset("data/geneList.mat")  # gene symbol
clusterID = load_dataset("data/clusterID.mat")  # cluster ID


numCluster = np.max(clusterID)  # get the number of clusters

print(numCluster)
