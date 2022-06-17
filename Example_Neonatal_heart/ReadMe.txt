Download mouse fetal lung dataset at: http://bis.zju.edu.cn/MCA/

I format the data into Matlab as follow
coordinate.mat: tsne coordinate
Expression.mat: gene expression matrix. Each gene corresponds to a column, each cell corresponds to a row
geneList.mat: gene list
clusterID.mat: cluster ID, from 1 to 9. -1 corresponds to non-clusterable cells.

- Run getAllRSMD.m to get the raw RSMD table for all genes in all clusters. Store in file matrixResult.txt
- Run OddRatioMarker.m to get differentially expressed gene. Save the result in file Stat.mat
- Run ComputePVal.m to get the p-value for the RSMD. Write the result into file RSMD_result.xlsx
- Run PrepareDataForSpartialID.m to prepare data for SpatialDE (DOI:10.1038/nmeth.4636)
Then run python code spatial_MOB_analysis.py in folder SpartialDE
- Run file PrepareDataForTrendsceek.m to prepare data for Trendsceek (DOI:10.1038/nmeth.4634)
Then run R code RunTrendsceek
- Run file SplitTrainTest.m to setup the data for cluster ID re-identification.
Then follow the ReadMe.txt file in folder train_test_exp to get the results