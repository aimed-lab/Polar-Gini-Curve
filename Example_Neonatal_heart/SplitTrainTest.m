allCoordinate = importdata('coordinate.mat');
allExpression = importdata('Expression.csv');
allClusterID = importdata('ClusterID.mat');

numCluster = max(allClusterID);

randIndex = randperm(length(allClusterID));

numTrain = round( length(randIndex) * 0.8 );
numTest = length(randIndex) - numTrain;

%% form the training set
coordinate = allCoordinate( randIndex(1:numTrain), : );
Expression = allExpression( randIndex(1:numTrain), : );
clusterID = allClusterID( randIndex(1:numTrain) );

save('train_test_exp/AllTrainData.mat', 'coordinate', 'Expression', 'clusterID');
% save the variables above and geneIndex into train_test_exp/AllTrainData.mat

%% form the test set
coordinate = allCoordinate( randIndex(numTrain + (1:numTest)), : );
Expression = allExpression( randIndex(numTrain + (1:numTest)), : );
clusterID = allClusterID( randIndex(numTrain + (1:numTest)) );

save('train_test_exp/AllTestData.mat', 'coordinate', 'Expression', 'clusterID');
% save the variables above and geneIndex into train_test_exp/AllTestData.mat

