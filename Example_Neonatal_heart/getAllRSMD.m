%% load the dataset
coordinate = importdata('coordinate.mat'); %2D visual coordinate
Expression = importdata('Expression.mat'); % gene expression
geneList = importdata('geneList.mat'); % gene symbol
clusterID = importdata('ClusterID.mat'); % cluster ID

numCluster = max(clusterID);

matrixResult = NaN(length(geneList), numCluster);

parfor i = 1 : length(geneList)
    for j = 1 : numCluster
        clusterIndex = find(clusterID == j);
        clusterCoor = coordinate(clusterIndex, :);
        clusterLbl = ones(length(clusterIndex), 1);
        
        geneIndex = find(clusterID == j & Expression(:, i) > 0);
        if (length(geneIndex) > 0)
            geneCoor = coordinate(geneIndex, :);
            geneLbl = 2*ones(length(geneIndex), 1);
            
            [angleList, allGini] = make2DGini([geneCoor; clusterCoor],...
                [geneLbl; clusterLbl], {geneList{i}, ['cluster ', num2str(j)]}, [1, 1]);
            
            y = abs( allGini{1} - allGini{2} ).^2;
            matrixResult(i, j) = sqrt(mean(y));
        end
    end
end

% save the result in file matrixResult.mat
save matrixResult.mat matrixResult -mat;