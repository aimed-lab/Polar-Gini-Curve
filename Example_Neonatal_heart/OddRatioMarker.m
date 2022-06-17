%% load the dataset

idx = importdata('ClusterID.mat'); %2D visual coordinate
expressMat = importdata('Expression.mat'); % gene expression
gene = importdata('geneList.mat'); % gene symbol
coordinate = importdata('coordinate.mat'); %2D visual coordinate


%% mean expression and percentage expression
expressMat = expressMat';

meanExp = zeros(length(gene), max(idx));
percenExp = zeros(length(gene), max(idx)); % percentage of cell expressing the specific gene in the cluster
pFisher = ones(length(gene), max(idx)); % p-value of fisher's exact test
oddRatio = ones(length(gene), max(idx)); % odd ratio
for i = 1 : length(gene)
    for j = 1 : max(idx)
        clusterIndex = find(idx == j);
        meanExp(i, j) = mean( expressMat(i, clusterIndex) );
        percenExp(i, j) = length( find( expressMat(i, clusterIndex) > 0 ) ) / length (clusterIndex);
        
        count1 = length ( find( expressMat(i, setdiff(1:length(idx), clusterIndex)) == 0 ) );
        count2 = length ( find( expressMat(i, clusterIndex) == 0 ) );
        count3 = length ( find( expressMat(i, setdiff(1:length(idx), clusterIndex)) > 0 ) );
        count4 = length ( find( expressMat(i, clusterIndex) > 0 ) );
        countNum = 1+[ count1, count2; count3, count4];
        
        [h,pFisher(i, j),stats] = fishertest(countNum,'Tail','right','Alpha',0.01);
        oddRatio(i, j) = stats.OddsRatio;
    end
end

save('Stat.mat','meanExp','percenExp', 'pFisher', 'oddRatio'); % save result into file Stat.mat

