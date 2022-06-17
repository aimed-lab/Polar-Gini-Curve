allData = importdata('ClusterID.mat');
idx = allData.clusterID;
expressMat = allData.Expression;
gene = allData.geneList;
coordinate = allData.coordinate;

mkdir trendsceek;
mkdir trendsceek/inputFile;
mkdir trendsceek/outputFile;
expressMat(:, 1) = [];

for clusterID = 1 : max(idx)

    index = find(clusterID == idx);
    
    subCoordinate = coordinate(index, :);
    subExpress = expressMat(index, :);
    
    %header = [{''}, gene'];
    
    inputTable = table(gene, subExpress');
    
    writetable(inputTable, ['trendsceek/inputFile/cluster', num2str(clusterID), '_expression.csv'],...
        'Delimiter', ',', 'WriteVariableNames', false);
    
    writetable(table(subCoordinate), ['trendsceek/inputFile/cluster', num2str(clusterID), '_coordinate.csv'],...
        'Delimiter', ',', 'WriteVariableNames', false);
end