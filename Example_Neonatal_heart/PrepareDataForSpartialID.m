%% load the dataset

idx = importdata('ClusterID.mat'); %2D visual coordinate
expressMat = importdata('Expression.mat'); % gene expression
gene = importdata('geneList.mat'); % gene symbol
coordinate = importdata('coordinate.mat'); %2D visual coordinate

mkdir SpartialDE;
mkdir SpartialDE/inputFile;
mkdir SpartialDE/outputFIle

%% format the data into the SpatiaDE readable file, put them into SpartialDE/inputFile folder.
% put the file 
for clusterID = 1 : max(idx)

    index = find(clusterID == idx);
    coorText = cell(length(index), 1);
    for i = 1 : length(index)
        coorText{i} = [ num2str(coordinate(index(i), 1)), 'x', num2str(coordinate(index(i), 2)) ];
    end
    
    subExpress = expressMat(index, :);
    
    header = [{''}, gene'];
    
    inputTable = table(coorText, subExpress);
    
    writetable(inputTable, ['SpartialDE/inputFile/cluster', num2str(clusterID), '.csv'],...
        'Delimiter', ',', 'WriteVariableNames', false);

end