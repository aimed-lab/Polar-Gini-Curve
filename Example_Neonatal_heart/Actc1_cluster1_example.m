clear all; close all; clc; % clear the Matlab workspace if needed.

% This example shows how to compute and plot the polar-gini curve for one
% gene (here I choose Actc1) in one cluster (here I choose cluster 1), in a
% real single-cell dataset
% For methodolgy, please refer to Nguyen et al. Polar Gini Curve: a Technique to Discover Single-cell Biomarker Using 2D Visual Information

% Installation:
% Matlab version 2015 or later
% Put two function-files into the Matlab workspace:
%    - computeGini.m. This function computes gini coefficient. Credit to Yvan Lengwiler (2020). Gini coefficient and the Lorentz curve (https://www.mathworks.com/matlabcentral/fileexchange/28080-gini-coefficient-and-the-lorentz-curve), MATLAB Central File Exchange. Retrieved March 4, 2020
%    - make2DGini.m This function computes the dsta structure to plot the
%    2D polar gini curve, and plot the curves.

% Prepare the dataset
% - I put the dataset into Matlab format (.mat) file. If the dataset is in
% .txt format, or in 10x, please refer to:
% - importdata (https://www.mathworks.com/help/matlab/ref/importdata.html).
% This function load the data from the matrix / vector format and
% automatically convert the data into the Matlab's numeric or text
% matrix/vector.
% - readtable(https://www.mathworks.com/help/matlab/ref/readtable.html)
%   readcsv(https://www.mathworks.com/help/matlab/ref/csvread.html)
%   readmatrix(https://www.mathworks.com/help/matlab/ref/readmatrix.html)
%   xlsread(https://www.mathworks.com/help/matlab/ref/xlsread.html)
%   table2cell(https://www.mathworks.com/help/matlab/ref/table2cell.html)
%   cell2mat(https://www.mathworks.com/help/matlab/ref/cell2mat.html)
% These function aboves can read the data from many formats and
% convert/preprocess the data into Matlab workspace matrix/vector.
% - spconvert (https://www.mathworks.com/help/matlab/ref/spconvert.html),
% for importing sparse matrix format (i.e. 10XGenomics expression file)

%% load the dataset of n = 5075 cells, m = 19493 genes
coordinate = importdata('coordinate.mat'); %2D x-y visual coordinate, in nx2 matrix format. 
Expression = importdata('Expression.mat'); % gene expression, in nxm matrix format
geneList = importdata('geneList.mat'); % gene symbol, in mx1 text vector
clusterID = importdata('ClusterID.mat'); % cluster/group ID, in nx1 vector

% visualize the dataset
figure, gscatter (coordinate(:, 1), coordinate(:, 2), clusterID); title('dataset visualization');

%% if the dataset only has the expression and gene list, a good pipeline to compute the coordinate and cluster ID is
% coordinate = tsne(Expression);
% clusterID = dbscan(coordinate, 0.5, 50);
% see https://www.mathworks.com/help/stats/dbscan-clustering.html to setup
% the proper parameters for dbscan clustering algorithm.

%% extract the dataset for Actc1 and cluster 1
geneName = 'Actc1';
clusterNum = 1;
[~, geneIndex] = ismember(geneName, geneList); % find where Actc1 gene is in the data
clusterIndex = find(clusterID == clusterNum); % find which cells are in cluster 1

clusterCoor = coordinate(clusterIndex, :); %2D coordinate of all cells in cluster 1
geneClusterExp = Expression(clusterIndex, geneIndex); % expression of Actc1 in all cluster 1 cell

% visualize cluster 1 (optional)
figure, gscatter(clusterCoor(:, 1), clusterCoor(:, 2)); title(['cluster ', num2str(clusterNum), 'visualization']);
% visualize Actc1 expression in cluster 1 (optional)
threshold = 0;
index1 = find(geneClusterExp <= threshold);
figure, scatter(clusterCoor(index1, 1), clusterCoor(index1, 2), 15, [192, 192, 192]/255, '.');
hold on
index1 = find(geneClusterExp > threshold );
scatter(clusterCoor(index1, 1), clusterCoor(index1, 2), 15, geneClusterExp(index1), '.');
colorbar;
colormap('jet');
title([geneName, ' expression in cluster ', num2str(clusterNum) ]);
hold off

%% compute and plot the polar gini curve

% find which cell in cluster 1 has Actc1 expression (>0)
threshold = 0;
geneExpressIndex = find(geneClusterExp > threshold );

% setup two groups. Group 1 includes all cells in cluster 1 with Actc
% expression > 0. Group 2 includes all ceels incluster 1.
plotGroup = [ones(size(geneExpressIndex)) ; 2*ones(size(clusterIndex))];

% setup the 2D coordinate matrix for the two groups above
plotCoor = [clusterCoor(geneExpressIndex, :); clusterCoor];

% label che group if needed
plotLbl = { [geneName, ' expresing cell'], ['cluster ', num2str(clusterNum)] };

% the polar gini curve
[angleList, allGini] = make2DGini(plotCoor, plotGroup, plotLbl);
title([geneName, ' polar gini curves']);

% to compute the RSMD as showed in the paper
diff = abs( allGini{1}- allGini{2} );
RSMD = sqrt(mean(diff.^2));
