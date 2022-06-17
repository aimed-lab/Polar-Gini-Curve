function [angleList, allGini] = make2DGini(x, clusterID, clusterName)
% make the 2D Gini curve plot on the 2d data
%
% input:
%    x: the 2D datapoints matrix (i.e., x-y coordinates after tSNE, the first
%       two PCs in PCA. each row corresponds to one datapoint
%    culsterID: cluster ID for each datapoint, in integer format: 1, 2, 3, ...
%       to the number of clusters. Prefered number of clusters: less than 7
%    clusterName: cluster name for each data cluster, correspond to clusterID,
%       in string format for figure legend
%
% output:
%    angleList: the list of angles, in radian, for reprodue in other polar
%       plot
%    allGini: all Gini index (radius) for each cluster in each angle
%
% example:
% x = [ mvnrnd([3 3], [1 0; 0 8],3000) ;  mvnrnd([-3 -3], [2 0.5; 0.5 1],3000)]; % randomize data with two clusters
% clusterID = [ones(3000, 1); 2*ones(3000, 1)]; % cluster ID
% [angleList, allGini] = make2DGini(x, clusterID); % make the gini 2D plot without cluster name
% clusterName = {'group 1', 'group 2'};
% [angleList, allGini] = make2DGini(x, clusterID, clusterName); % make the gini 2D plot with cluster name
% to reproduce one polar plot, for example in the first cluster:  polarplot( angleList, allGini{1});

assert(nargin >= 2, 'We need at least the 2D data and cluster ID.')

numCluster = length(unique(clusterID));

% automatically set the cluster name if needed
if nargin < 3
    clusterName = cell(numCluster, 1);
    for i = 1 : numCluster
        clusterName{i} = ['cluster ', num2str(i)];
    end
end

% get the angle list with resolution of 1000
resolution = 1000;
angleList = pi*(0:360/resolution:360) / 180;

ColorOdrDef = get(gca,'ColorOrder');
allGini = cell(1 : numCluster);

figure

% for each cluster, get the list of gini corresponding to each angle
for cluster = 1 : numCluster
    coordinate = x(find(clusterID == cluster), :);
    plotCoordinate = zeros(length(angleList), 2);
    gini = zeros(length(angleList), 1);
    for i = 1 : length(angleList)
        angle = angleList(i);
        value = coordinate * [cos(angle); sin(angle)];
        value = value - min(value);
        gini(i) = computeGini( ones(length(value), 1), value );
        plotCoordinate(i, :) = [ gini(i) * cos(angle), gini(i) * sin(angle) ];
        
    end
    allGini{cluster} = gini;
    
    polarplot( angleList, gini, 'color', ColorOdrDef(cluster, :) );
    hold on
    
end
legend(clusterName);
hold off

end