load ClusterID.mat;
load Stat.mat;
load('matrixResult.mat')

normalizedResult = ones(size(matrixResult));
pVal = ones(length(geneList), max(idx));

for j = 1 : max(idx)
    allClusterRSMD = matrixResult(:, j);
    index = find(percenExp(:, j) > 0.1);
    p_valueAble = allClusterRSMD(index);
    maxScore = max( p_valueAble );
    
    normalizedResult(index, j) = matrixResult( index, j ) / maxScore;
    p_valueAble = p_valueAble / maxScore;
    
    mu = mean(p_valueAble);
    sigma = std(p_valueAble);
    p = normcdf(p_valueAble,mu,sigma);
    pVal(index, j) = p;
end

save normalModelPVal.mat pVal -mat;

% write results into excel files
for i = 1 : max(idx)
    T = table (geneList, percenExp(:, i), matrixResult(:, i), normalizedResult(:, i), pVal(:, i));
    T.Properties.VariableNames = {'Gene', 'percentage_cell_expressing', 'raw_RSMD', 'normalized_RSMD', 'P_value'};
    writetable(T, ['RSMD_cluster', num2str(i), '.xlsx']);
end


