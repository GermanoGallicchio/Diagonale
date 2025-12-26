function results = di_analysis_parametricFeature_inferenceCluster_ttestInd(di_cfg, Y_orig, X_orig, rowIdx)
% perform independent sample t-tests on each column of X based on the
% splitting information in the column of Y. this analysis uses a parametric
% t test to derive p-values, but no inference is formally done.
% Instead the p values are used to form clusters. cluster metrics are
% extracted for observed and simulated data. Cluster metrics are used for either:
% - permutation null hypothesis testing (permutationH0testing):
%   - observed cluster metrics are compared against a null distribution of 
%   simulated-data cluster metrics
% - bootstrap stability (bootstrapStability):
%   - not yet implemented
%
% INPUT:
%   di_cfg        - analysis configuration structure
%   Y_orig        - group codes (m x 1) with exactly 2 unique values
%   X_orig        - data matrix (m x Nall) where Nall = product of all dimensions
%   rowIdx        - resampling row indices from di_reorderRowsGenerate
%
% OUTPUT:
%   results       - results structure:
%     .observed.statVal       (1 x Nall) observed t-statistics per feature
%     .observed.pVal          (1 x Nall) uncorrected p-values per feature
%     .observed.clusters      struct with:
%       .clusterMembership_obs  (1 x Nall) cluster ID for each feature (0=not in cluster)
%       .clustIDList_obs        vector of unique cluster IDs
%       .metrics_obs            struct with size and mass per cluster
%
%       .simulated      field referring to output from simulation-based
%                       analysis, depending on the objective: 
%     for permutation H0 testing
%       .simulated.permutationH0.clusterMetrics   array of struct, one per iteration
%         Each contains cluster metrics from that null iteration
%     
%     for bootstrap stability 
%       .simulated.bootstrapStability.values  (nIterations x Nall)
%         Bootstrap resamples of feature-level t-statistics
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% shortcuts

% number of iterations
nIterations = di_cfg.analysis.nIterations;

% Total number of features (product of all dimension sizes)
dimKeys  = fieldnames(di_cfg.dimensions);
dimSizes = cellfun(@(k) length(di_cfg.dimensions.(k).vec), dimKeys);
Nall     = prod(dimSizes);

%% sanity checks
% validate input dimensions

% --- general parameters ---

% num of dimensions
if Nall < 1
    error('Nall must be at least 1 (total number of features)');
end

% features in X matrix
if size(X_orig, 2) ~= Nall
    error('X_orig must have Nall columns, same number as features');
end

% num of iterations
if nIterations < 1
    error('nIterations must be at least 1');
end

% rowIdx has as many columns as iterations
if size(rowIdx, 2) ~= nIterations
    error('rowIdx must have nIterations columns');
end

% rows in Y and of rowIdx (resampling indices)
if size(rowIdx, 1) ~= size(Y_orig, 1)
    error('rowIdx must have same number of rows as Y_orig');
end

% X and Y same number of rows
if size(Y_orig, 1) ~= size(X_orig, 1)
    error('Y_orig and X_orig must have same number of rows');
end

% --- specific to this analysis ---

% Check that clustering parameters are configured
if ~isfield(di_cfg.analysis, 'clusterParams')
    error('di_cfg.clusterParams must be defined for theoretical L1 analysis');
end

% Y has only one column
if size(Y_orig, 2) ~= 1
    error('Y_orig must be a column vector (m x 1)');
end

% check that Y has exactly 2 unique values
unique_groups_orig = unique(Y_orig);
if numel(unique_groups_orig) ~= 2
    error('Y_orig must have exactly 2 unique group codes for independent t-test');
end

%% main analysis

for itIdx = 1:nIterations

    % apply row reordering based on indices in rowIdx
    [Y,X] = di_reorderRowsApply(di_cfg,Y_orig,X_orig,rowIdx,itIdx);

    % validate exactly two groups
    % not needed because Y has already been validated above
    % groupVals = unique(Y);
    % if numel(groupVals) ~= 2
    %     error('independent t-test requires exactly two group codes in Y');
    % end

    % perform independent sample t-test
    % note: category with larger Y compared vs category with lower Y
    varType = 'unequal';  % equal | unequal (for info see doc ttest2)
    [~,p,~,stats] = ttest2(X(Y==max(Y),:),X(Y==min(Y),:),'Vartype',varType);
    statVal = stats.tstat;
    pVal = p;

    % cluster forming
    [clusterMembership, clustIDList, metrics] = di_clusterForming(di_cfg, statVal, pVal);

    if itIdx==1
        statVal_obs           = statVal;
        pVal_obs              = pVal;
        clusterMembership_obs = clusterMembership ;
        clustIDList_obs       = clustIDList ;
        metrics_obs           = metrics;
    end

    switch di_cfg.analysis.objective
        case 'permutationH0testing'
            % initialize cluster-level metrics
            if itIdx==1
                simulatedMetrics = repmat(struct('id', [], 'size', [], 'mass', [], 'extremeVal', []), 1, nIterations);
            end

            simulatedMetrics(1,itIdx).id         = metrics.id;
            simulatedMetrics(1,itIdx).size       = metrics.size;
            simulatedMetrics(1,itIdx).mass       = metrics.mass;
            simulatedMetrics(1,itIdx).extremeVal = metrics.extremeVal;
            
        case 'bootstrapStability'
            error('not yet coded')
    end

    di_counter(itIdx,nIterations)  % iteration counter

end

%% collate results

% observed, feature level
results.observed.statVal = statVal_obs;     % [1 x Nall] t-statistic per feature
results.observed.pVal    = pVal_obs;        % [1 x Nall] p-value per feature

% observed, cluster level
results.observed.clusters.clusterMembership_obs = clusterMembership_obs;
results.observed.clusters.clustIDList_obs       = clustIDList_obs;
results.observed.clusters.metrics_obs           = metrics_obs;

% simulated, cluster level
switch di_cfg.analysis.objective
    case 'permutationH0testing'

        results.simulated.permutationH0.clusterMetrics = simulatedMetrics; % [1 x nIterations] array of struct
        % each element contains: .id, .size, .mass, .extremeVal for clusters in that null iteration
        % (these will be used by di_maxT.m to compute corrected p-values per cluster)
        % Used to estimate H0 distribution of cluster metrics
        
    case 'bootstrapStability'
        error('not coded yet')
end

%% validate output

% verify structure consistency
if size(results.observed.statVal, 2) ~= Nall
    error('observed.statVal does not have Nall columns');
end
if size(results.observed.pVal, 2) ~= Nall
    error('observed.pVal does not have Nall columns');
end
if size(results.observed.clusters.clusterMembership_obs, 2) ~= Nall
    error('clusters.clusterMembership_obs does not have Nall columns');
end

% check that clusters were actually found
if di_cfg.analysis.verbose
    if isempty(results.observed.clusters.clustIDList_obs)
        warning('Diagonale: no clusters were found in observed data');
    end
end

% verify resampling structure
switch di_cfg.analysis.objective
    case 'permutationH0testing'
        if length(results.simulated.permutationH0.clusterMetrics) ~= nIterations
            error('permutationH0.clusterMetrics does not have nIterations elements');
        end
    case 'bootstrapStability'
        error('not yet coded')
end

end
