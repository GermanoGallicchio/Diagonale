function results = di_analysis_parametricFeature_inferenceCluster_correlation(di_cfg, Y_orig, X_orig, rowIdx)
% perform correlation between each column of X and the one column of Y.
% Correlations can be Pearson (default), Spearman, or Kendall.
% this analysis uses a parametric (Pearson) or similar(Spearman, Kendall) correlation to derive p-values, but no inference is formally done.
% Instead the p values are used to form clusters. cluster metrics are
% extracted for observed and simulated data. Cluster metrics are used for either:
% - permutation null hypothesis testing (permutationH0testing):
%   - observed cluster metrics are compared against a null distribution of simulated-data cluster metrics
% - bootstrap stability (bootstrapStability):
%   - not yet implemented
%
% INPUT:
%   di_cfg        - validated analysis configuration struct
%   Y_orig        - group codes (m x 1) with exactly 2 unique values
%   X_orig        - data matrix (m x pX)
%   rowIdx        - resampling row indices from di_reorderRowsGenerate
%
% OUTPUT:
%   results       - results structure:
%     .observed.statVal       (1 x pX) observed r coefficient
%     .observed.pVal          (1 x pX) uncorrected p-values
%     .observed.clusters      struct() empty (no clustering done)
%
%     for permutation H0 testing
%       .simulated.permutationH0.statVal (nIterations x pX) null distribution
%
%     for bootstrap stability 
%       .simulated.bootstrapStability.values (nIterations x pX) bootstrap resamples of r coefficients
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)


% PHASING OUT
error('This function is being phased out. Already replaced by di_analysis_OLS and soon deleted entirely. currently here just for cross checking')

%% shortcuts

% number of iterations
nIterations = di_cfg.analysis.nIterations;

% total number of features (product of all dimension sizes)
dimKeys  = fieldnames(di_cfg.dimensions);
dimSizes = cellfun(@(k) length(di_cfg.dimensions.(k).vec), dimKeys);
pX       = prod(dimSizes);  % total number of X features

% analysis objective and type
analysisObjective = di_cfg.analysis.objective;
analysisType      = di_cfg.analysis.type;

%% sanity checks
% validate input dimensions

% --- general parameters ---

% num of dimensions must be at least 1
if pX < 1
    error('pX must be at least 1 (total number of features)');
end

% num of dimensions must be same as num of features in X matrix
if size(X_orig, 2) ~= pX
    error('X_orig must have pX columns, same number as features');
end

% num of iterations must be at least 1
if nIterations < 1
    error('nIterations must be at least 1');
end

% num of iterations must be same as num of columns of rowIdx (resampling indices)
if size(rowIdx, 2) ~= nIterations
    error('rowIdx must have nIterations columns');
end

% num of rows in Y must be same as num of rows of rowIdx (resampling indices)
if size(rowIdx, 1) ~= size(Y_orig, 1)
    error('rowIdx must have same number of rows as Y_orig');
end

% X and Y same number of rows
if size(Y_orig, 1) ~= size(X_orig, 1)
    error('Y_orig and X_orig must have same number of rows');
end


% --- specific to this analysis ---

% check that clustering parameters are configured
if ~isfield(di_cfg.analysis, 'clusterParams')
    error('di_cfg.clusterParams must be defined for %s + %s', analysisObjective, analysisType);
end

% Y has only one column for correlation: one column of Y for correlations
if size(Y_orig, 2) ~= 1
    error('Y_orig must be a column vector (m x 1)');
end

% Y is a continuous variable
% TO DO

%% choose correlation type

% choose analysis subtype
list = ["Pearson" "Spearman" "Kendall" "cylindrical"];
[idx,tf] = listdlg('ListString',list,'SelectionMode','single','ListSize',[160 100],'PromptString','choose correlation type');
if tf==0; idx=1; warning('You did not choose correlation type. I chose for you: Pearson'); end
corrType = list(idx);
if strcmp(corrType,'cylindrical')
    error('not coded yet') % it will require its own function to keep things tidy
end


%% main analysis

for itIdx = 1:nIterations

    % apply row reordering based on indices in rowIdx
    [Y,X] = di_reorderRowsApply(di_cfg, Y_orig, X_orig, rowIdx, itIdx);

    % perform correlation
    [statVal, pVal] = corr(Y, X, 'type', corrType);

    % cluster forming
    [clusterMembership, clustIDList, metrics] = di_clusterForming(di_cfg, statVal, pVal);

    % store observed values from first iteration
    if itIdx==1
        statVal_obs           = statVal;
        pVal_obs              = pVal;
        clusterMembership_obs = clusterMembership ;
        clustIDList_obs       = clustIDList ;
        metrics_obs           = metrics;
    end

    % store simulated values according to objective
    switch analysisObjective
        case 'permutationH0testing'
            % initialize on first iteration
            if itIdx==1
                simulatedMetrics = repmat(struct('id', [], 'size', [], 'mass', [], 'mostExtremeVal', []), 1, nIterations);
            end

            % fill structure
            simulatedMetrics(1,itIdx).id             = metrics.id;
            simulatedMetrics(1,itIdx).size           = metrics.size;
            simulatedMetrics(1,itIdx).mass           = metrics.mass;
            simulatedMetrics(1,itIdx).mostExtremeVal = metrics.mostExtremeVal;
            
        case 'bootstrapStability'
            error('not yet coded')
    end

    di_counter(itIdx,nIterations)  % iteration counter

end


%% collate results

% observed, feature level
results.observed.statVal = statVal_obs;     % [1 x pX] r coefficient per feature
results.observed.pVal    = pVal_obs;        % [1 x pX] p-value per feature

% observed, cluster level
results.observed.clusters.clusterMembership_obs = clusterMembership_obs;
results.observed.clusters.clustIDList_obs       = clustIDList_obs;
results.observed.clusters.metrics_obs           = metrics_obs;

% simulated, cluster level
switch analysisObjective
    case 'permutationH0testing'

        results.simulated.permutationH0.clusterMetrics = simulatedMetrics; % [1 x nIterations] array of struct
        % each element contains: .id, .size, .mass, .mostExtremeVal for clusters in that null iteration
        % (these will be used by di_maxT.m to compute corrected p-values per cluster)
        % Used to estimate H0 distribution of cluster metrics
        
    case 'bootstrapStability'
        error('not coded yet')
end

%% validate output

% verify structure consistency
if size(results.observed.statVal, 2) ~= pX
    error('observed.statVal does not have pX columns');
end
if size(results.observed.pVal, 2) ~= pX
    error('observed.pVal does not have pX columns');
end
if size(results.observed.clusters.clusterMembership_obs, 2) ~= pX
    error('clusters.clusterMembership_obs does not have pX columns');
end

% check that clusters were actually formed
if di_cfg.analysis.verbose
    if isempty(results.observed.clusters.clustIDList_obs)
        if di_cfg.clusterParams.clusterFormingThreshold > 0
            warning('Diagonale: no clusters were found in observed data');
        end
    end
end

% verify resampling structure
switch analysisObjective
    case 'permutationH0testing'
        if length(results.simulated.permutationH0.clusterMetrics) ~= nIterations
            error('permutationH0.clusterMetrics does not have nIterations elements');
        end
    case 'bootstrapStability'
        error('not yet coded')
end

end
