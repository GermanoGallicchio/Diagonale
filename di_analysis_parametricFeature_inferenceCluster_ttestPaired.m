function results = di_analysis_parametricFeature_inferenceCluster_ttestPaired(di_cfg, Y_orig, X_orig, rowIdx)
% perform paired sample t-tests on each column of X based on the
% splitting information in the column of Y. this analysis uses a parametric
% t test to derive p-values, but no inference is formally done.
% Instead the p values are used to form clusters. cluster metrics are
% extracted for observed and simulated data. Cluster metrics are used for either:
% - permutation null hypothesis testing (permutationH0testing):
%   - observed cluster metrics are compared against a null distribution of simulated-data cluster metrics
% - bootstrap stability (bootstrapStability):
%   - not yet implemented
%
% INPUT:
%   di_cfg        - validated analysis configuration struct
%   Y_orig        - condition codes (m x 1) with exactly 2 unique values
%   X_orig        - data matrix (m x pX) where pX = product of all dimensions
%   rowIdx        - resampling row indices from di_reorderRowsGenerate
%
% OUTPUT:
%   results       - results structure:
%     .observed.statVal       (1 x pX) observed t-statistics per feature
%     .observed.pVal          (1 x pX) uncorrected p-values per feature
%     .observed.clusters      struct with:
%       .clusterMembership_obs  (1 x pX) cluster ID for each feature (0=not in cluster)
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
%       .simulated.bootstrapStability.values  (nIterations x pX)
%         Bootstrap resamples of feature-level t-statistics
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

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
    error('pX must be at least 1 (total number of spatial features)');
end

% num of dimensions must be same as num of features in X matrix
if size(X_orig, 2) ~= pX
    error('X_orig must have pX columns, matching input parameter');
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

% Y has only one column for t-tests: one comparison only
if size(Y_orig, 2) ~= 1
    error('Y_orig must be a column vector (m x 1)');
end

% Y has exactly 2 unique values for t-tests: comparing two things
unique_groups_orig = unique(Y_orig);
if numel(unique_groups_orig) ~= 2
    error('Y_orig must have exactly 2 unique condition codes for paired t-test');
end

% check that dataStruct is a table and contains observationID values
if ~isfield(di_cfg.analysis, 'dataStruct') || ~istable(di_cfg.analysis.dataStruct)
    error('di_cfg.analysis.dataStruct must be a table for paired t-test');
end
if ~ismember('observationID', di_cfg.analysis.dataStruct.Properties.VariableNames)
    error('dataStruct must have observationID column for paired t-test');
end

%% shortcuts 2

obsIDVec = di_cfg.analysis.dataStruct.observationID;

%% main analysis

for itIdx = 1:nIterations

    % apply row reordering based on indices in rowIdx
    [Y,X] = di_reorderRowsApply(di_cfg, Y_orig, X_orig, rowIdx, itIdx);

    % --- find paired observations for paired t-test --- start
    % for each observation ID, find the two condition measurements and pair them.
    % this is critical if the order of data is not tidy
    
    % reorder observation IDs according to the current iteration's row permutation.
    % this preserves the pairing structure while shuffling row order for permutation testing.
    obs_perm = obsIDVec(rowIdx(:,itIdx)); % reordered observation ID vector
    
    % identify all unique observations in this (possibly permuted) dataset
    uniqueObs = unique(obs_perm);
    nObsPairs = numel(uniqueObs);
    
    % identify the two condition codes and determine which is max vs min.
    % (use min/max to ensure consistent ordering across iterations)
    condVals = unique(Y);
    condMax = max(condVals);
    condMin = min(condVals);
    
    % for each observation, extract all its measurements under both conditions
    % and pair them sequentially. This handles bootstrap resampling where the
    % same observation may appear multiple times with balanced conditions.

    % pre-allocate with case of needing most space (all pairs valid)
    % maximum possible pairs is half the total number of observations
    maxPairs = ceil(length(obs_perm) / 2);
    X_max_all = nan(maxPairs, pX);
    X_min_all = nan(maxPairs, pX);

    pairCount = 0;
    for obsIdx = 1:nObsPairs
        % logical mask for all rows belonging to this observation
        maskObs = obs_perm==uniqueObs(obsIdx);
        
        % find all row indices for this observation at each condition
        idxMax = find(maskObs & (Y==condMax));
        idxMin = find(maskObs & (Y==condMin));
        
        % each observation must appear the same
        % number of times in both conditions (pairs stay together)
        if ~isempty(idxMax) && ~isempty(idxMin)
            if numel(idxMax) ~= numel(idxMin)
                error('Observation %s has condition %d %d times but condition %d %d times . Pairs must be balanced.', ...
                    string(uniqueObs(obsIdx)), condMax, numel(idxMax), condMin, numel(idxMin));
            end
            nPairsForObs = numel(idxMax);
            for p = 1:nPairsForObs
                pairCount = pairCount + 1;
                X_max_all(pairCount,:) = X(idxMax(p),:);
                X_min_all(pairCount,:) = X(idxMin(p),:);
            end
        end
    end
    
    % verify that at least some observation pairs were formed
    if pairCount == 0
        error('no complete observation pairs found for paired t-test');
    end
    
    % trim to actual pairs used
    X_max = X_max_all(1:pairCount,:);
    X_min = X_min_all(1:pairCount,:);

    % --- find paired observations for paired t-test --- end
    
    % perform paired test on aligned per-observation rows
    % note: category with larger Y compared vs category with lower Y
    [~,p,~,stats] = ttest(X_max, X_min);
    statVal = stats.tstat;
    pVal = p;

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
results.observed.statVal = statVal_obs;     % [1 x pX] t-statistic per spatial feature
results.observed.pVal    = pVal_obs;        % [1 x pX] p-value per spatial feature

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
        if di_cfg.analysis.clusterParams.clusterFormingThreshold > 0
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
