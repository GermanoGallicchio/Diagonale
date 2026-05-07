function results = di_analysis_empiricalFeature_inferenceFeature_ttestPaired(di_cfg, Y_orig, X_orig, rowIdx)
% perform paired sample t-tests on each column of X based on the
% splitting information in the column of Y. Inference (permutation testing
% or bootstrap stability) is done at the feature level empirically, that is
% through simulated data:
% - permutation null hypothesis testing (permutationH0testing):
%   - empirical p-values are obtained for each feature by comparing against a null distribution
%   - FDR correction can be applied (no other correction is easily applied in this case)
% - bootstrap stability (bootstrapStability):
%   - bootstrap ratios and confidence intervals are obtained for each by using bootstrap samples
% Optionally, clusters are formed for descriptive purposes downstream
%
% INPUT:
%   di_cfg        - validated analysis configuration struct
%   Y_orig        - group codes (m x 1) with exactly 2 unique values
%   X_orig        - data matrix (m x pX)
%   rowIdx        - resampling row indices from di_reorderRowsGenerate
%
% OUTPUT:
%   results       - results structure:
%     .observed.statVal       (1 x pX) observed t-statistics
%     .observed.pVal          (1 x pX) uncorrected p-values
%     .observed.clusters      struct() empty (no clustering done)
%
%     for permutation H0 testing
%       .simulated.permutationH0.statVal (nIterations x pX) null distribution
%
%     for bootstrap stability 
%       .simulated.bootstrapStability.values (nIterations x pX) bootstrap resamples of t-statistics
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

% % check that clustering parameters are configured
% if ~isfield(di_cfg.analysis, 'clusterParams')
%     error('di_cfg.clusterParams must be defined for %s + %s', analysisObjective, analysisType);
% end
% not necessary for this analysis path // remove the lines above

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

    % store observed values from first iteration
    if itIdx == 1
        statVal_obs = statVal;
        pVal_obs    = pVal;
    end

    % store simulated values according to objective
    switch analysisObjective
        case 'permutationH0testing'
            % initialize on first iteration
            if itIdx == 1
                statVal_resamp = zeros(nIterations, pX);
            end

            % fill array
            statVal_resamp(itIdx,:) = statVal;

        case 'bootstrapStability'
            % initialize on first iteration
            if itIdx == 1
                statVal_boot = zeros(nIterations, pX);
            end

            % fill array
            statVal_boot(itIdx,:) = statVal;
    end

    di_counter(itIdx,nIterations)  % iteration counter
end

%% collate results

% observed, feature level
results.observed.statVal   = statVal_obs;     % [1 x pX] t-statistic per feature
results.observed.pVal      = pVal_obs;        % [1 x pX] p-value per feature

% observed, cluster level
results.observed.clusters  = struct();        % empty: no clustering for this analysis

% simulated, feature level
switch analysisObjective
    case 'permutationH0testing'

        results.simulated.permutationH0.statVal = statVal_resamp; % [nIterations x pX]
        % matrix storing the null distribution under H0 (exchangeability of
        % group labels) for each column of X

    case 'bootstrapStability'
        results.simulated.bootstrapStability.values = statVal_boot; % [nIterations x pX]
        % matrix storing the bootstrap resamples for each column of X
end

%% validate output

% verify structure consistency
if size(results.observed.statVal, 2) ~= pX
    error('observed.statVal does not have pX columns');
end
if size(results.observed.pVal, 2) ~= pX
    error('observed.pVal does not have pX columns');
end

% verify resampling structure
switch analysisObjective
    case 'permutationH0testing'
        if size(results.simulated.permutationH0.statVal, 1) ~= nIterations
            error('permutationH0.statVal does not have nIterations rows');
        end
        if size(results.simulated.permutationH0.statVal, 2) ~= pX
            error('permutationH0.statVal does not have pX columns');
        end
    case 'bootstrapStability'
        if size(results.simulated.bootstrapStability.values, 1) ~= nIterations
            error('bootstrapStability.values does not have nIterations rows');
        end
        if size(results.simulated.bootstrapStability.values, 2) ~= pX
            error('bootstrapStability.values does not have pX columns');
        end
end

end