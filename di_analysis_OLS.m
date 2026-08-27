function results = di_analysis_OLS(di_cfg, Y_orig, X_orig, rowIdx)
% Ordinary Least Squares (OLS) implementations of various statistical
% analyses, including correlations, independent and paired sample t-tests.
% The OLS implementation allows to group these analyses under a
% unified approach. it also allows extensions (e.g., by adding covariates
% to the above mentioned example analyses). 
% In classic statistical textbooks, it is presented as y = Xb where y is
% the data, X the dummy codes or more data, and b are the beta
% coefficients. In linear algebra, it's more often presented as b = Ax but
% it's identical. This script adopts a linear algebra approach in the
% computations but uses a statistics terminology (for now... perhaps I'll
% change terminology to linear algebra too)
%
% This function computes betas through the backslash operator b = X \ y;
% (or x = A \ b in linear algebra). If y is not a columns vector but a
% whole matrix, betas are computed separately for each column of y
% (features). X is either the design matrix or more data and sometimes
% optional terms (e.g., dummy codes for repeated measures). Betas are 
% computed always--in all cases. 
% Parametric p-values (via standard error and t-test scores) are computed 
% only when required by the analysis type. Parametric p values might be the
% end goal or the basis for cluster forming.
%
% betas (observed vs simulated) are always used for inference. They are put
% into a "statVal" variable, both observed and simulated (permutation or
% bootstrap). so, not using the more typical t or r. However, downstream 
% di_inference.m compares observed vs simulated betas. Because betas are
% proportional to t and r, the outcome of the comparison is the same. By
% comparing betas we gain in simplicity (same coefficients for all tests)
% and computational speed
%
%
% options (set via di_cfg.analysis.OLS struct):
%
%   varianceType : 'unequal' (default) | 'equal'
%                  'unequal' -> Welch t + Satterthwaite df
%                  'equal'   -> OLS pooled residual SE
%                  Only consulted for design [0 1] with parametric p-values.
%
%   ranked       - false (default) | true
%                  true applies a rank transform to X and each column of Y 
%                  before the solve. 
%                  Only consulted for design [0 0]. 
%
%   standardize  : false (default) | true
%                  z-score X and each column of Y before the solve.
%                  For design [0 0], beta(2,:) is on a correlation scale:
%                  Pearson r when ranked=false, Spearman rho when ranked=true.
%                  With extra regressors (e.g., subject dummies), beta(2,:)
%                  is a standardized partial slope (not a simple correlation).
% 
% Statistical equivalences:
%
%   Design=[0 0], ranked=false, standardize=true  =  Pearson correlation
%                                       statVal is Pearson r-scale
%                                       rVal reports Pearson r (regardless of standardize=true/false)
%
%   Design=[0 0], ranked=true,  standardize=true  =  Spearman correlation
%                                       statVal is Spearman rho-scale
%                                       rVal reports Spearman rho (regardless of standardize=true/false)
%
%   Design=[0 1], varianceType='equal'            =  Student's independent-samples t-test
%                                       statVal is the model slope (difference-scale effect)
%                                       the t values are in tVal           
%                                       rVal is a descriptive correlation with group coding
%
%   Design [0 1], varianceType='unequal'          =  Welch independent-samples t-test (with Satterthwaite df)
%                                       statVal is the model slope (difference-scale effect)
%                                       the t values are in tVal           
%                                       rVal is a descriptive correlation with group coding
%
%   Design [1 0]                                  =  paired-samples t-test (subjects as fixed effects)
%                                       statVal is the predictor partial slope (subject effects included)
%                                       the t values are in tVal           
%                                       rVal is descriptive partial r (controlling subject effects)
%
% INPUT:
%   di_cfg        - validated analysis configuration struct
%                   di_cfg.analysis.designCode must be set (done by di_analyze)
%   X_orig        - predictor / condition / group vector (m x 1)
%   Y_orig        - outcome/data matrix (m x pY)
%   rowIdx        - resampling row indices from di_reorderRowsGenerate (m x nIterations)
%
% OUTPUT:
%   results       - results structure with the following fields:
%
%   results.observed    subfield containing descriptive (not inferential) 
%                       values. it contains the following subfields:
%   results.observed.statVal         [1 x pY]  OLS beta
%   results.observed.pVal_parametric [1 x pY]  parametric p-values from t-statistics
%   results.observed.df_parametric   [1 x pY]  feature-wise degrees of freedom for pVal_parametric
%   results.observed.tVal            [1 x pY]  t-statistic, observed data only (descriptive)
%   results.observed.rVal            [1 x pY]  correlation coefficient, observed data only (descriptive)
%   results.observed.dVal            [1 x pY]  Cohen's d, observed data only (descriptive)
%   results.observed.clusters        struct()  empty for feature-level analyses
%                             struct with .clusterMembership_obs, .clustIDList_obs,
%                             .metrics_obs for cluster analyses
%
%   results.simulated   subfield containing the simulated values
%   results.simulated.permutationH0.statVal         [nIterations x pY]  null distribution of beta
%   results.simulated.permutationH0.clusterMetrics  [1 x nIterations] struct array
%   results.simulated.bootstrapStability.values     [nIterations x pY]  bootstrap betas
% TO DO: perhaps change "results.simulated.bootstrapStability.values" to "results.simulated.bootstrapStability.statVal" for consistency in terminology
%
% NOTES:
%   - pVal_parametric and df_parametric are observed-data companions from tVal.
%     They are not used for permutation/bootstrap inference.
%   - tVal and rVal are descriptive fields computed once from the observed data (itIdx==1)
%     at negligible cost. They are NOT stored in the null distribution.
%   - rVal is always computed via OLS on z-scored X_fit and Y_fit (column-wise), regardless
%     of the standardize option. This yields Pearson r (design [0 0], ranked=false)
%     or Spearman rho (design [0 0], ranked=true, because X_fit/Y_fit are ranks).
%   - For design [1 0], rVal is a partial r: subjectDummies are kept in D_std (unstandardized)
%     so subject effects are partialled out, matching the paired design logic.
%   - Welch SE (varianceType='unequal') is computed from the PERMUTED group
%     membership in each iteration, which is the correct permutation-test behavior.
%   - Subject fixed effects (design [1 0]) are built once outside the loop from
%     the original observationID; only the condition predictor X is permuted.
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% default options 

% ranked
ranked = false;
if isfield(di_cfg.analysis, 'OLS') && isfield(di_cfg.analysis.OLS, 'ranked')
    ranked = di_cfg.analysis.OLS.ranked;
end

% variance type
varianceType = 'unequal';
if isfield(di_cfg.analysis, 'OLS') && isfield(di_cfg.analysis.OLS, 'varianceType')
    varianceType = di_cfg.analysis.OLS.varianceType;
end

% standardize
standardize = false;
if isfield(di_cfg.analysis, 'OLS') && isfield(di_cfg.analysis.OLS, 'standardize')
    standardize = di_cfg.analysis.OLS.standardize;
end

%% shortcuts

nIterations       = di_cfg.analysis.nIterations;
inferenceLevel    = di_cfg.analysis.inferenceLevel;
if isfield(di_cfg.analysis, 'testingApproach')
    testingApproach = di_cfg.analysis.testingApproach;
else
    % bootstrap paths do not require explicit testingApproach
    testingApproach = 'empirical';
end
analysisObjective = di_cfg.analysis.objective;
designCode        = di_cfg.analysis.designCode;
ignore_col        = logical(di_cfg.analysis.ignore_col);

dimKeys  = fieldnames(di_cfg.featureDims);
dimSizes = cellfun(@(k) length(di_cfg.featureDims.(k).vec), dimKeys);
pY       = prod(dimSizes);

m = size(Y_orig, 1);

% flags controlling what gets computed (if just betas or also se/df/t/r)
needPVal    = ~strcmp(testingApproach, 'empirical');
needCluster = strcmp(inferenceLevel, 'cluster');

% sanity check: if we need cluster-based inference, we definitely also need parametric pvalues
if needCluster && ~needPVal
    error('Cluster-based inference requires parametric p-values. Set di_cfg.analysis.testingApproach = ''parametric''.');
end


%% sanity checks

if pY < 1
    error('pY must be at least 1');
end
if size(Y_orig, 2) ~= pY
    error('Y_orig must have pY columns');
end
if nIterations < 1
    error('nIterations must be at least 1');
end
if size(X_orig, 1) ~= size(Y_orig, 1)
    error('X_orig and Y_orig must have the same number of rows');
end
if size(X_orig, 2) ~= 1
    error('X_orig must be a column vector (m x 1)');
end
if size(rowIdx, 1) ~= size(X_orig, 1)
    error('rowIdx must have same number of rows as X_orig');
end
if size(rowIdx, 2) ~= nIterations
    error('rowIdx must have nIterations columns');
end
if numel(ignore_col) ~= pY
    error('di_cfg.analysis.ignore_col must have pY elements');
end
if ~any(~ignore_col)
    error('All features are ignored (ignore_col all true). Nothing to analyze.');
end

% design-specific checks
if isequal(designCode, [0 1])
    if numel(unique(X_orig)) ~= 2
        error('design [0 1]: X_orig must have exactly 2 unique group codes');
    end
end
if isequal(designCode, [1 0])
    if numel(unique(X_orig)) ~= 2
        error('design [1 0]: X_orig must have exactly 2 unique condition codes');
    end
    if ~isfield(di_cfg.analysis, 'dataStruct') || ~istable(di_cfg.analysis.dataStruct)
        error('design [1 0]: di_cfg.analysis.dataStruct must be a table');
    end
    if ~ismember('observationID', di_cfg.analysis.dataStruct.Properties.VariableNames)
        error('design [1 0]: dataStruct must contain observationID column');
    end
end
if needCluster
    if ~isfield(di_cfg.analysis, 'clusterParams')
        error('clusterParams must be defined for inferenceLevel = ''cluster''');
    end
end

% Subject centering is intentionally not supported in OLS.
% Repeated-measures effects are handled via subject fixed-effects dummies.
if isfield(di_cfg.analysis, 'subjectCentering') && strcmp(di_cfg.analysis.subjectCentering, 'center')
    error(['di_analysis_OLS: subjectCentering=''center'' is not supported for OLS. ' ...
        'OLS handles subject effects via subject dummies (fixed effects), not by subject-centering X.']);
end

% options consistency warnings
if ~isequal(designCode, [0 0]) && ranked
    warning('Diagonale:OLS', 'ranked is only used for design [0 0]; ignored here');
end
if ~isequal(designCode, [0 1]) && ~strcmp(varianceType, 'unequal')
    warning('Diagonale:OLS', 'varianceType is only used for design [0 1]; ignored here');
end

%% build subject dummies for matrix X 
% (from y = Xb)

% downstream we will augment X to become [intercept, predictor, subject_dummies]. 
% currently X contains only the predictor
predictor_colIdx = 2; % predictor column index
% TO DO (maybe): make the colIdx a parameter, because in the future there could
% be predictors in more than the 2nd column (if I allow OLS to compare more
% than 2 things)

% intercept is very simple (col of 1s) so it will be added later

% subject dummies
subjectDummies = []; % initialize to nothing
% subject dummies will be filled only for repeated-measures design [1 0]
% n-1 columns (i..e, drop last subject to keep full rank)
if designCode(1) == 1
    obsIDVec  = di_cfg.analysis.dataStruct.observationID;
    uniqueObs = unique(obsIDVec,'stable');
    nSubj     = numel(uniqueObs);
    subjectDummies = zeros(m, nSubj - 1);
    for sdcolIdx = 1:(nSubj-1)
        subjectDummies(:, sdcolIdx) = double(obsIDVec == uniqueObs(sdcolIdx));
    end
end

%%
% precompute rank transform of Y for rank-based correlation.
% this precomupation is only valid for permutations, because Y does not 
% change anyway (only X is permuted in the permutations)
if isequal(designCode, [0 0]) && ranked
    Y_orig_ranked = zeros(size(Y_orig));
    for colIdx = 1:pY
        Y_orig_ranked(:, colIdx) = tiedrank(Y_orig(:, colIdx));
    end
end

%% main analysis loop

for itIdx = 1:nIterations

    % apply row reordering:
    %   permutation : X rows shuffled, Y rows unchanged (breaks X-Y link)
    %   bootstrap   : both Y and X rows resampled with replacement
    [X, Y] = di_reorderRowsApply(di_cfg, X_orig, Y_orig, rowIdx, itIdx);

    % --- prepare X_fit and Y_fit (rank transform / standardize if requested) ---

    % Start from raw iteration data, then apply requested transforms in sequence.
    X_fit = X;
    Y_fit = Y;

    if isequal(designCode, [0 0]) && ranked
        % rank transform X at each iteration (both permutation and bootstrap)
        X_fit = tiedrank(X_fit);

        % rank-transform Y at each iteration (only for bootstrap). for
        % permutations not needed as Y does not get resampled and it's
        % already been ranked transformed outside this loop
        if strcmp(analysisObjective, 'bootstrapStability')
            Y_fit = zeros(size(Y_fit));
            for colIdx = 1:pY
                Y_fit(:, colIdx) = tiedrank(Y(:, colIdx));
            end
        else
            Y_fit = Y_orig_ranked;
        end
    end

    if standardize
        % z-score predictor and outcomes after optional ranking.
        % - ranked=false: slope equals Pearson r (single predictor).
        % - ranked=true : slope equals Spearman rho (single predictor).
        X_fit = (X_fit - mean(X_fit)) / std(X_fit);
        Y_fit = (Y_fit - mean(Y_fit, 1)) ./ std(Y_fit, 0, 1);
    end


    % --- build design matrix ---
    % X_fit_augmented = [intercept, predictor, subject_dummies (for paired design)]
    X_fit_augmented = [ones(m, 1), X_fit, subjectDummies];

    % --- OLS solve (vectorized across non-ignored features only) ---
    % beta_2use: [nPred x pY_2use]
    Y_fit_2use = Y_fit(:, ~ignore_col);
    beta_2use   = X_fit_augmented \ Y_fit_2use;
    statVal = nan(1, pY);
    statVal(1, ~ignore_col) = beta_2use(predictor_colIdx, :);

    % --- parametric p-values (only for parametric analysis types) ---
    if needPVal
        % compute parametric stats 
        % based on classical OLS or Welch approach
        pY_2use = nnz(~ignore_col);
        [tStat_2use, pVal_2use, df_iter_2use] = di_computeParametricStats( ...
            designCode, varianceType, Y_fit_2use, X_fit_augmented, beta_2use, predictor_colIdx, pY_2use);

        tStat = nan(1, pY);
        pVal = nan(1, pY);
        df_iter = nan(1, pY);
        tStat(1, ~ignore_col) = tStat_2use;
        pVal(1, ~ignore_col) = pVal_2use;
        df_iter(1, ~ignore_col) = df_iter_2use;

    end

    % --- cluster forming (inferenceCluster branch only) ---
    if needCluster
        if designCode(1)==0 && designCode(2)==0
            % beta for the current predictor (statVal) is already a correlation for [0 0] and standardize==true
            % otherwise (standardize==false) it's a slope, but still the best metric.
            clusterMetric = statVal; 
        else
            clusterMetric = tStat;
        end
        [clusterMembership, clustIDList, metrics] = di_clusterForming(di_cfg, clusterMetric, pVal);
    end

    % --- store observed (first iteration = identity permutation / original data) ---
    % for itIdx==1 (first iteration) under certain circumstances we have
    % only computed betas so far. however, for descriptive reasons, it
    % would be nice to have pvalue, df, tStat, rStat as well
    if itIdx == 1

        statVal_obs = statVal;  % beta: inference metric, echoed into null distribution each iteration

        % --- tStat_obs: t-statistic for observed data only (descriptive) ---
        % For parametric branch: tStat was already computed above; reuse it.
        % For empirical branch: compute SE here 
        if needPVal
            pVal_parametric_obs = pVal;
            df_parametric_obs   = df_iter;
            tStat_obs           = tStat;  % already computed in the needPVal block above
        else
            % compute parametric stats (only for iteration 1)
            % based on classical OLS or Welch approach
            pY_2use = nnz(~ignore_col);
            [tStat_obs_2use, pVal_parametric_obs_2use, df_parametric_obs_2use] = di_computeParametricStats( ...
                designCode, varianceType, Y_fit_2use, X_fit_augmented, beta_2use, predictor_colIdx, pY_2use);
            tStat_obs = nan(1, pY);
            pVal_parametric_obs = nan(1, pY);
            df_parametric_obs = nan(1, pY);
            tStat_obs(1, ~ignore_col) = tStat_obs_2use;
            pVal_parametric_obs(1, ~ignore_col) = pVal_parametric_obs_2use;
            df_parametric_obs(1, ~ignore_col) = df_parametric_obs_2use;
        end

        % --- rVal_obs: Pearson r or Spearman rho via OLS on z-scored inputs (descriptive) ---
        % Z-scoring is applied here regardless of the standardize option, ensuring
        % rVal is always on the [-1, 1] scale.
        % - Pearson case: z-score raw X_fit and Y_fit -> OLS slope = r exactly.
        % - Spearman case: X_fit and Y_fit are already ranks at this point;
        %   z-scoring ranks then running OLS gives rho exactly.
        % - Design [1 0]: subjectDummies remain unstandardized (they were precomputed outside the loop) in D_std, so rVal
        %   is partial r controlling for subject effects (consistent with the model).
        X_fit_z  = (X_fit - mean(X_fit)) / std(X_fit);                  % z-score predictor
        Y_fit_z  = (Y_fit - mean(Y_fit, 1)) ./ std(Y_fit, 0, 1);        % z-score each outcome
        Y_fit_z_2use = Y_fit_z(:, ~ignore_col);
        D_std    = [ones(m, 1), X_fit_z, subjectDummies];                % design matrix with z-scored predictor
        beta_std = D_std \ Y_fit_z_2use;                                 % [nPred x pY_2use]
        rVal_obs = nan(1, pY);
        rVal_obs(1, ~ignore_col) = beta_std(predictor_colIdx, :);     % [1 x pY]

        % --- dVal_obs: Cohen's d from observed raw data only (descriptive) ---
        % [0 1] -> pooled-SD Cohen's d (independent groups)
        % [1 0] -> Cohen's dz based on subject-level condition differences
        if isequal(designCode, [0 1])
            dVal_2use = di_computeCohensD_between(Y(:, ~ignore_col), X);
            dVal_obs = nan(1, pY);
            dVal_obs(1, ~ignore_col) = dVal_2use;
        elseif isequal(designCode, [1 0])
            obsID_for_d = di_cfg.analysis.dataStruct.observationID(rowIdx(:, itIdx));
            dVal_2use = di_computeCohensDz_within(Y(:, ~ignore_col), X, obsID_for_d);
            dVal_obs = nan(1, pY);
            dVal_obs(1, ~ignore_col) = dVal_2use;
        else
            dVal_obs = nan(1, pY);
        end

        if needCluster
            clusterMembership_obs = clusterMembership;
            clustIDList_obs       = clustIDList;
            metrics_obs           = metrics;
        end
    end

    % --- store resampled values ---
    switch analysisObjective

        case 'permutationH0testing'
            if itIdx == 1 % initialize
                if needCluster
                    simulatedMetrics = repmat( ...
                        struct('id', [], 'size', [], 'mass', [], 'mostExtremeVal', []), ...
                        1, nIterations);
                else
                    statVal_perm = nan(nIterations, pY);
                end
            end
            if needCluster
                simulatedMetrics(1, itIdx).id             = metrics.id;
                simulatedMetrics(1, itIdx).size           = metrics.size;
                simulatedMetrics(1, itIdx).mass           = metrics.mass;
                simulatedMetrics(1, itIdx).mostExtremeVal = metrics.mostExtremeVal;
            else
                statVal_perm(itIdx, :) = statVal;
            end

        case 'bootstrapStability'
            if itIdx == 1 % initialize
                statVal_boot = zeros(nIterations, pY);
            end
            statVal_boot(itIdx, :) = statVal;

    end

    di_counter(itIdx, nIterations);

end

%% collate results

% observed, feature level
results.observed.statVal = statVal_obs;  % [1 x pY] OLS beta — inference metric
results.observed.pVal_parametric = pVal_parametric_obs;  % [1 x pY] parametric p-value
results.observed.df_parametric   = df_parametric_obs;    % [1 x pY] parametric degrees of freedom

% descriptive statistics for observed data only (not stored in null distribution)
results.observed.tVal = tStat_obs;   % [1 x pY] t-statistic
results.observed.rVal = rVal_obs;    % [1 x pY] Pearson r or Spearman rho (standardized OLS)
results.observed.dVal = dVal_obs;    % [1 x pY] Cohen's d (pooled for [0 1], dz for [1 0])

% observed, cluster level
if needCluster
    results.observed.clusters.clusterMembership_obs = clusterMembership_obs;
    results.observed.clusters.clustIDList_obs       = clustIDList_obs;
    results.observed.clusters.metrics_obs           = metrics_obs;
else
    results.observed.clusters = struct();
end

% simulated
switch analysisObjective
    case 'permutationH0testing'
        if needCluster
            results.simulated.permutationH0.clusterMetrics = simulatedMetrics; % [1 x nIterations]
        else
            results.simulated.permutationH0.statVal = statVal_perm; % [nIterations x pY]
        end
    case 'bootstrapStability'
        results.simulated.bootstrapStability.values = statVal_boot; % [nIterations x pY]
end

%% validate output

if size(results.observed.statVal, 2) ~= pY
    error('observed.statVal does not have pY columns');
end
if size(results.observed.pVal_parametric, 2) ~= pY
    error('observed.pVal_parametric does not have pY columns');
end
if size(results.observed.df_parametric, 2) ~= pY
    error('observed.df_parametric does not have pY columns');
end
if size(results.observed.tVal, 2) ~= pY
    error('observed.tVal does not have pY columns');
end
if size(results.observed.rVal, 2) ~= pY
    error('observed.rVal does not have pY columns');
end
if size(results.observed.dVal, 2) ~= pY
    error('observed.dVal does not have pY columns');
end

if needCluster
    if size(results.observed.clusters.clusterMembership_obs, 2) ~= pY
        error('clusters.clusterMembership_obs does not have pY columns');
    end
    if di_cfg.analysis.verbose && isempty(results.observed.clusters.clustIDList_obs)
        if di_cfg.analysis.clusterParams.clusterFormingPvalThreshold > 0
            warning('Diagonale: no clusters found in observed data');
        end
    end
end

switch analysisObjective
    case 'permutationH0testing'
        if needCluster
            if length(results.simulated.permutationH0.clusterMetrics) ~= nIterations
                error('permutationH0.clusterMetrics does not have nIterations elements');
            end
        else
            if size(results.simulated.permutationH0.statVal, 1) ~= nIterations
                error('permutationH0.statVal does not have nIterations rows');
            end
            if size(results.simulated.permutationH0.statVal, 2) ~= pY
                error('permutationH0.statVal does not have pY columns');
            end
        end
    case 'bootstrapStability'
        if size(results.simulated.bootstrapStability.values, 1) ~= nIterations
            error('bootstrapStability.values does not have nIterations rows');
        end
        if size(results.simulated.bootstrapStability.values, 2) ~= pY
            error('bootstrapStability.values does not have pY columns');
        end
end

if all(isnan(results.observed.statVal))
    warning('Diagonale: observed.statVal is all NaN - possible computation failure');
end

end

function dVal = di_computeCohensD_between(Y, X)
% Cohen's d for independent groups using pooled SD.
% TO DO: in future, condider computing a Welch compatible version of Cohen's d that does not assume equal variance between groups
groupVals = unique(X);
if numel(groupVals) ~= 2
    error('Cohen''s d (between) requires exactly two group levels in X.');
end
mask1 = (X == max(groupVals));
mask0 = (X == min(groupVals));
n1 = sum(mask1);
n0 = sum(mask0);
if n1 <= 1 || n0 <= 1
    error('Cohen''s d (between) requires at least two observations in each group.');
end

mu1 = mean(Y(mask1, :), 1);
mu0 = mean(Y(mask0, :), 1);
s1sq = var(Y(mask1, :), 0, 1);
s0sq = var(Y(mask0, :), 0, 1);
sp = sqrt(((n1 - 1) * s1sq + (n0 - 1) * s0sq) / (n1 + n0 - 2));

dVal = (mu1 - mu0) ./ sp;
dVal(sp == 0) = NaN;
end

function dVal = di_computeCohensDz_within(Y, X, obsID)
% Cohen's dz for paired/repeated design: mean(subject differences) / SD(subject differences).
condVals = unique(X);
if numel(condVals) ~= 2
    error('Cohen''s dz (within) requires exactly two condition levels in X.');
end

condHi = max(condVals);
condLo = min(condVals);
uniqueObs = unique(obsID, 'stable');
nObs = numel(uniqueObs);

D = nan(nObs, size(Y, 2));
for obsIdx = 1:nObs
    id = uniqueObs(obsIdx);
    maskHi = (obsID == id) & (X == condHi);
    maskLo = (obsID == id) & (X == condLo);
    if ~any(maskHi) || ~any(maskLo)
        continue
    end
    yHi = mean(Y(maskHi, :), 1);
    yLo = mean(Y(maskLo, :), 1);
    D(obsIdx, :) = yHi - yLo;
end

validRows = all(~isnan(D), 2);
D = D(validRows, :);
if size(D, 1) <= 1
    error('Cohen''s dz (within) requires at least two complete paired observations.');
end

muD = mean(D, 1);
sdD = std(D, 0, 1);
dVal = muD ./ sdD;
dVal(sdD == 0) = NaN;
end

function [tStat, pVal, df_vec] = di_computeParametricStats( ...
    designCode, varianceType, Y_fit, X_fit_augmented, beta, predictor_colIdx, pY)
% Local function to keep the code above a bit tidier (since this function is used twice)
% Compute t, p, and df for OLS predictor coefficients.
%
% MODE 1 (Welch): design [0 1] with varianceType='unequal'.
% MODE 2 (Classical OLS): all other supported designs.

if isequal(designCode, [0 1]) && strcmp(varianceType, 'unequal')
    modeIdx = 1;
else
    modeIdx = 2;
end

b = beta(predictor_colIdx, :);

switch modeIdx
    case 1
    
    % here we assume X has only two-level group codes, so we can estimate
    % the SE for the coefficient of the only predictor (groups) directly, without needing to compute the full covariance matrix as in classical OLS (as in mode 2, where we are not assuming one two-level predictor only).
    % Use the actual predictor column from the fitted design matrix (e.g., in case it gets standardized).
    predictorVals = X_fit_augmented(:, predictor_colIdx);
    % the OLS slope is the mean difference divided by the spacing
    % between the two predictor levels, the Welch SE must use that same spacing.
    groupVals = unique(predictorVals);
    mask1 = (predictorVals == max(groupVals));
    mask0 = (predictorVals == min(groupVals));
    n1    = sum(mask1);
    n0    = sum(mask0);
    % sanity checks for number of observations per group
    if n1 <= 1 || n0 <= 1
        error('Welch mode requires at least two observations in each group.');
    end
    % Compyte spacing between the two predictor levels used in the model.
    codeSpacing = abs(max(groupVals) - min(groupVals));
    if codeSpacing == 0
        error('Welch mode requires two distinct predictor levels.');
    end
    s1sq  = var(Y_fit(mask1, :), 0); % use n-1 in denominator for sample variance (unbiased estimator), as per Welch's method
    s0sq  = var(Y_fit(mask0, :), 0); 
    SE = sqrt(s1sq/n1 + s0sq/n0) / codeSpacing; % Welch SE for the OLS slope, per feature

    % Step 2: compute t-statistic for each feature.
    tStat = b ./ SE;

    % Step 3: compute degrees of freedom (Satterthwaite, feature-wise).
    num_df = (s1sq / n1 + s0sq / n0) .^ 2;
    den_df = (s1sq / n1) .^ 2 / (n1 - 1) + (s0sq / n0) .^ 2 / (n0 - 1);
    df_vec = num_df ./ den_df;


case 2
    
    % while in this case we are assuming equal variances across levels, we are also being more general and allowing multiple predictors
    % (e.g., subject dummies in the repeated measure design), so we need an approach that generalized to any predictor, differently from mode 1.
    % such general approach would be headache for mode 1. but mode 1 is restricted to one predictor with two-levels.

    % e.g., p 158 of Fox (2009) A mathematical primer for social statistics
    % e.g., p.134 of Kutner et al. (2004) Applied Linear Statistical Models, Ed. 5 

    % Step 1: compute sigma squared (residual variance) from the fitted model, per each feature of Y_fit.
    E       = Y_fit - X_fit_augmented * beta; % residuals (per feature) = difference between the original data and predicted (beta based linear combination of predictors)
    RSS     = sum(E .^ 2, 1); % residual sum of squares (per feature)
    % residuals degrees of freedom = n-(k+1), where
    % n = number of observations
    n = size(X_fit_augmented, 1);
    % k + 1 = number of linearly independent parameters in the design matrix,
    % including the intercept. Using rank() is more robust than counting
    % columns size(X_fit_augmented,2) because it handles collinearity.
    kPlus1 = rank(X_fit_augmented);
    df_resid = n - kPlus1;
    % sanity check on df_resid to avoid division by zero or negative df
    if df_resid <= 0
        error('Residual degrees of freedom is zero or negative. Check design matrix for full rank and ensure n > k+1.');
    end
    % residual variance (unbiased estimate) (feature-wise):
    sigma2  = RSS / df_resid;

    % Step 2: compute SE for the predictor coefficient 
    % using  sigma2 and the design matrix X_fit_augmented.
    % covariance matrix for a single feature and for beta_hat would be: sigma2 * (X'X)^(-1)
    % the standard errors for the regression coefficients would be the squared roots of the diagonal entries of the covariance matrix:
    % but in this case, sigma2 can be a vector (per feature of Y), so 
    % let's compute the inverse of X'X first (as that is dependent only on X and not Y)
    % then let's use only the diagonal entry corresponding to that predictor,
    % and then multiply by sigma2 (which is per feature) to get the SE associated to that predictor for each Y feature.
    XtX_inv = inv(X_fit_augmented' * X_fit_augmented); % [nPred x nPred]
    XtX_inv_jj = XtX_inv(predictor_colIdx, predictor_colIdx); % the j,j entry of XtX_inv corresponding to the predictor coefficient (a scalar)
    SE = sqrt(sigma2 * XtX_inv_jj); % SE for the predictor coefficient, per feature (1 x pY)
    
    % Step 3: compute t-statistic for each feature.
    tStat = b ./ SE; % technically, it is beta_hat minus the null value (zero here), divided by SE

    % Step 4: collate degrees of freedom for all features into a vector 
    % (same df for all features in classical OLS, so just repeat the same value)
    df_vec = repmat(df_resid, 1, pY);

end

% compute two-tailed p-values using t and df.
pVal = 2 * tcdf(-abs(tStat), df_vec); % 1 x pY

end
