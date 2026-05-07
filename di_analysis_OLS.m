function results = di_analysis_OLS(di_cfg, X_orig, Y_orig, rowIdx)
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


% TO DO: the parametric computations can be done in a separate function (at the bottom) and be reused internally
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
analysisType      = di_cfg.analysis.type;
analysisObjective = di_cfg.analysis.objective;
designCode        = di_cfg.analysis.designCode;

dimKeys  = fieldnames(di_cfg.dimensions);
dimSizes = cellfun(@(k) length(di_cfg.dimensions.(k).vec), dimKeys);
pY       = prod(dimSizes);

m = size(Y_orig, 1);

% flags controlling what gets computed (if just betas or also se/df/t/r)
needPVal    = ~strcmp(analysisType, 'empiricalFeature_inferenceFeature');
needCluster = strcmp(analysisType, 'parametricFeature_inferenceCluster');

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
        error('clusterParams must be defined for parametricFeature_inferenceCluster');
    end
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

    % --- OLS solve (vectorized across all pY outcomes) ---
    % beta: [nPred x pY]
    beta   = X_fit_augmented \ Y_fit;
    statVal = beta(predictor_colIdx, :);  % [1 x pY]

    % --- parametric p-values (only for parametric analysis types) ---
    if needPVal
        keyboard % double check this
        % TO DO: convert the code below to a function at the bottom that can be reused in the two parts // easier for maintenance and tidier
        if isequal(designCode, [0 1]) && strcmp(varianceType, 'unequal')
            % MODE: Independent groups, unequal variances (Welch).
            % Uses group-specific variance and feature-wise Satterthwaite df.
            % Group masks are computed from the current X (permuted in H0 runs).
            condVals = unique(X);
            mask1 = (X == max(condVals));
            mask0 = (X == min(condVals));
            n1    = sum(mask1);
            n0    = sum(mask0);
            % Step 1: build variance ingredients per feature (group-wise variances).
            s1sq  = var(Y_fit(mask1, :));        % [1 x pY]
            s0sq  = var(Y_fit(mask0, :));        % [1 x pY]
            % Step 2: compute coefficient SE from those variance ingredients.
            SE_W  = sqrt(s1sq / n1 + s0sq / n0); % [1 x pY]
            % Step 3: compute t-statistic for each feature.
            tStat = statVal ./ SE_W;
            % Step 4: compute degrees of freedom (Satterthwaite, feature-wise).
            num_df  = (s1sq / n1 + s0sq / n0) .^ 2;
            den_df  = (s1sq / n1) .^ 2 / (n1 - 1) + (s0sq / n0) .^ 2 / (n0 - 1);
            df_satt = num_df ./ den_df;           % [1 x pY]
            % Step 5: compute two-tailed p-values from t and df.
            pVal    = 2 * tcdf(-abs(tStat), df_satt);
            df_iter = df_satt;

        else
            % MODE: Classical OLS residual-variance t-test.
            % Used for equal-variance groups, paired design, and correlation design.
            % Degrees of freedom are model residual df from the design matrix rank.
            % Step 1: build variance ingredients per feature (RSS from residuals).
            E       = Y_fit - X_fit_augmented * beta;                          % [m x pY] residuals
            RSS     = sum(E .^ 2, 1);                            % [1 x pY]
            % Step 2: compute coefficient SE from those variance ingredients.
            df_resid = m - rank(X_fit_augmented);                              % scalar, robust to collinearity
            sigma2  = RSS / df_resid;                            % [1 x pY]
            DtD_inv = inv(X_fit_augmented' * X_fit_augmented);                               % [nPred x nPred]
            SE      = sqrt(sigma2 * DtD_inv(predictor_colIdx, predictor_colIdx)); % [1 x pY]
            % Step 3: compute t-statistic for each feature.
            tStat   = statVal ./ SE;                             % [1 x pY]
            % Step 4: compute degrees of freedom (residual df, shared across features).
            % Step 5: compute two-tailed p-values from t and df.
            pVal    = 2 * tcdf(-abs(tStat), df_resid);           % [1 x pY]
            df_iter = repmat(df_resid, 1, pY);
        end

    end

    % --- cluster forming (inferenceCluster branch only) ---
    if needCluster
        keyboard % double check this
        [clusterMembership, clustIDList, metrics] = di_clusterForming(di_cfg, statVal, pVal);
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
            % empirical branch: one-time SE computation for descriptive tVal
            if isequal(designCode, [0 1]) && strcmp(varianceType, 'unequal')
                % MODE: Observed-data Welch branch (unequal variances).
                % Uses original (unshuffled) group membership at itIdx==1.
                condVals_obs = unique(X);
                mask1_obs    = (X == max(condVals_obs));
                mask0_obs    = (X == min(condVals_obs));
                n1_obs       = sum(mask1_obs);
                n0_obs       = sum(mask0_obs);
                % Step 1: build variance ingredients per feature (group-wise variances).
                s1sq_obs     = var(Y_fit(mask1_obs, :));                    % [1 x pY]
                s0sq_obs     = var(Y_fit(mask0_obs, :));                    % [1 x pY]
                % Step 2: compute coefficient SE from those variance ingredients.
                SE_obs       = sqrt(s1sq_obs / n1_obs + s0sq_obs / n0_obs); % [1 x pY]
                % Step 4: compute degrees of freedom (Satterthwaite, feature-wise).
                num_df_obs   = (s1sq_obs / n1_obs + s0sq_obs / n0_obs) .^ 2;
                den_df_obs   = (s1sq_obs / n1_obs) .^ 2 / (n1_obs - 1) + (s0sq_obs / n0_obs) .^ 2 / (n0_obs - 1);
                df_parametric_obs = num_df_obs ./ den_df_obs;
            else
                % MODE: Observed-data classical OLS residual-variance branch.
                % Applies to equal-variance groups, paired design, and correlation design.
                % Step 1: build variance ingredients per feature (RSS from residuals).
                E_obs        = Y_fit - X_fit_augmented * beta;                              % [m x pY]
                RSS_obs      = sum(E_obs .^ 2, 1);                            % [1 x pY]
                % Step 2: compute coefficient SE from those variance ingredients.
                df_resid_obs = m - rank(X_fit_augmented);                                   % scalar, robust to collinearity
                sigma2_obs   = RSS_obs / df_resid_obs;                        % [1 x pY]
                DtD_inv_obs  = inv(X_fit_augmented' * X_fit_augmented);                                   % [nPred x nPred]
                SE_obs       = sqrt(sigma2_obs * DtD_inv_obs(predictor_colIdx, predictor_colIdx)); % [1 x pY]
                % Step 4: compute degrees of freedom (residual df, shared across features).
                df_parametric_obs = repmat(df_resid_obs, 1, pY);
            end
            % Step 3: compute t-statistic for each feature.
            tStat_obs = statVal ./ SE_obs;  % [1 x pY]
            % Step 5: compute two-tailed p-values from t and df.
            pVal_parametric_obs = 2 * tcdf(-abs(tStat_obs), df_parametric_obs);
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
        D_std    = [ones(m, 1), X_fit_z, subjectDummies];                % design matrix with z-scored predictor
        beta_std = D_std \ Y_fit_z;                                      % [nPred x pY]
        rVal_obs = beta_std(predictor_colIdx, :);                           % [1 x pY]

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
                    statVal_perm = zeros(nIterations, pY);
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

if needCluster
    if size(results.observed.clusters.clusterMembership_obs, 2) ~= pY
        error('clusters.clusterMembership_obs does not have pY columns');
    end
    if di_cfg.analysis.verbose && isempty(results.observed.clusters.clustIDList_obs)
        if di_cfg.analysis.clusterParams.clusterFormingThreshold > 0
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
