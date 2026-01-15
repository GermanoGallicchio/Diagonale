function results = di_inference(di_cfg,results)
% make statistical inference based on results (observed and simulated). the
% approach depends on analysis objective (permutation vs. bootstrap) and
% analysis type:
%
% - permutation H0 testing: compare observed statistics against null distribution 
%   of those values to compute p-values.
%   * empiricalFeature_inferenceFeature:
%   * empiricalFeature_inferenceCluster:
%   * parametricFeature_inferenceFeature:
%   * parametricFeature_inferenceCluster:
%   * PLS_SVD
%   * AJIVE
%
% - bootstrap stability inference: use resampled statistics to compute
%   bootstrap ratio and confidence intervals.
%   note 1: no multiple comparison correction: none (exploratory). 
%   note 2: purpose is not testing, but stability estimation (especially
%       useful after testing has already been done)
%
%   * empiricalFeature_inferenceFeature:
%   * empiricalFeature_inferenceCluster:
%   * parametricFeature_inferenceFeature:
%   * parametricFeature_inferenceCluster:
%   * PLS_SVD:
%   * AJIVE:
%
% INPUT:
%   di_cfg      - validated analysis configuration struct
%   results     - results structure from di_analyze():
%                 .observed 
%                   .statVal            (1 x pX) observed statistics
%                   .pVal               (1 x pX) p-values (parametric analyses only)
%                   .clusters           cluster info (cluster-based analyses only)
%                     .clusterMembership_obs
%                     .clustIDList_obs
%                     .metrics_obs
%                 .simulated
%                   .permutationH0  
%                       .statVal        (nIterations x pX) null distribution
%                       .clusterMetrics (1 x nIterations) struct array
%                       .modes          (for PLS_SVD) mode-level null metrics
%                   .bootstrapStability
%                       .values         (nIterations x pX) bootstrap samples
%                       .loadings       (for PLS_SVD) bootstrap loading samples
%                 .PLS_SVD            (for PLS_SVD analyses only)
%                   .loadings.U_obs, .V_obs, .YU_obs, .XV_obs
%                   .modes.s_obs, .p_obs, .r_obs, .nModes
%
% OUTPUT:
%   results     - same as input but augmented with inference fields:
%       .inference
%           for permutation, empiricalFeature_inferenceFeature
%                   .feature.pVal_emp           - empirical pvalue per feature
%                   .feature.pVal_emp_FDR       - empirical pvalue per feature (FDR corrected)
%           for permutation, parametricFeature_inferenceFeature
%                   .feature.pVal_maxT          - pvalue per feature (maxT corrected)
%                   .feature.thresholds.*       - maxT thresholds per metric
%           for permutation, parametricFeature_inferenceCluster
%                   .cluster.pVal_maxT          - pvalue per cluster (maxT corrected)
%                   .cluster.thresholds.*       - maxT thresholds per metric
%           for permutation, PLS_SVD
%                   .mode.pVal_maxT             - pvalue per mode (maxT corrected)
%                   .mode.thresholds.*          - maxT thresholds per metric
%           for permutation, AJIVE
%                   .mode       not yet implemented
%           for bootstrap, empiricalFeature_inferenceFeature
%                   .feature.BR_rob             - bootstrap ratio (robust)
%                   .feature.CIlo               - lower 95% confidence interval
%                   .feature.CIup               - upper 95% confidence interval
%           for bootstrap, PLS_SVD
%                   stored in .simulated.bootstrapStability.loadings.inference:
%                   .U_sd, .V_sd                - standard deviation of loadings
%                   .U_ci, .V_ci                - 95% confidence intervals [lower, upper]
%                   .U_aligned, .V_aligned      - Procrustes-aligned bootstrap samples
%           for bootstrap, AJIVE
%                   .mode       not yet implemented
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)


%% shortcuts

% p value for inference
p_crit = di_cfg.analysis.p_crit;

% num of iterations
nIterations = di_cfg.analysis.nIterations;

% extract dimensions in d# order (fieldnames order)
dimKeys  = fieldnames(di_cfg.dimensions);
dimTypes = cellfun(@(k) di_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);
dimSizes = cellfun(@(k) length(di_cfg.dimensions.(k).vec), dimKeys);
nDims    = numel(dimKeys);
pX       = prod(dimSizes);  % total number of X features across all dimensions

% ignore col
ignore_col = di_cfg.analysis.ignore_col;

%% sanity checks

% only one spherical dimenion is allowed
if nnz(strcmp(dimTypes, 'spherical')) > 1
    error('Only one spherical dimension is supported');
end

% ignore_col and dimensions must have same numerosity
if numel(ignore_col) ~= pX
    error('ignore_col must have one entry per feature across all dimensions');
end

% validate that results structure exists and has required fields
if ~isstruct(results)
    error('results must be a struct');
end
if ~isfield(results, 'observed')
    error('results must have .observed field');
end
if ~isfield(results.observed, 'statVal')
    error('results.observed must have .statVal field');
end


% check that di_cfg has objective and analysis fields
if ~isfield(di_cfg, 'analysis')
    error('di_cfg must have .analysis field');
end
if ~isfield(di_cfg.analysis, 'objective')
    error('di_cfg.analysis must have .objective field (permutationH0testing or bootstrapStability)');
end
if ~isfield(di_cfg.analysis, 'type')
    error('di_cfg.analysis must have .type field');
end


%% implementation

% build key using analysis objective and analysis type 
key = sprintf('%s + %s', di_cfg.analysis.objective, di_cfg.analysis.type);

switch key
    
    case 'permutationH0testing + empiricalFeature_inferenceFeature'
        % - per each feature, compares observed statistical scores to a null distribution
        % - derive empirical p-values 
        % - applies FDR correction across features belonging to chosen dimensions (eg, it could be all, some, or none)
        % - optional, form clusters for descriptive purposes (done elsewhere, downstream function)

        % get observed statistics and null distribution from unified structure
        statVal_obs = results.observed.statVal;  % [1 x pX]
        null_values = results.simulated.permutationH0.statVal;  % [nIterations x pX]
        
        % Compute empirical p-values: fraction of null values as extreme as observed
        % Obs: minimum p-value granularity is 1/nIterations. So, for p < 0.05, need nIterations >= 20; for p < 0.01, need nIterations >= 100
        if nIterations < 100
            warning(['nIterations = ' num2str(nIterations) ' is small; p-value resolution = 1/' num2str(nIterations)]);
        end
        

        pval_emp = nan(1, size(statVal_obs, 2)); % initialize
        for colIdx = 1:size(statVal_obs, 2)
            if di_cfg.analysis.ignore_col(colIdx)
                continue
            end
            vals_perm = null_values(:, colIdx);
            val_obs = statVal_obs(colIdx);
            % p-value (two tailed): proportion of |null values| >= |observed|
            pval_emp(1, colIdx) = sum(abs(vals_perm) >= abs(val_obs)) / nIterations;
        end
        
        % --- start FDR correction --- 
        % TO DO: move to its own function
        % TO DO: make my own fdr correction (Benjamini Hochberg) without Bioinfomatics Toolbox's mafdr

        % Description: restructure p values (and ignore_col) so that p values are in columns 
        % where each column is an FDR pool. For exmample, if FDR_dimensions are all 1s, there is just one big pool
        % the p values will be reshaped to the original after the FDR
        % correction.

        % initialize
        pval_emp_FDR = nan(1, size(statVal_obs, 2));

        % FDR pooling selection across dimensions (logical vector over all dims)
        % shortcut
        poolDims = di_cfg.analysis.FDR_dimensions;
        % sanity check: num of FDR instructions match num of dimensions
        if numel(poolDims) ~= nDims
            error('analysis.FDR_dimensions must match the number of dimensions');
        end
        % find idx of dimensions to FDR correct and not correct
        FDRdim_idx = find(poolDims); % dimensions to pool for the correction
        otherdim_idx = setdiff(1:nDims, FDRdim_idx, 'stable'); % the other dimensions

        % permute to bring upfront the dimensions over which pvalues will be pooled
        perm  = [FDRdim_idx, otherdim_idx];
        pval_emp_permuted = permute(reshape(pval_emp, dimSizes'), perm);
        R_ignore_permuted = permute(reshape(ignore_col, dimSizes'), perm); % same thing to ignore_col

        % reshape to matrix to have columns of p values upon which the correction is applied (one iteration per column)
        sz = size(pval_emp_permuted);
        L  = prod(sz(1:numel(FDRdim_idx)));
        M  = prod(sz(numel(FDRdim_idx)+1:end));
        pval_emp_matrix = reshape(pval_emp_permuted, L, M);
        R_ignore_matrix = reshape(R_ignore_permuted, L, M); % same thing to ignore_col

        % apply FDR-BH along each column
        pval_emp_FDR_matrix = nan(size(pval_emp_matrix));
        for colIdx = 1:M
            pvalVec    = pval_emp_matrix(:,colIdx);
            RignoreVec = R_ignore_matrix(:,colIdx);
            pvalVec2use = pvalVec(~RignoreVec); % remove the p values corresponding with features that are not of interest in this analysis
            % if pvalVec2use is empty (e.g., because this feature is totally ignored), move on
            if isempty(pvalVec2use)
                continue
            end
            pvalVec2use_FDR = mafdr(pvalVec2use , 'BHFDR', true);  % FDR BH correction
            pval_emp_FDR_matrix(~RignoreVec,colIdx) = pvalVec2use_FDR; % collocate pvalues in the longer matrix where they belong
        end

        % reshape back
        pval_emp_FDR_permuted = reshape(pval_emp_FDR_matrix, sz);
        pval_emp_FDR = reshape(ipermute(pval_emp_FDR_permuted, perm),[1 pX]);

        % --- end of FDR correction --

        % figure
        if di_cfg.analysis.figFlag
            figure(); clf
            f = gcf; f.Units = 'normalized'; f.Position = [0.2    0    0.4    0.9];
            vertJitterVec = randn(1,length(pval_emp(~di_cfg.analysis.ignore_col)));
            lp1 = semilogx(pval_emp(~di_cfg.analysis.ignore_col),vertJitterVec.*ones(1,length(pval_emp(~di_cfg.analysis.ignore_col))),'o');
            lp1.Parent.XLim = [0 1];
            lp1.Parent.XTick = [0 0.01 0.05 0.1 0.2 0.5 1];
            lp1.Parent.XAxis.Label.String = 'p value';
            lp1.Parent.XMinorTick = 'off';
            lp1.Parent.YAxis.Visible = 'off';
            hold on
            lp2 = semilogx(pval_emp_FDR(~di_cfg.analysis.ignore_col),vertJitterVec.*ones(1,length(pval_emp(~di_cfg.analysis.ignore_col))),'x');
            ln = xline(di_cfg.analysis.p_crit);
            legend([lp1 lp2 ln],["uncorrected" "FDR corrected" "p_{crit}"],'Location','southoutside')
            title('pvalues before vs after FDR correction')
        end

        % Store inference results in unified structure
        results.inference.feature.pVal_emp     = pval_emp;
        results.inference.feature.pVal_emp_FDR = pval_emp_FDR;

    case 'permutationH0testing + empiricalFeature_inferenceCluster'
        % not implemented, not difficult to code, but difficult for a
        % computer to run: it requires nIterations per each iteration, so nIterations^2...
        error('not coded. probably I will never code this one')

    case 'permutationH0testing + parametricFeature_inferenceFeature'
        % - per feature: compare observed statistic against null distribution (permutation)
        % - apply maxT correction on feature-level statistics
        % Null distribution must be provided upstream as results.simulated.permutationH0.statVal
        results = di_maxT(di_cfg, results);

    case 'permutationH0testing + parametricFeature_inferenceCluster'
        % - (per each feature, use parametric statistics to get a p-value) (done upstream)
        % - (form clusters and compute cluster metrics for observed and simulated data) (done upstream)
        % - compares each observed cluster metric against null distribution
        % - (applies maxT correction on cluster-level statistics) (done downstream)
      
        results = di_maxT(di_cfg, results);

        %keyboard; % check this analysis works after the changes

    case 'permutationH0testing + PLS_SVD'
        % - per each mode, compares observed singular values against its null distribution
        % two types of correction (on singular values, inertia, correlation): 
        % - sequential correction
        % - maxT correction on mode-level statistics
        %    OBS    - Not yet fully implemented
        warning('not yet done');
        % TODO: Implement mode-wise maxT inference similar to di_maxT.m but for modes
    
    case 'permutationH0testing + AJIVE'
        error('not coded yet')

    case 'bootstrapStability + empiricalFeature_inferenceFeature'
        % - per each feature, use resampled distribution to compute statistical scores stability
        % - optional, form clusters for descriptive purposes (done elsewhere, downstream)

        keyboard; % check this works after the changes

        % Get observed statistics and bootstrap distribution from unified structure
        statVal_obs = results.observed.statVal;  % [1 x pX]
        boot_values = results.simulated.bootstrapStability.values;  % [nIterations x pX]
        
        % Compute bias-corrected robust estimates and confidence intervals
        % Bootstrap Robust (BR): ratio of bootstrap median to observed
        BR_rob = nan(1, size(statVal_obs, 2));
        CIlo = nan(1, size(statVal_obs, 2));
        CIup = nan(1, size(statVal_obs, 2));
        
        for colIdx = 1:size(statVal_obs, 2)
            boot_dist = boot_values(:, colIdx);
            val_obs = statVal_obs(colIdx);
            
            % Robust estimate: bootstrap median relative to observed
            boot_median = median(boot_dist);
            % Safeguard: if observed value is very close to zero, ratio becomes unstable
            if abs(val_obs) > 1e-10  % threshold for "effectively nonzero"
                BR_rob(1, colIdx) = boot_median / val_obs;  % ratio: if ~1, stable under resampling
            else
                % For near-zero observed values, use difference instead
                BR_rob(1, colIdx) = boot_median - val_obs;
                warning(['Feature ' num2str(colIdx) ' has near-zero observed value; using difference for BR_rob']);
            end
            
            % 95% confidence interval from bootstrap distribution
            CIlo(1, colIdx) = quantile(boot_dist, 0.025);
            CIup(1, colIdx) = quantile(boot_dist, 0.975);
        end
        
        % Store inference results in unified structure
        results.inference.feature.BR_rob = BR_rob;
        results.inference.feature.CIlo = CIlo;
        results.inference.feature.CIup = CIup;

    case 'bootstrapStability + empiricalFeature_inferenceCluster'
        % not implemented // unsure how to implement this
        error('not coded. not gonna happen any time soon')

    case 'bootstrapStability + parametricFeature_inferenceFeature'
        % not implemented because not meaningful (or I could use both observed values and resampling distribution) 
        error('not coded. it might happen but low priority')

    case 'bootstrapStability + parametricFeature_inferenceCluster'
        % - not implemented // unsure how to implement this
        error('not coded yet. it would be nice but low priority');

    case 'bootstrapStability + PLS_SVD'
        %       - per each mode, uses resampled distribution to compute singular
%       vectors stability
%       - optional, form clusters for descriptive purposes
%    OBS    - Not yet fully implemented
% 
        keyboard; % check it works after the changes

% PLS-SVD + BOOTSTRAP: Assess stability of PLS loadings via bootstrap resampling
        %
        % STATISTICAL RATIONALE:
        %   Bootstrap resampling estimates sampling variability of PLS loadings.
        %   Narrow confidence intervals → stable loadings (reliable interpretation)
        %   Wide confidence intervals → unstable loadings (interpret with caution)
        %
        % WORKFLOW:
        %   1. Extract observed loadings and bootstrap samples
        %   2. Align bootstrap samples using Procrustes (remove rotational ambiguity)
        %   3. Compute stability metrics (SD, confidence intervals)
        %
        % OUTPUT:
        %   .inference.U_sd, .V_sd        - Standard deviation of loadings
        %   .inference.U_ci, .V_ci        - 95% confidence intervals
        %   .inference.U_aligned, .V_aligned - Aligned bootstrap samples
        
        % STEP 1: Extract observed loadings (reference for alignment)
        U_obs = results.PLS_SVD.loadings.U_obs;      % (pY x nModes) observed Y loadings
        V_obs = results.PLS_SVD.loadings.V_obs;      % (pX x nModes) observed X loadings
        
        % STEP 2: Extract bootstrap loading samples (from di_analysis_plsSVD.m)
        % These already have sign convention applied (aligned to observed via di_signConvention)
        % but may still have rotational ambiguity in multi-mode subspace
        U_boot = results.simulated.bootstrapStability.loadings.U_boot;   % (pY x nModes x nIterations)
        V_boot = results.simulated.bootstrapStability.loadings.V_boot;   % (pX x nModes x nIterations)
        nIterations = di_cfg.analysis.nIterations;
        
        % STEP 3: Procrustes alignment
        % Remove rotational ambiguity by finding optimal orthogonal rotation
        % that aligns each bootstrap sample to observed loadings.
        % Why needed: SVD has two ambiguities:
        %   - Sign ambiguity (±): handled by di_signConvention in di_analysis_plsSVD.m
        %   - Rotation (multi-mode): handled here via Procrustes
        % Example: With 2 modes, bootstrap samples can rotate in 2D plane
        %          even after sign correction. Procrustes finds the rotation
        %          that best matches observed orientation.
        [U_aligned, V_aligned] = di_procrustesAlign(U_obs, V_obs, U_boot, V_boot, nIterations);
        
        % STEP 4: Compute loading variability metrics
        % Standard deviation quantifies element-wise variability across bootstrap samples
        % Dimension 3 is the iteration dimension (nIterations)
        results.simulated.bootstrapStability.loadings.inference.U_sd = std(U_aligned, 0, 3);  % (pY x nModes)
        results.simulated.bootstrapStability.loadings.inference.V_sd = std(V_aligned, 0, 3);  % (pX x nModes)
        
        % STEP 5: Compute bootstrap confidence intervals (percentile method)
        % 95% CI = [2.5th percentile, 97.5th percentile] across bootstrap distribution
        % Output dimensions: (pY x nModes x 2) for U_ci, (pX x nModes x 2) for V_ci
        % where (:,:,1) is lower bound and (:,:,2) is upper bound
        results.simulated.bootstrapStability.loadings.inference.U_ci = prctile(U_aligned, [2.5 97.5], 3);
        results.simulated.bootstrapStability.loadings.inference.V_ci = prctile(V_aligned, [2.5 97.5], 3);
        
        % STEP 6: Store aligned loadings for diagnostic purposes
        % Users may want to inspect the full aligned bootstrap distribution
        % (e.g., check for multimodality, non-normality, outliers)
        results.simulated.bootstrapStability.loadings.inference.U_aligned = U_aligned;  % (pY x nModes x nIterations)
        results.simulated.bootstrapStability.loadings.inference.V_aligned = V_aligned;  % (pX x nModes x nIterations)
    
    case 'bootstrapStability + AJIVE'
        error('not yet coded')

end

%% output validation

% check that inference has been applied (at least one new field added)
inference_added = false; % initialize to false

% check for empirical-specific fields
if isfield(results, 'inference') && isfield(results.inference, 'feature') && isfield(results.inference.feature, 'pVal_emp_FDR')
    if ~all(isnan(results.inference.feature.pVal_emp_FDR)) || ~isempty(results.inference.feature.pVal_emp_FDR)
        inference_added = true;
    end
end

% check for bootstrap-specific fields
if isfield(results, 'inference') && isfield(results.inference, 'feature')
    if isfield(results.inference.feature, 'BR_rob')
        inference_added = true;
    end
end

% check for cluster-specific fields
if isfield(results, 'inference') && isfield(results.inference, 'cluster')
    inference_added = true;
end

% warn if no inference fields were detected... it might indicate unimplemented path
if ~inference_added
    switch [di_cfg.analysis.objective ' & ' di_cfg.analysis.type]
        case {'permutationH0testing & empiricalFeature_inferenceFeature', 'permutationH0testing & parametricFeature_inferenceCluster'}
            warning('No inference fields detected in results - check that inference was actually computed');
    end
end

end
