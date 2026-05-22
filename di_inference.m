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
%                   .U_BR, .V_BR                - bootstrap ratio (observed / bootstrap SD)
%                   .U_95CIlo, .U_95CIup        - lower/upper 95% confidence interval for U loadings
%                   .V_95CIlo, .V_95CIup        - lower/upper 95% confidence interval for V loadings
%                   .U_ci, .V_ci                - 95% confidence intervals [lower, upper] (stacked)
%                   .U_aligned, .V_aligned      - Procrustes-aligned bootstrap samples
%           for bootstrap, AJIVE
%                   .mode       not yet implemented
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% shortcuts

% p value for inference (only used by permutation testing paths; not set for bootstrapStability)
if isfield(di_cfg.analysis, 'p_crit')
    p_crit = di_cfg.analysis.p_crit;
else
    p_crit = [];  % not needed for bootstrap; will error naturally if mistakenly used
end

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
    error('\\ Only one spherical dimension is supported');
end

% ignore_col and dimensions must have same numerosity
if numel(ignore_col) ~= pX
    error('\\ ignore_col must have one entry per feature across all dimensions');
end

% validate that results structure exists and has required fields
if ~isstruct(results)
    error('\\ results must be a struct');
end
if ~isfield(results, 'observed')
    error('\\ results must have .observed field');
end
if ~isfield(results.observed, 'statVal')
    error('\\ results.observed must have .statVal field');
end


% check that di_cfg has objective and analysis fields
if ~isfield(di_cfg, 'analysis')
    error('\\ di_cfg must have .analysis field');
end
if ~isfield(di_cfg.analysis, 'objective')
    error('\\ di_cfg.analysis must have .objective field (permutationH0testing or bootstrapStability)');
end
if ~isfield(di_cfg.analysis, 'type')
    error('\\ di_cfg.analysis must have .type field');
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
        
        % FDR correction 
        pval_emp_FDR = di_fdrPoolCorrect(pval_emp, ignore_col, dimSizes, di_cfg.analysis.FDR_dimensions);

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
        error('\\ not coded. probably I will never code this one')

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
        % Mode-level permutation inference of SV metrics with maxT correction
        % - Per each mode, compares observed SV-based metrics against their null distributions
        % - Applies maxT correction across modes for multiple metrics:
        %   * singular values (s): covariance captured by each mode
        %   * total inertia: sum of all singular values
        %   * Wilks lambda: cumulative variance (from end)
        %   * sequential variance: cumulative variance (from start)
        % - maxT controls FWER across modes within each metric
        
        % Compute uncorrected empirical p-values per mode
        % Each mode compared against its own null distribution
        results = di_empiricalP(di_cfg, results);
        
        % Compute maxT-corrected p-values (FWER control)
        % Each mode compared against max across all modes per iteration
        results = di_maxT(di_cfg, results);
    
    case 'permutationH0testing + AJIVE'
        error('\\ not coded yet')

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

% PLS-SVD + BOOTSTRAP: Assess stability of PLS loadings via bootstrap resampling
        %
        % STATISTICAL RATIONALE:
        %   Bootstrap resampling estimates sampling variability of PLS loadings.
        %   stable loadings (reliable interpretation) have greater
        %   bootstrap ratios and narrower confidence intervals

        % STEP 1: Extract observed loadings (reference for alignment)
        U_obs = results.PLS_SVD.loadings.U_obs;      % (pY x nModes) observed Y loadings
        V_obs = results.PLS_SVD.loadings.V_obs;      % (pX x nModes) observed X loadings

        nIterations = di_cfg.analysis.nIterations;
        nModes = size(U_obs,2);
        
        % STEP 2: Extract bootstrapped loadings (these are already
        % sign-aligned and rotated)
        U_boot = results.simulated.bootstrapStability.loadings.U_boot;   % (pY x nModes x nIterations)
        V_boot = results.simulated.bootstrapStability.loadings.V_boot;   % (pX x nModes x nIterations)
        
        % STEP 3: Compute loading variability metrics
        % Bootstrap ratio (observed / bootstrap SD), conventional
        % set tolerance
        eps0 = 1e-10;
        % initialize
        % TO DO: I might set these to zero, so if below toleratnce, they
        % get zero rather than NaN
        U_BR = nan(size(U_obs)); 
        V_BR = nan(size(V_obs));
        % compute SD across iterations
        U_sd = std(U_boot, 0, 3);  % (pY x nModes)
        V_sd = std(V_boot, 0, 3);  % (pX x nModes)
        % find idx across iterations of loadings with SD above tolerance
        % (meaningful variation across iterations) and to avoid diving by a
        % tiny number (potentially zero)
        U_valid = abs(U_sd) > eps0;
        V_valid = abs(V_sd) > eps0;
        % compute BR only for valid loadings
        U_BR(U_valid) = U_obs(U_valid) ./ U_sd(U_valid);
        V_BR(V_valid) = V_obs(V_valid) ./ V_sd(V_valid);
        results.simulated.bootstrapStability.loadings.inference.U_BR = U_BR;
        results.simulated.bootstrapStability.loadings.inference.V_BR = V_BR;
        
        % STEP 4: Compute bootstrap confidence intervals (percentile method)
        % 95% CI = [2.5th percentile, 97.5th percentile] across bootstrap distribution
        U_95CIlo = prctile(U_boot, 2.5, 3);
        U_95CIup = prctile(U_boot, 97.5, 3);
        V_95CIlo = prctile(V_boot, 2.5, 3);
        V_95CIup = prctile(V_boot, 97.5, 3);
        results.simulated.bootstrapStability.loadings.inference.U_95CIlo = U_95CIlo;
        results.simulated.bootstrapStability.loadings.inference.U_95CIup = U_95CIup;
        results.simulated.bootstrapStability.loadings.inference.V_95CIlo = V_95CIlo;
        results.simulated.bootstrapStability.loadings.inference.V_95CIup = V_95CIup;
        
    case 'bootstrapStability + AJIVE'
        error('not yet coded')

end

%% output validation
% TO DO: think of something to include, but very low priority


end
