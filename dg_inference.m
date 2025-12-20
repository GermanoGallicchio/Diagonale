function results = dg_inference(dg_cfg,results)

% Descr
% permutation (hypothesis-testing inference) => pvalues (1 tail)
% bootstrap (estimation-based inference) => BR, 2-tail 95%CI (2 tails)
%
% empiricalL1_FDR
%   permutation - feature-wise inference // multiple-comparison corrected: FDR over selected dimensions
%   bootstrap   - feature-wise inference // multiple-comparison uncorrected 
%
% cluster
%   permutation - cluster-level inference // multiple-comparison corrected: maxT
%   bootstrap   - cluster-level inference (mass or size) // multiple-comparison uncorrected (but fewer comparisons than feature-wise) 
%
% pls_svd
%   permutation - mode-level inference // multiple-comparison corrected: maxT and perMode
%   bootstrap   - loading-level inference // multiple comparison uncorrected
%
% INPUT:
%
%   dg_cfg
%
%
% OUTPUT:
%
% rowIdx
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)


%% shortcuts

p_crit = dg_cfg.p_crit;
nIterations = dg_cfg.nIterations;

ny1 = dg_cfg.dimensions.y1_num;
nx2 = dg_cfg.dimensions.x2_num;
nz3 = dg_cfg.dimensions.z3_num;

R_ignore = dg_cfg.R_ignore;

%% implementation

switch [dg_cfg.objective ' & ' dg_cfg.analysis]
    case 'permutationH0testing & empiricalL1_FDR'

        % find how "extreme" each statVal is compared with its distribution of statVals under H0
        pval_emp = nan(1,size(results.statVal_obs,2)); % initialize
        for colIdx = 1:size(results.statVal_obs,2)
            if dg_cfg.R_ignore(colIdx)
                continue
            end
            vals_perm = results.resampling.statVal_resamp(:,colIdx);
            val_obs = results.statVal_obs(colIdx);
            pval_emp(1,colIdx) = sum(abs(vals_perm)>=abs(val_obs)) / nIterations;
        end

        % --- TO DO: put in its own function ---
        % FDR correction
        pval_emp_FDR = nan(1,size(results.statVal_obs,2)); % initialize
        FDRdim_idx = find(dg_cfg.FDR_params.dimensions); % dimensions to pool for the correction
        otherdim_idx = setdiff(1:3, FDRdim_idx, 'stable'); % the other dimensions

        % permute to bring upfront the dimensions over which pvalues will be pooled
        perm  = [FDRdim_idx, otherdim_idx];
        pval_emp_permuted = permute(reshape(pval_emp,[ny1 nx2 nz3]),perm);
        R_ignore_permuted = permute(reshape(R_ignore,[ny1 nx2 nz3]),perm);

        % reshape to matrix to have columns of p values upon which the correction is applied (one iteration per column)
        sz = size(pval_emp_permuted);
        L  = prod(sz(1:numel(FDRdim_idx)));
        M  = prod(sz(numel(FDRdim_idx)+1:end));
        pval_emp_matrix = reshape(pval_emp_permuted, L, M);
        R_ignore_matrix = reshape(R_ignore_permuted, L, M);

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
        pval_emp_FDR = reshape(ipermute(pval_emp_FDR_permuted, perm),[1 ny1*nx2*nz3]);

        % figure
        if dg_cfg.figFlag
            figure(); clf
            f = gcf; f.Units = 'normalized'; f.Position = [0.2    0    0.4    0.9];
            vertJitterVec = randn(1,length(pval_emp(~dg_cfg.R_ignore)));
            lp1 = semilogx(pval_emp(~dg_cfg.R_ignore),vertJitterVec.*ones(1,length(pval_emp(~dg_cfg.R_ignore))),'o');
            lp1.Parent.XLim = [0 1];
            lp1.Parent.XTick = [0 0.01 0.05 0.1 0.2 0.5 1];
            lp1.Parent.XAxis.Label.String = 'p value';
            lp1.Parent.XMinorTick = 'off';
            lp1.Parent.YAxis.Visible = 'off';
            hold on
            lp2 = semilogx(pval_emp_FDR(~dg_cfg.R_ignore),vertJitterVec.*ones(1,length(pval_emp(~dg_cfg.R_ignore))),'x');
            ln = xline(dg_cfg.p_crit);
            legend([lp1 lp2 ln],["uncorrected" "FDR corrected" "p_{crit}"],'Location','southoutside')
            title('pvalues before vs after FDR correction')
        end

        % add to results
        results.resampling.pVal_emp     = pval_emp;
        results.resampling.pVal_emp_FDR = pval_emp_FDR;
        % --- ---

    case 'permutationH0testing & theoreticalL1_clusterMaxT'
        results = dg_maxT(dg_cfg,results); % maxT

    case 'permutationH0testing & PLS_SVD'
        % not implemented here
end

end
