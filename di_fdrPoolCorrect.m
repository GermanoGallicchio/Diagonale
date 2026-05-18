function pval_emp_FDR = di_fdrPoolCorrect(pval_emp, ignore_col, dimSizes, poolDims)
% Apply FDR-BH correction within pools defined across selected dimensions.
%
% Description: restructure p values (and ignore_col) so that p values are in columns
% where each column is an FDR pool. For example, if FDR_dimensions are all 1s,
% there is just one big pool. The p values are reshaped back after correction.
%
% INPUT:
%   pval_emp    - 1 x pX vector of empirical p-values
%   ignore_col  - 1 x pX logical/numeric mask (true/1 means ignore feature)
%   dimSizes    - 1 x nDims vector with feature-grid sizes per dimension
%   poolDims    - 1 x nDims logical/numeric vector selecting pooled dimensions
%
% OUTPUT:
%   pval_emp_FDR - 1 x pX vector of BH-FDR corrected p-values
%
% notes: it previously used mafdr from the Bioinformatics toolbox.
% it now used a custom implementation of the same FDR-BH procedure via di_FDR_BH.m that is base MATLAB only.
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

% initialize
pval_emp_FDR = nan(1, numel(pval_emp));

% FDR pooling selection across dimensions (logical vector over all dims)
% sanity check: num of FDR instructions match num of dimensions
nDims = numel(dimSizes);
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
    pvalVec2use_FDR = di_FDR_BH(pvalVec2use);  % FDR BH correction (base MATLAB)
    pval_emp_FDR_matrix(~RignoreVec,colIdx) = pvalVec2use_FDR; % collocate pvalues in the longer matrix where they belong
end

% reshape back
pX = prod(dimSizes);
pval_emp_FDR_permuted = reshape(pval_emp_FDR_matrix, sz);
pval_emp_FDR = reshape(ipermute(pval_emp_FDR_permuted, perm),[1 pX]);

end
