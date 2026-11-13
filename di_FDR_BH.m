function pAdj = di_FDR_BH(p)
% Benjamini-Hochberg FDR correction
%
% INPUT:
%   p     - vector of raw p-values (row or column)
%
% OUTPUT:
%   pAdj  - vector of BH-adjusted p-values, same size as p
%
% Notes:
%   - NaN values are preserved as NaN.
%   - Valid p-values are expected in [0,1].
%
% note: alternative to the mafdr in Bioinformatics toolbox but with a much smaller scope (this only does BH procedure)
%
% RATIONALE: BH is a STEP-UP procedure on the ordered p-values. The false
% discovery rate (the expected proportion of false positives among the
% features declared significant) is a less strict criterion than the
% family-wise error rate: instead of guarding against ANY false positive,
% it tolerates a controlled fraction of them, buying power in large
% families. Each sorted p(i) is scaled by m/i (the BH line rising with
% rank); the adjusted value reported here is the smallest FDR level at
% which that feature would be declared significant. The cummin step below
% enforces monotonicity of the adjusted p-values (a feature can never be
% "more significant" than one with a smaller raw p).
%
% Reference: Benjamini, Y., & Hochberg, Y. (1995). Controlling the false
% discovery rate: a practical and powerful approach to multiple testing.
% Journal of the Royal Statistical Society B, 57(1), 289-300.
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

if ~isvector(p)
    error('\\ di_FDR_BH: InputNotVector Input p must be a vector.');
end

pAdj = nan(size(p));

% work on non-NaN values only
validMask = ~isnan(p);
if ~any(validMask)
    return
end

pValid = p(validMask);

if any(pValid < 0 | pValid > 1)
    error('\\ di_FDR_BH: InputOutOfRange Valid p-values must be in the interval [0,1].');
end

m = numel(pValid);
[pSorted, sortIdx] = sort(pValid, 'ascend');
ranks = (1:m)';

% BH raw adjusted values
pAdjSorted = pSorted(:) .* (m ./ ranks);

% enforce monotonicity from high rank to low rank
pAdjSorted = flipud(cummin(flipud(pAdjSorted)));

% cap at 1
pAdjSorted = min(pAdjSorted, 1);

% map back to original (valid-only) order
pAdjValid = nan(size(pValid));
pAdjValid(sortIdx) = pAdjSorted;
pAdj(validMask) = pAdjValid;

end
