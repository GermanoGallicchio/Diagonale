function [U_aligned, V_aligned] = di_procrustesRotate(U_obs, V_obs, U_boot, V_boot)
% Align bootstrap loadings to observed loadings using Procrustes rotation.
%
% Applies orthogonal rotation (Procrustes transformation) to align bootstrap loading
% matrices to observed loadings, removing rotational ambiguity in the
% subspace spanned by multiple modes.
%
% INPUT:
%   U_obs       - (pX x nModes) observed X-side loadings (reference)
%   V_obs       - (pY x nModes) observed Y-side loadings (reference)
%   U_boot      - (pX x nModes) bootstrap X-side loadings
%   V_boot      - (pY x nModes) bootstrap Y-side loadings
%
% OUTPUT:
%   U_aligned   - (pX x nModes) Procrustes-aligned X-side loadings
%   V_aligned   - (pY x nModes) Procrustes-aligned Y-side loadings
%
% ALGORITHM:
%     1. Build alignment matrix M from both blocks at once, between the
%     observed (reference) and bootstrapped (to align) loadings
%       M = U_boot' * U_obs  +  V_boot' * V_obs
%     2. Compute SVD of M = L * S * R' to get orthogonal singular vectors
%     3. Compute transformation matrix as the product T of L * R' (orthogonal rotation: Procrustes)
%     4. Apply transformation to each bootstrapped matrix:
%       U_aligned = U_boot * T
%       V_aligned = V_boot * T
%
% WHY THIS WORKS: T = L*R' is the closed-form solution of the orthogonal
% Procrustes problem -- among all orthogonal matrices, it is the one that
% minimises the misalignment ||[U_boot;V_boot]*T - [U_obs;V_obs]||
% (Frobenius norm), i.e. the best possible re-orientation of the bootstrap
% mode subspace onto the observed one.
%
% why both blocks: the same T is applied to U and V, therefore it must be
% the one that minimises the misalignment of both at the same time.
% Computing it from a single block would be a big problem. The number of
% modes is min(pX,pY), so the loading matrix of the smaller block is always
% square, and a square matrix of orthonormal columns can be rotated exactly
% onto another one. A T obtained from that block alone would map every
% bootstrap iteration exactly onto the observed loadings (i.e., zero
% variation across iterations), therefore the resulting stability estimates
% would describe this alignment step and not the data. Stacking the two
% blocks avoids this, because the stack has pX+pY rows for only min(pX,pY)
% columns and is therefore never square.
%
% note on reflections: T can be a reflection as well as a rotation.
% Resampling can flip the sign of a mode, and correcting those flips is part
% of what this alignment is for, therefore constraining T to rotations only
% would leave some of them uncorrected.
%
% note on other implementations: McIntosh's PLS toolbox (pls_analysis.m,
% through rri_bootprocrust) computes the rotation from one block and applies
% it to the other, both for permutation and for bootstrap. This is not a
% problem there, because the ratios it reports are for the brain block,
% which is large and keeps varying across iterations no matter which
% rotation is applied, and because it resamples saliences already scaled by
% the singular values rather than unit singular vectors. Here either block
% can be the
% smaller one and the loadings are unit singular vectors, therefore the
% rotation is computed from both. A later variant (neuropm-lab/svd_pls) does
% use both blocks, but averages a rotation obtained from each: the average
% of two orthogonal matrices is not orthogonal, therefore what it applies is
% not a rotation anymore. Summing the cross-products before the SVD, as
% below, keeps T orthogonal.
%
% References:
%   Schonemann, P. H. (1966). A generalized solution of the orthogonal
%     Procrustes problem. Psychometrika, 31(1), 1-10.
%   Milan, L., & Whittaker, J. (1995). Application of the parametric
%     bootstrap to models that incorporate a singular value decomposition.
%     Applied Statistics, 44(1), 31-49.
%   Krishnan, A., Williams, L. J., McIntosh, A. R., & Abdi, H. (2011).
%     Partial Least Squares (PLS) methods for neuroimaging: a tutorial and
%     review. NeuroImage, 56(2), 455-475.
%
% TO DO: the name procrustes is terrible and I would prefer a different
% one based on geometry... maybe something like orthogonal rotation
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% SANITY CHECKS

[pX, nModes] = size(U_obs);
pY = size(V_obs, 1);

if size(V_obs, 2) ~= nModes
    error('U_obs and V_obs must have same number of modes');
end
if ~isequal(size(U_boot), size(U_obs))
    error('U_boot dimensions must match U_obs');
end
if ~isequal(size(V_boot), size(V_obs))
    error('V_boot dimensions must match V_obs');
end

%% implementation

% 1. cross-product matrix of observed and bootstrapped loadings, from both
% blocks. summing the two cross-products is the same as stacking the blocks
% (i.e., [U; V]) and taking one cross-product of the stack, which gives the
% single rotation that minimises the misalignment of X-side and Y-side
% loadings at the same time
M = U_boot' * U_obs  +  V_boot' * V_obs;

% 2. singular vectors via SVD
[L, ~, R] = svd(M, 'econ');

% 3. transformation matrix
T = L * R';

% 4. apply transformation to X-side and Y-side bootstrapped loadings
U_aligned = U_boot * T;
V_aligned = V_boot * T;
% note: identical transformation matrix for both U and V to preserve the
% geometric relationship between X and Y loading spaces. After this transformation, both U_aligned and V_aligned
% are now in the same orientation as the observed loadings.

end
