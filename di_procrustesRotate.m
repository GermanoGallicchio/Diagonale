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
%     1. Build alignment matrix M = U_boot' * U_obs between the observed
%     (reference) and bootstrapped (to align) loadings
%     2. Compute SVD of M = L * S * R' to get orthogonal singular vectors
%     3. Compute transformation matrix as the product T of L * R' (orthogonal rotation: Procrustes)
%     4. Enforce a proper rotation (det(T)=+1), removing reflections
%     5. Apply transformation to each bootstrapped matrix: 
%       U_aligned = U_boot * T
%       V_aligned = V_boot * T
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

% initialize output arrays to store aligned loadings
U_aligned = zeros(size(U_boot));
V_aligned = zeros(size(V_boot));
    
% 1. cross-product matrix of observed and bootstrapped loading matrices
M = U_boot' * U_obs;

% 2. singular vectors via SVD
[L, ~, R] = svd(M, 'econ');    
    
% 3. transformation matrix    
T = L * R';

% 4. Enforce a proper rotation by excluding reflections
% det(T)=+1
if det(T) < 0
    L(:, end) = -L(:, end);
    T = L * R';
end
% NOTE: this step is needed to avoid improper rotations (reflections).
    
% 5. Apply transformation (rotation) to X-side and Y-side bootstrapped loadings 
U_aligned = U_boot * T;
V_aligned = V_boot * T;
% note: identical transformation matrix for both U and V to preserve the
% geometric relationship between X and Y loading spaces. After this transformation, both U_aligned and V_aligned
% are now in the same orientation as the observed loadings.

end
