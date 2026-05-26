function [U, V] = di_signAlign(U, V, U_ref, V_ref)
% di_signAlign - Align signs of PLS loadings for reproducibility across bootstrap iterations
%
% Anchors signs to reference loadings via dot-product alignment to ensure
% consistent sign convention across resampling iterations.
%
% Usage:
%   [U, V] = di_signAlign(U, V, U_ref, V_ref)
%
% Inputs:
%   U, V           - Current PLS loadings (nObs × nModes for U, nFeatures × nModes for V)
%   U_ref, V_ref   - Reference loadings (e.g., from itIdx 1)
%
% Outputs:
%   U, V           - Aligned loadings
%
% Details:
%   Each mode is checked: if dot(U_iter, U_ref) + dot(V_iter, V_ref) < 0,
%   both U and V for that mode are flipped to align with reference.
%   This ensures reproducible sign convention based on the chosen reference.

nModes = size(U, 2);
for modeIdx = 1:nModes
    % Align score combines U and V mode alignments
    alignScore = dot(U(:, modeIdx), U_ref(:, modeIdx)) + ...
                 dot(V(:, modeIdx), V_ref(:, modeIdx));
    % If score is negative, the mode is opposed; flip both U and V for that mode
    if alignScore < 0
        U(:, modeIdx) = -U(:, modeIdx);
        V(:, modeIdx) = -V(:, modeIdx);
    end
end
end
