function [Y, X] = di_reorderRowsApply(di_cfg,Y_orig,X_orig,rowIdx,itIdx)
% reorder rows of matrices according to indices (for iteration itIdx of the nIterations included in the cols of matrix rowIdx)
%
% Applies row reordering based on indices from di_reorderRowsGenerate. Rows
% are reordered in equal way across all columns respecting whatever
% structure the columns represent.
% This function is agnostic to how indices were generated. It simply uses them 
% to reorder rows.
% All reordering logic (permutation strategy, structure preservation, bootstrap
% behavior) is determined by di_reorderRows, not this function.
%
% INPUT:
%   di_cfg      - struct (used only for reference, not logic)
%   Y_orig      - Original response/design/data matrix (m x pY)
%   X_orig      - Original data matrix (m x pX)
%   rowIdx      - Row index matrix (m x nIterations) from di_reorderRows
%   itIdx       - Current iteration index (scalar, 1 to nIterations)
%
% OUTPUT:
%   Y           - Reordered Y matrix using rowIdx(:,itIdx)
%   X           - Reordered X matrix using rowIdx(:,itIdx)
%
% NOTES:
%   - This function performs simple row indexing: no decision logic
%   - Structure preservation, permutation constraints, and bootstrap logic
%     all originate from di_reorderRowsGenerate, not this function.
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% sort rows as appropriate
switch di_cfg.analysis.objective
    case 'permutationH0testing'
        % idea: break the link betweeen X and Y but permuting one but not
        % the other. Permute Y's rows (typically fewer columns than X so faster computation)
        Y(:,:) = Y_orig(rowIdx(:,itIdx),:);
        X(:,:) = X_orig;

        
    case 'bootstrapStability'
        % idea: sample whole observations with all their values from all
        % their matrices
        Y(:,:) = Y_orig(rowIdx(:,itIdx),:); % resample rows of Y
        X(:,:) = X_orig(rowIdx(:,itIdx),:); % resample rows of X
end

