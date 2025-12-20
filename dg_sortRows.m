function [L, R] = dg_sortRows(dg_cfg,L_orig,R_orig,rowIdx,itIdx)

% INPUT:
%   L           - 
%
%   R           - 
%
%   dg_cfg
%
% OUTPUT:
%
%
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)


%% sort rows as appropriate
switch dg_cfg.objective
    case 'permutationH0testing'
        % either permute L or R
        % it's faster to permute L's rows (fewer columns) but with
        % rmFactors we need to permute R's rows
        switch num2str(dg_cfg.designCode)
            case num2str([0 0 1])
                L(:,:) = L_orig;
                R(:,:) = R_orig(rowIdx(:,itIdx),:); % permute rows of R
            otherwise
                L(:,:) = L_orig(rowIdx(:,itIdx),:); % permute rows of L
                R(:,:) = R_orig;
        end

        
    case 'bootstrapStability'
        L(:,:) = L_orig(rowIdx(:,itIdx),:); % permute rows of L
        R(:,:) = R_orig(rowIdx(:,itIdx),:); % permute rows of R
end
