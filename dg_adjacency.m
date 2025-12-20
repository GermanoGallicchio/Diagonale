function [adjMatrix] = dg_adjacency(dg_cfg, figFlag)
% Creates an adjacency vector for 3D physiological data
%
% generates a sparse adjacency vector describing neighborhood relationships
% across three vectorized dimensions: y1, x2, and z3.
% 
% INPUT: 
%       dg_cfg
%
% OUTPUT:
%       adjMatrix - Sparse N×N adjacency matrix where N = ny1*nx2*nz3
%                  Elements are 1 where points are considered "adjacent"
%                  and might join the same cluster, 0 otherwise

%% input sanity checks

fieldNeeded = 'clusterParams';
if ~any(contains(fieldnames(dg_cfg),fieldNeeded))
    error(['dg_cfg needs field: ' fieldNeeded])
end

%% shortcuts 

ny1 = dg_cfg.dimensions.y1_num;
nx2 = dg_cfg.dimensions.x2_num;
nz3 = dg_cfg.dimensions.z3_num;
Nall = ny1*nx2*nz3;

distance_y1x2_euclidean = dg_cfg.clusterParams.distance_y1x2_euclidean^2; % squared for speed
distance_z3_angular = dg_cfg.clusterParams.distance_z3_angular;
z3_distanceMatrix = dg_cfg.clusterParams.z3_distanceMatrix;

%% additional data

[y1Matrix, x2Matrix, z3Matrix] = ndgrid(1:ny1,1:nx2,1:nz3);

%% choose approach to implementation

if Nall < 100000
    sparseApproach = 'fromLogical'; 
else
    sparseApproach = 'fromIdx';
end

%% implementation

switch sparseApproach
    case 'fromLogical'

        adjMatrix_logical = false(Nall,Nall); 

        for pIdx = 1:Nall
            [pnt_y1, pnt_x2, pnt_z3] = ind2sub([ny1 nx2 nz3], pIdx);

            criterion1 = (y1Matrix-pnt_y1).^2 + (x2Matrix-pnt_x2).^2 <= distance_y1x2_euclidean;
            chanNeighIdx = find(z3_distanceMatrix(pnt_z3,:) <= distance_z3_angular);
            criterion2 = ismember(z3Matrix,chanNeighIdx);

            adjMatrix_logical(pIdx,:) = reshape(criterion1 & criterion2,[1 ny1*nx2*nz3]);

            if pIdx==1
                disp('computing the adjacency matrix')
            end
            dg_counter(pIdx,Nall);
        end
        
        adjMatrix = sparse(adjMatrix_logical);
        clear adjMatrix_logical

    case 'fromIdx'

        a_idx = [];
        b_idx = [];

        for pIdx = 1:Nall
            [pnt_y1, pnt_x2, pnt_z3] = ind2sub([ny1 nx2 nz3], pIdx);

            criterion1 = (y1Matrix-pnt_y1).^2 + (x2Matrix-pnt_x2).^2 <= distance_y1x2_euclidean;
            chanNeighIdx = find(z3_distanceMatrix(pnt_z3,:) <= distance_z3_angular);
            criterion2 = ismember(z3Matrix,chanNeighIdx);

            criteria = criterion1 & criterion2;
            criteria_idx = find(criteria);

            a_idx = [a_idx; repmat(pIdx, numel(criteria_idx), 1)];
            b_idx = [b_idx; criteria_idx(:)];

            if pIdx==1
                disp('computing the adjacency matrix')
            end
            dg_counter(pIdx,Nall);
        end

        value = 1;
        adjMatrix = sparse(a_idx, b_idx, value, Nall, Nall);
        
end

%% adjacency matrix figure
if dg_cfg.figFlag
    figure()
    spy(adjMatrix)
    xlabel('all dimensions ny1*nx2*nz3')
    ylabel('all dimensions ny1*nx2*nz3')
    title('adjacency matrix in y1x2z3 space')
    f = gcf;
    f.Units = 'normalized';
    f.Position = [0.6 0.6 0.3 0.3];
end

%% fig proximity

if dg_cfg.figFlag
    distMatrix_halfExtent = 3; 
    if ny1>1
        if nx2>1
            [distMatrix_y1, distMatrix_x2] = ndgrid(-distMatrix_halfExtent:distMatrix_halfExtent,-distMatrix_halfExtent:distMatrix_halfExtent);
        else
            [distMatrix_y1, distMatrix_x2] = ndgrid(-distMatrix_halfExtent:distMatrix_halfExtent,0);
        end
    else
        if nx2>1
            [distMatrix_y1, distMatrix_x2] = ndgrid(0,-distMatrix_halfExtent:distMatrix_halfExtent);
        else
            [distMatrix_y1, distMatrix_x2] = 0;
        end
    end
    distMatrix_val = sqrt(distMatrix_y1.^2 + distMatrix_x2.^2);
    distMatrix_val_thresholded = zeros(size(distMatrix_val));
    distMatrix_val_thresholded(distMatrix_val <= distance_y1x2_euclidean) = 1;
    [~, minIdx] = min(distMatrix_val(:));
    distMatrix_val_thresholded(minIdx) = 2;
    figure(); clf
    f = gcf; f.Units = 'normalized'; f.Position = [0.6    0.05   0.24    0.4];
    imagesc(distMatrix_val_thresholded)
    xlabel('dimension 1');
    set(gca, 'XTick', 0.5:1:size(distMatrix_val, 2)+0.5)
    set(gca, 'XTickLabel', [])
    ylabel('dimension 2');
    set(gca, 'YTick', 0.5:1:size(distMatrix_val, 1)+0.5)
    set(gca, 'YTickLabel', [])
    grid on; 
    set(gca, 'GridColor', [0 0 0], 'GridAlpha', 1, 'LineWidth', 0.5);
    axis equal; axis tight;
    colormap([0.8 0.8 0.8;     0.8500    0.3250    0.0980;      0    0.4470    0.7410]);
    for rowIdx = 1:size(distMatrix_val,1)
        for colIdx = 1:size(distMatrix_val,2)
            switch mod(distMatrix_val(rowIdx,colIdx), 1)
                case 0
                    textStr = sprintf('%d', distMatrix_val(rowIdx, colIdx));
                otherwise
                    textStr = sprintf('√%d', distMatrix_y1(rowIdx, colIdx)^2 + distMatrix_x2(rowIdx, colIdx)^2);
            end
            switch distMatrix_val_thresholded(rowIdx,colIdx)
                case 0
                    textCol = [0 0 0];
                case 1
                    textCol = [1 1 1];
            end
            text(colIdx, rowIdx, textStr, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'Color', textCol, 'FontSize', 10, 'FontWeight', 'normal');
        end
    end
    sgtitle(sprintf('cluster search area (threshold ≤ %d)', distance_y1x2_euclidean), 'FontWeight', 'normal');
    f2 = gcf; f2.Units = 'normalized'; f2.Position = [0.5 0.05 0.4 0.4];
end

end
