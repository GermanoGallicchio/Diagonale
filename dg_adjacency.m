function [adjMatrix] = dg_adjacency(dg_cfg)
% Creates an adjacency matrix for multi-dimensional physiological data
%
% Generates a sparse adjacency matrix describing neighborhood relationships
% across any number of continuous dimensions and optionally one spherical dimension.
% 
% Criterion 1: Euclidean distance over all continuous dimensions
% Criterion 2: Angular distance over spherical dimension (if present)
% Categorical dimensions are not supported for adjacency computation.
% 
% INPUT: 
%       dg_cfg - Configuration structure with dimensions and clusterParams
%
% OUTPUT:
%       adjMatrix - Sparse N×N adjacency matrix where N = product of all dimension sizes
%                  Elements are 1 where points are considered "adjacent"
%                  and might join the same cluster, 0 otherwise

%% input sanity checks

if ~isfield(dg_cfg.analysis,'clusterParams')
    error('dg_cfg.analysis.clusterParams needed')
end


%% shortcuts 

% Extract dimensions in d# order (fieldnames order)
dimKeys = fieldnames(dg_cfg.dimensions);
dimTypes = cellfun(@(k) dg_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);
dimSizes = cellfun(@(k) length(dg_cfg.dimensions.(k).vec), dimKeys);
nDims = numel(dimKeys);

% Identify dimension types
contIdx = find(strcmp(dimTypes, 'continuous'));
sphIdx = find(strcmp(dimTypes, 'spherical'));
catIdx = find(strcmp(dimTypes, 'categorical'));

% Check for categorical dimensions (not supported in adjacency)
if ~isempty(catIdx)
    error('Adjacency cannot be computed for categorical dimensions');
end

% Total number of points
Nall = prod(dimSizes);

% Distance thresholds
distance_cont_euclid_sq = dg_cfg.analysis.clusterParams.distance_continuous_index^2;

% Handle spherical dimension if present
if ~isempty(sphIdx)
    distance_chan_angular = dg_cfg.analysis.clusterParams.distance_spherical_radians;
    % Channel angular distance matrix (compute if absent)
    if isfield(dg_cfg.analysis.clusterParams,'sphericalDistanceMatrix')
        sphericalDistanceMatrix = dg_cfg.analysis.clusterParams.sphericalDistanceMatrix;
    else
        sphericalDistanceMatrix = dg_sphericalDistance(dg_cfg);
        dg_cfg.analysis.clusterParams.sphericalDistanceMatrix = sphericalDistanceMatrix;
    end
else
    distance_chan_angular = [];
    sphericalDistanceMatrix = [];
end

%% additional data

% Create coordinate grids for all dimensions
gridArgs = cell(1, nDims);
for dimIdx = 1:nDims
    gridArgs{dimIdx} = 1:dimSizes(dimIdx);
end
dimGrids = cell(1, nDims);
[dimGrids{:}] = ndgrid(gridArgs{:});

% Prepare spherical grid if present
if ~isempty(sphIdx)
    sphGrid = dimGrids{sphIdx(1)};
else
    sphGrid = [];
end

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
            subs = cell(1, nDims);
            [subs{:}] = ind2sub(dimSizes, pIdx);

            % Euclidean criterion across continuous dimensions
            if isempty(contIdx)
                criterion1 = true(size(dimGrids{1}));
            else
                distSq = zeros(size(dimGrids{1}));
                for c = contIdx
                    distSq = distSq + (dimGrids{c} - subs{c}).^2;
                end
                criterion1 = distSq <= distance_cont_euclid_sq;
            end

            % Angular criterion for spherical dimension (if present)
            if ~isempty(sphIdx)
                chanNeighIdx = find(sphericalDistanceMatrix(subs{sphIdx(1)}, :) <= distance_chan_angular);
                criterion2 = ismember(sphGrid, chanNeighIdx);
            else
                criterion2 = true(size(dimGrids{1}));
            end

            adjMatrix_logical(pIdx,:) = reshape(criterion1 & criterion2,[1, Nall]);

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
            subs = cell(1, nDims);
            [subs{:}] = ind2sub(dimSizes, pIdx);

            % Euclidean criterion across continuous dimensions
            if isempty(contIdx)
                criterion1 = true(size(dimGrids{1}));
            else
                distSq = zeros(size(dimGrids{1}));
                for c = contIdx
                    distSq = distSq + (dimGrids{c} - subs{c}).^2;
                end
                criterion1 = distSq <= distance_cont_euclid_sq;
            end

            % Angular criterion for spherical dimension (if present)
            if ~isempty(sphIdx)
                chanNeighIdx = find(sphericalDistanceMatrix(subs{sphIdx(1)}, :) <= distance_chan_angular);
                criterion2 = ismember(sphGrid, chanNeighIdx);
            else
                criterion2 = true(size(dimGrids{1}));
            end

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
if dg_cfg.analysis.figFlag
    figure()
    spy(adjMatrix)
    xlabel('vectorized points')
    ylabel('vectorized points')
    title('adjacency matrix over all dimensions')
    f = gcf;
    f.Units = 'normalized';
    f.Position = [0.6 0.6 0.3 0.3];
end

%% fig proximity (first two continuous dimensions)

if dg_cfg.analysis.figFlag && ~isempty(contIdx)
    distMatrix_halfExtent = dg_cfg.analysis.clusterParams.distance_continuous_index + 1;

    if numel(contIdx) >= 2
        [distMatrix_dim1, distMatrix_dim2] = ndgrid(-distMatrix_halfExtent:distMatrix_halfExtent, -distMatrix_halfExtent:distMatrix_halfExtent);
        dimLabel1 = dimKeys{contIdx(2)};
        dimLabel2 = dimKeys{contIdx(1)};
        
    else
        [distMatrix_dim1, distMatrix_dim2] = ndgrid(0,-distMatrix_halfExtent:distMatrix_halfExtent);
        dimLabel1 = 'constant';
        dimLabel2 = dimKeys{contIdx(1)};
        
    end

    distMatrix_val = sqrt(distMatrix_dim1.^2 + distMatrix_dim2.^2);
    distMatrix_val_thresholded = zeros(size(distMatrix_val));
    distMatrix_val_thresholded(distMatrix_val <= sqrt(distance_cont_euclid_sq)) = 1;
    [~, minIdx] = min(distMatrix_val(:));
    distMatrix_val_thresholded(minIdx) = 2;

    figure(); clf
    f = gcf; f.Units = 'normalized'; f.Position = [0.6    0.05   0.24    0.4];
    imagesc(distMatrix_val_thresholded)
    xlabel(dimLabel2);
    set(gca, 'XTick', 0.5:1:size(distMatrix_val, 2)+0.5)
    set(gca, 'XTickLabel', [])
    ylabel(dimLabel1);
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
                    textStr = sprintf('√%d', distMatrix_dim1(rowIdx, colIdx)^2 + distMatrix_dim2(rowIdx, colIdx)^2);
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
    sgtitle(sprintf('cluster search area (threshold ≤ %.3g)', sqrt(distance_cont_euclid_sq)), 'FontWeight', 'normal');
    f2 = gcf; f2.Units = 'normalized'; f2.Position = [0.5 0.05 0.4 0.4];
end

end
