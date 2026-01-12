function [adjMatrix] = di_adjacency(di_cfg)
% Create an adjacency matrix for multi-dimensional data based on
% continuity/neighborhood along dimensions in the data, which can only be
% defined for continuous dimensions (e.g., time) or spherical (e.g., EEG
% channels). Using information about the dimensions (see
% di_validateDimensions), it creates a (sparse) adjacency matrix describing
% continuity/neighborhood across a vectorized pool of all levels of all
% dimensions (i.e., all levels vs all levels). 1=neighbor, 0=not neighbor.
% Having an adjacency matrix is essential for the cluster forming
% algorithm, which will not waste a huge amount of time checking whether
% far away points are in the same cluster.
% 
% Criterion 1: Euclidean distance over all continuous dimensions (if present)
% Criterion 2: Angular distance over spherical dimension (if present)
% Categorical dimensions are not supported for adjacency computation.
% 
% INPUT: 
%       di_cfg - Configuration structure with dimensions and clusterParams
%
% OUTPUT:
%       adjMatrix - Sparse N×N adjacency matrix where N = product of all dimension sizes
%                  Elements are 1 where points are considered "adjacent"
%                  and might join the same cluster (they are candidate in
%                  the cluster forming algorithm), 0 otherwise
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% input sanity checks

if ~isfield(di_cfg,'dimensions')
    error('di_cfg.dimensions needed')
end

if ~isfield(di_cfg.analysis,'clusterParams')
    error('di_cfg.analysis.clusterParams needed')
end

%% shortcuts 

% Extract dimensions in d# order (fieldnames order)
dimKeys = fieldnames(di_cfg.dimensions);
dimTypes = cellfun(@(k) di_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);
dimSizes = cellfun(@(k) length(di_cfg.dimensions.(k).vec), dimKeys);
nDims = numel(dimKeys);

% Identify dimension types
contIdx = find(strcmp(dimTypes, 'continuous'));
sphIdx = find(strcmp(dimTypes, 'spherical'));
catIdx = find(strcmp(dimTypes, 'categorical'));

% sanity check categorical dimensions (not supported in adjacency)
% TO DO (maybe.. or not) in future allow adjacency for the other dimensions (to work
% for each level of the categorical dimension(s))
if ~isempty(catIdx)
    error('Adjacency cannot be computed for categorical dimensions');
end

% Total number of points
pX = prod(dimSizes);

% distance thresholds
% continuous dimensios, if present
if ~isempty(contIdx)
    distance_cont_euclid_sq = di_cfg.analysis.clusterParams.distance_continuous_index^2; % it's faster because it avoid the sqrt() later
else
    distance_cont_euclid_sq = [];
end
% spherical dimension, if present
if ~isempty(sphIdx)
    distance_spherical_radians = di_cfg.analysis.clusterParams.distance_spherical_radians;
    sphericalDistanceMatrix    = di_cfg.analysis.clusterParams.sphericalDistanceMatrix;
else
    distance_spherical_radians = [];
    sphericalDistanceMatrix = [];
end

%% additional data

% coordinate grids for all dimensions
dimVecIdxCell = cell(1, nDims); % initialize cell with sizes per dimension
for dimIdx = 1:nDims
    dimVecIdxCell{dimIdx} = 1:dimSizes(dimIdx);
end
dimGrids = cell(1, nDims);
[dimGrids{:}] = ndgrid(dimVecIdxCell{:});

% % prepare spherical grid if present  // DELETE
% if ~isempty(sphIdx)
%     sphGrid = dimGrids{sphIdx(1)};
% else
%     sphGrid = [];
% end

%% choose approach to implementation

if pX < 100000
    sparseApproach = 'fromLogical';
    % this was the first-ish approach
    % it creates first a logical N x N matrix, fills it, then convert it to sparse
    % it fails for large N at the first step, but otherwise it's faster for small N only

else
    sparseApproach = 'fromIdx';
    % it creates sparse matrix from indices
    % idx approach is less memory intensive, but slower
end

%% implementation

switch sparseApproach
    case 'fromLogical'

        adjMatrix_logical = false(pX,pX); % initialize

        % loop through each point and find its adjacent across all dimensions
        for pIdx = 1:pX

            % find the coordinates in the multidimensional grid for point pIdx
            pIdx_subs = cell(1, nDims);
            [pIdx_subs{:}] = ind2sub(dimSizes, pIdx);

            % adjacency along all continuous dimensions, if present
            % Euclidean distance between indices (not physical units)
            if ~isempty(contIdx)
                distSq = zeros(size(dimGrids{1})); % initialize
                for contDimIdx = contIdx
                    distSq = distSq + (dimGrids{contDimIdx} - pIdx_subs{contDimIdx}).^2;
                end
                criterion1 = distSq <= distance_cont_euclid_sq;
            else
                criterion1 = true(size(dimGrids{1}));
            end

            % adjacency along the pherical dimension, if present
            % angular distance using angular distance matrix
            if ~isempty(sphIdx)
                seedChanIdx = pIdx_subs{sphIdx(1)};
                chanNeighIdx = find(sphericalDistanceMatrix(seedChanIdx, :) <= distance_spherical_radians);
                %chanNeighLbl =
                %di_cfg.dimensions.(dimKeys{sphIdx}).vec(chanNeighIdx);  % just for curiosity
                criterion2 = ismember(dimGrids{sphIdx}, chanNeighIdx);
            else
                criterion2 = true(size(dimGrids{1}));
            end

            % nnz(criterion1 & criterion2); % just for curiosity
            adjMatrix_logical(pIdx,:) = reshape(criterion1 & criterion2,[1, pX]);
            
            if pIdx==1
                disp('computing the adjacency matrix')
            end
            di_counter(pIdx,pX);
        end
        
        adjMatrix = sparse(adjMatrix_logical);
        clear adjMatrix_logical % to save memory

    case 'fromIdx'

% TO DO: double check this still works after the diagonale change with
% multiple dimensions

        a_idx = [];
        b_idx = [];

        for pIdx = 1:pX
            pIdx_subs = cell(1, nDims);
            [pIdx_subs{:}] = ind2sub(dimSizes, pIdx);

            % Euclidean criterion across continuous dimensions
            if isempty(contIdx)
                criterion1 = true(size(dimGrids{1}));
            else
                distSq = zeros(size(dimGrids{1}));
                for contDimIdx = contIdx
                    distSq = distSq + (dimGrids{contDimIdx} - pIdx_subs{contDimIdx}).^2;
                end
                criterion1 = distSq <= distance_cont_euclid_sq;
            end

            % angular criterion for spherical dimension (if present)
            if ~isempty(sphIdx)
                chanNeighIdx = find(sphericalDistanceMatrix(pIdx_subs{sphIdx(1)}, :) <= distance_spherical_radians);
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
            di_counter(pIdx,pX);
        end

        value = 1;
        adjMatrix = sparse(a_idx, b_idx, value, pX, pX);
        
end

%% adjacency matrix figure
if di_cfg.analysis.figFlag
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

if di_cfg.analysis.figFlag && ~isempty(contIdx)
    distMatrix_halfExtent = di_cfg.analysis.clusterParams.distance_continuous_index + 1;

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
