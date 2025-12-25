function [clusterMembership, clustIDList, clusterMetrics] = ...
    dg_clusterForming(dg_cfg, statVal, pVal)
% this algorithm identifies clusters along any number of
% continuous dimensions and/or one spherical dimension
%
% INPUT:
%   dg_cfg   - configuration with analysis.clusterParams.clusterFormingThreshold
%   statVal  - (1 x Nall) test statistics
%   pVal     - (1 x Nall) p-values or pseudo-significance mask
%
% OUTPUT:
%   clusterMembership - (1 x Nall) cluster ID per feature (0 = not in cluster)
%   clustIDList       - unique cluster IDs
%   clusterMetrics    - struct with descriptive metrics per cluster:
%       .id            - cluster IDs matching clustIDList
%       .pointIdx      - {1 x nClusters} cell array of linear indices of cluster members
%       .size          - number of features in each cluster
%       .mass          - sum of statistics in each cluster
%       .medianVal     - median statistic value in each cluster
%       .extremeIdx    - linear index of peak/trough (max absolute value)
%       .extremeVal    - statistic value at extreme (preserves sign)
%       .leastExtremeIdx - linear index of weakest point (min absolute value)
%       .leastExtremeVal - statistic value at weakest point (preserves sign)
%       .extent_cont_dim  - struct with subfields for each continuous dimension
%           .(dimKey).minIdx/maxIdx  - dimensional indices of extent bounds
%           .(dimKey).minVal/maxVal  - dimensional values at extent bounds
%       .extent_sph_dim   - struct with subfields for each spherical dimension
%           .(dimKey).channelCount   - number of unique channels in cluster
%           .(dimKey).channelIdx     - array of specific channel indices
%
% Note: dual role
%   - inferential
%   clusters are created, metrics are extracted, and inference is done on
%   those cluster metrics. question: are my observed cluster metrics
%   significantly larger than simulated cluster metrics? clusters are 
%   localized but inference is not done on the localization, just on their
%   existance (see a couple of papers: Maris and Oostenveltd; Groppe;
%   somethig else)
%
%   - descriptive
%   inference is done before forming clusters. clusters are formed post hoc merely
%   to describe the pre-cluster results, but it's strictly not needed. this
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% shortcuts and input checks

% Extract dimensions in d# order (fieldnames order)
dimKeys  = fieldnames(dg_cfg.dimensions);
dimTypes = cellfun(@(k) dg_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);
dimSizes = cellfun(@(k) length(dg_cfg.dimensions.(k).vec), dimKeys);
nDims    = numel(dimKeys);
Nall     = prod(dimSizes);

% categorical dimensions are not supported for clustering
if any(strcmp(dimTypes, 'categorical'))
    error('Categorical dimensions are not supported for cluster analyses');
end
if nnz(strcmp(dimTypes, 'spherical')) > 1
    error('Only one spherical dimension is supported');
end

if isfield(dg_cfg.analysis,'clusterParams') && isfield(dg_cfg.analysis.clusterParams,'adjacencyMatrix')
    adjacencyMatrix = dg_cfg.analysis.clusterParams.adjacencyMatrix;
else
    error('dg_cfg.analysis.clusterParams.adjacencyMatrix is needed')
end

pThreshold = dg_cfg.analysis.clusterParams.clusterFormingPvalThreshold;
ignoreMask = dg_cfg.analysis.ignore_col;

if size(statVal,2) ~= Nall
    warning('the dimensions of statMatrix do not agree with the configured dimensions');
    keyboard
end
if numel(ignoreMask) ~= Nall
    warning('ignore mask length does not match the configured dimensions');
    keyboard
end
if ~isequal(size(adjacencyMatrix,1), Nall) || ~isequal(size(adjacencyMatrix,2), Nall)
    warning('adjacency matrix size does not match the configured dimensions');
    keyboard
end

%% cluster forming 
% by thresholding pvalues (or anything similar)
% algorithm based on connected components analysis (REFERENCE)
% starting point: 
%   - a vector of (1 x Nall) values (eg, statistical values), 
%   - adjacencyMatrix (Nall x Nall) telling which values are physically adjacent to which values

% step 1: full graph from the adjacency matrix (Nall x Nall)
% graph() builds a graph made of nodes and edges using a matrix
% - nodes are all vlaues 1 to Nall
% - edges reèresent nonzero value at the intersection between two nodes
%   - note: if the matrix is symmetric, the graph is nondirectional
%   - note2: edge's magnitude does not matter

G = graph(adjacencyMatrix);

% step 2: keep only a subset of the graph G including only the nodes 
% respecting these criteria:
criterion1 = pVal < pThreshold; % significant points (eg, p-value below threshold) 
criterion2 = ~ignoreMask;       % not to be ignored 
sigMask = criterion1 & criterion2; % binary mask: 1 if point respects criteria, 0 otherwise
sigIdx = find(sigMask); % idx in full graph space for all points to be considered
% subset of graph G
% note: downside is that subgraph does not retain the same indices as the
% original graph G. the indices are in 1..nnz(sigMask)
subG = subgraph(G, sigMask);

% step 3: connect points into clusters
% conncomp() labels each node of the subgraph with a component ID such
% that nodes with the same ID are mutually reachable via the subgraph's
% edges. If there are K connected components, IDs are 1..K.
bins = conncomp(subG,'Type','weak');  % cluster IDs per significant node (in subgraph order)

% step 4: 
clusterMembership = zeros(1, Nall); % initialize cluster membership array (all zeros = not in any cluster)
clustIDList = unique(bins);
% edge case: no significant points (= no clusters)
if nnz(sigIdx)==0 || isempty(bins)
    % return empty metrics
    clusterMetrics.id = clustIDList;
    clusterMetrics.size = zeros(1, numel(clustIDList));
    clusterMetrics.mass = zeros(1, numel(clustIDList));
    clusterMetrics.medianVal = [];
    clusterMetrics.extremeIdx = [];
    clusterMetrics.extremeVal = [];
    clusterMetrics.leastExtremeIdx = [];
    clusterMetrics.leastExtremeVal = [];
    clusterMetrics.extent_cont_dim = struct();
    clusterMetrics.extent_sph_dim = struct();
    clusterMetrics.pointIdx = cell(1, 0);
    return
end
% assign cluster membership
% map cluster IDs from subgraph back to full graph space
% for each cluster k of the subgraph, assign all its idx members to cluster k in full graph space
for k = 1:max(bins)
    members = sigIdx(bins==k);  % get all full graph space points belonging to this cluster
    clusterMembership(members) = k;  % assign these points to cluster k
end

% imporatant TO DO!!: this implementation confounds positive and negative values if they
% are adjacent and significant, so do this process twice by adding a
% criterion 3 to step 2 and looping from there to here. differently from
% before use natural numbers only and maybe another sign vector

%keyboard % UNTIL HERE OK


%% sanity check: statVal points within a cluster have the same sign
% all points within each cluster must have the same sign of statVal

for k = 1:numel(clustIDList)
    clusterIdx = clusterMembership==clustIDList(k);
    clusterStats = statVal(clusterIdx);
    clusterSigns = sign(clusterStats);  % +1, -1, or 0 for each value
    if sum(clusterSigns.^2) ~= length(clusterSigns)
        error('Diagonale: likely a bug: cluster %d contains both positive and negative values.', clustIDList(k));
    end
end

%% compute cluster descriptors
% size
% mass
% extreme values
% median
% dimensional extent

% initialize 
clusterMetrics.id = clustIDList;
clusterMetrics.pointIdx = cell(1, numel(clustIDList));  % linear indices of cluster members
clusterMetrics.size = zeros(1, numel(clustIDList));
clusterMetrics.mass = zeros(1, numel(clustIDList));
clusterMetrics.medianVal = zeros(1, numel(clustIDList));
clusterMetrics.extremeIdx = zeros(1, numel(clustIDList));  % Peak/trough (max absolute value)
clusterMetrics.extremeVal = zeros(1, numel(clustIDList));
clusterMetrics.leastExtremeIdx = zeros(1, numel(clustIDList));  % Weakest point (min absolute value)
clusterMetrics.leastExtremeVal = zeros(1, numel(clustIDList));

% initialize extent structures for continuous and spherical dimensions
clusterMetrics.extent_cont_dim = struct();
clusterMetrics.extent_sph_dim = struct();

% get continuous and spherical dimension keys
contDimIdx = find(strcmp(dimTypes, 'continuous'));
sphDimIdx = find(strcmp(dimTypes, 'spherical'));

% initialize subfields for each continuous dimension (if it exists)
if ~isempty(contDimIdx)
for dimIdx = contDimIdx
    dimKey = dimKeys{dimIdx};
    clusterMetrics.extent_cont_dim.(dimKey).minIdx = zeros(1, numel(clustIDList));
    clusterMetrics.extent_cont_dim.(dimKey).maxIdx = zeros(1, numel(clustIDList));
    clusterMetrics.extent_cont_dim.(dimKey).minVal = zeros(1, numel(clustIDList));
    clusterMetrics.extent_cont_dim.(dimKey).maxVal = zeros(1, numel(clustIDList));
end
end

% initialize subfields for spherical dimension (if exists)
if ~isempty(sphDimIdx)
    sphKey = dimKeys{sphDimIdx};
    clusterMetrics.extent_sph_dim.(sphKey).channelCount = zeros(1, numel(clustIDList));
    clusterMetrics.extent_sph_dim.(sphKey).channelIdx = cell(1, numel(clustIDList));
end

% fill the structure
for k = 1:numel(clustIDList)
    % cluster membership mask: logical array where 1=member of cluster k, 0=non-member of that cluster
    idx = clusterMembership == clustIDList(k);
    clusterMetrics.pointIdx{k} = find(idx);
        
    % size (measured in points, not physical units)
    clusterMetrics.size(k) = nnz(idx);

    % mass: sum of statistic values in cluster (preserves sign, reflects magnitude)
    clusterMetrics.mass(k) = sum(statVal(idx));

    % median of statistics in cluster
    clusterMetrics.medianVal(k) = median(statVal(idx));
    
    % most extreme: max absolute value (peak or trough, depending on sign)
    % most extreme point in the cluster
    clusterStats = statVal(idx);
    [~, localExtIdx] = max(abs(clusterStats));  % idx of max absolute value
    globalExtIdx = find(idx);
    globalExtIdx = globalExtIdx(localExtIdx);  % convert back to global index
    
    clusterMetrics.extremeIdx(k) = globalExtIdx;  % idx of extreme
    clusterMetrics.extremeVal(k) = statVal(globalExtIdx);  % statVal of extreme
    
    % least extreme: same as most extreme but min absolute value
    % the weakest point in the cluster
    [~, localLeastIdx] = min(abs(clusterStats));  % idx of min absolute value
    globalLeastIdx = find(idx);
    globalLeastIdx = globalLeastIdx(localLeastIdx);  % convert back to global index
    
    clusterMetrics.leastExtremeIdx(k) = globalLeastIdx;  % idx or least extreme
    clusterMetrics.leastExtremeVal(k) = statVal(globalLeastIdx);  % statVal of least extreme
    
    % continuous dimensions: extent (min/max) per dimension
    % For each continuous dimension, find the range of the cluster
    if ~isempty(contDimIdx)
        for dimIdx = contDimIdx
            dimKey = dimKeys{dimIdx};
            dimVec = dg_cfg.dimensions.(dimKey).vec;

            % in2sub to extract indices (per dimension) of cluster members
            [subs{1:nDims}] = ind2sub(dimSizes, find(idx));
            dimIndicesInCluster = subs{dimIdx};  % indices of cluster members in this dimension
            dimValsInCluster = dimVec(dimIndicesInCluster);  % actual values based on dimVec (eg, frequencies, times)

            % find extent of cluster along this dimension
            [minVal, minLocalIdx] = min(dimValsInCluster);  % min value (ie, lowest therefore earliest)
            [maxVal, maxLocalIdx] = max(dimValsInCluster);  % max value (ie, largest therefore latest)


            % Store both the dimensional indices and values separatele per dimension
            clusterMetrics.extent_cont_dim.(dimKey).minVal(k) = minVal;
            clusterMetrics.extent_cont_dim.(dimKey).maxVal(k) = maxVal;
            clusterMetrics.extent_cont_dim.(dimKey).minIdx(k) = dimIndicesInCluster(minLocalIdx);
            clusterMetrics.extent_cont_dim.(dimKey).maxIdx(k) = dimIndicesInCluster(maxLocalIdx);

        end
    end
    
    % spherical dimension: unique channels in cluster
    % identify which spatial sensors/channels participate in the cluster
    if ~isempty(sphDimIdx)

        sphKey = dimKeys{sphDimIdx};
        [subs{1:nDims}] = ind2sub(dimSizes, find(idx));
        sphIndicesInCluster = subs{sphDimIdx};  % channel idxof cluster members
        uniqueSphIdx = unique(sphIndicesInCluster);  % remove duplicates
        
        clusterMetrics.extent_sph_dim.(sphKey).channelCount(k) = numel(uniqueSphIdx); % num of unique channels
        clusterMetrics.extent_sph_dim.(sphKey).channelIdx{k} = uniqueSphIdx;  % cell array for flexibility
    end
end
