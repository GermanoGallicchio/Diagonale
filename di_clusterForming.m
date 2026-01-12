function [clusterMembership, clustIDList, clusterMetrics] = ...
    di_clusterForming(di_cfg, statVal, pVal)
% this algorithm identifies clusters along any number of
% continuous dimensions and/or one spherical dimension
%
% INPUT:
%   di_cfg   - configuration with analysis.clusterParams.clusterFormingThreshold
%   statVal  - (1 x pX) test statistics
%   pVal     - (1 x pX) p-values or pseudo-significance mask
%
% OUTPUT:
%   clusterMembership - (1 x pX) cluster ID per feature (0 = not in cluster)
%   clustIDList       - unique cluster IDs
%   clusterMetrics    - struct with descriptive metrics per cluster:
%       .id            - cluster IDs matching clustIDList (natural numbers: 1, 2, 3, ...)
%       .sign          - sign of each cluster (+1 for positive, -1 for negative)
%       .pointIdx      - {1 x nClusters} cell array of linear indices of cluster members
%       .size          - number of features in each cluster
%       .mass          - sum of statistics in each cluster
%       .medianVal     - median statistic value in each cluster
%       .mostExtremeIdx    - linear index of peak/trough (max absolute value)
%       .mostExtremeVal    - statistic value at extreme (preserves sign)
%       .leastExtremeIdx - linear index of weakest point (min absolute value)
%       .leastExtremeVal - statistic value at weakest point (preserves sign)
%       .mostExtreme_coord - struct with per-dimension indices/values for extreme point
%           .(dimKey).idx  - dimensional index of extreme in this dimension
%           .(dimKey).val  - dimensional value of extreme in this dimension
%       .leastExtreme_coord - struct with per-dimension indices/values for leastExtreme point
%           .(dimKey).idx  - dimensional index of leastExtreme in this dimension
%           .(dimKey).val  - dimensional value of leastExtreme in this dimension
%       .extent_cont_dim  - struct with subfields for each continuous dimension
%           .(dimKey).minIdx/maxIdx  - dimensional indices of extent bounds
%           .(dimKey).minVal/maxVal  - dimensional values at extent bounds
%       .extent_sph_dim   - struct with subfields for each spherical dimension
%           .(dimKey).channelCount   - number of unique channels in cluster
%           .(dimKey).channelIdx     - cell array of specific channel indices
%           .(dimKey).channelLabel   - cell array of channel labels/names
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
dimKeys  = fieldnames(di_cfg.dimensions);
dimTypes = cellfun(@(k) di_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);
dimSizes = cellfun(@(k) length(di_cfg.dimensions.(k).vec), dimKeys);
nDims    = numel(dimKeys);
pX       = prod(dimSizes);  % total number of X features across all dimensions

% categorical dimensions are not supported for clustering
if any(strcmp(dimTypes, 'categorical'))
    error('Categorical dimensions are not supported for cluster analyses');
end
if nnz(strcmp(dimTypes, 'spherical')) > 1
    error('Only one spherical dimension is supported');
end

if isfield(di_cfg.analysis,'clusterParams') && isfield(di_cfg.analysis.clusterParams,'adjacencyMatrix')
    adjacencyMatrix = di_cfg.analysis.clusterParams.adjacencyMatrix;
else
    error('di_cfg.analysis.clusterParams.adjacencyMatrix is needed')
end

pThreshold = di_cfg.analysis.clusterParams.clusterFormingPvalThreshold;
ignoreMask = di_cfg.analysis.ignore_col;

if size(statVal,2) ~= pX
    warning('the dimensions of statMatrix do not agree with the configured dimensions');
    keyboard
end
if numel(ignoreMask) ~= pX
    warning('ignore mask length does not match the configured dimensions');
    keyboard
end
if ~isequal(size(adjacencyMatrix,1), pX) || ~isequal(size(adjacencyMatrix,2), pX)
    warning('adjacency matrix size does not match the configured dimensions');
    keyboard
end

%% cluster forming 
% by thresholding pvalues (or anything similar)
% algorithm based on connected components analysis (REFERENCE)
% starting point: 
%   - a vector of (1 x pX) values (eg, statistical values), 
%   - adjacencyMatrix (pX x pX) telling which values are physically adjacent to which values

% step 1: full graph from the adjacency matrix (pX x pX)
% graph() builds a graph made of nodes and edges using a matrix
% - nodes are all vlaues 1 to pX
% - edges reèresent nonzero value at the intersection between two nodes
%   - note: if the matrix is symmetric, the graph is nondirectional
%   - note2: edge's magnitude does not matter

G = graph(adjacencyMatrix);

% step 2: keep only a subset of the graph G including only the nodes 
% respecting these criteria:
criterion1 = pVal < pThreshold; % significant points (eg, p-value below threshold) 
criterion2 = ~ignoreMask;       % not to be ignored 
% criterion3 will separate positive and negative values to process them independently

% Initialize cluster membership array (all zeros = not in any cluster)
clusterMembership = zeros(1, pX);

% Initialize cluster sign tracker (maps cluster ID to its sign: +1 or -1)
clusterSign = [];

% Process positive and negative significant points separately to avoid sign mixing
% All clusters use natural number IDs (1, 2, 3, ...), with sign stored separately

% POSITIVE CLUSTERS
criterion3_pos = statVal > 0;  % positive values only
sigMask_pos = criterion1 & criterion2 & criterion3_pos;
sigIdx_pos = find(sigMask_pos);

nPosClusters = 0;  % track number of positive clusters
if nnz(sigIdx_pos) > 0
    % Create subgraph containing only positive significant points
    subG_pos = subgraph(G, sigMask_pos);
    
    % step 3: connect points into clusters
    % conncomp() labels each node of the subgraph with a component ID such
    % that nodes with the same ID are mutually reachable via the subgraph's
    % edges. If there are K connected components, IDs are 1..K.
    bins_pos = conncomp(subG_pos);
    nPosClusters = max(bins_pos);
    
    % step 4: assign cluster membership with positive (regardless of sign) integers IDs
    % map cluster IDs from subgraph back to full graph space
    for k = 1:nPosClusters
        members = sigIdx_pos(bins_pos == k);  % get all full graph space points belonging to this cluster
        clusterMembership(members) = k;  % assign cluster ID
        clusterSign(k) = +1;  % store positive sign for this cluster
    end
end

% NEGATIVE CLUSTERS
criterion3_neg = statVal < 0;  % negative values only
sigMask_neg = criterion1 & criterion2 & criterion3_neg;
sigIdx_neg = find(sigMask_neg);

if nnz(sigIdx_neg) > 0
    % Create subgraph containing only negative significant points
    subG_neg = subgraph(G, sigMask_neg);
    
    % step 3: connect points into clusters
    bins_neg = conncomp(subG_neg);
    nNegClusters = max(bins_neg);
    
    % step 4: assign cluster membership with positive (regardless of sign) integers IDs
    % map cluster IDs from subgraph back to full graph space
    for k = 1:nNegClusters
        clustID = nPosClusters + k;  % continue numbering after positive clusters
        members = sigIdx_neg(bins_neg == k);  % get all full graph space points belonging to this cluster
        clusterMembership(members) = clustID;  % assign cluster ID
        clusterSign(clustID) = -1;  % store negative sign for this cluster
    end
end

% get list of all cluster IDs (natural numbers only)
clustIDList = unique(clusterMembership(clusterMembership ~= 0));

% edge case: no significant points (= no clusters)
if isempty(clustIDList)
    % return empty metrics
    clusterMetrics.id = clustIDList;
    clusterMetrics.sign = [];
    clusterMetrics.size = zeros(1, 0);
    clusterMetrics.mass = zeros(1, 0);
    clusterMetrics.medianVal = [];
    clusterMetrics.mostExtremeIdx = [];
    clusterMetrics.mostExtremeVal = [];
    clusterMetrics.leastExtremeIdx = [];
    clusterMetrics.leastExtremeVal = [];
    clusterMetrics.extent_cont_dim = struct();
    clusterMetrics.extent_sph_dim = struct();
    clusterMetrics.pointIdx = cell(1, 0);
    return
end

%% compute cluster descriptors
% size
% mass
% extreme values
% median
% dimensional extent

% initialize 
clusterMetrics.id = clustIDList;
clusterMetrics.sign = clusterSign(clustIDList);  % sign (+1 or -1) for each cluster
clusterMetrics.pointIdx = cell(1, numel(clustIDList));  % linear indices of cluster members
clusterMetrics.size = zeros(1, numel(clustIDList));
clusterMetrics.mass = zeros(1, numel(clustIDList));
clusterMetrics.medianVal = zeros(1, numel(clustIDList));
clusterMetrics.mostExtremeIdx = zeros(1, numel(clustIDList));  % Peak/trough (max absolute value)
clusterMetrics.mostExtremeVal = zeros(1, numel(clustIDList));
clusterMetrics.leastExtremeIdx = zeros(1, numel(clustIDList));  % Weakest point (min absolute value)
clusterMetrics.leastExtremeVal = zeros(1, numel(clustIDList));

% initialize extent and coordinate structures
clusterMetrics.extent_cont_dim = struct();
clusterMetrics.extent_sph_dim = struct();
clusterMetrics.mostExtreme_coord = struct();
clusterMetrics.leastExtreme_coord = struct();

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
    % coordinate structures for extreme points
    clusterMetrics.mostExtreme_coord.(dimKey).idx = zeros(1, numel(clustIDList));
    clusterMetrics.mostExtreme_coord.(dimKey).val = zeros(1, numel(clustIDList));
    clusterMetrics.leastExtreme_coord.(dimKey).idx = zeros(1, numel(clustIDList));
    clusterMetrics.leastExtreme_coord.(dimKey).val = zeros(1, numel(clustIDList));
end
end

% initialize subfields for spherical dimension (if exists)
if ~isempty(sphDimIdx)
    sphKey = dimKeys{sphDimIdx};
    clusterMetrics.extent_sph_dim.(sphKey).channelCount = zeros(1, numel(clustIDList));
    clusterMetrics.extent_sph_dim.(sphKey).channelIdx = cell(1, numel(clustIDList));
    clusterMetrics.extent_sph_dim.(sphKey).channelLabel = cell(1, numel(clustIDList));
    % coordinate structures for extreme points
    clusterMetrics.mostExtreme_coord.(sphKey).idx = zeros(1, numel(clustIDList));
    clusterMetrics.mostExtreme_coord.(sphKey).label = cell(1, numel(clustIDList));
    clusterMetrics.leastExtreme_coord.(sphKey).idx = zeros(1, numel(clustIDList));
    clusterMetrics.leastExtreme_coord.(sphKey).label = cell(1, numel(clustIDList));
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
    
    clusterMetrics.mostExtremeIdx(k) = globalExtIdx;  % idx of extreme
    clusterMetrics.mostExtremeVal(k) = statVal(globalExtIdx);  % statVal of extreme
    
    % least extreme: same as most extreme but min absolute value
    % the weakest point in the cluster
    [~, localLeastIdx] = min(abs(clusterStats));  % idx of min absolute value
    globalLeastIdx = find(idx);
    globalLeastIdx = globalLeastIdx(localLeastIdx);  % convert back to global index
    
    clusterMetrics.leastExtremeIdx(k) = globalLeastIdx;  % idx or least extreme
    clusterMetrics.leastExtremeVal(k) = statVal(globalLeastIdx);  % statVal of least extreme
    
    % compute per-dimension coordinates for extreme points
    [subs{1:nDims}] = ind2sub(dimSizes, globalExtIdx);  % extreme point subscripts
    [subs_least{1:nDims}] = ind2sub(dimSizes, globalLeastIdx);  % leastExtreme point subscripts
    
    if ~isempty(contDimIdx)
        for dimIdx = contDimIdx
            dimKey = dimKeys{dimIdx};
            dimVec = di_cfg.dimensions.(dimKey).vec;
            
            % Most extreme point coordinates
            dimIdx_extreme = subs{dimIdx};
            clusterMetrics.mostExtreme_coord.(dimKey).idx(k) = dimIdx_extreme;
            clusterMetrics.mostExtreme_coord.(dimKey).val(k) = dimVec(dimIdx_extreme);
            
            % Least extreme point coordinates
            dimIdx_least = subs_least{dimIdx};
            clusterMetrics.leastExtreme_coord.(dimKey).idx(k) = dimIdx_least;
            clusterMetrics.leastExtreme_coord.(dimKey).val(k) = dimVec(dimIdx_least);
        end
    end
    
    % Spherical dimension extreme point coordinates
    if ~isempty(sphDimIdx)
        sphKey = dimKeys{sphDimIdx};
        sphVec = di_cfg.dimensions.(sphKey).vec;
        
        % Most extreme point channel
        sphIdx_extreme = subs{sphDimIdx};
        clusterMetrics.mostExtreme_coord.(sphKey).idx(k) = sphIdx_extreme;
        clusterMetrics.mostExtreme_coord.(sphKey).label{k} = sphVec{sphIdx_extreme};
        
        % Least extreme point channel
        sphIdx_least = subs_least{sphDimIdx};
        clusterMetrics.leastExtreme_coord.(sphKey).idx(k) = sphIdx_least;
        clusterMetrics.leastExtreme_coord.(sphKey).label{k} = sphVec{sphIdx_least};
    end
    
    % continuous dimensions: extent (min/max) per dimension
    % For each continuous dimension, find the range of the cluster
    if ~isempty(contDimIdx)
        for dimIdx = contDimIdx
            dimKey = dimKeys{dimIdx};
            dimVec = di_cfg.dimensions.(dimKey).vec;

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
        sphVec = di_cfg.dimensions.(sphKey).vec;  % channel labels/names
        [subs{1:nDims}] = ind2sub(dimSizes, find(idx));
        sphIndicesInCluster = subs{sphDimIdx};  % channel idx of cluster members
        uniqueSphIdx = unique(sphIndicesInCluster);  % remove duplicates
        
        clusterMetrics.extent_sph_dim.(sphKey).channelCount(k) = numel(uniqueSphIdx); % num of unique channels
        clusterMetrics.extent_sph_dim.(sphKey).channelIdx{k} = uniqueSphIdx;  % cell array for flexibility
        clusterMetrics.extent_sph_dim.(sphKey).channelLabel{k} = sphVec(uniqueSphIdx);  % actual channel labels
    end
end
