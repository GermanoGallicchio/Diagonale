function [clusterMembership, clustIDList, clusterMetrics] = ...
    dg_clusterForming(dg_cfg, statVal, pVal)
% Cluster forming algorithm(s). This function identifies clusters in a 3D
% statistical matrix (y1, x2, z3) using methods (described below):
% - threshold

%% shortcuts

ny1 = dg_cfg.dimensions.y1_num;
nx2 = dg_cfg.dimensions.x2_num;
nz3 = dg_cfg.dimensions.z3_num;

adjacencyMatrix = dg_cfg.clusterParams.adjacencyMatrix;
pThreshold = dg_cfg.clusterParams.clusterFormingThreshold;
ignoreMask = dg_cfg.R_ignore;

%% sanity checks (abbreviated)
if ~isequal(size(statVal,2),[ny1*nx2*nz3])
    warning('the dimensions of statMatrix do not agree with the y1, x2, z3 dg_cfg.dimensions structure');
    keyboard
end

%% cluster forming by threshold on p-values

% binarize significant points (1=significant, 0=otherwise)
sigMask = pVal < pThreshold;
sigMask = sigMask .* ~ignoreMask;

% build graph of adjacency among significant points
G = graph(adjacencyMatrix);
sigIdx = find(sigMask);
subG = subgraph(G, sigIdx);
bins = conncomp(subG);

clusterMembership = zeros(1, ny1*nx2*nz3);
clustIDList = unique(bins);
for k = 1:max(bins)
    members = sigIdx(bins==k);
    clusterMembership(members) = k;
end

% metrics: size and mass
clusterMetrics.id = clustIDList;
clusterMetrics.size = zeros(1, numel(clustIDList));
clusterMetrics.mass = zeros(1, numel(clustIDList));
for k = 1:numel(clustIDList)
    idx = clusterMembership==clustIDList(k);
    clusterMetrics.size(k) = nnz(idx);
    clusterMetrics.mass(k) = sum(statVal(idx));
end
