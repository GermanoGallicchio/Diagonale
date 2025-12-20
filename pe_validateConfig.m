function pe_validateConfig(pe_cfg, required)
% Validate presence of required fields in pe_cfg and their basic shapes.
% Usage:
%   pe_validateConfig(pe_cfg, required)
% where 'required' is a struct with logical flags for sub-structs, e.g.,
%   required.dimensions = true; required.clusterParams = true;
%
% This function throws informative errors if fields are missing or invalid.
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

if nargin<2; required = struct(); end

% top-level checks
if isfield(required,'dimensions') && required.dimensions
    if ~isfield(pe_cfg,'dimensions')
        error('pe_cfg needs field: dimensions');
    end
    dims = pe_cfg.dimensions;
    f = {'y1_num','x2_num','z3_num'};
    for k=1:numel(f)
        if ~isfield(dims, f{k}) || ~isscalar(dims.(f{k})) || ~isnumeric(dims.(f{k}))
            error('pe_cfg.dimensions.%s must be a numeric scalar', f{k});
        end
    end
end

if isfield(required,'clusterParams') && required.clusterParams
    if ~isfield(pe_cfg,'clusterParams')
        error('pe_cfg needs field: clusterParams');
    end
    cp = pe_cfg.clusterParams;
    if ~isfield(cp,'distance_y1x2_euclidean') || ~isscalar(cp.distance_y1x2_euclidean)
        error('pe_cfg.clusterParams.distance_y1x2_euclidean must be a scalar');
    end
    if ~isfield(cp,'distance_z3_angular') || ~isscalar(cp.distance_z3_angular)
        error('pe_cfg.clusterParams.distance_z3_angular must be a scalar');
    end
    if ~isfield(cp,'z3_distanceMatrix') || ~ismatrix(cp.z3_distanceMatrix)
        error('pe_cfg.clusterParams.z3_distanceMatrix must be a square matrix');
    end
end

end
