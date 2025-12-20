function dg_validateConfig(dg_cfg, required)
% Validate presence of required fields in dg_cfg and their basic shapes.

if nargin<2; required = struct(); end

if isfield(required,'dimensions') && required.dimensions
    if ~isfield(dg_cfg,'dimensions')
        error('dg_cfg needs field: dimensions');
    end
    dims = dg_cfg.dimensions;
    f = {'y1_num','x2_num','z3_num'};
    for k=1:numel(f)
        if ~isfield(dims, f{k}) || ~isscalar(dims.(f{k})) || ~isnumeric(dims.(f{k}))
            error('dg_cfg.dimensions.%s must be a numeric scalar', f{k});
        end
    end
end

if isfield(required,'clusterParams') && required.clusterParams
    if ~isfield(dg_cfg,'clusterParams')
        error('dg_cfg needs field: clusterParams');
    end
    cp = dg_cfg.clusterParams;
    if ~isfield(cp,'distance_y1x2_euclidean') || ~isscalar(cp.distance_y1x2_euclidean)
        error('dg_cfg.clusterParams.distance_y1x2_euclidean must be a scalar');
    end
    if ~isfield(cp,'distance_z3_angular') || ~isscalar(cp.distance_z3_angular)
        error('dg_cfg.clusterParams.distance_z3_angular must be a scalar');
    end
    if ~isfield(cp,'z3_distanceMatrix') || ~ismatrix(cp.z3_distanceMatrix)
        error('dg_cfg.clusterParams.z3_distanceMatrix must be a square matrix');
    end
end

end
