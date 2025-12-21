function coord_2d = dg_sphericalProjection(dg_cfg, projectionType)
% Project spherical coordinates to 2D plane
% Inputs:
%   dg_cfg - configuration struct with validated dimensions containing spherical coordinates
%   projectionType - 'orthographic', 'azimuthalEquidistant', or 'azimuthalConformal'
% Output:
%   coord_2d - [nChan x 2] matrix of 2D coordinates

% Extract dimension metadata
dimKeys = fieldnames(dg_cfg.dimensions);
dimTypes = cellfun(@(k) dg_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);

% Find spherical dimension
sphIdx = find(strcmp(dimTypes, 'spherical'));
if isempty(sphIdx)
    error('dg_z3Projection requires at least one spherical dimension in dg_cfg.dimensions');
end
if length(sphIdx) > 1
    warning('Multiple spherical dimensions found. Using the first one: %s', dimKeys{sphIdx(1)});
end

% Extract spherical coordinates
sphKey = dimKeys{sphIdx(1)};
sphCoord = dg_cfg.dimensions.(sphKey).coord;
nChan = dg_cfg.dimensions.(sphKey).num;

% Add nasion reference point for projection
sphTheta = [sphCoord.sphTheta(:); 90];
sphPhi = [sphCoord.sphPhi(:); 0];

radius = 1;

projectionOptions = ["orthographic" "azimuthalEquidistant" "azimuthalConformal"];
if ~ismember(projectionType, projectionOptions)
    disp(projectionOptions)
    error('projection option must be one of the above')
end

az  = deg2rad(sphTheta(:)');
el  = deg2rad(sphPhi(:)');
colat = (pi/2) - el;

switch projectionType
    case 'orthographic'
        [x3, y3, z3] = sph2cart(az, el, radius);
        xCoord = x3; yCoord = y3;
    case 'azimuthalEquidistant'
        [xTopo, yTopo] = pol2cart(az, radius * colat);
        xCoord = xTopo; yCoord = yTopo;
    case 'azimuthalConformal'
        rho = tan(colat ./ 2);
        [xStereo, yStereo] = pol2cart(az, radius * rho);
        xCoord = xStereo; yCoord = yStereo;
end

coord_2d = [xCoord' yCoord'];
