function coord_2d = di_sphericalProjection(di_cfg, projectionType)
% Project spherical coordinates to 2D plane
% Inputs:
%   di_cfg - configuration struct with validated dimensions containing spherical coordinates
%   projectionType - one of the following ones:
%           'orthographic',
%           'azimuthalEquidistant',
%           'azimuthalConformal'
%           'azimuthalEqualArea' (Lambert)
%
% Output:
%   coord_2d - [nChan x 2] matrix of 2D coordinates

% Extract dimension metadata
dimKeys = fieldnames(di_cfg.dimensions);
dimTypes = cellfun(@(k) di_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);

% Find spherical dimension
sphIdx = find(strcmp(dimTypes, 'spherical'));
if isempty(sphIdx)
    error('\\ di_sphericalDistance requires a spherical dimension in di_cfg.dimensions');
end
if length(sphIdx) > 1
    warning('\\ Multiple spherical dimensions found. Using the first one: %s', dimKeys{sphIdx(1)});
end

% Extract spherical coordinates
sphKey = dimKeys{sphIdx(1)};
sphCoord = di_cfg.dimensions.(sphKey).coord;
nChan = di_cfg.dimensions.(sphKey).num;

% Add nasion reference point for projection
sphTheta = [sphCoord.sphTheta(:); 90];
sphPhi = [sphCoord.sphPhi(:); 0];

radius = 1;

projectionOptions = ["orthographic" "azimuthalEquidistant" "azimuthalConformal" "azimuthalEqualArea"];
if ~ismember(projectionType, projectionOptions)
    disp(projectionOptions)
    error('\\ projection option must be one of the above')
end

az  = deg2rad(sphTheta(:)');
el  = deg2rad(sphPhi(:)');
colat = (pi/2) - el;

switch projectionType
    case 'orthographic',         rho = sin(colat);
    case 'azimuthalEquidistant', rho = colat;
    case 'azimuthalConformal',   rho = tan(colat ./ 2);
    case 'azimuthalEqualArea',   rho = sqrt(2) * sin(colat ./ 2);
    otherwise, error('\\ unhandled projection: %s', projectionType)
end
[xCoord, yCoord] = pol2cart(az, radius * rho);

if strcmp(projectionType,'orthographic') && any(colat > pi/2)
    warning('\\ orthographic is not injective past the equator; %d channel(s) with negative elevation will fold inward', sum(colat > pi/2))
end

coord_2d = [xCoord' yCoord'];
