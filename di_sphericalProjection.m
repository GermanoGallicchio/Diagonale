function [coord_2d, rhoFun] = di_sphericalProjection(di_cfg, projectionType)
% Project spherical coordinates to 2D plane
% Inputs:
%   di_cfg - configuration struct with validated dimensions containing spherical coordinates
%            coord.sphTheta = azimuth   in degrees (0 = right, 90 = front)
%            coord.sphPhi   = elevation in degrees (90 = vertex, 0 = equator, <0 = below)
%            same convention as EEGLAB chanlocs sph_theta / sph_phi
%
%   projectionType - one of:
%       'orthographic'          rho = sin(colat)            head seen from above; NOT injective past the equator
%       'azimuthalEquidistant'  rho = colat                 distance from the vertex is exact (topoplot default)
%       'azimuthalConformal'    rho = tan(colat/2)          stereographic; preserves angles, inflates the rim ~3.5x
%       'azimuthalEqualArea'    rho = sqrt(2)*sin(colat/2)  Lambert; area on the page = area on the scalp
%
% Output:
%   coord_2d - [(nChan+1) x 2] matrix of 2D coordinates
%              rows 1:nChan   channels, in di_cfg dimension order
%              row  nChan+1   nasion reference; carries no data (see below)
%   rhoFun   - function handle colat -> native radial distance. Callers that
%              rasterise (di_view_Sph_voronoi) need rhoFun(pi/2) to convert
%              between native and normalised radii.
%
% NOTE: each projection has its own radial scale -- at the equator rho is
% 1, pi/2, 1 and sqrt(2) respectively. coord_2d is therefore NOT comparable
% across projectionType values. Callers normalise against the nasion row.
%
% Author: germano.gallicchio@gmail.com

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

% reference coordinate
% Add nasion reference point for projection
% appended as an extra "channel" at azimuth 90, elevation 0 (equator, front).
% it is not associated with data -- it gives callers a fixed point on the
% equator to normalise against, so the head circumference lands consistently
% regardless of projectionType. callers must drop row nChan+1 before use.
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
    case 'orthographic',         rhoFun = @(c) sin(c);
    case 'azimuthalEquidistant', rhoFun = @(c) c;
    case 'azimuthalConformal',   rhoFun = @(c) tan(c ./ 2);
    case 'azimuthalEqualArea',   rhoFun = @(c) sqrt(2) * sin(c ./ 2);
    otherwise, error('\\ unhandled projection: %s', projectionType)
end
rho = rhoFun(colat);
[xCoord, yCoord] = pol2cart(az, radius * rho);

% equidistant, conformal and equal-area are monotonic in colat, so they
% handle sub-equatorial electrodes correctly. orthographic is not.
foldTol = deg2rad(5);
if strcmp(projectionType,'orthographic') && any(colat > pi/2 + foldTol)
    warning('\\ orthographic is not injective past the equator; %d channel(s) more than 5 deg below it will fold inward', sum(colat > pi/2 + foldTol))
end

coord_2d = [xCoord' yCoord'];
