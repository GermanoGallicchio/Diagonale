function colat = di_sphericalProjectionInverse(rho, projectionType)
% Inverse azimuthal projection: page radius -> colatitude.
%
% Inputs:
%   rho            - NATIVE radial distance, i.e. the same scale returned by
%                    di_sphericalProjection. NOT the normalised coordinates
%                    that di_sphericalPlot produces. Multiply normalised
%                    radii by rhoFun(pi/2) before calling.
%   projectionType - same four options as di_sphericalProjection
%
% Output:
%   colat - colatitude in radians (0 = vertex, pi/2 = equator)
%
% Used by di_view_Sph_voronoi to map each pixel back onto the sphere, so
% that the nearest-electrode decision is taken with the geodesic metric.
% This is what makes the tessellation invariant to projectionType.
%
% Author: germano.gallicchio@gmail.com

switch projectionType
    case 'orthographic',         colat = asin(min(rho,1));
    case 'azimuthalEquidistant', colat = rho;
    case 'azimuthalConformal',   colat = 2*atan(rho);
    case 'azimuthalEqualArea',   colat = 2*asin(min(rho./sqrt(2),1));
    otherwise, error('\\ unhandled projection: %s', projectionType)
end

end
