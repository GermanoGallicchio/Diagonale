function colat = di_sphericalProjectionInverse(rho, projectionType)
switch projectionType
    case 'orthographic',         colat = asin(min(rho,1));
    case 'azimuthalEquidistant', colat = rho;
    case 'azimuthalConformal',   colat = 2*atan(rho);
    case 'azimuthalEqualArea',   colat = 2*asin(min(rho./sqrt(2),1));
    otherwise, error('\\ unhandled projection: %s', projectionType)
end
