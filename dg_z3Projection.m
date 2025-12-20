function coord_2d = dg_z3Projection(coord_3d, projectionType)

radius = 1;

projectionOptions = ["orthographic" "azimuthalEquidistant" "azimuthalConformal"];
if ~ismember(projectionType,projectionOptions)
    disp(projectionOptions)
    error('projection option must be one of the above')
end

az  = deg2rad([coord_3d.sph_theta]);
el  = deg2rad([coord_3d.sph_phi]);
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
