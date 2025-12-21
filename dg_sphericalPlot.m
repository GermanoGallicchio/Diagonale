function dg_sphericalPlot(dg_cfg, sphValues, params)
% Plot spherical data projected to 2D
% Inputs:
%   dg_cfg - configuration struct with validated dimensions containing spherical coordinates
%   sphValues - (1 x nChan) values to categorize channels: 
%       1=seed channel
%       2=a neighbor
%       0=not a neighbor
%
%   params - struct with optional fields:
%     projectionType, drawLines, lineWidth, lineCol, chanMarkerSize,
%     chanMarkerChar, chanLbl, colBar, colMap, cLim



% Extract dimension metadata
dimKeys = fieldnames(dg_cfg.dimensions);
dimTypes = cellfun(@(k) dg_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);

% Find spherical dimension
sphIdx = find(strcmp(dimTypes, 'spherical'));
if isempty(sphIdx)
    error('dg_z3Plot requires at least one spherical dimension in dg_cfg.dimensions');
end
if length(sphIdx) > 1
    warning('Multiple spherical dimensions found. Using the first one: %s', dimKeys{sphIdx(1)});
end

% Extract spherical coordinates and labels
sphKey = dimKeys{sphIdx(1)};
sphCoord = dg_cfg.dimensions.(sphKey).coord;
nChan = dg_cfg.dimensions.(sphKey).num;

% Build coord_sph struct for projection
coord_sph = repmat(struct('sph_theta', [], 'sph_phi', [], 'labels', ''), nChan, 1);
for i = 1:nChan
    coord_sph(i).sph_theta = sphCoord.sphTheta(i);
    coord_sph(i).sph_phi = sphCoord.sphPhi(i);
    if isfield(sphCoord, 'labels')
        coord_sph(i).labels = sphCoord.labels{i};
    else
        coord_sph(i).labels = sprintf('Chan%d', i);
    end
end
fieldNames = fieldnames(params);

if any(strcmp(fieldNames,'projectionType'))
    projectionType = params.projectionType;
else
    projectionType = 'azimuthalEquidistant';
end

if any(strcmp(fieldNames,'drawLines'))
    drawLines = params.drawLines;
else
    drawLines = true;
end

if any(strcmp(fieldNames,'lineWidth'))
    lineWidth = params.lineWidth;
else
    lineWidth = 1;
end

if any(strcmp(fieldNames,'lineCol'))
    lineCol = params.lineCol;
else
    lineCol = [0 0 0 1];
end

if any(strcmp(fieldNames,'chanMarkerSize'))
    chanMarkerSize  = params.chanMarkerSize;
else
    chanMarkerSize = 5;
end

if any(strcmp(fieldNames,'chanMarkerChar'))
    chanMarkerChar = params.chanMarkerChar;
else
    chanMarkerChar = 'o';
end

if any(strcmp(fieldNames,'chanLbl'))
    chanLbl = params.chanLbl;
else
    chanLbl = false;
end

if any(strcmp(fieldNames,'colBar'))
colBar = params.colBar;
else
    colBar = false;
end

if any(strcmp(fieldNames,'colMap'))
    colMap = params.colMap;
else
    colMap = turbo;
end

if any(strcmp(fieldNames,'cLim'))
    cLim = params.cLim;
else
    cLim = prctile(sphValues, [0 100]);
end

%% reference coordinates
% for the very top and the very right of the circumference
% very top = nasion
% include 3d coordinates of nasion as extra channel
% these coordinates are not associated with data
% they simply work as reference to scale the units and make the circumference pass by

% nasion
coord_sph(nChan+1).sph_theta = 90;
coord_sph(nChan+1).sph_phi   = 0;
coord_sph(nChan+1).labels    = 'nasionRef';

%% project electrode 3D locations to 2D
coord_2d = dg_sphericalProjection(dg_cfg, projectionType);

% scale 2D coordinates based on the top and right reference
% so that they represent respectively height of 1 and length of 1 
% in a [0, 1] axis space
coord_2d(:,1) = (coord_2d(:,1)-min(coord_2d(:,1)));
coord_2d(:,2) = (coord_2d(:,2)-min(coord_2d(:,2)));
%references length to normalize
yRef = coord_2d(nChan+1,2);
xRef = coord_2d(nChan+1,1)*2;
coord_2d(:,1) = coord_2d(:,1) / (xRef - min(coord_2d(:,1)));
coord_2d(:,2) = coord_2d(:,2) / (yRef - min(coord_2d(:,2)));

%% value normalization for colors

sphValues_norm = (sphValues - cLim(1)) / (cLim(2) - cLim(1));
idx = round(1 + sphValues_norm * (size(colMap, 1) - 1));
idx = max(1, min(size(colMap, 1), idx));  % clip to colormap range

%% figure

% channel markers
for chanIdx = 1:nChan
    plot(coord_2d(chanIdx,1), coord_2d(chanIdx,2), ...
        chanMarkerChar, ...
        'MarkerSize', chanMarkerSize, ...
        'MarkerFaceColor',colMap(idx(chanIdx),:), ...
        'MarkerEdgeColor',lineCol(1,1:3));
    hold on
end

% nasion 'x'
wantNasion = false;
if wantNasion
    plot(coord_2d(chanIdx+1,1), coord_2d(chanIdx+1,2), 'X', 'MarkerSize', chanMarkerSize*2); % nasion
end

% channel labels
if chanLbl
    text(coord_2d(1:nChan, 1), coord_2d(1:nChan, 2), {coord_sph(1:nChan).labels}, ...
        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontSize', 8)
end

% color bar
if colBar
    cb = colorbar;
    cb.Ticks = [0 1];
    cb.TickLabels = [ cLim ];
end

% anatomical lines 

if drawLines
    % circumference
    hold on
    theta = linspace(0, 2*pi, 200);           % 200 points around the circle
    r     = coord_2d(nChan+1,2)/2;                   % radius based on nasion's coordinate
    cx    = 0.5; cy = 0.5;                       % center coordinates
    x     = cx + r*cos(theta);
    y     = cy + r*sin(theta);
    circumf = plot(x, y, 'LineWidth', lineWidth, 'Color',lineCol);
    uistack(circumf,'bottom')

    axis equal

    % add left and right ears
    leln = line([-0.03  -0.025],[0.38 0.56],'LineWidth', lineWidth, 'Color',lineCol);
    reln = line([ 1.03   1.025],[0.38 0.56],'LineWidth', lineWidth, 'Color',lineCol);
    % add nose
    lnln = line([0.4 0.5],[0.95 1]+0.05,'LineWidth', lineWidth, 'Color',lineCol);
    rnln = line([0.6 0.5],[0.95 1]+0.05,'LineWidth', lineWidth, 'Color',lineCol);

    uistack(leln,"bottom")
    uistack(reln,"bottom")
    uistack(lnln,"bottom")
    uistack(rnln,"bottom")

end


ax = gca;
ax.Color = [1 1 1 0];
ax.XAxis.Color = [1 1 1 0];
ax.YAxis.Color = [1 1 1 0];
ax.Colormap = colMap;
axis square; 
axis equal; 
%title(projectionType);

hold off
