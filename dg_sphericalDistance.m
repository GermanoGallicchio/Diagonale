function sphericalDistanceMatrix = dg_sphericalDistance(dg_cfg)
% Compute spherical distance matrix from dg_cfg.dimensions
% Input: dg_cfg - configuration struct with validated dimensions
% Output: sphericalDistanceMatrix - angular distance in radians between all channel pairs


% Extract dimension metadata
dimKeys = fieldnames(dg_cfg.dimensions);
dimTypes = cellfun(@(k) dg_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);

% Find spherical dimension
sphIdx = find(strcmp(dimTypes, 'spherical'));
if isempty(sphIdx)
    error('dg_z3distance requires at least one spherical dimension in dg_cfg.dimensions');
end
if length(sphIdx) > 1
    warning('Multiple spherical dimensions found. Using the first one: %s', dimKeys{sphIdx(1)});
end

% Extract spherical coordinates
sphKey = dimKeys{sphIdx(1)};
sphCoord = dg_cfg.dimensions.(sphKey).coord;
sphCoord.labels = dg_cfg.dimensions.(sphKey).vec;
nChan = dg_cfg.dimensions.(sphKey).num;

theta = deg2rad(sphCoord.sphTheta(:)');
phi   = deg2rad(sphCoord.sphPhi(:)');

% Compute angular distance matrix
sphericalDistanceMatrix = zeros(nChan, nChan);
for chanAIdx = 1:nChan
    for chanBIdx = 1:nChan
        if chanAIdx == chanBIdx
            angDist = 0;
        else
            angDist = acos( sin(phi(chanAIdx))*sin(phi(chanBIdx)) + ...
                           cos(phi(chanAIdx))*cos(phi(chanBIdx))*cos(theta(chanAIdx)-theta(chanBIdx)) );
        end
        sphericalDistanceMatrix(chanAIdx, chanBIdx) = angDist;
    end
end

%% figure

if nChan > 2
    % get distance threshold for adjacency
    distance_spherical_angular = dg_cfg.analysis.clusterParams.distance_spherical_radians;
    
    % Bin distances for neighbor determination
    angularEdges = 0:distance_spherical_angular:max(sphericalDistanceMatrix(:));
    if angularEdges(end) ~= max(sphericalDistanceMatrix(:))
        angularEdges(end+1) = max(sphericalDistanceMatrix(:));
    end
    Dbinned = -999*ones(size(sphericalDistanceMatrix));
    dIdx = sphericalDistanceMatrix == angularEdges(1);
    Dbinned(dIdx) = angularEdges(1);
    for eIdx = 1:length(angularEdges)-1
        dIdx = sphericalDistanceMatrix > angularEdges(eIdx) & sphericalDistanceMatrix <= angularEdges(eIdx+1);
        Dbinned(dIdx) = eIdx;
    end
    neighborMatrix = Dbinned == 1;
    
    % Create figure
    figure(); clf
    f = gcf; f.Units = 'normalized'; f.Position = [0.1 0.01 0.8 0.9];
    panelA = uipanel('Position', [0  0.4  1  0.6]);
    panelB = uipanel('Position', [0  0    0.5  0.4]);
    
    % Panel B: Show neighborhood examples for random channels
    tldB = tiledlayout(panelB, 'flow');
    tldB.Padding = 'tight';
    tldB.TileSpacing = 'tight';
    topoPlotFig_chanIdx = sort(randperm(nChan, max(8, ceil(nChan/4))));
    
    for cIdx = 1:length(topoPlotFig_chanIdx)
        nexttile(tldB)
        chanIdx = topoPlotFig_chanIdx(cIdx);
        sphValues = zeros(1, nChan);
        sphValues(chanIdx) = 1;
        sphValues(neighborMatrix(chanIdx, :)) = 2;
        params.projectionType = 'azimuthalEquidistant';
        params.drawLines = true;
        params.lineWidth = 1;
        params.lineCol = [0 0 0 1];
        params.chanMarkerSize = 5;
        params.chanMarkerChar = 'o';
        params.chanLbl = false;
        params.colBar = false;
        params.colMap = [1 1 1; lines(2)];
        dg_sphericalPlot(dg_cfg, sphValues, params)
    end
    
    % Panel A: Distance matrix and histogram
    tldA = tiledlayout(panelA, 'flow');
    
    % Distance matrix heatmap
    nexttile(tldA)
    if isfield(sphCoord, 'labels')
        chanLabels = sphCoord.labels;
    else
        chanLabels = arrayfun(@(i) sprintf('Ch%d', i), 1:nChan, 'UniformOutput', false);
    end
    mt = imagesc(1:nChan, 1:nChan, sphericalDistanceMatrix);
    mt.Parent.YTick = 1:nChan;
    mt.Parent.YTickLabel = string(chanLabels);
    mt.Parent.YLabel.String = 'channel labels';
    mt.Parent.XTick = 1:nChan;
    mt.Parent.XTickLabel = string(chanLabels);
    mt.Parent.XLabel.String = 'channel labels';
    mt.Parent.Colormap = flipud(copper);
    cb = colorbar;
    cb.Label.String = 'angular distance [rad]';
    
    % Histogram of pairwise distances
    nexttile(tldA)
    counts = tril(sphericalDistanceMatrix);
    counts(counts == 0) = NaN;
    hst = histogram(counts(:), 'BinWidth', 0.1);
    hst.Parent.YLabel.String = {'num of channel pairs' 'at a certain distance'};
    hst.Parent.XLabel.String = 'angular distance [rad]';
end

end
