function sphericalDistanceMatrix = di_sphericalDistance(di_cfg)
% Compute distance matrix between all sensors in the spherical dimension of di_cfg.dimensions
% Input: di_cfg - configuration struct with validated dimensions
% Output: sphericalDistanceMatrix - angular distance in radians between all channel pairs

%% check dimensions have been validated

% make sure the dimensions are validated prior to this function
dimensionValidation = false; % initialize it false
if isfield(di_cfg,'validation')
    if isfield(di_cfg.validation,'dimensions')
        dimensionValidation = true;
    end
end

if dimensionValidation==false
    error('\\ di_cfg.dimensions not validated. use: di_cfg = di_validateDimensions(di_cfg)')
end

%% spherical dimension checks

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
sphCoord.labels = di_cfg.dimensions.(sphKey).vec;
nChan = di_cfg.dimensions.(sphKey).num;

%% computations

% convert to radians
% TO DO: in future, require radians (not degrees) as input
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
if isfield(di_cfg, 'analysis') && isfield(di_cfg.analysis, 'figFlag') && di_cfg.analysis.figFlag

    if nChan > 2
        % get distance threshold for adjacency
        distance_spherical_angular = di_cfg.analysis.clusterParams.distance_spherical_radians;

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
        panelA = uipanel('Position', [0    0.4  1    0.6]);
        panelB = uipanel('Position', [0    0    0.5  0.4]);
        panelC = uipanel('Position', [0.5  0    0.5  0.4]);

        % Panel B:

        % choose a random set of channels to plot
        if nChan >= 8
            nRandChans = 8;
        else
            nRandChans = ceil(nChan/4);
        end
        topoPlotFig_chanIdx = sort(randperm(nChan, nRandChans));

        % topographic plot of spherical sites, highlighting neighbors for random seeds
        if nRandChans==8
            nRow = 2;
            nCol = 4;
            tldB = tiledlayout(panelB, nRow, nCol);
        else
            tldB = tiledlayout(panelB, 'flow');
        end
        tldB.Padding = 'tight';
        tldB.TileSpacing = 'tight';

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
            di_sphericalPlot(di_cfg, sphValues, params)

            %ttl = title({['seed channel: ' chanLbl{cIdx}] [num2str(sum(chan2plot)) ' neighbors' ]});
            ttl = title(['seed: ' sphCoord.labels{chanIdx}]);
            ttl.FontWeight = 'normal';
            ttl.VerticalAlignment = 'bottom';
            ttl.FontSize = 10;
            %ttl.Position(2) = 0.55;
        end

        %% Panel A: Distance matrix and histogram
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

        %% Panel C:
        % which threshold yields the most similar number of neighbors across channels

        tldC = tiledlayout(panelC,'flow');
        dVal = logspace(log10(0.25+1),log10(3+1),400)-1;
        MinNumNeigh = nan(size(dVal));
        MaxNumNeigh = nan(size(dVal));
        StdNumNeigh = nan(size(dVal));
        for nIdx = 1:length(dVal)
            MinNumNeigh(nIdx) = min(sum(sphericalDistanceMatrix<=dVal(nIdx),2)-1);
            MaxNumNeigh(nIdx) = max(sum(sphericalDistanceMatrix<=dVal(nIdx),2)-1);
            StdNumNeigh(nIdx) = std(sum(sphericalDistanceMatrix<=dVal(nIdx),2)-1);
        end
        StdNumNeigh = sqrt(StdNumNeigh);

        nexttile(tldC)
        plot(dVal,MinNumNeigh); hold on
        plot(dVal,MaxNumNeigh)

        xTick = logspace(log10(0.5+1),log10(3+1),6)-1;
        set(gca,'XScale','log','XTick',xTick,'XLim',[dVal(1) dVal(end)])
        sgtitle('effect of distance threshold on neighbor count')
        ylabel('number of neighbor per channel')
        xlabel('angular distance threshold [rad]')

        yyaxis right
        plot(dVal,StdNumNeigh)

        xline(di_cfg.analysis.clusterParams.distance_spherical_radians,'Color',0.5*[1 1 1],'LineStyle','--','LineWidth',0.25)
        legend(["minimum" "maximum" "sqrt of st. deviation" num2str(di_cfg.analysis.clusterParams.distance_spherical_radians)])
    end

end

end
