function sphericalDistanceMatrix = dg_z3distance(dg_cfg)

chanlocs = dg_cfg.dimensions.z3_chanLocs;
nz3 = dg_cfg.dimensions.z3_num;

theta = deg2rad([chanlocs(:).sph_theta]);
phi   = deg2rad([chanlocs(:).sph_phi]);

distance_z3_angular = dg_cfg.clusterParams.distance_z3_angular;

sphericalDistanceMatrix = zeros(nz3, nz3);
for chanAIdx = 1:nz3
    for chanBIdx = 1:nz3
        if chanAIdx==chanBIdx
            angDist = 0;
        else
            angDist = acos( sin(phi(chanAIdx))*sin(phi(chanBIdx))  + cos(phi(chanAIdx))*cos(phi(chanBIdx))*cos(theta(chanAIdx)-theta(chanBIdx)) );
        end
        sphericalDistanceMatrix(chanAIdx,chanBIdx) = angDist;
    end
end

angularEdges = 0:distance_z3_angular:max(sphericalDistanceMatrix(:));
if angularEdges(end)~=max(sphericalDistanceMatrix(:));     angularEdges(end+1) = max(sphericalDistanceMatrix(:));  end
Dbinned = -999*ones(size(sphericalDistanceMatrix));
dIdx = sphericalDistanceMatrix==angularEdges(1);
Dbinned(dIdx) = angularEdges(1);
for eIdx = 1:length(angularEdges)-1
    dIdx = sphericalDistanceMatrix>angularEdges(eIdx)  &  sphericalDistanceMatrix<=angularEdges(eIdx+1);
    Dbinned(dIdx) = eIdx;
end
neighborMatrix = Dbinned==1;
neighborNum = sum(neighborMatrix,1);

if dg_cfg.figFlag==1  &&  dg_cfg.dimensions.z3_num>2
    figure(); clf
    f = gcf; f.Units = 'normalized'; f.Position = [0.1 0.01 0.8 0.9];
    panelA = uipanel('Position',[0  0.4  1  0.6]);
    panelB = uipanel('Position',[0  0    0.5  0.4]);
    panelC = uipanel('Position',[0.5  0    0.5  0.4]);

    tldB = tiledlayout(panelB,'flow');
    tldB.Padding = 'tight';
    tldB.TileSpacing = 'tight';
    topoPlotFig_chanIdx = sort(randperm(nz3,max(8,ceil(nz3/4))));
    chanIdx = topoPlotFig_chanIdx;

    plotCounter = 0;
    for cIdx = 1:length(chanIdx)
        plotCounter = plotCounter + 1;
        nexttile(tldB)
        coord_3d = chanlocs;
        z3Values = zeros(1,nz3);
        z3Values(chanIdx(cIdx)) = 1;
        z3Values(neighborMatrix(chanIdx(cIdx),:)) = 2;
        params.projectionType = 'azimuthalEquidistant';
        params.drawLines = true; params.lineWidth = 1; params.lineCol = [0 0 0 1];
        params.chanMarkerSize = 5; params.chanMarkerChar = 'o'; params.chanLbl = false; params.colBar = false;
        params.colMap = [1 1 1; lines(2)];
        dg_z3Plot(coord_3d,z3Values,params)
    end

    tldA = tiledlayout(panelA,'flow');
    nexttile(tldA)
    mt = imagesc(1:nz3, 1:nz3, sphericalDistanceMatrix);
    mt.Parent.YTick = 1:nz3; mt.Parent.YTickLabel = string({chanlocs.labels}); mt.Parent.YLabel.String = 'channel labels';
    mt.Parent.XTick = 1:nz3; mt.Parent.XTickLabel = string({chanlocs.labels}); mt.Parent.XLabel.String = 'channel labels';
    mt.Parent.Colormap = flipud(copper);
    cb = colorbar; cb.Label.String = 'angular distance [rad]';

    nexttile(tldA)
    counts = tril(sphericalDistanceMatrix); counts(counts==0)=NaN;
    hst = histogram(counts(:), 'BinWidth',0.1);
    hst.Parent.YLabel.String = {'num of channel pairs' 'at a certain distance'};
    hst.Parent.XLabel.String = 'angular distance [rad]';
end

end
