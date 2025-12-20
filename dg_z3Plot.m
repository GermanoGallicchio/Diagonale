function dg_z3Plot(coord_3d, z3Values, params)

nChan = size(coord_3d,2);
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
    cLim = prctile(z3Values,[0 100]);
end

coord_3d(nChan+1).sph_theta = 90;
coord_3d(nChan+1).sph_phi   = 0;
coord_3d(nChan+1).labels    = 'nasionRef';

coord_2d = dg_z3Projection(coord_3d, projectionType);

coord_2d(:,1) = (coord_2d(:,1)-min(coord_2d(:,1)));
coord_2d(:,2) = (coord_2d(:,2)-min(coord_2d(:,2)));
yRef = coord_2d(nChan+1,2);
xRef = coord_2d(nChan+1,1)*2;
coord_2d(:,1) = coord_2d(:,1) / (xRef-min(coord_2d(:,1)));
coord_2d(:,2) = coord_2d(:,2) / (yRef-min(coord_2d(:,2)));

z3Values_norm = (z3Values - cLim(1)) / (cLim(2) - cLim(1));
idx = round(1 + z3Values_norm * (size(colMap, 1) - 1));
idx = max(1, min(size(colMap,1), idx));

for chanIdx = 1:nChan
    plot(coord_2d(chanIdx,1), coord_2d(chanIdx,2), ...
        chanMarkerChar, ...
        'MarkerSize', chanMarkerSize, ...
        'MarkerFaceColor',colMap(idx(chanIdx),:), ...
        'MarkerEdgeColor',lineCol(1,1:3));
    hold on
end

if chanLbl
    text(coord_2d(1:nChan,1), coord_2d(1:nChan,2),{coord_3d(1:nChan).labels},'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',8)
end

if colBar
    cb = colorbar;
    cb.Ticks = [0 1];
    cb.TickLabels = [ cLim ];
end
