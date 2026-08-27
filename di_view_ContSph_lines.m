function usedViewParams = di_view_ContSph_lines(di_cfg, data, highlightMasks, viewParams)
% viewer for dimensions: continuous x spherical 
% rendered as lines
%
% USAGE:
%   di_view_ContSph_lines(di_cfg, data, highlightMasks, viewParams)
%   usedViewParams = di_view_ContSph_lines(di_cfg, data, highlightMasks, viewParams)
%
% REQUIRED INPUTS:
%   di_cfg.featureDims   - dimension metadata (same schema used in Diagonale)
%   data                - struct array with one element per line and fields:
%                         .name      label used in the legend
%                         .values    [nCont x nSph] values
%                         .color     [1x3] RGB (optional)
%                         .lineWidth scalar (optional)
%   highlightMasks      - optional struct array with fields:
%                         .name   (char/string)
%                         .mask   logical [nTime x nChan]
%                         .color  [1x3] RGB (optional)
%                         .alpha  scalar in [0,1] (optional)
%   viewParams          - optional struct
%
% OUTPUT:
%   usedViewParams      - the viewParams after defaults were applied. 

%% input check

% sanity check: di_cfg is a structure
if ~isstruct(di_cfg)
    error('\\ di_cfg must be a structure')
end

% sanity check: di_cfg.featureDims exists
if ~isfield(di_cfg,'dimensions')
    error('\\ di_cfg.featureDims must exist')
end

% sanity check: data is a structure
if ~isstruct(data)
    error('\\ data must be a structure')
end

% sanity check: highlightMasks is a structure
if ~isstruct(highlightMasks)
    error('\\ highlightMasks must be a structure')
end

nLines = numel(data);
if nLines == 0
    error('\\ data must contain at least one line specification')
end

% sanity check: viewParams is a structure
if ~isstruct(viewParams)
    error('\\ viewParams must be a structure')
end

%% shortcuts (d# order)

dimKeys  = fieldnames(di_cfg.featureDims);
dimVecs  = cellfun(@(k) di_cfg.featureDims.(k).vec, dimKeys, 'UniformOutput', false);
dimLbls  = cellfun(@(k) di_cfg.featureDims.(k).lbl, dimKeys, 'UniformOutput', false);
dimUnits = cellfun(@(k) di_cfg.featureDims.(k).units, dimKeys, 'UniformOutput', false);
dimTypes = cellfun(@(k) di_cfg.featureDims.(k).type, dimKeys, 'UniformOutput', false);
dimSizes = cellfun(@(k) length(di_cfg.featureDims.(k).vec), dimKeys);

contDimIdx = find(strcmp(dimTypes, 'continuous'), 1);
sphIdx = find(strcmp(dimTypes, 'spherical'), 1);

if isempty(contDimIdx)
    error('\\ di_cfg.featureDims requires one continuous dimension')
end
if isempty(sphIdx)
    error('\\ di_cfg.featureDims requires one spherical dimension')
end

xVals = dimVecs{contDimIdx};
nx2 = dimSizes(contDimIdx);
nz3 = dimSizes(sphIdx);
chanLabels = string(di_cfg.featureDims.(dimKeys{sphIdx}).vec);

%% validate line specifications 
% from data structure

for lIdx = 1:nLines

    % name
    if ~isfield(data(lIdx), 'name') || isempty(data(lIdx).name)
        error('\\ data(%d).name is required', lIdx)
    end

    % values
    if ~isfield(data(lIdx), 'values') || isempty(data(lIdx).values)
        error('\\ data(%d).values is required', lIdx)
    end
    lineMat = data(lIdx).values;
    if ~isnumeric(lineMat)
        error('\\ data(%d).values must be numeric', lIdx)
    end
    if ~isequal(size(lineMat), [nx2 nz3])
        error('\\ data(%d).values must be [%d x %d] (continuous x spherical)', lIdx, nx2, nz3)
    end

    % color
    if isfield(data(lIdx), 'color') && ~isempty(data(lIdx).color)
        if ~(isnumeric(data(lIdx).color) && numel(data(lIdx).color) == 3)
            error('\\ data(%d).color must be a numeric RGB triplet', lIdx)
        end
    end

    % lineWidth
    if isfield(data(lIdx), 'lineWidth') && ~isempty(data(lIdx).lineWidth)
        if ~(isnumeric(data(lIdx).lineWidth) && isscalar(data(lIdx).lineWidth))
            error('\\ data(%d).lineWidth must be a numeric scalar', lIdx)
        end
    end

    % lineStyle
    if isfield(data(lIdx), 'lineStyle') && ~isempty(data(lIdx).lineStyle)
        if ~ischar(data(lIdx).lineStyle)
            error('\\ data(%d).lineStyle must be a string / see doc plot for more info', lIdx)
        end
    end
end

%% validate masks
nMasks = numel(highlightMasks);
for mIdx = 1:nMasks
    if ~isfield(highlightMasks(mIdx), 'mask')
        error('\\ highlightMasks(%d) must include field .mask', mIdx)
    end
    if ~islogical(highlightMasks(mIdx).mask)
        error('\\ highlightMasks(%d).mask must be logical', mIdx)
    end
    if ~isequal(size(highlightMasks(mIdx).mask), [nx2 nz3])
        error('\\ highlightMasks(%d).mask must be [%d x %d]', mIdx, nx2, nz3)
    end
end

%% defaults
fieldNames = fieldnames(viewParams);
autoFilledFields = {};

if ~any(strcmp(fieldNames, 'xLabel'))
    viewParams.xLabel = [dimLbls{contDimIdx} ' [' dimUnits{contDimIdx} ']'];
    autoFilledFields{end+1} = 'xLabel';
end
if ~any(strcmp(fieldNames, 'hLim'))
    viewParams.hLim = prctile(xVals, [0 100]);
    autoFilledFields{end+1} = 'hLim';
end
if ~any(strcmp(fieldNames, 'hStep'))
    viewParams.hStep = range(xVals);
    autoFilledFields{end+1} = 'hStep';
end
if ~any(strcmp(fieldNames, 'xAxis4z3Lbl'))
    viewParams.xAxis4z3Lbl = chanLabels([1 end]);
    autoFilledFields{end+1} = 'xAxis4z3Lbl';
end

if ~any(strcmp(fieldNames, 'yLabel'))
    viewParams.yLabel = 'Value';
    autoFilledFields{end+1} = 'yLabel';
end
if ~any(strcmp(fieldNames, 'vLim'))
    allVals = [];
    for lIdx = 1:nLines
        allVals = [allVals; data(lIdx).values(:)];
    end
    vLim = prctile(allVals, [0 100]);
    if any(~isfinite(vLim)) || diff(vLim) == 0
        maxAbs = max(abs(allVals(~isnan(allVals))));
        if isempty(maxAbs) || maxAbs == 0
            maxAbs = 1;
        end
        vLim = [-maxAbs maxAbs];
    end
    viewParams.vLim = vLim;
    autoFilledFields{end+1} = 'vLim';
end
if ~any(strcmp(fieldNames, 'vStep'))
    viewParams.vStep = range(viewParams.vLim);
    autoFilledFields{end+1} = 'vStep';
end
if ~any(strcmp(fieldNames, 'hTickLabelRotation'))
    viewParams.hTickLabelRotation = 90;
    autoFilledFields{end+1} = 'hTickLabelRotation';
end
if ~any(strcmp(fieldNames, 'vTickLabelRotation'))
    viewParams.vTickLabelRotation = 0;
    autoFilledFields{end+1} = 'vTickLabelRotation';
end
if ~any(strcmp(fieldNames, 'yAxis4z3Lbl'))
    viewParams.yAxis4z3Lbl = chanLabels([1 end]);
    autoFilledFields{end+1} = 'yAxis4z3Lbl';
end
if ~any(strcmp(fieldNames, 'title4z3Lbl'))
    viewParams.title4z3Lbl = strings(1, 0);
    autoFilledFields{end+1} = 'title4z3Lbl';
end

if ~any(strcmp(fieldNames, 'xyLength'))
    viewParams.xyLength = [0.155 0.09];
    autoFilledFields{end+1} = 'xyLength';
end
if ~any(strcmp(fieldNames, 'xyPadding'))
    viewParams.xyPadding = [0.030 0.010];
    autoFilledFields{end+1} = 'xyPadding';
end
if ~any(strcmp(fieldNames, 'panelPosition'))
    viewParams.panelPosition = [0 0 1 1];
    autoFilledFields{end+1} = 'panelPosition';
end
if ~any(strcmp(fieldNames, 'projectionType'))
    viewParams.projectionType = 'azimuthalConformal';
    autoFilledFields{end+1} = 'projectionType';
end
if ~any(strcmp(fieldNames, 'projectionWarping'))
    viewParams.projectionWarping = [0.75 1];
    autoFilledFields{end+1} = 'projectionWarping';
end
if ~any(strcmp(fieldNames, 'defaultLineWidth'))
    viewParams.defaultLineWidth = 0.9;
    autoFilledFields{end+1} = 'defaultLineWidth';
end
if ~any(strcmp(fieldNames, 'defaultLineStyle'))
    viewParams.defaultLineStyle = '-';
    autoFilledFields{end+1} = 'defaultLineStyle';
end
if any(strcmp(fieldNames, 'lineColorPalette')) && size(viewParams.lineColorPalette, 1) < nLines
    error('viewParams.lineColorPalette must have at least one row per line')
end
if ~any(strcmp(fieldNames, 'legendVisible'))
    viewParams.legendVisible = true;
    autoFilledFields{end+1} = 'legendVisible';
end
if ~any(strcmp(fieldNames, 'legendPosition'))
    panel = viewParams.panelPosition;
    viewParams.legendPosition = [panel(1) + 0.05 * panel(3) panel(2) + 0.88 * panel(4) 0.2 * panel(3) 0.08 * panel(4)];
    autoFilledFields{end+1} = 'legendPosition';
end
if ~any(strcmp(fieldNames, 'verboseDefaults'))
    viewParams.verboseDefaults = false;
    autoFilledFields{end+1} = 'verboseDefaults';
end

%% defaults for mask style
if any(strcmp(fieldNames, 'lineColorPalette'))
    defaultLineColors = viewParams.lineColorPalette;
else
    defaultLineColors = lines(max(nLines, 2));
end
for lIdx = 1:nLines
    if ~isfield(data(lIdx), 'color') || isempty(data(lIdx).color)
        data(lIdx).color = defaultLineColors(lIdx, :);
    end
    if ~isfield(data(lIdx), 'lineWidth') || isempty(data(lIdx).lineWidth)
        data(lIdx).lineWidth = viewParams.defaultLineWidth;
    end
    if ~isfield(data(lIdx), 'lineStyle') || isempty(data(lIdx).lineStyle)
        data(lIdx).lineStyle = viewParams.defaultLineStyle;
    end
end

for mIdx = 1:nMasks
    if ~isfield(highlightMasks(mIdx), 'name') || isempty(highlightMasks(mIdx).name)
        highlightMasks(mIdx).name = ['mask' num2str(mIdx)];
    end
    if ~isfield(highlightMasks(mIdx), 'color') || isempty(highlightMasks(mIdx).color)
        tmp = lines(max(2, nMasks));
        highlightMasks(mIdx).color = tmp(mIdx, :);
    end
    if ~isfield(highlightMasks(mIdx), 'alpha') || isempty(highlightMasks(mIdx).alpha)
        highlightMasks(mIdx).alpha = 0.25;
    end
end

%% axis (subplot) locations
coord_2d = di_sphericalProjection(di_cfg, viewParams.projectionType);
coord_2d = coord_2d(1:nz3, :);
coord_2d(:,1) = sign(coord_2d(:,1)) .* (abs(coord_2d(:,1)).^(viewParams.projectionWarping(1)));
coord_2d(:,2) = sign(coord_2d(:,2)) .* (abs(coord_2d(:,2)).^(viewParams.projectionWarping(2)));

Position = coord_2d;
panel = viewParams.panelPosition;
posX_min = panel(1) + viewParams.xyPadding(1) * panel(3);
posX_max = panel(1) + panel(3) - viewParams.xyLength(1) * panel(3) - viewParams.xyPadding(1) * panel(3);
posY_min = panel(2) + viewParams.xyPadding(2) * panel(4);
posY_max = panel(2) + panel(4) - viewParams.xyLength(2) * panel(4) - viewParams.xyPadding(2) * panel(4);
Position(:,1) = normalize(Position(:,1), 'range', [posX_min posX_max]);
Position(:,2) = normalize(Position(:,2), 'range', [posY_min posY_max]);
Position = [Position repmat(viewParams.xyLength(1) * panel(3), nz3, 1) repmat(viewParams.xyLength(2) * panel(4), nz3, 1)];

%% draw
legendHandles = gobjects(1, nLines);
legendLabels = strings(1, nLines);

for chanIdx = 1:nz3
    ax(chanIdx) = axes();
    ax(chanIdx).Units = 'normalized';
    ax(chanIdx).Position = Position(chanIdx,:);
    hold(ax(chanIdx), 'on')

    % mask areas (global across lines)
    for mIdx = 1:nMasks
        currentMask = highlightMasks(mIdx).mask(:, chanIdx);
        d_sigIdx = diff([0 currentMask(:)' 0]);
        area_startIdx = find(d_sigIdx == 1);
        area_endIdx = find(d_sigIdx == -1) - 1;
        if numel(xVals) > 1
            step = mean(diff(xVals)) / 2;
        else
            step = 0;
        end
        for islandIdx = 1:length(area_startIdx)
            a = xVals(area_startIdx(islandIdx)) - step;
            b = xVals(area_endIdx(islandIdx)) + step;
            ar = area(ax(chanIdx), [a b], range(viewParams.vLim) * ones(1, 2), viewParams.vLim(1), ...
                'FaceColor', highlightMasks(mIdx).color, ...
                'FaceAlpha', highlightMasks(mIdx).alpha, ...
                'EdgeColor', 'none');
            uistack(ar, 'bottom')
        end
    end

    % data lines
    for lIdx = 1:nLines
        lineData = data(lIdx).values;
        lh = plot(ax(chanIdx), xVals, lineData(:, chanIdx), ...
            'Color', data(lIdx).color, ...
            'LineStyle', data(lIdx).lineStyle, ...
            'LineWidth', data(lIdx).lineWidth);
        if chanIdx == 1
            legendHandles(lIdx) = lh;
            legendLabels(lIdx) = string(data(lIdx).name);
        end
    end

    ax(chanIdx).XLim = viewParams.hLim;
    xTick = unique([fliplr(0:-viewParams.hStep:min(xVals)) 0:viewParams.hStep:max(xVals)], 'stable');
    ax(chanIdx).XTick = xTick;
    ax(chanIdx).XTickLabelRotation = viewParams.hTickLabelRotation;

    ax(chanIdx).YLim = viewParams.vLim;
    yTick = unique([fliplr(0:-viewParams.vStep:min(viewParams.vLim)) 0:viewParams.vStep:max(viewParams.vLim)], 'stable');
    ax(chanIdx).YTick = yTick;
    ax(chanIdx).YTickLabelRotation = viewParams.vTickLabelRotation;

    if ismember(chanLabels(chanIdx), viewParams.xAxis4z3Lbl)
        ax(chanIdx).XLabel.String = viewParams.xLabel;
    else
        ax(chanIdx).XTickLabel = [];
    end
    if ismember(chanLabels(chanIdx), viewParams.yAxis4z3Lbl)
        ax(chanIdx).YLabel.String = viewParams.yLabel;
    else
        ax(chanIdx).YTickLabel = [];
    end
    if ismember(chanLabels(chanIdx), viewParams.title4z3Lbl)
        title(chanLabels(chanIdx), 'VerticalAlignment', 'bottom')
    end

    box(ax(chanIdx), 'off')
    ax(chanIdx).Color = [0 0 0 0];
    ax(chanIdx).TickLength = [0.05 0.025];
end

if viewParams.legendVisible
    ax2 = axes();
    ax2.Position = viewParams.panelPosition;
    ax2.Visible = 'off';
    legend(ax2, legendHandles, legendLabels, 'Position', viewParams.legendPosition, 'Box', 'off')
end



if isfield(viewParams, 'verboseDefaults') && viewParams.verboseDefaults && ~isempty(autoFilledFields)
    disp('di_view_ContSph_lines auto-filled these viewParams fields:')
    disp(autoFilledFields)
end

usedViewParams = viewParams;

end
