function dg_x2z3View(dg_cfg, results, L, R, viewParams)
% ... Description ...
%
% INPUT:
%
%
% OUTPUT:
%
% Author:

%% shortcuts

ny1 = dg_cfg.dimensions.y1_num;
nx2 = dg_cfg.dimensions.x2_num;
nz3 = dg_cfg.dimensions.z3_num;
ny1x2 = [ny1 nx2];

chanlocs = dg_cfg.dimensions.z3_chanLocs;

%% default parameters (if not entered)

fieldNames = fieldnames(viewParams);

% keep track of what is implicitly used by default and print it is dg_cfg.verbose = true
enteredByDefault = repmat("",0,0); % initialize

% line color palette
if ~any(strcmp(fieldNames,'lineColorPalette'))
	enteredByDefault = [enteredByDefault "lineColorPalette"];
	tmpLines = lines(5);
	viewParams.lineColorPalette = tmpLines([4 5],:); % lines/waveforms
end


% x axis
if ~any(strcmp(fieldNames,'xLabel'))
	enteredByDefault = [enteredByDefault "xLabel"];
	viewParams.xLabel = [dg_cfg.dimensions.x2_lbl ' ' '[' dg_cfg.dimensions.x2_units ']'];
end
if ~any(strcmp(fieldNames,'hLim'))
	enteredByDefault = [enteredByDefault "hLim"];
	viewParams.hLim = prctile(dg_cfg.dimensions.x2_vec,[0 100]);
end
if ~any(strcmp(fieldNames,'hStep'))
	enteredByDefault = [enteredByDefault "hStep"];
	viewParams.hStep = range(dg_cfg.dimensions.x2_vec);
end
if ~any(strcmp(fieldNames,'xAxis4z3Lbl'))
	enteredByDefault = [enteredByDefault "xAxis4z3Lbl"];
	viewParams.xAxis4z3Lbl = string({dg_cfg.dimensions.z3_chanLocs([1 end]).labels});
end




% y axis
if ~any(strcmp(fieldNames,'yLabel'))
	enteredByDefault = [enteredByDefault "yLabel"];
	viewParams.yLabel = [dg_cfg.dimensions.y1_lbl ' ' '[' dg_cfg.dimensions.y1_units ']'];
end
if ~any(strcmp(fieldNames,'vLim'))
	enteredByDefault = [enteredByDefault "vLim"];
	viewParams.vLim = prctile(R(:),[0 100]);
end
if ~any(strcmp(fieldNames,'yAxis4z3Lbl'))
	enteredByDefault = [enteredByDefault "yAxis4z3Lbl"];
	viewParams.yAxis4z3Lbl = string({dg_cfg.dimensions.z3_chanLocs([1 end]).labels});
end

% title
if ~any(strcmp(fieldNames,'title4z3Lbl'))
	enteredByDefault = [enteredByDefault "title4z3Lbl"];
	viewParams.title4z3Lbl = string({dg_cfg.dimensions.z3_chanLocs([1 end]).labels});
end

% xyLength and xyPadding
if ~any(strcmp(fieldNames,'xyLength'))
	enteredByDefault = [enteredByDefault "xyLength"];
	viewParams.xyLength  = [0.10  0.10];
end
if ~any(strcmp(fieldNames,'xyPadding'))
	enteredByDefault = [enteredByDefault "xyPadding"];
	viewParams.xyPadding = [0.05  0.05];
end

% projection
if ~any(strcmp(fieldNames,'projectionType'))
	enteredByDefault = [enteredByDefault "projectionType"];
	viewParams.projectionType = 'azimuthalEquidistant';
end
if ~any(strcmp(fieldNames,'projectionWarping'))
	enteredByDefault = [enteredByDefault "projectionWarping"];
	viewParams.projectionWarping = [1 1];
end

% highlight mask
if ~any(strcmp(fieldNames,'highlightMask'))
	enteredByDefault = [enteredByDefault "highlightMask"];
	viewParams.highlightMask = false;
end

% area color palette
if viewParams.highlightMask
if ~any(strcmp(fieldNames,'areaColorPalette'))
	enteredByDefault = [enteredByDefault "areaColorPalette"];
	tmpLines = lines(5);
	viewParams.areaColorPalette = tmpLines([2 1],:); % highlight area
end
end



if ~any(strcmp(fieldNames,'statValAxis'))
	enteredByDefault = [enteredByDefault "statValAxis"];
	viewParams.statValAxis = false;
end
if viewParams.statValAxis
if ~any(strcmp(fieldNames,'statValLim'))
	enteredByDefault = [enteredByDefault "statValLim"];
	viewParams.statValLim = prctile(results.statVal_obs,[0 100]);
end
if ~any(strcmp(fieldNames,'y2Axis4z3Lbl'))
	enteredByDefault = [enteredByDefault "y2Axis4z3Lbl"];
	viewParams.y2Axis4z3Lbl = viewParams.yAxis4z3Lbl;
end
end


if dg_cfg.verbose  &&  ~isempty(enteredByDefault)
	disp('Note: the following fields of viewParams were added by default')
	disp(enteredByDefault)
	disp('Check below the "viewParams" used by the code to see what fields can be changed')
	disp(viewParams)
end

%% decide content of horizontal and vertical axes

% validate dimensions and suggest alternative functions to the user
switch num2str(ny1x2>1)
	case num2str([1   0])
		error('not coded yet / probably I will never enable this feature. simply swap dimensions y1 and x2 earlier and the job is done')

	case num2str([0  1])
		hAxisDim = 'x2';
		xVals = dg_cfg.dimensions.([hAxisDim '_vec']);

	case num2str([1 1])
		error('not coded yet / dimensions y1 and x2 are both non-null. perhaps you want to use y1x2z3View?')

	otherwise
		error('dimensions y1 and x2 are both null. perhaps you want to use z3Plot?')
end



%% highlight area mask
switch [dg_cfg.objective ' & ' dg_cfg.analysis]
	case 'permutationH0testing & empiricalL1_FDR'
		mask(:,:,:,1) = reshape(results.resampling.pVal_emp_FDR<dg_cfg.p_crit & results.statVal_obs>0, ny1, nx2, nz3); % positive
		mask(:,:,:,2) = reshape(results.resampling.pVal_emp_FDR<dg_cfg.p_crit & results.statVal_obs<0, ny1, nx2, nz3); % negative

	case 'permutationH0testing & theoreticalL1_clusterMaxT'
		mask(:,:,:,1) = reshape(results.clusters.clusterMembership_mass_obs_Corrected~=0 & results.statVal_obs>0, ny1, nx2, nz3); % positive
		mask(:,:,:,2) = reshape(results.clusters.clusterMembership_mass_obs_Corrected~=0 & results.statVal_obs<0, ny1, nx2, nz3); % negative

	case 'bootstrapStability & empiricalL1_FDR'
		mask(:,:,:,1) = reshape(results.resampling.inference.BR_rob>2, ny1, nx2, nz3); % positive
		mask(:,:,:,2) = reshape(results.resampling.inference.BR_rob<-2, ny1, nx2, nz3); % negative

	case 'bootstrapStability & theoreticalL1_clusterMaxT'
		error('not yet coded')

	otherwise
		error('not yet coded')
end

%% axis (subplot) locations 

% 3d --> 2d projections
coord_3d = chanlocs;
coord_2d = dg_z3Projection(coord_3d, viewParams.projectionType);

% coordinate warping (nonlinear scaling) to fit the drawing area better
xSign = sign(coord_2d(:,1));
ySign = sign(coord_2d(:,2));
coord_2d(:,1) = xSign.*(abs(coord_2d(:,1)).^(viewParams.projectionWarping(1))); % 0.85
coord_2d(:,2) = ySign.*(abs(coord_2d(:,2)).^(viewParams.projectionWarping(2))); % 1.11

% create position vector
Position = coord_2d;
Position(:,1) = normalize(Position(:,1),'range',[0+viewParams.xyPadding(1) 1-viewParams.xyLength(1)-viewParams.xyPadding(1)]);
Position(:,2) = normalize(Position(:,2),'range',[0+viewParams.xyPadding(2) 1-viewParams.xyLength(2)-viewParams.xyPadding(2)]);
Position = [Position(:,:) repmat(viewParams.xyLength(1),nz3,1) repmat(viewParams.xyLength(2),nz3,1) ];

%% decide values to plot
% depending on designCode (comparing groups/condition vs association between variables)


switch num2str(results.designCode)
	case num2str([1 0 0])
		% idx corresponding with 1st and 4th quartiles
		L_perObs = L(~dg_cfg.row_ignore);
		qIdx = 1*(L<=prctile(L_perObs,25)) + ...
			2*(L<=prctile(L_perObs,50)&L>prctile(L_perObs,25)) + ...
			3*(L<=prctile(L_perObs,75)&L>prctile(L_perObs,50)) + ...
			4*(L>prctile(L_perObs,75));
		qIdx = qIdx.*~dg_cfg.row_ignore; % remove rows to ignore
		% row indices
		lp1_rowIdx = qIdx==1;
		lp2_rowIdx = qIdx==4;
		% associated labels
		lp1_lbl    = "Q1";
		lp2_lbl    = "Q4";
        
	case num2str([0 1 0])
		% idx corresponding with the two groups
		% row indices
		lp1_rowIdx = L==max(L);
		lp2_rowIdx = L==min(L);
		% associated labels
		lp1_lbl = unique(dg_cfg.designTbl.groupID(L==max(L)));
		lp2_lbl = unique(dg_cfg.designTbl.groupID(L==min(L)));
        
	case num2str([0 0 1])
		% idx corresponding with the two conditions
		% row indices
		lp1_rowIdx = L==max(L);
		lp2_rowIdx = L==min(L);
		% associated labels 
		rmFactorColIdx = contains(fieldnames(dg_cfg.designTbl),'rmFactor');
		lp1_lbl = unique(join(dg_cfg.designTbl{L==max(L),rmFactorColIdx},', ',2));
		lp2_lbl = unique(join(dg_cfg.designTbl{L==min(L),rmFactorColIdx},', ',2));
  
	otherwise
		error('not coded yet')
end

%% draw

for chanIdx = 1:length(chanlocs)

	% create axis set (subplot)
	ax(chanIdx) = axes();
	ax(chanIdx).Units = 'normalized';
	ax(chanIdx).Position = Position(chanIdx,:);
    
	% plot data series/waveform of data corresponding with lp1
	rowIdx = lp1_rowIdx;
	val2plot_tmp = reshape(R(rowIdx,:),[sum(rowIdx) ny1 nx2 nz3]);
	val2plot = shiftdim(mean(val2plot_tmp(:,:,:,chanIdx),1),1);
	lp1 = plot(ax(chanIdx), xVals, val2plot); hold on
	lp1.LineWidth = 1;

	% plot data series/waveform of data corresponding with lp2
	rowIdx = lp2_rowIdx;
	val2plot_tmp = reshape(R(rowIdx,:),[sum(rowIdx) ny1 nx2 nz3]);
	val2plot = shiftdim(mean(val2plot_tmp(:,:,:,chanIdx),1),1);
	lp2 = plot(ax(chanIdx), xVals, val2plot);
	lp2.LineWidth = 0.8;

	% colors
	lp1.Color = viewParams.lineColorPalette(1,:);
	lp2.Color = viewParams.lineColorPalette(2,:);


	% x-axis
	ax(chanIdx).XLim = viewParams.hLim;
	xTick = unique([fliplr(0:-viewParams.hStep:min(xVals))  0:viewParams.hStep:max(xVals)],'stable');
	ax(chanIdx).XTick = xTick;
	if ismember(chanlocs(chanIdx).labels,viewParams.xAxis4z3Lbl)
		ax(chanIdx).XLabel.String = viewParams.xLabel;
	else
		ax(chanIdx).XTickLabel = [];
	end

	% y-axis
	ax(chanIdx).YLim = viewParams.vLim;
	if ismember(chanlocs(chanIdx).labels,viewParams.yAxis4z3Lbl)
		ax(chanIdx).YLabel.String = viewParams.yLabel;
	else
		ax(chanIdx).YTickLabel = [];
	end

	% title channel label
	if ismember(chanlocs(chanIdx).labels,viewParams.title4z3Lbl)
		title(chanlocs(chanIdx).labels,'VerticalAlignment','bottom')
	end

	% significance area
	if viewParams.highlightMask 
	for signIdx = 1:size(mask,4)
		mask2plot = mask(1,:,chanIdx,signIdx);

		% semi-transparent areas
		d_sigIdx = diff([0 mask2plot 0]);  % highlight starting and ending points (i.e., moving from 0 to 1 gives "1" (starting), moving from 1 to 0 gives "-1" (ending), everything else is 0 and indicates no changes)
		area_startIdx = find(d_sigIdx==1);
		area_endIdx = find(d_sigIdx==-1)-1;
		step = mean(diff(xVals))/2;
		for islandIdx = 1:length(area_startIdx)
			a = xVals(area_startIdx(islandIdx));
			b = xVals(area_endIdx(islandIdx));
			% add step to allow an area instead of a line
			a = a - step;
			b = b + step;
			% draw
			ar(signIdx) = area(ax(chanIdx),[a b], range(viewParams.vLim)*ones(1,2), viewParams.vLim(1), ...
				'FaceColor' , viewParams.areaColorPalette(signIdx,:), ...
				'FaceAlpha' , 0.25, ...        % 25% opacity
				'EdgeColor' , 'none');         % no edge line
			uistack(ar(signIdx),'bottom');
		end
	end
	end

	% plot (statVal)
	if viewParams.statValAxis 
	yyaxis right
	val2plot_tmp = reshape(results.statVal_obs,[ny1 nx2 nz3]);
	val2plot = shiftdim(mean(val2plot_tmp(:,:,chanIdx),1),1);
	lp3 = plot(ax(chanIdx), xVals, val2plot, 'LineWidth',0.8); hold on
	% colors
	lp3.Color = 0.5*[1 1 1];
	lp3.Color(4) = 0.5;
	lp3.LineWidth = 0.1;
	lp3.LineStyle = '-';
	ax(chanIdx).YAxis(2).Limits = viewParams.statValLim;
	ax(chanIdx).YAxis(2).Color  = lp3.Color;
	% y-axis (right)
	ax(chanIdx).YAxis(2).TickValues = [viewParams.statValLim(1) viewParams.statValLim(2)];
	if ismember(chanlocs(chanIdx).labels,viewParams.y2Axis4z3Lbl)
		ax(chanIdx).YLabel.String = 'statistic score';
	else
		ax(chanIdx).YTickLabel = [];
	end
	end


	box off
	ax(chanIdx).Color = [0 0 0   0];
	ax(chanIdx).TickLength = [0.05 0.025];
    
end

% legend
ax2 = axes();
ax2.Position = [0 0 1 1];
ax2.Visible = 'off';
legend(ax2,[lp1 lp2],[ lp1_lbl lp2_lbl ],'Position',[0.05 0.88 0.125 0.08],'Box','off')

end
