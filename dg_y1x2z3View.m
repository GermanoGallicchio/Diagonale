function dg_y1x2z3View(dg_cfg, results, L, R, viewParams)
% ... Description ...
%
% INPUT:
%           results ... structure variable containing either 
%                       1- results of the dg_stats
%                       or
%                       2- an empty structure (e.g., struct(repmat('ciao',0,0))). 
%                       if empty, plot no statistical values but the
%                       row-averages of R according to a contrast column L
%                       for all y1x2z3 variables in the columns of R
%
% OUTPUT:
%
% Author:


% TO DO: allow any aggregating function to work not just the mean (e.g., std, median)

%% shortcuts

ny1 = dg_cfg.dimensions.y1_num;
nx2 = dg_cfg.dimensions.x2_num;
nz3 = dg_cfg.dimensions.z3_num;
ny1x2 = [ny1 nx2];

chanlocs = dg_cfg.dimensions.z3_chanLocs;

%% decide values to plot
% depending on results content 

xVals = dg_cfg.dimensions.x2_vec;
yVals = dg_cfg.dimensions.y1_vec;



if isempty(results)
	rowIdx = logical(L);
	values2plot = reshape(mean(R(rowIdx,:),1), ny1, nx2, nz3);

else
	values2plot = reshape(results.statVal_obs, ny1, nx2, nz3);
end




%% default parameters (if not entered)

fieldNames = fieldnames(viewParams);

% keep track of what is implicitly used by default and print it is dg_cfg.verbose = true
enteredByDefault = repmat("",0,0); % initialize

% colorMap
if ~any(strcmp(fieldNames,'colorMap'))
	enteredByDefault = [enteredByDefault "colorMap"];
	viewParams.colorMap = turbo;
end

% color limits
if ~any(strcmp(fieldNames,'cLim'))
	enteredByDefault = [enteredByDefault "cLim"];
	if isempty(results)
		viewParams.cLim = prctile(R(:),[1 99]);
	else
		viewParams.cLim = prctile(results.statVal_obs(:),[1 99]);
	end
end

% color bar visibility
if ~any(strcmp(fieldNames,'colorBarVisibility'))
	enteredByDefault = [enteredByDefault "colorBarVisibility"];
	viewParams.colorBarVisibility = true;
end

% color bar position
if viewParams.colorBarVisibility
	if ~any(strcmp(fieldNames,'colorBarPosition'))
		enteredByDefault = [enteredByDefault "colorBarPosition"];
		viewParams.colorBarPosition = [0.750         0.2    0.05    0.6];
	end
end




% color bar tickLength
if ~any(strcmp(fieldNames,'tickLength'))
	enteredByDefault = [enteredByDefault "tickLength"];
	viewParams.tickLength = [0.03    0.03];
end

% contour mask
if ~any(strcmp(fieldNames,'contourMask'))
	enteredByDefault = [enteredByDefault "contourMask"];
	viewParams.contourMask = false;
end

% contour parameters
if viewParams.contourMask
    
	% contour color
	if ~any(strcmp(fieldNames,'contourColor'))
		enteredByDefault = [enteredByDefault "contourColor"];
		%viewParams.contourColor = [0 0 0];
		viewParams.contourColor = viewParams.colorMap(round(prctile(1:size(viewParams.colorMap,1),[75 25])),:);
	end

	% contour line width
	if ~any(strcmp(fieldNames,'contourLineWidth'))
		enteredByDefault = [enteredByDefault "contourLineWidth"];
		viewParams.contourLineWidth = 1.5;
	end
end




% x axis
if ~any(strcmp(fieldNames,'hLim'))
	enteredByDefault = [enteredByDefault "hLim"];
	viewParams.hLim = prctile(dg_cfg.dimensions.x2_vec,[0 100]);
end

if ~any(strcmp(fieldNames,'xLabel'))
	enteredByDefault = [enteredByDefault "xLabel"];
	viewParams.xLabel = [dg_cfg.dimensions.x2_lbl ' ' '[' dg_cfg.dimensions.x2_units ']'];
end

if ~any(strcmp(fieldNames,'xAxis4z3Lbl'))
	enteredByDefault = [enteredByDefault "xAxis4z3Lbl"];
	viewParams.xAxis4z3Lbl = string({dg_cfg.dimensions.z3_chanLocs([1 end]).labels});
end




% y axis

if ~any(strcmp(fieldNames,'vLim'))
	enteredByDefault = [enteredByDefault "vLim"];
	viewParams.vLim = prctile(dg_cfg.dimensions.y1_vec,[0 100]);
end

if ~any(strcmp(fieldNames,'yScale'))
	enteredByDefault = [enteredByDefault "yScale"];
	viewParams.yScale = 'linear';
end

if ~any(strcmp(fieldNames,'yLabel'))
	enteredByDefault = [enteredByDefault "yLabel"];
	viewParams.yLabel = [dg_cfg.dimensions.y1_lbl ' ' '[' dg_cfg.dimensions.y1_units ']'];
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




if dg_cfg.verbose  &&  ~isempty(enteredByDefault)
	disp('Note: the following fields of viewParams were added by default')
	disp(enteredByDefault)
	disp('Check below the "viewParams" used by the code to see what fields can be changed')
	disp(viewParams)
end


%% contour mask
if viewParams.contourMask
	switch [dg_cfg.objective ' & ' dg_cfg.analysis]
		case 'permutationH0testing & empiricalL1_FDR'
			mask(:,:,:,1) = reshape(results.resampling.pVal_emp_FDR<dg_cfg.p_crit & results.statVal_obs>0, ny1, nx2, nz3); % positive
			mask(:,:,:,2) = reshape(results.resampling.pVal_emp_FDR<dg_cfg.p_crit & results.statVal_obs<0, ny1, nx2, nz3); % negative

		case 'permutationH0testing & theoreticalL1_clusterMaxT'
			mask(:,:,:,1) = reshape(results.clusters.clusterMembership_mass_obs_Corrected~=0 & results.statVal_obs>0, ny1, nx2, nz3); % positive
			mask(:,:,:,2) = reshape(results.clusters.clusterMembership_mass_obs_Corrected~=0 & results.statVal_obs<0, ny1, nx2, nz3); % negative


			% --new section-- experimental does not work
			%         statVal_obsMasked = reshape(results.statVal_obs, ny1, nx2, nz3);
			%         masked(:,:,:,1) = statVal_obsMasked .* ~mask(:,:,:,1);
			%         masked(:,:,:,2) = statVal_obsMasked .* ~mask(:,:,:,2);
			%         mask = masked;
			% ---

		case 'bootstrapStability & empiricalL1_FDR'
			mask(:,:,:,1) = reshape(results.resampling.inference.BR_rob>2, ny1, nx2, nz3); % positive
			mask(:,:,:,2) = reshape(results.resampling.inference.BR_rob<-2, ny1, nx2, nz3); % negative

		case 'bootstrapStability & theoreticalL1_clusterMaxT'
			error('not yet coded')

		otherwise
			error('not yet coded')
	end
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


%% figure

for chanIdx = 1:length(chanlocs)

	% create axis set (subplot)
	ax(chanIdx) = axes();
	ax(chanIdx).Units = 'normalized';
	ax(chanIdx).Position = Position(chanIdx,:);
    
	% draw
	im = imagesc(xVals,yVals,values2plot(:,:,chanIdx));

	% colormap
	ax(chanIdx).Colormap = viewParams.colorMap;

	% color limits
	ax(chanIdx).CLim = viewParams.cLim;



	% x-axis
	ax(chanIdx).XLim = viewParams.hLim;
	if any(strcmp(fieldNames,'xTick'))
		xTick = viewParams.xTick;
		ax(chanIdx).XTick = xTick;
	end

	if ismember(chanlocs(chanIdx).labels,viewParams.xAxis4z3Lbl)
		ax(chanIdx).XLabel.String = viewParams.xLabel;
	else
		ax(chanIdx).XTickLabel = [];
	end
	ax(chanIdx).XMinorTick = "off";


	% y-axis
	ax(chanIdx).YLim = viewParams.vLim;
	if any(strcmp(fieldNames,'yTick'))
		yTick = viewParams.yTick;
		ax(chanIdx).YTick = yTick;
	end
	ax(chanIdx).YDir = 'normal';
	ax(chanIdx).YScale = viewParams.yScale;
	if ismember(chanlocs(chanIdx).labels,viewParams.yAxis4z3Lbl)
		ax(chanIdx).YLabel.String = viewParams.yLabel;
	else
		ax(chanIdx).YTickLabel = [];
	end
	ax(chanIdx).YMinorTick = "off";

	% title channel label
	if ismember(chanlocs(chanIdx).labels,viewParams.title4z3Lbl)
		title(chanlocs(chanIdx).labels,'VerticalAlignment','bottom')
	end

	% cone of influence
	if any(strcmp(fieldNames,'coiLeft'))
		hold on
		plot(viewParams.coiLeft, yVals, '-', 'Color',[0 0 0]);
	end
	if any(strcmp(fieldNames,'coiRight'))
		hold on
		plot(viewParams.coiRight, yVals, '-', 'Color',[0 0 0]);
	end




	% contour
	if viewParams.contourMask
		hold on
		for mIdx = 1:size(mask,4)
			if nnz(mask(:,:,chanIdx,mIdx)) > 0
				contour2plot = sum(mask(:,:,chanIdx,mIdx),4);
				col = viewParams.contourColor(mIdx,:);
				col = col.^(1/4);
				[~,ct] = contour(xVals,yVals,contour2plot,1, ...
					'LineWidth',viewParams.contourLineWidth);
				ct.EdgeColor = col;
				ct.ShowText = 'off';

			end
		end
	end

	% misc
	ax(chanIdx).TickLength = viewParams.tickLength;
    
end

% colorbar
if viewParams.colorBarVisibility

	ax2 = axes();
	ax2.Position = [0 0 1 1];
	ax2.Visible = 'off';

	cb = colorbar;
	cb.Limits = [0 1]; 
	cb.Ticks = [0 1];
	cb.TickLabels = viewParams.cLim;

	cb.Position = viewParams.colorBarPosition;

	cb.Parent.Colormap = viewParams.colorMap;


	if any(strcmp(fieldNames,'colorBarLabel'))
		cbTxt = viewParams.colorBarLabel;
		cb.Label.String = cbTxt;
	end

end


%% design label

switch num2str(results.designCode)
	case num2str([1 0 0])
		% lbl corresponding with the group/condition used
		lbl = unique(join(dg_cfg.designTbl{~dg_cfg.row_ignore,~ismember(dg_cfg.designTbl.Properties.VariableNames,'subjID')},', '));
        
	case num2str([0 1 0])
		% lbls corresponding with the two groups
		lp1_lbl = unique(dg_cfg.designTbl.groupID(L==max(L)));
		lp2_lbl = unique(dg_cfg.designTbl.groupID(L==min(L)));
        
	case num2str([0 0 1])
		% lbls corresponding with the two conditions
		rmFactorColIdx = contains(fieldnames(dg_cfg.designTbl),'rmFactor');
		lp1_lbl = unique(join(dg_cfg.designTbl{L==max(L),rmFactorColIdx},', '));
		lp2_lbl = unique(join(dg_cfg.designTbl{L==min(L),rmFactorColIdx},', '));
  
	otherwise
		error('not coded yet')
end

% add panel title
ax3 = axes();
ax3.Position = [0 0 1 1];
ax3.Visible = 'off';
switch num2str(results.designCode)
	case num2str([1 0 0])
		txtLbl = split(lbl);
	case {num2str([0 0 1]) num2str([0 1 0])}
		txtLbl = {lp1_lbl 'vs' lp2_lbl};
end
text(0.01,0.95,txtLbl,'VerticalAlignment','top','HorizontalAlignment','left','FontWeight','normal','FontSize',9)

end
