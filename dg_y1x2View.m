function dg_y1x2View(dg_cfg, results, L, R, viewParams)

ny1 = dg_cfg.dimensions.y1_num;
nx2 = dg_cfg.dimensions.x2_num;
nz3 = dg_cfg.dimensions.z3_num;
ny1x2 = [ny1 nx2];

chanlocs = dg_cfg.dimensions.z3_chanLocs;

fieldNames = fieldnames(viewParams);
enteredByDefault = repmat("",0,0);

if ~any(strcmp(fieldNames,'colorMap'))
	enteredByDefault = [enteredByDefault "colorMap"];
	viewParams.colorMap = turbo;
end

if ~any(strcmp(fieldNames,'cLim'))
	enteredByDefault = [enteredByDefault "cLim"];
	viewParams.cLim = prctile(results.statVal_obs(:),[1 99]);
end

if ~any(strcmp(fieldNames,'colorBarVisibility'))
	enteredByDefault = [enteredByDefault "colorBarVisibility"];
	viewParams.colorBarVisibility = true;
end

if viewParams.colorBarVisibility
	cTicks = unique([0 viewParams.cLim],'sorted');
	viewParams.cTicks = cTicks;
end

if ~any(strcmp(fieldNames,'contourMask'))
	enteredByDefault = [enteredByDefault "contourMask"];
	viewParams.contourMask = false;
end

if viewParams.contourMask
	if ~any(strcmp(fieldNames,'contourColor'))
		enteredByDefault = [enteredByDefault "contourColor"];
		viewParams.contourColor = viewParams.colorMap(round(prctile(1:size(viewParams.colorMap,1),[75 25])),:);
	end
	if ~any(strcmp(fieldNames,'contourLineWidth'))
		enteredByDefault = [enteredByDefault "contourLineWidth"];
		viewParams.contourLineWidth = 1.5;
	end
end

if ~any(strcmp(fieldNames,'hLim'))
	enteredByDefault = [enteredByDefault "hLim"];
	viewParams.hLim = prctile(dg_cfg.dimensions.x2_vec,[0 100]);
end

if ~any(strcmp(fieldNames,'xAxisVisibility'))
	enteredByDefault = [enteredByDefault "xAxisVisibility"];
	viewParams.xAxisVisibility = 'true';
end

if ~any(strcmp(fieldNames,'xLabel'))
	enteredByDefault = [enteredByDefault "xLabel"];
	viewParams.xLabel = [dg_cfg.dimensions.x2_lbl ' ' '[' dg_cfg.dimensions.x2_units ']'];
end

if ~any(strcmp(fieldNames,'vLim'))
	enteredByDefault = [enteredByDefault "vLim"];
	viewParams.vLim = prctile(dg_cfg.dimensions.y1_vec,[0 100]);
end

if ~any(strcmp(fieldNames,'yScale'))
	enteredByDefault = [enteredByDefault "yScale"];
	viewParams.yScale = 'linear';
end

if ~any(strcmp(fieldNames,'yAxisVisibility'))
	enteredByDefault = [enteredByDefault "yAxisVisibility"];
	viewParams.yAxisVisibility = 'true';
end

if ~any(strcmp(fieldNames,'yLabel'))
	enteredByDefault = [enteredByDefault "yLabel"];
	viewParams.yLabel = [dg_cfg.dimensions.y1_lbl ' ' '[' dg_cfg.dimensions.y1_units ']'];
end

if dg_cfg.verbose  &&  ~isempty(enteredByDefault)
	disp('Note: the following fields of viewParams were added by default')
	disp(enteredByDefault)
	disp('Check below the "viewParams" used by the code to see what fields can be changed')
	disp(viewParams)
end

if viewParams.contourMask
	switch [dg_cfg.objective ' & ' dg_cfg.analysis]
		case 'permutationH0testing & empiricalL1_FDR'
			mask(:,:,:,1) = reshape(results.resampling.pVal_emp_FDR<dg_cfg.p_crit & results.statVal_obs>0, ny1, nx2, nz3);
			mask(:,:,:,2) = reshape(results.resampling.pVal_emp_FDR<dg_cfg.p_crit & results.statVal_obs<0, ny1, nx2, nz3);
		case 'permutationH0testing & theoreticalL1_clusterMaxT'
			mask(:,:,:,1) = reshape(results.clusters.clusterMembership_mass_obs_Corrected~=0 & results.statVal_obs>0, ny1, nx2, nz3);
			mask(:,:,:,2) = reshape(results.clusters.clusterMembership_mass_obs_Corrected~=0 & results.statVal_obs<0, ny1, nx2, nz3);
		case 'bootstrapStability & empiricalL1_FDR'
			mask(:,:,:,1) = reshape(results.resampling.inference.BR_rob>2, ny1, nx2, nz3);
			mask(:,:,:,2) = reshape(results.resampling.inference.BR_rob<-2, ny1, nx2, nz3);
		case 'bootstrapStability & theoreticalL1_clusterMaxT'
			error('not yet coded')
		otherwise
			error('not yet coded')
	end
end

xVals = dg_cfg.dimensions.x2_vec;
yVals = dg_cfg.dimensions.y1_vec;
values2plot = reshape(results.statVal_obs, ny1, nx2, nz3);

im = imagesc(xVals,yVals,values2plot);
im.Parent.Colormap = viewParams.colorMap;
im.Parent.CLim = viewParams.cLim;

im.Parent.XLim = viewParams.hLim;
if any(strcmp(fieldNames,'xTick'))
	xTick = viewParams.xTick;
	im.Parent.XTick = xTick;
end
if viewParams.xAxisVisibility
	im.Parent.XAxis.Label.String = viewParams.xLabel;
else
	im.Parent.XAxis.TickLabels = [];
end
im.Parent.XMinorTick = "off";

im.Parent.YLim = viewParams.vLim;
if any(strcmp(fieldNames,'yTick'))
	yTick = viewParams.yTick;
	im.Parent.YTick = yTick;
end
im.Parent.YDir = 'normal';
im.Parent.YScale = viewParams.yScale;
if viewParams.yAxisVisibility
	im.Parent.YAxis.Label.String = viewParams.yLabel;
else
	im.Parent.YAxis.TickLabels = [];
end
im.Parent.YMinorTick = "off";

if any(strcmp(fieldNames,'title'))
	titleTxt = viewParams.title;
	title(titleTxt);
end

if any(strcmp(fieldNames,'coiLeft'))
	hold on
	plot(viewParams.coiLeft, yVals, '-', 'Color',[0 0 0]);
end
if any(strcmp(fieldNames,'coiRight'))
	hold on
	plot(viewParams.coiRight, yVals, '-', 'Color',[0 0 0]);
end

if viewParams.contourMask
	hold on
	for mIdx = 1:size(mask,4)
		if nnz(mask(:,:,:,mIdx)) > 0
			contour2plot = sum(mask(:,:,:,mIdx),4);
			col = viewParams.contourColor(mIdx,:);
			col = col.^(1/4);
			[~,ct] = contour(xVals,yVals,contour2plot,1, ...
				'LineWidth',viewParams.contourLineWidth);
			ct.EdgeColor = col;
			ct.ShowText = 'off';
		end
	end
end

if viewParams.colorBarVisibility 
	cb = colorbar;
	cb.Limits = viewParams.cLim;
	cb.Ticks = viewParams.cTicks;
	if any(strcmp(fieldNames,'colorBarLabel'))
		cbTxt = viewParams.colorBarLabel;
		cb.Label.String = cbTxt;
	end
end

end
