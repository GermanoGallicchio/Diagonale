function results = dg_maxT(dg_cfg,results)

nIterations = dg_cfg.analysis.nIterations;
metrics = results.resampling.metrics;
p_crit = dg_cfg.p_crit;

metrics_lbl = fieldnames(results.resampling.metrics);
metrics_lbl = metrics_lbl(~strcmp(metrics_lbl,'id'));
metrics_num = length(metrics_lbl);
switch dg_cfg.analysis
    case 'theoreticalL1_clusterMaxT'
        nClust = length(metrics(end).(metrics_lbl{1}));
        r = nClust;
    case 'PLS_SVD'
        nModes = results.lvl2.nModes;
        r = nModes;
end

metrics_maxT = metrics;
for itIdx = 1:nIterations
    for metIdx = 1:metrics_num
        metricVal = [metrics(1,itIdx).(metrics_lbl{metIdx})]; 
        [maxVal, maxIdx] = max(abs(metricVal));
        if maxVal==0; keyboard; end
        if ~isempty(maxVal)
            metrics_maxT(1,itIdx).(metrics_lbl{metIdx}) = metricVal(maxIdx);
            if any(strcmp(fieldnames(metrics_maxT),'id'))  &&  metIdx==1
                metrics_maxT(1,itIdx).id = metrics_maxT(1,itIdx).id(maxIdx);
            end
        else
            metrics_maxT(1,itIdx).(metrics_lbl{metIdx}) = 0;
            if any(strcmp(fieldnames(metrics_maxT),'id'))  &&  metIdx==1
                metrics_maxT(1,itIdx).id = NaN;
            end
        end
    end
end

for metIdx = 1:metrics_num
    varLbl = metrics_lbl{metIdx};
    inference_maxT.thresholds.(varLbl) = prctile(abs([metrics_maxT.(metrics_lbl{metIdx})]),100-p_crit *100);
end

for msIdx = 1:metrics_num
    H0distribution = [metrics_maxT.(metrics_lbl{msIdx})];
    obsVal = metrics(1).(metrics_lbl{msIdx});
    pval_maxT = nan(1,length(obsVal));
    for nIdx = 1:length(obsVal)
        pval_maxT(1,nIdx) = sum(abs(H0distribution)>=abs(obsVal(nIdx))) / length(H0distribution);
    end
    inference_maxT.pval.(metrics_lbl{msIdx}) = pval_maxT;
end

if dg_cfg.figFlag  &&  nIterations>2
    figure(); clf;
    fig = gcf; fig.Units = "normalized"; fig.Position = [0.05 0.2 0.9 0.6];
    tld = tiledlayout('flow');
    yLbl = ['num of simulated samples (out of ' num2str(nIterations) ' iterations)'];
    tld.YLabel.String = yLbl;
    tld.Title.String = 'maxT null distribution';
    for metIdx = 1:metrics_num
        nexttile(tld);
        histogram(abs([metrics_maxT.(metrics_lbl{metIdx})]), 'FaceColor', [1 0 0],'FaceAlpha',0.5);
        xlabel(metrics_lbl{metIdx})
        hold on
        xline(inference_maxT.thresholds.([metrics_lbl{metIdx}]), 'Color', [1 0 0], 'LineStyle', '-', 'LineWidth', 2);
        if ~isempty(metrics(1).id)
            xline([metrics(1).(metrics_lbl{metIdx})],'LineStyle','-.','Color',[0.5 0.5 0],'LineWidth',0.5);
        else
            xline(0,'LineStyle','-.','Color',[0.5 0.5 0],'LineWidth',2);
        end
        yyaxis right
        nVals = length([metrics_maxT.(metrics_lbl{metIdx})]);
        plot([metrics_maxT.(metrics_lbl{metIdx})],randn(1,nVals)+ones(1,nVals),'x','Color',[0 0 1]);
        set(gca,'YTick',[])
        legend({metrics_lbl{metIdx} [metrics_lbl{metIdx} '_{threshold}'] [metrics_lbl{metIdx} '_{observed}']},'Location','Best')
    end
end

results.clusters.inference_maxT = inference_maxT;

end
