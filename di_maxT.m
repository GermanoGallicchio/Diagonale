function results = di_maxT(di_cfg,results)
% Apply maxT correction for cluster-level or mode-level permutation inference
%
% Implements family-wise error rate (FWER) correction using the maxT approach.
% For each permutation iteration, identifies the most extreme value per each metric (e.g.,
% a cluster metric such as its mass or size or a mode's singular value based metric)
% This builds a null distribution of "worst-case" largest values under the null hypothesis of exchangeability,
% and therefore controlling for multiple comparisons at the cluster/mode level.
% 1. For each null hypothesis iteration, find max(abs(metric)) across all clusters/modes
% 2. Build H0 distribution from these maxima
% 3. Compute p-value: proportion of H0 maxima >= observed metric
% 4. Compute threshold: percentile of H0 distribution at p_crit (not crucial, a bit usesell given the pvalue but effortless)
%
% INPUT:
%   di_cfg  - configuration with:
%             .analysis.type (supported: empiricalFeature_inferenceFeature,
%                              parametricFeature_inferenceFeature,
%                              parametricFeature_inferenceCluster, PLS_SVD)
%             .analysis.nIterations and .analysis.p_crit
%   results - contains permutation null values for the selected analysis type:
%             - feature level: .simulated.permutationH0.statVal
%             - cluster level: .simulated.permutationH0.clusterMetrics
%             - mode level:    .simulated.permutationH0.modes
%
% OUTPUT:
%   results - augmented with one of:
%             .inference.feature.*
%             .inference.cluster.*
%             .inference.mode.*
%             Feature level stores canonical outputs in:
%               .feature.pVal_corr
%               .feature.thresholds.*
%             Cluster/mode levels store metric-specific outputs in:
%               .pVal_maxT.(metric_name)
%               .thresholds.(metric_name)
%
% the metrics depend on the analysis type and are detected from the results structure.
% Clusters
% - mass: sum of statistic values (spatial extent × intensity)
% - size: number of significant points (spatial extent only)
% - extremeVal: peak/trough value (focal intensity only) // work in progress so don't rely on this metric
% Modes:
% - s
% - wilk
% Features:
% - statVal
%
% Note that the observed data is considered like a possible permutation (itIdx==1 is the observed/original)
% so it is included in the null distribution.
%        
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% shortcuts

nIterations = di_cfg.analysis.nIterations;

p_crit = di_cfg.analysis.p_crit;


inferenceLevel = di_cfg.analysis.inferenceLevel;

switch inferenceLevel
    case 'feature'
        % Feature-level maxT: build metrics from permutation null of feature stats
        if ~isfield(results, 'simulated') || ~isfield(results.simulated, 'permutationH0') || ~isfield(results.simulated.permutationH0, 'statVal')
            error('\\ feature-level maxT: missing results.simulated.permutationH0.statVal')
        end
        statVal_null = results.simulated.permutationH0.statVal; % [nIterations x pX]
        pX = size(statVal_null, 2);
        if ~isfield(di_cfg.analysis, 'ignore_col')
            error('\\ feature-level maxT requires di_cfg.analysis.ignore_col')
        end
        ignore_col = logical(di_cfg.analysis.ignore_col);
        if numel(ignore_col) ~= pX
            error('\\ ignore_col length must match number of feature columns in permutationH0.statVal')
        end
        feature_2use = ~ignore_col;
        if ~any(feature_2use)
            error('\\ all features are ignored (ignore_col all true)')
        end
        if size(statVal_null,1) ~= nIterations
            error('\\ statVal null rows must equal nIterations')
        end
        metrics = repmat(struct('statVal', []), 1, nIterations);
        for itIdx = 1:nIterations
            metrics(1,itIdx).statVal = statVal_null(itIdx, :);
        end

    case 'cluster'
        
        % cluster metrics available (mass, size, extremeVal, etc.)
        metrics = results.simulated.permutationH0.clusterMetrics;  % [1 x nIterations] struct array

    case 'latent'
        % Build metrics struct array from mode-level permutation distribution
        % Structure: results.simulated.permutationH0.modes.s (nIterations x nModes)
        if ~isfield(results, 'simulated') || ~isfield(results.simulated, 'permutationH0') || ~isfield(results.simulated.permutationH0, 'modes')
            error('PLS_SVD: missing results.simulated.permutationH0.modes')
        end
        
        modes = results.simulated.permutationH0.modes;
        % sanity check we have the required fields
        if ~isfield(modes, 's') || ~isfield(modes, 'wilk')
            error('PLS_SVD: missing required permutation mode fields (s and/or wilk)')
        end
        
        % Create metrics struct array from ALL iterations (observed + permutations)
        metrics = repmat(struct('s', [], 'wilk', []), 1, nIterations);
        
        for itIdx = 1:nIterations
            metrics(1,itIdx).s    = modes.s(itIdx, :);
            metrics(1,itIdx).wilk = modes.wilk(itIdx, :);
        end

    otherwise
        error('not yet coded')
end


metrics_lbl = fieldnames(metrics);
metrics_lbl = metrics_lbl(~strcmp(metrics_lbl,'id'));  % exclude id field (only existing for clusters)
metrics_num = length(metrics_lbl);

%% build H0 of maximum values (maxT)
% for each iteration (simulated data) select the most extreme 
% value (feature level), cluster metric (cluster level), or singular value
% (mode level). If there are multiple metrics (eg, cluster and mode
% levels), do it separately per metric

metrics_maxT = metrics;  % initialize with full metrics, will be reduced to maxima
for metIdx = 1:metrics_num
    for itIdx = 1:nIterations
        
        % get all feature/cluster/mode level values for this metric in this iteration
        metricVal = [metrics(1,itIdx).(metrics_lbl{metIdx})];  % e.g., [cluster1_mass, cluster2_mass, ...]

        % For feature-level analyses, enforce ignore_col in maxT family definition.
% TO DO: add support for di_cfg.analysis.maxT_dimensions (analogous to di_cfg.analysis.FDR_dimensions)
% so maxT can be applied within user-defined dimension pools instead of as done currently with one global feature family.
        if strcmp(inferenceLevel, 'feature')
            metricVal = metricVal(feature_2use);
        end

        % guard against NaNs in null values (e.g., ignored features from upstream).
        metricVal = metricVal(~isnan(metricVal));
        
        % find the most extreme value
        [maxVal, maxIdx] = max(abs(metricVal));
        
        if ~isempty(maxVal)
            % keep only the most extreme value (preserving sign)
            metrics_maxT(1,itIdx).(metrics_lbl{metIdx}) = metricVal(maxIdx);
            % if cluster IDs exist, store which cluster was the maximum
            if any(strcmp(fieldnames(metrics_maxT),'id'))  &&  metIdx==1
                metrics_maxT(1,itIdx).id = metrics_maxT(1,itIdx).id(maxIdx);
            end
        else % no clusters found in this iteration
            metrics_maxT(1,itIdx).(metrics_lbl{metIdx}) = 0;
            if any(strcmp(fieldnames(metrics_maxT),'id'))  &&  metIdx==1
                metrics_maxT(1,itIdx).id = NaN;
            end
        end
    end
end

%% compute significance thresholds and p-values

thresholds = struct();
pval_maxT = struct();
for metIdx = 1:metrics_num
    varLbl = metrics_lbl{metIdx};
    
    % observed feature/cluster/mode values (iteration 1)
    itIdx = 1;
    obsVal = metrics(itIdx).(varLbl);
    
    % extract H0 distribution: vector of maximum values (abs) across each
    % iteration
    H0distribution = abs([metrics_maxT.(varLbl)]);

    % threshold: (1 - p_crit) percentile of H0
    thresholds.(varLbl) = prctile(H0distribution, 100*(1 - p_crit));

    % p-values: proportion of H0 maxima >= each observed |value|
    if strcmp(inferenceLevel, 'feature')
        pvals = nan(size(obsVal));
        obsActive = obsVal(feature_2use);
        validActive = ~isnan(obsActive);
        pvalsActive = nan(size(obsActive));
        pvalsActive(validActive) = arrayfun(@(v) sum(H0distribution >= abs(v)) / numel(H0distribution), obsActive(validActive));
        pvals(feature_2use) = pvalsActive;
    else
        pvals = arrayfun(@(v) sum(H0distribution >= abs(v)) / numel(H0distribution), obsVal);
    end
    pval_maxT.(varLbl) = pvals;
end

%% store values

switch inferenceLevel
    case 'feature'
        results.inference.feature.pVal_corr = pval_maxT.statVal;
        results.inference.feature.thresholds = thresholds;
    case 'cluster'
        results.inference.cluster.pVal_maxT = pval_maxT;
        results.inference.cluster.thresholds = thresholds;
    case 'latent'
        results.inference.mode.pVal_maxT = pval_maxT;
        results.inference.mode.thresholds = thresholds;
    otherwise
        error('\\ likely a bug. unclear analysis path')
end


%% sanity check figure

if di_cfg.analysis.figFlag  &&  nIterations>2

    figure(); clf;
    fig = gcf; fig.Units = "normalized"; fig.Position = [0.05 0.2 0.9 0.6];
    tld = tiledlayout('flow');
    yLbl = ['num of simulated samples (out of ' num2str(nIterations) ' iterations)'];
    tld.YLabel.String = yLbl;
    tld.Title.String = 'maxT null distribution';
    for metIdx = 1:metrics_num
        nexttile(tld);

        % null distribution (abs value)
        histogram(abs([metrics_maxT.(metrics_lbl{metIdx})]), ...
            'EdgeColor',[0 0 1],'FaceColor', 0.75*[1 1 1],'FaceAlpha',1);
        xlabel(metrics_lbl{metIdx})
        hold on

        % threshold
        xline(thresholds.(metrics_lbl{metIdx}), 'Color', [0 0 1], 'LineStyle', '--', 'LineWidth', 2);

        % observed values (abs value): draw if metric has observed entries
        if ~isempty(metrics(1).(metrics_lbl{metIdx}))
            xline(abs([metrics(1).(metrics_lbl{metIdx})]),'LineStyle','-','Color',[1 0.5 0.5],'LineWidth',0.5);
        else
            xline(0,'LineStyle','-','Color',[0.5 0.5 0.5],'LineWidth',0.5);
        end

        % each single H0 value (abs value)
        yyaxis right
        nVals = length([metrics_maxT.(metrics_lbl{metIdx})]);
        plot(abs([metrics_maxT.(metrics_lbl{metIdx})]),randn(1,nVals)+ones(1,nVals),'x','Color',[0 0 1]);
        set(gca,'YTick',[])

        % legend
        legend({metrics_lbl{metIdx} [metrics_lbl{metIdx} '_{threshold}'] [metrics_lbl{metIdx} '_{observed}']},'Location','Best')
    end
end


end
