function results = di_clusterReport(di_cfg, results)
% Generate cluster descriptive statistics and summary report
% Takes the raw cluster metrics from di_clusterForming and synthesizes them
% into a report with dimensional coordinates, variability metrics, and summary table.
%
% INPUT:
%   di_cfg   - validated analysis configuration
%   results  - struct with .clusters.metrics_obs (from di_clusterForming)
%
% OUTPUT:
%   results  - augmented with .clusters.report containing:
%       .statLabel           - label of descriptive metric used in the report ('r', 't', or generic 'statVal' most likely OLS's beta)
%       .statVal_median       - median statistic per cluster (of the metric descriebd in the previous field)
%       .statVal_iqr          - interquartile range per cluster
%       .statVal_mad          - median absolute deviation per cluster
%       .mostExtreme_idx      - index of peak/trough per cluster
%       .mostExtreme_val      - value at peak/trough per cluster
%       .mostExtreme_coord    - dimensional coordinates of peak/trough (struct with dim names)
%       .leastExtreme_idx     - index of weakest point per cluster
%       .leastExtreme_val     - value at weakest point per cluster
%       .leastExtreme_coord   - dimensional coordinates of weakest point (struct with dim names)
%       .extent_summary       - min/max ranges per dimension (struct with dim names)
%       .summaryTable        - table with one row per cluster, all metrics
%       .pVal_*               - maxT p-values per metric when results.inference.cluster is present
%       .thresholds           - maxT thresholds per metric when results.inference.cluster is present
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

% TO DO: include (at descriptive level) measures of effect size fo some sort (e.g., Cohen's d, Pearson's r) (medium-low priority)
% TO DO: include empirical pvalues and thresholds (low priority // just for curiosity)

%% sanity checks

if ~isfield(results, 'observed') || ~isfield(results.observed, 'clusters')
    error('results.observed.clusters not found. Call di_clusterDescribe first.');
end
if ~isfield(results.observed.clusters, 'metrics_obs')
    error('results.observed.clusters.metrics_obs not found. Call di_clusterDescribe first.');
end

%% shortcuts

dimKeys = fieldnames(di_cfg.dimensions);
dimTypes = cellfun(@(k) di_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);
dimLbls = cellfun(@(k) di_cfg.dimensions.(k).lbl, dimKeys, 'UniformOutput', false);
dimUnits = cellfun(@(k) di_cfg.dimensions.(k).units, dimKeys, 'UniformOutput', false);

% get cluster information
clusterMetrics = results.observed.clusters.metrics_obs;
if isfield(results.observed.clusters, 'statVal_display')
    statVal = results.observed.clusters.statVal_display;
else
    statVal = results.observed.statVal;
end
if isfield(results.observed.clusters, 'statLabel')
    statLabel = results.observed.clusters.statLabel;
else
    statLabel = 'statVal';
end
nClusters = numel(clusterMetrics.id);

% identify continuous and spherical dimensions
contDimIdx = find(strcmp(dimTypes, 'continuous'));
sphDimIdx = find(strcmp(dimTypes, 'spherical'));

%% report structure

% initialization
report = struct();
report.statLabel = statLabel;

%% extract and enhance metrics

% some basic metrics are already in clusterMetrics
report.statVal_median = clusterMetrics.medianVal;
report.statVal_mass   = clusterMetrics.mass;
report.size           = clusterMetrics.size;
report.mostExtreme_val = clusterMetrics.mostExtremeVal;
report.mostExtreme_idx = clusterMetrics.mostExtremeIdx;
report.leastExtreme_val = clusterMetrics.leastExtremeVal;
report.leastExtreme_idx = clusterMetrics.leastExtremeIdx;


% coordinates already available in clusterMetrics
if nClusters > 0  % Only process if clusters exist
    if ~isempty(contDimIdx)
        for dimIdx = contDimIdx
            dimKey = dimKeys{dimIdx};
            report.mostExtreme_coord.(dimKey) = clusterMetrics.mostExtreme_coord.(dimKey);
            report.leastExtreme_coord.(dimKey) = clusterMetrics.leastExtreme_coord.(dimKey);
            
            % Extent (already in clusterMetrics)
            report.extent_summary.(dimKey) = clusterMetrics.extent_cont_dim.(dimKey);
        end
    end
end

% channel info already available in clusterMetrics (if spherical dimension exists)
if nClusters > 0 && ~isempty(sphDimIdx)
    sphKey = dimKeys{sphDimIdx};
    report.spherical_channels.count = clusterMetrics.extent_sph_dim.(sphKey).channelCount;
    report.spherical_channels.idx = clusterMetrics.extent_sph_dim.(sphKey).channelIdx;
    report.spherical_channels.label = clusterMetrics.extent_sph_dim.(sphKey).channelLabel;
end

% extract p-values if available (from di_maxT.m for parametricFeature_inferenceCluster analyses)
% p-values are stored per metric (mass, size, extremeVal) in inference results
if isfield(results, 'inference') && isfield(results.inference, 'cluster')
    inference_maxT.pval       = results.inference.cluster.pVal_maxT;
    inference_maxT.thresholds = results.inference.cluster.thresholds;
    % store all available p-values as struct for flexibility
    report.pVal_all = inference_maxT.pval;
    report.thresholds = inference_maxT.thresholds;
end

% Compute IQR and MAD (median absolute deviation) per cluster
report.statVal_iqr = zeros(1, nClusters);
report.statVal_mad = zeros(1, nClusters);
for clIdx = 1:nClusters
    pointIdx = clusterMetrics.pointIdx{clIdx};
    clusterStats = statVal(pointIdx);
    report.statVal_iqr(clIdx) = iqr(clusterStats);
    report.statVal_mad(clIdx) = mad(clusterStats, 1); % flag=1 uses median (not mean)
end

%% create summary table

% Build a struct array then convert to table for easy viewing
summaryStruct = struct();
for clIdx = 1:nClusters
    
    % cluster ID
    summaryStruct(clIdx).cluster_id = clusterMetrics.id(clIdx);

    % p-values (if available from maxT inference) - add all available p-values dynamically
    if isfield(report, 'pVal_all')
        pval_metrics = fieldnames(report.pVal_all);
        for pIdx = 1:length(pval_metrics)
            metricName = pval_metrics{pIdx};
            summaryStruct(clIdx).(sprintf('pVal_%s', metricName)) = report.pVal_all.(metricName)(clIdx);
        end
    end
    
    % thresholds (if available from maxT inference) - add all available thresholds dynamically
    if isfield(report, 'thresholds')
        threshold_metrics = fieldnames(report.thresholds);
        for tIdx = 1:length(threshold_metrics)
            metricName = threshold_metrics{tIdx};
            summaryStruct(clIdx).(sprintf('threshold_%s', metricName)) = report.thresholds.(metricName);
        end
    end

    summaryStruct(clIdx).size = clusterMetrics.size(clIdx);
    summaryStruct(clIdx).mass = clusterMetrics.mass(clIdx);
    summaryStruct(clIdx).mostExtreme_val = report.mostExtreme_val(clIdx);
    summaryStruct(clIdx).leastExtreme_val = report.leastExtreme_val(clIdx);

    % based on the descriptive cluster metric chosen upstream
    summaryStruct(clIdx).median = report.statVal_median(clIdx);
    summaryStruct(clIdx).iqr = report.statVal_iqr(clIdx);
    summaryStruct(clIdx).mad = report.statVal_mad(clIdx);
    
    % Dimensional coordinates (construct field names dynamically)
    for dimIdx = contDimIdx
        dimKey = dimKeys{dimIdx};
        dimLbl = dimLbls{dimIdx};
        dimUnit = dimUnits{dimIdx};
        
        % Most extreme location
        fieldname_mostExtreme = sprintf('%s_mostExtreme_%s', dimLbl, dimUnit);
        summaryStruct(clIdx).(fieldname_mostExtreme) = report.mostExtreme_coord.(dimKey).val(clIdx);
        
        % Least extreme location
        fieldname_leastExtreme = sprintf('%s_leastExtreme_%s', dimLbl, dimUnit);
        summaryStruct(clIdx).(fieldname_leastExtreme) = report.leastExtreme_coord.(dimKey).val(clIdx);
        
        % Extent range
        fieldname_range = sprintf('%s_range_%s', dimLbl, dimUnit);
        summaryStruct(clIdx).(fieldname_range) = ...
            sprintf('[%.2f, %.2f]', report.extent_summary.(dimKey).minVal(clIdx), ...
                                    report.extent_summary.(dimKey).maxVal(clIdx));
    end
    
    % Spherical dimension
    if ~isempty(sphDimIdx)
        sphKey = dimKeys{sphDimIdx};
        sphLbl = dimLbls{sphDimIdx};
        summaryStruct(clIdx).(sprintf('%s_count', sphLbl)) = report.spherical_channels.count(clIdx);
        summaryStruct(clIdx).(sprintf('%s_idx', sphLbl)) = join(string(report.spherical_channels.idx{clIdx}));
        summaryStruct(clIdx).(sprintf('%s_label', sphLbl)) = join(report.spherical_channels.label{clIdx});
        % Most extreme channel (if available)
        if isfield(clusterMetrics.mostExtreme_coord, sphKey)
            summaryStruct(clIdx).(sprintf('%s_mostExtreme_idx', sphLbl)) = clusterMetrics.mostExtreme_coord.(sphKey).idx(clIdx);
            if isfield(clusterMetrics.mostExtreme_coord.(sphKey), 'label')
                summaryStruct(clIdx).(sprintf('%s_mostExtreme_label', sphLbl)) = clusterMetrics.mostExtreme_coord.(sphKey).label{clIdx};
            end
        end
        % Least extreme channel (if available)
        if isfield(clusterMetrics.leastExtreme_coord, sphKey)
            summaryStruct(clIdx).(sprintf('%s_leastExtreme_idx', sphLbl)) = clusterMetrics.leastExtreme_coord.(sphKey).idx(clIdx);
            if isfield(clusterMetrics.leastExtreme_coord.(sphKey), 'label')
                summaryStruct(clIdx).(sprintf('%s_leastExtreme_label', sphLbl)) = clusterMetrics.leastExtreme_coord.(sphKey).label{clIdx};
            end
        end
    end
end

% Convert to table and sort by mass (descending)
report.summaryTable = struct2table(summaryStruct);

if ~isempty(summaryStruct)
    if height(report.summaryTable)~=0
        %[~, sortIdx] = sort(clusterMetrics.mass, 'descend');
        [~, sortIdx] = sort(abs(report.summaryTable.mass), 'descend');
        report.summaryTable = report.summaryTable(sortIdx, :);
    end
end

%% store in results

results.observed.clusters.report = report;

end
