function results = di_empiricalP(di_cfg, results)
% Compute empirical p-values from permutation distributions (uncorrected)
% Computes uncorrected empirical p-values by comparing each observed value
% against its own null distribution (without maxT correction).
%   - Feature-level: each feature vs. its null distribution
%   - Cluster-level: each cluster metric vs. null distribution of all clusters and all iterations 
%   - Mode-level: each mode vs. its null distribution
%
% ALGORITHM:
%   For each element (feature/cluster/mode):
%     p-value = proportion of |null values| >= |observed value|
%
% INPUT:
%   di_cfg  - configuration with .analysis.nIterations
%   results - contains observed values and null distributions
%
% OUTPUT:
%   results - augmented with .inference.{feature|cluster|mode}.pVal_emp
%
% NOTE: these pvalues make sense probably only for the PLS-SVD case when
% used as testing sequential "significan" (i.e., if mode 1 is significant,
% test mode 2 and so on. it is not valid to check mode n if mode n-1 is not
% significant). it is more of a curiosity for the other analysis types, as this does not have a correction
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% shortcuts

nIterations = di_cfg.analysis.nIterations;

%% implementation

switch di_cfg.analysis.type
    
    case {'empiricalFeature_inferenceFeature'}
keyboard % until here ok but double check this section
        % Feature-level empirical p-values
        if ~isfield(results, 'simulated') || ~isfield(results.simulated, 'permutationH0') || ~isfield(results.simulated.permutationH0, 'statVal')
            error('empiricalFeature_inferenceFeature: missing results.simulated.permutationH0.statVal')
        end
        
        statVal_obs = results.observed.statVal;  % (1 x pX)
        statVal_null = results.simulated.permutationH0.statVal;  % (nIterations x pX)
        pX = size(statVal_obs, 2);
        
        pVal_emp = nan(1, pX);
        for featIdx = 1:pX
            if di_cfg.analysis.ignore_col(featIdx)
                continue
            end
            % Empirical p-value: proportion of |null| >= |observed|
            pVal_emp(featIdx) = sum(abs(statVal_null(:, featIdx)) >= abs(statVal_obs(featIdx))) / nIterations;
        end
        
        results.inference.feature.pVal_emp = pVal_emp;
        
    case {'parametricFeature_inferenceFeature'}
keyboard % for debugging  // this path is not yet validated
        % Feature-level empirical p-values (from permutation null)
        if ~isfield(results, 'simulated') || ~isfield(results.simulated, 'permutationH0') || ~isfield(results.simulated.permutationH0, 'statVal')
            error('parametricFeature_inferenceFeature: missing results.simulated.permutationH0.statVal')
        end
        
        statVal_obs = results.observed.statVal;  % (1 x pX)
        statVal_null = results.simulated.permutationH0.statVal;  % (nIterations x pX)
        pX = size(statVal_obs, 2);
        
        pVal_emp = nan(1, pX);
        for featIdx = 1:pX
            if di_cfg.analysis.ignore_col(featIdx)
                continue
            end
            pVal_emp(featIdx) = sum(abs(statVal_null(:, featIdx)) >= abs(statVal_obs(featIdx))) / nIterations;
        end
        
        results.inference.feature.pVal_emp = pVal_emp;
        
    case {'parametricFeature_inferenceCluster'}
keyboard % for debugging  // this path is not yet validated -- in fact I don't see any reason for using it... I might integrate it for completeness and curiosity
        % TO DO: integrate this branch (low priority)
        % Cluster-level empirical p-values
        if ~isfield(results, 'simulated') || ~isfield(results.simulated, 'permutationH0') || ~isfield(results.simulated.permutationH0, 'clusterMetrics')
            error('parametricFeature_inferenceCluster: missing results.simulated.permutationH0.clusterMetrics')
        end
        
        metrics_obs = results.observed.clusters.metrics_obs;
        metrics_null = results.simulated.permutationH0.clusterMetrics;  % (1 x nIterations) struct array
        
        % Get metric names (exclude 'id' field)
        metrics_lbl = fieldnames(metrics_obs);
        metrics_lbl = metrics_lbl(~strcmp(metrics_lbl, 'id'));
        
        % Compute empirical p-values for each cluster metric
        pVal_emp = struct();
        for metIdx = 1:length(metrics_lbl)
            metricName = metrics_lbl{metIdx};
            obsVals = [metrics_obs.(metricName)];  % observed cluster metrics
            nClusters = length(obsVals);
            
            % For each observed cluster, compare against null distribution
            pVals = nan(1, nClusters);
            for cIdx = 1:nClusters
                % Build null distribution for this metric (all clusters, all iterations)
                nullVals = [];
                for itIdx = 1:nIterations
                    if ~isempty(metrics_null(itIdx).(metricName))
                        nullVals = [nullVals, metrics_null(itIdx).(metricName)];
                    end
                end
                % Empirical p-value
                if ~isempty(nullVals)
                    pVals(cIdx) = sum(abs(nullVals) >= abs(obsVals(cIdx))) / length(nullVals);
                end
            end
            pVal_emp.(metricName) = pVals;
        end
        
        results.inference.cluster.pVal_emp = pVal_emp;
        
    case {'PLS_SVD'}
        % Mode-level empirical p-values (uncorrected, per-mode inference)
        % SVD returns singular values sorted in descending order by magnitude.
        % the nth position in the sorted SV sequence represents the  nth-strongest 
        % latent structure across all permutations. We compare magnitudes (singular values) only

        % check simulated metrics exist
        if ~isfield(results, 'simulated') || ~isfield(results.simulated, 'permutationH0') || ~isfield(results.simulated.permutationH0, 'modes')
            error('\\ PLS_SVD: missing results.simulated.permutationH0.modes')
        end
        % check observed metrics exist
        if ~isfield(results, 'PLS_SVD') || ~isfield(results.PLS_SVD, 'modes')
            error('PLS_SVD: missing results.PLS_SVD.modes')
        end
        
        modes_obs = results.PLS_SVD.modes;
        modes_null = results.simulated.permutationH0.modes;
        nModes = modes_obs.nModes;

        % Only singular values are used for confirmatory sequential mode inference.
        if ~isfield(modes_null, 's') || ~isfield(modes_obs, 's_obs')
            error('PLS_SVD: missing singular value fields (modes.s and/or modes.s_obs)')
        end

        obsVals = modes_obs.s_obs;      % (1 x nModes)
        nullVals = modes_null.s;        % (nIterations x nModes)

        % compute empirical p-values per mode for singular values
        pVal_emp = struct();
        pVals = nan(1, nModes);
        for modeIdx = 1:nModes
            % Empirical p-value: proportion of |null| >= |observed| for
            % this mode; however, all values are strictly positive anyway
            % Because SVs are sorted descending, modeIdx=1 corresponds to
            % "strongest latent structure", modeIdx=2 to "2nd strongest", etc.
            pVals(modeIdx) = sum(abs(nullVals(:, modeIdx)) >= abs(obsVals(modeIdx))) / nIterations;
        end
        pVal_emp.s = pVals;
        
        results.inference.mode.pVal_emp = pVal_emp;
        
    case {'AJIVE'}
        error('\\ not yet implemented')
        
    otherwise
        error('\\ Analysis type not supported: %s', di_cfg.analysis.type)
end

%% sanity check figure 
% for cluster
% TO DO (low priority)

%% sanity check figure 
% for PLS_SVD mode-level inference

if strcmpi(di_cfg.analysis.type, 'PLS_SVD') && di_cfg.analysis.figFlag && nIterations > 2

    % Get metric names and modes info
    modes_obs = results.PLS_SVD.modes;
    modes_null = results.simulated.permutationH0.modes;
    nModes = modes_obs.nModes;
    metrics_lbl = fieldnames(pVal_emp);

    % Cap the number of modes displayed per metric
    capModes = 3;  % maximum modes to show in figure
    nModesToShow = min(nModes, capModes);

    figure(); clf;
    fig = gcf; fig.Units = "normalized"; fig.Position = [0.1 0.1 0.8 0.7];
    nMetrics = length(metrics_lbl);
    tld = tiledlayout(nModesToShow, nMetrics);
    yLbl = ['num of simulated samples (out of ' num2str(nIterations) ' iterations)'];
    tld.YLabel.String = yLbl;
    tld.Title.String = ['Empirical p-values per mode (max ' num2str(capModes) ' modes)'];

    for metIdx = 1:nMetrics
        metricName = metrics_lbl{metIdx};

        % Get observed and null values
        obsName = [metricName '_obs'];
        obsVals = modes_obs.(obsName);  % (1 x nModes)
        nullVals = modes_null.(metricName);  % (nIterations x nModes)

        % Shared x-axis range across modes (based on shown modes only)
        maxVal = 0;
        for modeIdx = 1:nModesToShow
            maxVal = max(maxVal, max(abs(nullVals(:, modeIdx))));
            maxVal = max(maxVal, abs(obsVals(modeIdx)));
        end
        if maxVal <= 0
            maxVal = 1;
        end

        for modeIdx = 1:nModesToShow
            tileIdx = (modeIdx - 1) * nMetrics + metIdx;
            nexttile(tld, tileIdx);

            % Mode-specific null and threshold (per-mode inference)
            nullAbs = abs(nullVals(:, modeIdx));
            hNull = histogram(nullAbs, ...
                'EdgeColor',[0 0 1],'FaceColor',0.75*[1 1 1],'FaceAlpha',1);
            hold on

            hThresh = [];
            if isfield(di_cfg, 'analysis') && isfield(di_cfg.analysis, 'p_crit')
                thrVal = prctile(nullAbs, 100*(1 - di_cfg.analysis.p_crit));
                hThresh = xline(thrVal, 'Color', [0 0 1], 'LineStyle', '--', 'LineWidth', 2);
            end

            hObs = xline(abs(obsVals(modeIdx)), 'LineStyle','-', 'Color',[1 0.5 0.5], 'LineWidth', 0.8);

            xlim([0 maxVal]);

            xlabel([metricName ' | mode ' num2str(modeIdx)])

            % Legend matching di_maxT style
            legHandles = hNull;
            legEntries = {metricName};
            if ~isempty(hThresh)
                legHandles(end+1) = hThresh;
                legEntries{end+1} = [metricName '_{threshold}'];
            end
            if ~isempty(hObs)
                legHandles(end+1) = hObs(1);
                legEntries{end+1} = [metricName '_{observed}'];
            end
            legend(legHandles, legEntries, 'Location', 'Best')
        end
    end
end

end

