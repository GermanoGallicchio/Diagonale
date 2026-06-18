function results = di_clusterDescribe(di_cfg,results)
% Form descriptive clusters after feature-level inference 
% has been completed
%
%   Identifies contiguous regions of significant features.
%   This function is called after statistical inference has already taken place(FDR, bootstrap) to
%   group significant points into interpretable clusters. But not a cluster inference itself 
%
% INPUT:
%   di_cfg  - configuration with .analysis.objective and .analysis.inferenceLevel
%   results - contains observed descriptive metrics and inference results
%
% OUTPUT:
%   results - augmented with .observed.clusters fields:
%             .clusterMembership, .clustIDList, .metrics_obs,
%             .statVal_display, .statLabel
%
% NOTE:
%   For *_inferenceCluster modes, inferential clusters are formed upstream,
%   so this function does not form additional descriptive clusters.
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% shortcuts (flexible dimensions in d# order)

% num of dimensions and their numerosity 
% in d# order (fieldnames order)
dimKeys  = fieldnames(di_cfg.dimensions);
dimTypes = cellfun(@(k) di_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);

% TO DO: if this check is done upstream in the validation of dimensions, this can be deleted
if nnz(strcmp(dimTypes, 'spherical')) > 1
    error('only one spherical dimension is supported');
end


% Choose the descriptive metric used to characterize post-hoc clusters.
% Inference remains based on results.observed.statVal upstream/downstream.
if isequal(di_cfg.analysis.designCode, [0 0]) && isfield(results.observed, 'rVal')
    clusterStatVal = results.observed.rVal;
    clusterStatLabel = 'r';
elseif isfield(results.observed, 'tVal')
    clusterStatVal = results.observed.tVal;
    clusterStatLabel = 't';
else
    clusterStatVal = results.observed.statVal;
    clusterStatLabel = 'statVal';
end

results.observed.clusters.statVal_display = clusterStatVal;
results.observed.clusters.statLabel = clusterStatLabel;



%% descriptive clusters 

key = sprintf('%s + %s', di_cfg.analysis.objective, di_cfg.analysis.inferenceLevel);

switch key
    
    case 'permutationH0testing + feature'
        % descriptive cluster forming 
        % Inference was already done at feature level. Clustering is post-hoc for easier description only.
        if ~isfield(results, 'inference') || ~isfield(results.inference, 'feature') || ~isfield(results.inference.feature, 'pVal_corr')
            error('missing results.inference.feature.pVal_corr for descriptive cluster forming')
        end
        pVals2use = results.inference.feature.pVal_corr;
        pVals2use = pVals2use .* ~di_cfg.analysis.ignore_col; % apply ignore mask
        [clusterMembership, clustIDList, metrics] = di_clusterForming(di_cfg, clusterStatVal, pVals2use);
        results.observed.clusters.clusterMembership = clusterMembership;
        results.observed.clusters.clustIDList       = clustIDList;
        results.observed.clusters.metrics_obs       = metrics;

    
    case 'permutationH0testing + cluster'
        % Inferential clusters already formed upstream. Nothing else to do here.

    case 'permutationH0testing + latent'
        % No descriptive cluster-forming step implemented for this mode.
        



    case 'bootstrapStability + feature'
        % Descriptive cluster forming
        % Form clusters from stable features in bootstrap analysis.
        mask_BRrob = logical(abs(results.inference.feature.BR_rob) > 2);
        % apply ignore_col
        mask_BRrob = mask_BRrob .* ~di_cfg.analysis.ignore_col;
        % di_clusterForming expects p-like values and uses p < threshold;
        % invert stability mask so stable points behave like "significant" points.
        [clusterMembership, clustIDList, metrics] = di_clusterForming(di_cfg, clusterStatVal, ~mask_BRrob);
        results.observed.clusters.clusterMembership = clusterMembership;
        results.observed.clusters.clustIDList       = clustIDList;
        results.observed.clusters.metrics_obs       = metrics;
    
    case 'bootstrapStability + cluster'
        % this path does not exist
    
    case 'bootstrapStability + latent'
        % TO DO: if needed, define descriptive clustering on loadings.
end


end
