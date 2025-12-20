function results = dg_clusterPruning(dg_cfg,results)

metrics_lbl = fieldnames(results.resampling.metrics);
metrics_lbl = metrics_lbl(~strcmp(metrics_lbl,'id'));
metrics_num = length(metrics_lbl);

for fIdx = 1:metrics_num
    results.clusters.(['clusterMembership_' metrics_lbl{fIdx} '_obs_Corrected']) = results.clusters.clusterMembership_obs;
    results.clusters.(['clustIDList_' metrics_lbl{fIdx} '_obs_Corrected'])       = results.clusters.clustIDList_obs;

    clusters2remove_idx  = abs(results.clusters.metrics_obs.(metrics_lbl{fIdx})) <= results.clusters.inference_maxT.thresholds.(metrics_lbl{fIdx});
    for clIdx = find(clusters2remove_idx)
        idx = results.clusters.clusterMembership_obs==results.clusters.clustIDList_obs(clIdx);
        results.clusters.(['clusterMembership_' metrics_lbl{fIdx} '_obs_Corrected'])(idx) = 0;
        results.clusters.(['clustIDList_' metrics_lbl{fIdx} '_obs_Corrected'])(clIdx) = NaN;
    end
end
