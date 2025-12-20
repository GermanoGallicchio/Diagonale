function results = dg_describeClusters(dg_cfg,results)

%% shortcuts

[y1Matrix_units, x2Matrix_units, z3Matrix_lbl ] = ndgrid(dg_cfg.dimensions.y1_vec,  dg_cfg.dimensions.x2_vec,  dg_cfg.dimensions.z3_chanLbl);

y1Lbl = dg_cfg.dimensions.y1_lbl;
y1Units = dg_cfg.dimensions.y1_units;
x2Lbl = dg_cfg.dimensions.x2_lbl;
x2Units = dg_cfg.dimensions.x2_units;
z3Lbl = dg_cfg.dimensions.z3_lbl;

%% descriptive clusters 
switch [dg_cfg.objective ' & ' dg_cfg.analysis]
    case 'permutationH0testing & empiricalL1_FDR'
        mask = results.resampling.pVal_emp_FDR<dg_cfg.p_crit;  % binarize the pvalue
        mask = mask .* ~dg_cfg.R_ignore;  % apply R_ignore mask

        [clusterMembership, clustIDList, metrics] = dg_clusterForming(dg_cfg, results.statVal_obs , ~mask);
        results.clusters.clusterMembership_obs = clusterMembership;
        results.clusters.clustIDList_obs       = clustIDList;
        results.clusters.metrics_obs           = metrics;

    case 'bootstrapStability & empiricalL1_FDR'
        mask_BRrob = logical(abs(results.resampling.inference.BR_rob) > 2);
        mask_CI    = logical(abs(sum(sign([results.resampling.inference.CIlo; results.resampling.inference.CIup]),1))>0);
        mask_BRrob = mask_BRrob .* ~dg_cfg.R_ignore;
        mask_CI    = mask_CI    .* ~dg_cfg.R_ignore;

        [clusterMembership, clustIDList, metrics] = dg_clusterForming(dg_cfg, results.statVal_obs , ~mask_BRrob);
        results.clusters.clusterMembership_obs = clusterMembership;
        results.clusters.clustIDList_obs       = clustIDList;
        results.clusters.metrics_obs           = metrics;

    case 'permutationH0testing & theoreticalL1_clusterMaxT'
        % clusters already created

    case 'bootstrapStability & theoreticalL1_clusterMaxT'
        warning('double-check this works'); keyboard
        mask_BRrob = logical(abs(results.lvl1.inference_maxT.BR_rob) > 2);
        mask_CI    = logical(abs(sum(sign([results.lvl1.inference_maxT.CIlo; results.lvl1.inference_maxT.CIup]),1))/2);
        mask_BRrob = mask_BRrob .* ~dg_cfg.R_ignore;
        mask_CI    = mask_CI    .* ~dg_cfg.R_ignore;
        [clusterMembership, clustIDList, metrics] = dg_clusterForming(dg_cfg, results.lvl1.statVal , ~mask_BRrob);
        results.clusters.clusterMembership = clusterMembership;
        results.clusters.clustIDList       = clustIDList;
        results.clusters.metrics           = metrics;

    case 'permutationH0testing & PLS_SVD'
        % nothing to do
    case 'bootstrapStability & PLS_SVD'
        error('not coded yet')
end

end
