function [summaryLabel, testLabel] = di_analysisSummary(di_cfg)
% Build the human-readable identity of the analysis about to run.
% Called by di_analyze AFTER validation (designCode and resolved OLS options
% must exist), printed once per run and stamped into results.analysis.
%
% OUTPUT:
%   summaryLabel - one line: test | design | inference machinery
%   testLabel    - short test name, e.g. 'Welch two-sample t-test
%                  (Satterthwaite df)' - quote it in a methods section
%
% The test name is DERIVED from the data structure and options, so this line
% is the ground truth of what Diagonale is doing - if it does not match the
% analysis you intended, check dataStruct (see di_testHelp).
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

a = di_cfg.analysis;
isLatent = strcmp(a.inferenceLevel, 'latent');

%% test identity
if isLatent
    testLabel = 'PLS-SVD (latent modes)';
else
    ranked = isfield(a, 'OLS') && isfield(a.OLS, 'ranked') && a.OLS.ranked;
    if isequal(a.designCode, [0 0])
        if ranked
            testLabel = 'Spearman correlation';
        else
            testLabel = 'Pearson correlation';
        end
    elseif isequal(a.designCode, [0 1])
        if strcmp(a.OLS.varianceType, 'equal')
            testLabel = 'Student two-sample t-test';
        else
            testLabel = 'Welch two-sample t-test (Satterthwaite df)';
        end
    elseif isequal(a.designCode, [1 0])
        testLabel = 'paired-samples t-test (subject fixed effects)';
    else
        testLabel = 'mixed design (not supported)';
    end
end

%% resampled statistic
if isLatent
    if strcmp(a.objective, 'permutationH0testing')
        statStr = 'singular values (s, wilk)';
    else
        statStr = 'loadings (U, V)';
    end
elseif isfield(a, 'OLS') && isfield(a.OLS, 'inferenceStatistic_resolved')
    statStr = a.OLS.inferenceStatistic_resolved;
else
    statStr = 'unresolved';
end

%% objective block
switch a.objective
    case 'permutationH0testing'
        corrStr = a.correction;
        if strcmp(a.correction, 'maxT')
            if isfield(a, 'maxT') && isfield(a.maxT, 'stepDown') && a.maxT.stepDown
                corrStr = 'maxT (step-down)';
            else
                corrStr = 'maxT (single-step)';
            end
        end
        objBlock = sprintf('permutation H0 test (%s, %d iterations) | statistic: %s | correction: %s | alpha = %g', ...
            a.permuteUnit, a.nIterations, statStr, corrStr, a.p_crit);
    case 'bootstrapStability'
        nGroups = numel(unique(a.dataStruct.indFactor1));
        if nGroups > 1
            stratStr = sprintf('stratified by group (%d groups)', nGroups);
        else
            stratStr = 'whole-sample';
        end
        objBlock = sprintf(['bootstrap stability estimation (%d resamples, %s) | ' ...
            'resampled statistic: %s | outputs: BR, BR_rob, 95%%CI (descriptive, not a test)'], ...
            a.nIterations, stratStr, statStr);
    otherwise
        objBlock = sprintf('objective: %s', a.objective);
end

%% cluster forming (when applicable)
extraStr = '';
if strcmp(a.inferenceLevel, 'cluster') && isfield(a, 'clusterParams') ...
        && isfield(a.clusterParams, 'clusterFormingPvalThreshold')
    extraStr = sprintf(' | cluster-forming p < %g', a.clusterParams.clusterFormingPvalThreshold);
end

%% feature count
dimKeys = fieldnames(di_cfg.featureDims);
pY = prod(cellfun(@(k) length(di_cfg.featureDims.(k).vec), dimKeys));
if isfield(a, 'ignore_col') && any(a.ignore_col)
    featStr = sprintf('%d features (%d ignored)', pY, nnz(a.ignore_col));
else
    featStr = sprintf('%d features', pY);
end

summaryLabel = sprintf('%s | design [%d %d] | %s | %s%s', ...
    testLabel, a.designCode(1), a.designCode(2), featStr, objBlock, extraStr);

end
