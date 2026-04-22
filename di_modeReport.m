function results = di_modeReport(di_cfg, results)
% Generate comprehensive mode-level summary report
%
% INPUT:
%   di_cfg   - validated analysis configuration
%   results  - struct with 
%                 .PLS_SVD.modes with s_obs, p_obs, r_obs, nModes
%                 .inference.mode with pVal_maxT (for permutation testing)
%                 .simulated.permutationH0.modes with s for null comparison
%
% OUTPUT:
%   results     - augmented with .PLS_SVD.modes.report.summaryTable
%
%   Mode                      - mode number (1, 2, 3, ...)
%   SingularValue             - observed singular value
%   VarianceExplained_pct     - proportion of total covariance (%,  individual mode)
%   CumulativeVariance_pct    - cumulative % explained (all modes up to current)
%   CorrelationXY             - correlation between X and Y latent scores
%   pVal_maxT                 - maxT p-values per metric when results.inference.cluster is present
%   pVal_emp                  - empirical p-values (specific for each mode) for sequential evaluation only
%
%
% AUTHOR: Germano Gallicchio (germano.gallicchio@gmail.com)

%% sanity checks

if ~isfield(results, 'PLS_SVD')
    error('results must contain .PLS_SVD from a PLS_SVD analysis');
end

if ~isfield(results.PLS_SVD, 'modes')
    error('results.PLS_SVD must contain .modes');
end

if ~isfield(results.PLS_SVD.modes, 's_obs')
    error('results.PLS_SVD.modes must contain .s_obs (singular values)');
end

if ~isfield(results.PLS_SVD.modes, 'p_obs')
    error('results.PLS_SVD.modes must contain .p_obs (proportion variance explained)');
end

if ~isfield(results.PLS_SVD.modes, 'r_obs')
    error('results.PLS_SVD.modes must contain .r_obs (XY correlations)');
end

if ~isfield(results.PLS_SVD.modes, 'nModes')
    error('results.PLS_SVD.modes must contain .nModes');
end

%% shortcuts

nModes = results.PLS_SVD.modes.nModes;
s_obs = results.PLS_SVD.modes.s_obs;  % singular values: [1 x nModes]
p_obs = results.PLS_SVD.modes.p_obs;  % proportion variance: [1 x nModes]
p_obs_cumsum = cumsum(p_obs);
r_obs = results.PLS_SVD.modes.r_obs;  % XY correlations: [1 x nModes]
pVal_s_maxT = results.inference.mode.pVal_maxT.s;
pVal_s_emp = results.inference.mode.pVal_emp.s;
pVal_wilk_maxT = results.inference.mode.pVal_maxT.wilk;

%% 

report = struct(); % initialization

for modeIdx = 1:nModes
    report(modeIdx).modeId = modeIdx;
    report(modeIdx).s_obs = s_obs(modeIdx);
    report(modeIdx).p_obs = p_obs(modeIdx);
    report(modeIdx).p_obs_cumsum = p_obs_cumsum(modeIdx);
    report(modeIdx).r_obs = r_obs(modeIdx);
    report(modeIdx).pVal_s_maxT = pVal_s_maxT(modeIdx);
    report(modeIdx).pVal_s_emp = pVal_s_emp(modeIdx);
    report(modeIdx).pVal_wilk_maxT = pVal_wilk_maxT(modeIdx);
end

summaryTable = struct2table(report);

%% store in results

results.PLS_SVD.modes.report.summaryTable = summaryTable;

end
