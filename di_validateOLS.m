function di_cfg = di_validateOLS(di_cfg)
% Validate di_cfg.analysis.OLS
% Runs for permutationH0testing and bootstrapStability (options of di_cfg.analysis.objective), 
% immediately after di_cfg.analysis.designCode is set.
% DIfferently from di_validateH0 (which runs only for permutationH0testing), 
% this validator must always run for OLS-based analyses.
%
% 1. di_cfg.analysis.OLS.varianceType
%    default and design-mismatch warning
% 2. di_cfg.analysis.OLS.inferenceStatistic
%    resolves 'auto' to an objective-dependent statistic label, without
%    mutating the user-facing 'auto' setting (so a saved di_cfg stays
%    portable across objectives). The resolved value is exposed via
%    di_cfg.analysis.OLS.inferenceStatistic_resolved.
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% varianceType
% DELETE THIS COMMENT AT SOME POINT: (this validation was previosly in di_analysis_OLS.m)

if ~isfield(di_cfg.analysis, 'OLS')
    di_cfg.analysis.OLS = struct();
end

if ~isfield(di_cfg.analysis.OLS, 'varianceType')
    di_cfg.analysis.OLS.varianceType = 'unequal';
end

% design-mismatch warning: varianceType only affects design [0 1] (Welch vs Student)
if ~isequal(di_cfg.analysis.designCode, [0 1]) && ~strcmp(di_cfg.analysis.OLS.varianceType, 'unequal')
    warning('Diagonale:OLS', 'varianceType is only used for design [0 1]; ignored here');
end

%% inferenceStatistic

% options available
inferenceStatisticTypes = {'auto', 't', 'r', 'beta'};

if isfield(di_cfg.analysis.OLS, 'inferenceStatistic')
    inferenceStatistic_raw = di_cfg.analysis.OLS.inferenceStatistic;
else
    inferenceStatistic_raw = 'auto';
end

if ~any(strcmp(inferenceStatistic_raw, inferenceStatisticTypes))
    disp('\\ allowed di_cfg.analysis.OLS.inferenceStatistic values: ')
    disp(inferenceStatisticTypes)
    error('\\ di_cfg.analysis.OLS.inferenceStatistic must be one of the allowed values mentioned above');
end

% resolve 'auto' per objective. The resolved value is NOT written back into
% di_cfg.analysis.OLS.inferenceStatistic: 'auto' must stay portable if the
% same di_cfg is reused for both objectives.
if strcmp(inferenceStatistic_raw, 'auto')
    switch di_cfg.analysis.objective
        case 'permutationH0testing'
            inferenceStatistic_resolved = 't';
        case 'bootstrapStability'
            inferenceStatistic_resolved = 'beta';
    end
else
    % explicit user choice is honoured as-is, even where it is a poor fit
    % for the chosen objective (see warnings below)
    inferenceStatistic_resolved = inferenceStatistic_raw;
end

% set of warnings for poor inferenceStatistic choices

% (A) 'r' + design [0 1] + varianceType 'unequal' (Welch):
% under Satterthwaite each feature has its own df, so r is not related to t
% by a map shared across features; the maxT family is not equivalent.
if strcmp(inferenceStatistic_resolved, 'r') && isequal(di_cfg.analysis.designCode, [0 1]) ...
        && strcmp(di_cfg.analysis.OLS.varianceType, 'unequal')
    warning('Diagonale:OLS', ['inferenceStatistic = ''r'' with design [0 1] and varianceType = ''unequal'': ' ...
        'under Satterthwaite each feature has its own df, so r is not related to t by a map shared ' ...
        'across features; the maxT family is not equivalent. Better choice: ''t''']);
end

% (B) 'beta' + permutationH0testing + correction 'maxT':
% Bad one! Raw betas are not scale-free -- their permutation null SD is proportional
% to each Y column's SD, so the maxT threshold is dominated by high-variance
% features and results are not invariant to per-feature rescaling.
% di_cfg.analysis.correction does not exist under bootstrapStability
% (di_validateAnalysis.m errors if it is set there), hence the isfield guard.
if strcmp(inferenceStatistic_resolved, 'beta') && strcmp(di_cfg.analysis.objective, 'permutationH0testing') ...
        && isfield(di_cfg.analysis, 'correction') && strcmp(di_cfg.analysis.correction, 'maxT')
    warning('Diagonale:OLS', ['inferenceStatistic = ''beta'' with objective = ''permutationH0testing'' and correction = ''maxT'': ' ...
        'raw betas are not scale-free! -- their permutation null SD is proportional to each Y column''s SD, ' ...
        'so the maxT threshold is dominated by high-variance features and results are not invariant to ' ...
        'per-feature rescaling. Better choice: ''t''']);
end

% (C) 't' or 'r' + bootstrapStability:
% bootstrap ratios already studentise by the resampling SD, so BR on t/r is
% a double normalisation, approaching t/r itself with added estimation
% noise. 'beta' is the appropriate input for bootstrapStability. The value
% is honoured as given, not overridden.
if any(strcmp(inferenceStatistic_resolved, {'t', 'r'})) && strcmp(di_cfg.analysis.objective, 'bootstrapStability')
    warning('Diagonale:OLS', ['inferenceStatistic = ''' inferenceStatistic_resolved ''' with objective = ''bootstrapStability'': ' ...
        'bootstrap ratios already studentise by the resampling SD, so BR on ' inferenceStatistic_resolved ...
        ' is a double normalisation, approaching ' inferenceStatistic_resolved ' itself with added estimation noise. ' ...
        'Better choice: ''beta'' is the appropriate input for bootstrapStability.']);
end

% write the resolved choice to the OLS field
di_cfg.analysis.OLS.inferenceStatistic_resolved = inferenceStatistic_resolved;

%% validation done

di_cfg.validation.OLS = true;

disp('\\ OLS parameters validated')

end
