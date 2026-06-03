function di_cfg = di_validateH0(di_cfg)
% Validate di_cfg.analysis.permuteUnit
% only in case of permutationH0testing
% 1. require explicit permuteUnit from the user
% 2. check it is compatible with the design code
% 3. add a "validated" flag
%
% IMPORTANT:
% - permuteUnit controls permutation mechanics only (row reordering in
%   di_reorderRowsGenerate).
% - Subject centering is configured separately via
%   di_cfg.analysis.subjectCentering and applied inside analysis engines
%   (e.g., di_analysis_plsSVD). It is NOT encoded in permuteUnit.
%
% Allowed permuteUnit labels:
%   wholeObservation   : permute whole observation blocks
%   withinObservation  : permute repeated-measure rows within each observation

%% require explicit permuteUnit

permuteUnitTypes = {'wholeObservation' 'withinObservation'};

if ~isfield(di_cfg.analysis, 'permuteUnit')
    disp('\\ permutationH0testing requires di_cfg.analysis.permuteUnit to be set explicitly. Allowed values:')
    disp(permuteUnitTypes)
    error('\\ Set di_cfg.analysis.permuteUnit to one of the options mentioned above.')
end

%% check compatibility with design

key = sprintf('%d  %d', di_cfg.analysis.designCode(1), di_cfg.analysis.designCode(2));
switch key
    case {'0  0', '0  1'} % Correlation / independent observations: whole-observation only
        if ~strcmp(di_cfg.analysis.permuteUnit, 'wholeObservation')
            error(['\\ Design [' key '] only supports permuteUnit = ''wholeObservation''.'])
        end

    case '1  0' % Repeated measures: within-observation only
        if ~strcmp(di_cfg.analysis.permuteUnit, 'withinObservation')
            error(['\\ Design [' key '] only supports permuteUnit = ''withinObservation''.'])
        end

    case '1  1' % Mixed design: both types are allowed — no further check needed
end


%% validation done

% add a validated flag
di_cfg.validation.permuteUnit = true;

disp('\\ Permutation unit validated')
end
