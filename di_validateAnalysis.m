function di_cfg = di_validateAnalysis(di_cfg)
% Validate di_cfg.analysis
% 0. check dimensions have been validated upstream
% 1. sanity checking required fields (and potentially add some as defaults)
% 2. add a "validated" flag

%% analysis options and objectives supported by Diagonale

analysisObjectives = {'permutationH0testing' 'bootstrapStability'};
inferenceLevels = {'feature' 'cluster' 'latent'};
compatibilityRules = struct( ...
    'objective', {'permutationH0testing', 'permutationH0testing', 'permutationH0testing', 'bootstrapStability', 'bootstrapStability', 'bootstrapStability'}, ...
    'inferenceLevel', {'feature', 'cluster', 'latent', 'feature', 'cluster', 'latent'}, ...
    'supported', {true, true, true, true, false, true}, ...
    'testingApproaches', {{'parametric', 'empirical'}, {'parametric'}, {'empirical'}, {}, {}, {}}, ...
    'corrections', {{'FDR', 'maxT'}, {'maxT'}, {'maxT'}, {}, {}, {}} ...
    );

%% check dimensions have been validated

% make sure the dimensions are validated prior to this function
dimensionValidation = false; % initialize it false
if isfield(di_cfg,'validation')
    if isfield(di_cfg.validation,'dimensions')
        dimensionValidation = true;
    end
end

if dimensionValidation==false
    error('di_cfg.dimensions not validated. use: di_cfg = di_validateDimensions(di_cfg)')
end

%% shortcuts

% dimensions
dimKeys  = fieldnames(di_cfg.dimensions);
dimTypes = cellfun(@(k) di_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);
dimSizes = cellfun(@(k) length(di_cfg.dimensions.(k).vec), dimKeys);

%% validation checks

% sanity check: di_cfg.analysis exists
if ~isfield(di_cfg,'analysis')
    error('di_cfg needs field: analysis');
end

% sanity check: di_cfg.analysis.inferenceLevel exists
if ~isfield(di_cfg.analysis,'inferenceLevel')
    error('di_cfg.analysis needs field: inferenceLevel');
end

% sanity check: di_cfg.analysis.objective exists
if ~isfield(di_cfg.analysis,'objective')
    error('di_cfg.analysis needs field: objective');
end

% sanity check: analysis.inferenceLevel must be one of the allowed values
if ~any(strcmp(di_cfg.analysis.inferenceLevel,inferenceLevels))
    disp("\\ allowed inferenceLevel values: ")
    disp(inferenceLevels')
    error('\\ di_cfg.analysis.inferenceLevel must be one of the allowed values');
end

% sanity check: analysis.objective must be of the allowed objectives
if ~any(strcmp(di_cfg.analysis.objective,analysisObjectives))
    disp("\\ allowed objectives: ")
    disp(analysisObjectives)
    error('\\ di_cfg.analysis.objective must be of any of the allowed objectives');
end

% subject-centering mode
% accepted values:
%   'center'   -> remove subject means before model decomposition/fit
%   'noCenter' -> do not remove subject means
% subject centering is applied after ignoring rows (di_cfg.analysis.ignore_row) 
if isfield(di_cfg.analysis, 'subjectCentering')
    subjectCenteringTypes = {'center', 'noCenter'};
    if ~any(strcmp(di_cfg.analysis.subjectCentering, subjectCenteringTypes))
        disp("\\ allowed subjectCentering values: ")
        disp(subjectCenteringTypes)
        error('\\ di_cfg.analysis.subjectCentering must be one of the allowed values mentioned above');
    end
else
    di_cfg.analysis.subjectCentering = 'noCenter';
    warning('\\ di_cfg.analysis.subjectCentering not provided by the user. Using default: ''noCenter'' (no subject centering). Set to ''center'' to remove subject means before model decomposition/fit')
end

% pre-check for permuteUnit validity (just its existance with one of supported options)
% design-dependent validations are done later in di_validateH0 (after designCode is known).
if strcmp(di_cfg.analysis.objective,'permutationH0testing')
    permuteUnitTypes = {'wholeObservation', 'withinObservation'};
    if ~isfield(di_cfg.analysis,'permuteUnit')
        error('\\ permutation testing requires to speficy a type of permutation in di_cfg.analysis.permuteUnit')
    else
        if ~any(strcmp(di_cfg.analysis.permuteUnit, permuteUnitTypes))
            disp('allowed permuteUnit values:')
            disp(permuteUnitTypes)
            error('\\ di_cfg.analysis.permuteUnit must be one of the allowed values, mentioned above');
        end
    end
end



% --- permutation-specific compatibility checks ---
isPermutation = strcmp(di_cfg.analysis.objective,'permutationH0testing');
isFeatureLevel = strcmp(di_cfg.analysis.inferenceLevel,'feature');
isClusterLevel = strcmp(di_cfg.analysis.inferenceLevel,'cluster');
isLatentLevel  = strcmp(di_cfg.analysis.inferenceLevel,'latent');

% find the specific compatibility rule for the chosen objective and inference level,
% to check if they are compatible and to know which settings are required for the chosen 
% combination of objective and inference level.
ruleIdx = find(strcmp(di_cfg.analysis.objective, {compatibilityRules.objective}) & ...
    strcmp(di_cfg.analysis.inferenceLevel, {compatibilityRules.inferenceLevel}), 1);
if isempty(ruleIdx)
    error('Internal validation error: no compatibility rule found for objective = ''%s'' and inferenceLevel = ''%s''', ...
        di_cfg.analysis.objective, di_cfg.analysis.inferenceLevel)
end
selectedRule = compatibilityRules(ruleIdx);

if ~selectedRule.supported
    error('di_cfg.analysis.objective = ''%s'' is not supported with di_cfg.analysis.inferenceLevel = ''%s''', ...
        di_cfg.analysis.objective, di_cfg.analysis.inferenceLevel)
end

if isPermutation
    if ~isfield(di_cfg.analysis,'testingApproach')
        error('permutationH0testing requires di_cfg.analysis.testingApproach')
    end
    if ~any(strcmp(di_cfg.analysis.testingApproach, selectedRule.testingApproaches))
        disp('allowed testingApproach values:')
        disp(selectedRule.testingApproaches)
        error('\\ di_cfg.analysis.testingApproach must be one of the allowed values, mentioned above');
    end

    if ~isfield(di_cfg.analysis,'correction')
        error('permutationH0testing requires di_cfg.analysis.correction')
    end
    if ~any(strcmp(di_cfg.analysis.correction, selectedRule.corrections))
        disp('allowed correction values:')
        disp(selectedRule.corrections)
        error('\\ di_cfg.analysis.correction must be one of the allowed values, mentioned above');
    end

else
    if isfield(di_cfg.analysis,'testingApproach')
        error('di_cfg.analysis.testingApproach is only valid for di_cfg.analysis.objective = ''permutationH0testing''')
    end
    if isfield(di_cfg.analysis,'correction')
        error('di_cfg.analysis.correction is only valid for di_cfg.analysis.objective = ''permutationH0testing''')
    end
end

% --- cluster specific ---
% cluster-based analysis: validate required cluster parameters
if isClusterLevel
    if ~isfield(di_cfg.analysis,'clusterParams')
        error('di_cfg.analysis.inferenceLevel = ''cluster'' requires di_cfg.analysis.clusterParams')
    else
        if any(strcmp(dimTypes,'continuous'))
            if ~isfield(di_cfg.analysis.clusterParams,'distance_continuous_index')
                error('cluster inference with continuous dimensions requires clusterParams.distance_continuous_index (e.g., = 1)')
            end
        end
        if any(strcmp(dimTypes,'spherical'))
            if ~isfield(di_cfg.analysis.clusterParams,'distance_spherical_radians')
                % Try to infer default based on number of channels
                spherical_idx = find(strcmp(dimTypes,'spherical'));
                nChannels = dimSizes(spherical_idx);
                if nChannels == 32
                    di_cfg.analysis.clusterParams.distance_spherical_radians = 0.63;
                    warning('clusterParams.distance_spherical_radians not provided for 32-channel data. Using default: 0.63')
                elseif nChannels == 128
                    di_cfg.analysis.clusterParams.distance_spherical_radians = 0.36;
                    warning('clusterParams.distance_spherical_radians not provided for 128-channel data. Using default: 0.36')
                else
                    error('cluster inference with spherical dimensions requires clusterParams.distance_spherical_radians (e.g., = 0.63 for 32-channels or 0.36 for 128-channels)')
                end
            end
        end
        % adjacency matrix is necessary for cluster forming
        if ~isfield(di_cfg.analysis.clusterParams,'adjacencyMatrix')
            error('cluster inference requires clusterParams.adjacencyMatrix')
        end
    end
elseif isfield(di_cfg.analysis,'clusterParams')
    warning('\\ di_cfg.analysis.clusterParams provided but inferenceLevel is not ''cluster''. I will ignore di_cfg.analysis.clusterParams')
end

% --- feature-level multiple comparison correction specific ---
% Feature-level + permutation analyses may require additional settings
if isFeatureLevel && isPermutation

    % FDR settings are required only when correction='FDR'.
    if strcmp(di_cfg.analysis.correction,'FDR')
        if ~isfield(di_cfg.analysis,'FDR_dimensions')
            di_cfg.analysis.FDR_dimensions = ones(1,length(fieldnames(di_cfg.dimensions)));
            warning(['di_cfg.analysis.FDR_dimensions not provided. I am using ' num2str(di_cfg.analysis.FDR_dimensions) ' by default'])
        else
            if length(di_cfg.analysis.FDR_dimensions) ~= length(fieldnames(di_cfg.dimensions))
                error('analysis.FDR_dimensions must match the number of dimensions');
            end
        end
    elseif isfield(di_cfg.analysis,'FDR_dimensions')
        error('di_cfg.analysis.FDR_dimensions is only valid when di_cfg.analysis.correction = ''FDR''')
    end
else
    if isfield(di_cfg.analysis,'FDR_dimensions')
        error('di_cfg.analysis.FDR_dimensions is only valid for feature-level inference with permutation analyses with di_cfg.analysis.correction = ''FDR''')
    end
end
  

% --- PLS SVD specific ---
% PLS SVD parameters 
if isLatentLevel
    if ~isfield(di_cfg.analysis,'plssvdParams')
        di_cfg.analysis.plssvdParams = struct();
    end
    % Default zscoring vector: [Y_zscored, X_zscored]
    % Default [false, true]: suitable for Task PLS-SVD (Y is task codes/contrast, X needs z-scoring)
    % For continuous associations (PLSCorrelation), manually set to [true, true] to z-score both X and Y
    if ~isfield(di_cfg.analysis.plssvdParams,'zscoringVec')
        di_cfg.analysis.plssvdParams.zscoringVec = [false, true];
        warning('\\ analysis.plssvdParams.zscoringVec not provided. Using default: [false, true] (X zscored, Y not zscored). For continuous associations, set to [true, true].')
    end
elseif isfield(di_cfg.analysis,'plssvdParams')
    % provided but not applicable for the chosen inference level/objective
    warning('\\ di_cfg.analysis.plssvdParams provided but inferenceLevel/objective do not use latent inference. I will ignore di_cfg.analysis.plssvdParams')
end

%% dataStruct
% analysis.dataStruct table explains how the data are structured, what each row represents)

% sanity check: analysis.dataStruct exists and is a table 
if ~isfield(di_cfg.analysis,'dataStruct')
    error('di_cfg.analysis.dataStruct is needed informing on the structure of the data matrix')
else
    if ~istable(di_cfg.analysis.dataStruct)
        error('di_cfg.analysis.dataStruct must be a table');
    end
end

% normalise missing columns so downstream code (di_parseDesign,
% di_reorderRowsGenerate) always sees a complete dataStruct.
% only adds columns that are not already present; never overwrites user-supplied columns.
varNames = di_cfg.analysis.dataStruct.Properties.VariableNames;
nRows    = height(di_cfg.analysis.dataStruct);
if ~any(strcmp(varNames, 'observationID'))
    di_cfg.analysis.dataStruct.observationID = (1:nRows)';
    warning('\\ dataStruct has no observationID column. Adding observationID = (1:%d)'' (assuming each row is a unique independent observation).', nRows)
    disp('\\ note: if my assumption is not correct, create a di_cfg.analysis.dataStruct.observationID variable)')
end
if ~any(startsWith(varNames, 'indFactor'))
    di_cfg.analysis.dataStruct.indFactor1 = ones(nRows, 1);
    warning('\\ dataStruct has no indFactor# column. Adding indFactor1 = ones(%d,1) (treating all observations as one group).', nRows)
    disp('\\ I am assuming you do not want to compare between independent observations.')
    disp('\\ note: to achieve the same and not see the warning above, create a di_cfg.analysis.dataStruct.indFactor1 variable and use all 1s)')
end
if ~any(startsWith(varNames, 'repFactor'))
    di_cfg.analysis.dataStruct.repFactor1 = ones(nRows, 1);
    warning('\\ dataStruct has no repFactor# column. Adding repFactor1 = ones(%d,1) (treating all observations as one condition).', nRows)
    disp('\\ I am assuming you do not want to compare between repeated observations.')
    disp('\\ note: to achieve the same and not see the warning above, create a di_cfg.analysis.dataStruct.repFactor1 variable and use all 1s)')
end

%% defaults

if ~isfield(di_cfg.analysis,'randomSeed')
    di_cfg.analysis.randomSeed = 42;
    warning(['\\ di_cfg.analysis.randomSeed not provided. I am using ' num2str(di_cfg.analysis.randomSeed) ' by default'])
end

% random seed validation (primary validation point)
if ~isscalar(di_cfg.analysis.randomSeed) || ~isnumeric(di_cfg.analysis.randomSeed) || ...
        ~isfinite(di_cfg.analysis.randomSeed)
    error('\\ analysis.randomSeed must be a finite numeric scalar')
end
if di_cfg.analysis.randomSeed ~= round(di_cfg.analysis.randomSeed)
    error('\\ analysis.randomSeed must be integer-valued (e.g., 42)')
end


% p value critical (only if doing permutation testing)
if strcmp(di_cfg.analysis.objective,'permutationH0testing')
    if ~isfield(di_cfg.analysis,'p_crit')
        di_cfg.analysis.p_crit = 0.05;
        warning(['\\ di_cfg.analysis.p_crit not provided. I am using ' num2str(di_cfg.analysis.p_crit) ' by default'])
    end
end


% num of iterations 
nIterations_lowerEdge = 5000; % note: a bit arbitrary
if ~isfield(di_cfg.analysis,'nIterations')
    di_cfg.analysis.nIterations = nIterations_lowerEdge;
    warning(['\\ di_cfg.analysis.nIterations not provided. I am using ' num2str(di_cfg.analysis.nIterations) ' by default'])
else
    if di_cfg.analysis.nIterations < nIterations_lowerEdge
        warning('\\ analysis might be unstable if di_cfg.analysis.nIterations are fewer than %d', nIterations_lowerEdge)
    end
end

% Default cluster forming p-threshold when clustering is used
if isClusterLevel
    if ~isfield(di_cfg.analysis,'clusterParams') || ~isfield(di_cfg.analysis.clusterParams,'clusterFormingPvalThreshold')
        if ~isfield(di_cfg.analysis,'clusterParams')
            di_cfg.analysis.clusterParams = struct();
        end
        di_cfg.analysis.clusterParams.clusterFormingPvalThreshold = 0.05;
        warning(['\\ di_cfg.analysis.clusterParams.clusterFormingPvalThreshold not provided. I am using ' num2str(di_cfg.analysis.clusterParams.clusterFormingPvalThreshold) ' by default'])
    end
end


if ~isfield(di_cfg.analysis,'verbose')
    di_cfg.analysis.verbose = true;
    warning(['\\ di_cfg.analysis.verbose not provided. I am using logical ' num2str(di_cfg.analysis.verbose) ' by default'])
end


if ~isfield(di_cfg.analysis,'figFlag')
    di_cfg.analysis.figFlag = true;
    warning(['\\ di_cfg.analysis.figFlag not provided. I am using logical ' num2str(di_cfg.analysis.figFlag) ' by default'])
end

%% validation done
% add a validated flag
di_cfg.validation.analysis = true;

disp('\\ analysis validated')
end
