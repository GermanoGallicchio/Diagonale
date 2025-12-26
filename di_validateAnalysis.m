function di_cfg = di_validateAnalysis(di_cfg)
% Validate presence of required fields in di_cfg.analysis

% make sure the dimensions are validated prior to this function
if ~isfield(di_cfg,'validation')
    dimensionUnvalidated = true;
else
    dimensionUnvalidated = false;
    if ~isfield(di_cfg.validation,'dimensions')
        dimensionUnvalidated = true;
    else
        if ~di_cfg.validation.dimensions
            dimensionUnvalidated = true;
        end
    end
end
if dimensionUnvalidated
    error('di_cfg.dimensions not validated. use: di_cfg = pe_validateDimensions(di_cfg)')
else
    dimKeys  = fieldnames(di_cfg.dimensions);
    dimTypes = cellfun(@(k) di_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);
end


analysisTypes = {'empiricalFeature_inferenceFeature' 'parametricFeature_inferenceFeature' 'parametricFeature_inferenceCluster' 'PLSSVD' 'AJIVE'};
analysisObjectives = {'permutationH0testing' 'bootstrapStability'};

% di_cfg.analysis exists
if ~isfield(di_cfg,'analysis')
    error('di_cfg needs field: analysis');
end

% di_cfg.analysis.type exists
if ~isfield(di_cfg.analysis,'type')
    error('di_cfg.analysis needs field: type');
end

% di_cfg.analysis.objective exists
if ~isfield(di_cfg.analysis,'objective')
    error('di_cfg.analysis needs field: objective');
end

% analysis.type must be of the allowed kind
if ~any(strcmp(di_cfg.analysis.type,analysisTypes))
    disp("allowed types: ")
    disp(analysisTypes')
    error('di_cfg.analysis.type must be of any of the allowed types');
end

% analysis.objective must be of the allowed objectives
if ~any(strcmp(di_cfg.analysis.objective,analysisObjectives))
    disp("allowed objectives: ")
    disp(analysisObjectives)
    error('di_cfg.analysis.objective must be of any of the allowed objectives');
end



if strcmp(di_cfg.analysis.type,'theoretical_cluster')
    if ~isfield(di_cfg.analysis,'clusterParams')
        error('di_cfg.analysis.type = theoretical_cluster requires di_cfg.analysis.clusterParams')
    else
        if any(strcmp(dimTypes,'continuous'))
            if ~isfield(di_cfg.analysis.clusterParams,'distance_continuous_index')
                error('di_cfg.analysis.type = parametricFeature_inferenceCluster with continuous dimensions requires di_cfg.analysis.clusterParams.distance_continuous_index (e.g., = 1)')
            end
        end
        if any(strcmp(dimTypes,'spherical'))
            if ~isfield(di_cfg.analysis.clusterParams,'distance_spherical_radians')
                error('di_cfg.analysis.type = parametricFeature_inferenceCluster with spherical dimensions requires di_cfg.analysis.clusterParams.distance_spherical_radians (e.g., = 0.63)')
            end
        end
    end
    
end

if strcmp(di_cfg.analysis.type,'empirical_FDR')
    if ~isfield(di_cfg.analysis,'FDR_dimensions')
        error('di_cfg.analysis.type = empirical_FDR requires di_cfg.analysis.FDR_dimensions')
        % apply FDR correction to which dimensions? 
        % example: 
        % [1 1 1] corrects pvalues from all three dimensions at once; 
        % [0 1 0] corrects pvalues for dimension 2 separately for each level of the other dimensions
    end
end
  
% must have a dataStruct table (this explains how the data are structured, what each row represents)
if ~isfield(di_cfg.analysis,'dataStruct')
    error('di_cfg.analysis.dataStruct is needed informing on the structure of the data matrix')
else
    if ~istable(di_cfg.analysis.dataStruct)
        error('di_cfg.analysis.dataStruct must be a table');
    end
end

%% defaults

if ~isfield(di_cfg.analysis,'randomSeed')
    di_cfg.analysis.randomSeed = 42;
    warning(['di_cfg.analysis.randomSeed not provided. I am using ' num2str(di_cfg.analysis.randomSeed) ' by default'])
end


% p value critical (only if doing permutation testing)
if strcmp(di_cfg.analysis.objective,'permutationH0testing')
    if ~isfield(di_cfg.analysis,'p_crit')
        di_cfg.analysis.p_crit = 0.05;
        warning(['di_cfg.analysis.p_crit not provided. I am using ' num2str(di_cfg.analysis.p_crit) ' by default'])
    end
end


% num of iterations 
nIterations_lowerEdge = 5000; % note: a bit arbitrary
if ~isfield(di_cfg.analysis,'nIterations')
    di_cfg.analysis.nIterations = nIterations_lowerEdge;
    warning(['di_cfg.analysis.nIterations not provided. I am using ' num2str(di_cfg.analysis.nIterations) ' by default'])
else
    if di_cfg.analysis.nIterations < nIterations_lowerEdge
        warning('analysis might be unstable if di_cfg.analysis.nIterations are fewer than %d', nIterations_lowerEdge)
    end
end


if strcmp(di_cfg.analysis.type,'theoreticalL1_clusterMaxT')
    if ~isfield(di_cfg.analysis.clusterParams,'clusterFormingPvalThreshold')
        di_cfg.analysis.clusterParams.clusterFormingPvalThreshold = 0.05;
        warning(['di_cfg.analysis.clusterParams.clusterFormingPvalThreshold not provided. I am using ' num2str(di_cfg.analysis.clusterParams.clusterFormingPvalThreshold) ' by default'])
    end
end

if isfield(di_cfg.analysis,'FDR_dimensions')
    if strcmp(di_cfg.analysis.type,'empiricalL1_FDR')
        di_cfg.analysis.FDR_dimensions = ones(1,length(fieldnames(di_cfg.dimensions)));
        warning(['di_cfg.analysis.FDR_dimensions not provided. I am using ' num2str(di_cfg.analysis.FDR_dimensions) ' by default'])
    else
        warning('you have defined FDR dimension but you have chosen di_cfg.analysis.type that does not use FDR correction. I will ignore di_cfg.analysis.FDR_dimensions')
    end
end





if ~isfield(di_cfg.analysis,'verbose')
    di_cfg.analysis.verbose = true;
    warning(['di_cfg.analysis.verbose not provided. I am using logical ' num2str(di_cfg.analysis.verbose) ' by default'])
end


if ~isfield(di_cfg.analysis,'figFlag')
    di_cfg.analysis.figFlag = true;
    warning(['di_cfg.analysis.figFlag not provided. I am using logical ' num2str(di_cfg.analysis.figFlag) ' by default'])
end


% add a validated flag
di_cfg.validation.analysis = true;

%%
disp('Diagonale: analysis validated')
end
