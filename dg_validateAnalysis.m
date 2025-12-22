function dg_cfg = dg_validateAnalysis(dg_cfg)
% Validate presence of required fields in dg_cfg.analysis

% make sure the dimensions are validated prior to this function
if ~isfield(dg_cfg,'validation')
    dimensionUnvalidated = true;
else
    dimensionUnvalidated = false;
    if ~isfield(dg_cfg.validation,'dimensions')
        dimensionUnvalidated = true;
    else
        if ~dg_cfg.validation.dimensions
            dimensionUnvalidated = true;
        end
    end
end
if dimensionUnvalidated
    error('dg_cfg.dimensions not validated. use: dg_cfg = pe_validateDimensions(dg_cfg)')
else
    dimKeys  = fieldnames(dg_cfg.dimensions);
    dimTypes = cellfun(@(k) dg_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);
end


analysisTypes = {'empiricalL1_FDR' 'theoreticalL1_clusterMaxT' 'PLSSVD' 'AJIVE'};
analysisObjectives = {'permutationH0testing' 'bootstrapStability'};

% dg_cfg.analysis exists
if ~isfield(dg_cfg,'analysis')
    error('dg_cfg needs field: analysis');
end

% dg_cfg.analysis.type exists
if ~isfield(dg_cfg.analysis,'type')
    error('dg_cfg.analysis needs field: type');
end

% dg_cfg.analysis.objective exists
if ~isfield(dg_cfg.analysis,'objective')
    error('dg_cfg.analysis needs field: objective');
end

% analysis.type must be of the allowed kind
if ~any(strcmp(dg_cfg.analysis.type,analysisTypes))
    disp("allowed types: ")
    disp(analysisTypes)
    error('dg_cfg.analysis.type must be of any of the allowed types');
end

% analysis.objective must be of the allowed objectives
if ~any(strcmp(dg_cfg.analysis.objective,analysisObjectives))
    disp("allowed objectives: ")
    disp(analysisObjectives)
    error('dg_cfg.analysis.objective must be of any of the allowed objectives');
end



if strcmp(dg_cfg.analysis.type,'theoreticalL1_clusterMaxT')
    if ~isfield(dg_cfg.analysis,'clusterParams')
        error('dg_cfg.analysis.type = theoreticalL1_clusterMaxT requires dg_cfg.analysis.clusterParams')
    else
        if any(strcmp(dimTypes,'continuous'))
            if ~isfield(dg_cfg.analysis.clusterParams,'distance_continuous_index')
                error('dg_cfg.analysis.type = theoreticalL1_clusterMaxT with continuous dimensions requires dg_cfg.analysis.clusterParams.distance_continuous_index (e.g., = 1)')
            end
        end
        if any(strcmp(dimTypes,'spherical'))
            if ~isfield(dg_cfg.analysis.clusterParams,'distance_spherical_radians')
                error('dg_cfg.analysis.type = theoreticalL1_clusterMaxT with spherical dimensions requires dg_cfg.analysis.clusterParams.distance_spherical_radians (e.g., = 0.63)')
            end
        end
    end
    
end

if strcmp(dg_cfg.analysis.type,'empiricalL1_FDR')
    if ~isfield(dg_cfg.analysis,'FDR_dimensions')
        error('dg_cfg.analysis.type = empiricalL1_FDR requires dg_cfg.analysis.FDR_dimensions')
        % apply FDR correction to which dimensions? 
        % example: 
        % [1 1 1] corrects pvalues from all three dimensions at once; 
        % [0 1 0] corrects pvalues for dimension 2 separately for each level of the other dimensions
    end
end
  
% must have a dataStruct table (this explains how the data are structured, what each row represents)
if ~isfield(dg_cfg.analysis,'dataStruct')
    error('dg_cfg.analysis.dataStruct is needed informing on the structure of the data matrix')
else
    if ~istable(dg_cfg.analysis.dataStruct)
        error('dg_cfg.analysis.dataStruct must be a table');
    end
end

%% defaults

if ~isfield(dg_cfg.analysis,'randomSeed')
    dg_cfg.analysis.randomSeed = 42;
    warning(['dg_cfg.analysis.randomSeed not provided. I am using ' num2str(dg_cfg.analysis.randomSeed) ' by default'])
end


% p value critical (only if doing permutation testing)
if strcmp(dg_cfg.analysis.objective,'permutationH0testing')
    if ~isfield(dg_cfg.analysis,'p_crit')
        dg_cfg.analysis.p_crit = 0.05;
        warning(['dg_cfg.analysis.p_crit not provided. I am using ' num2str(dg_cfg.analysis.p_crit) ' by default'])
    end
end


% num of iterations 
if ~isfield(dg_cfg.analysis,'nIterations')
    dg_cfg.analysis.nIterations = 5000;
    warning(['dg_cfg.analysis.nIterations not provided. I am using ' num2str(dg_cfg.analysis.nIterations) ' by default'])
else
    if dg_cfg.analysis.nIterations < 5000
        warning('analysis should be more stable if dg_cfg.analysis.nIterations was at least 5000')
    end
end


if strcmp(dg_cfg.analysis.type,'theoreticalL1_clusterMaxT')
    if ~isfield(dg_cfg.analysis.clusterParams,'clusterFormingPvalThreshold')
        dg_cfg.analysis.clusterParams.clusterFormingPvalThreshold = 0.05;
        warning(['dg_cfg.analysis.clusterParams.clusterFormingPvalThreshold not provided. I am using ' num2str(dg_cfg.analysis.clusterParams.clusterFormingPvalThreshold) ' by default'])
    end
end

if isfield(dg_cfg.analysis,'FDR_dimensions')
    if strcmp(dg_cfg.analysis.type,'empiricalL1_FDR')
        dg_cfg.analysis.FDR_dimensions = ones(1,length(fieldnames(dg_cfg.dimensions)));
        warning(['dg_cfg.analysis.FDR_dimensions not provided. I am using ' num2str(dg_cfg.analysis.FDR_dimensions) ' by default'])
    else
        warning('you have defined FDR dimension but you have chosen dg_cfg.analysis.type that does not use FDR correction. I will ignore dg_cfg.analysis.FDR_dimensions')
    end
end





if ~isfield(dg_cfg.analysis,'verbose')
    dg_cfg.analysis.verbose = true;
    warning(['dg_cfg.analysis.verbose not provided. I am using logical ' num2str(dg_cfg.analysis.verbose) ' by default'])
end


if ~isfield(dg_cfg.analysis,'figFlag')
    dg_cfg.analysis.figFlag = true;
    warning(['dg_cfg.analysis.figFlag not provided. I am using logical ' num2str(dg_cfg.analysis.figFlag) ' by default'])
end


% add a validated flag
dg_cfg.validation.analysis = true;

%%
disp('Diagonale: analysis validated')
end
