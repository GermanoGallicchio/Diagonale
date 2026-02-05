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
    error('di_cfg.dimensions not validated. use: di_cfg = di_validateDimensions(di_cfg)')
else
    dimKeys  = fieldnames(di_cfg.dimensions);
    dimTypes = cellfun(@(k) di_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);
end


analysisTypes = {'empiricalFeature_inferenceFeature' 'parametricFeature_inferenceFeature' 'parametricFeature_inferenceCluster' 'PLS_SVD' 'AJIVE'};
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


% Cluster-based analysis: validate required cluster parameters
if strcmp(di_cfg.analysis.type,'parametricFeature_inferenceCluster')
    if ~isfield(di_cfg.analysis,'clusterParams')
        error('di_cfg.analysis.type = parametricFeature_inferenceCluster requires di_cfg.analysis.clusterParams')
    else
        if any(strcmp(dimTypes,'continuous'))
            if ~isfield(di_cfg.analysis.clusterParams,'distance_continuous_index')
                error('parametricFeature_inferenceCluster with continuous dimensions requires clusterParams.distance_continuous_index (e.g., = 1)')
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
                    error('parametricFeature_inferenceCluster with spherical dimensions requires clusterParams.distance_spherical_radians (e.g., = 0.63 for 32-channels or 0.36 for 128-channels)')
                end
            end
        end
        % adjacency matrix is necessary for cluster forming
        if ~isfield(di_cfg.analysis.clusterParams,'adjacencyMatrix')
            error('parametricFeature_inferenceCluster requires clusterParams.adjacencyMatrix')
        end
    end
end

% Empirical feature-level analysis with permutation: FDR dimensions default
if strcmp(di_cfg.analysis.type,'empiricalFeature_inferenceFeature') && strcmp(di_cfg.analysis.objective,'permutationH0testing')
    if ~isfield(di_cfg.analysis,'FDR_dimensions')
        di_cfg.analysis.FDR_dimensions = ones(1,length(fieldnames(di_cfg.dimensions)));
        warning(['di_cfg.analysis.FDR_dimensions not provided. I am using ' num2str(di_cfg.analysis.FDR_dimensions) ' by default'])
    else
        if length(di_cfg.analysis.FDR_dimensions) ~= length(fieldnames(di_cfg.dimensions))
            error('analysis.FDR_dimensions must match the number of dimensions');
        end
    end
elseif isfield(di_cfg.analysis,'FDR_dimensions')
    % Provided but not applicable for the chosen type/objective
    warning('di_cfg.analysis.FDR_dimensions provided but type/objective do not use FDR correction. I will ignore di_cfg.analysis.FDR_dimensions')
end
  


% PLS SVD parameters 
if strcmp(di_cfg.analysis.type,'PLS_SVD')
    if ~isfield(di_cfg.analysis,'plssvdParams')
        di_cfg.analysis.plssvdParams = struct();
    end
    % Default zscoring vector: [Y_zscored, X_zscored]
    % Default [false, true]: suitable for Task PLS-SVD (Y is task codes/contrast, X needs z-scoring)
    % For continuous associations (PLSCorrelation), manually set to [true, true] to z-score both X and Y
    if ~isfield(di_cfg.analysis.plssvdParams,'zscoringVec')
        di_cfg.analysis.plssvdParams.zscoringVec = [false, true];
        warning('plssvdParams.zscoringVec not provided. Using default: [false, true] (X zscored, Y not zscored). For continuous associations, set to [true, true].')
    end
elseif isfield(di_cfg.analysis,'plssvdParams')
    % provided but not applicable for the chosen type/objective
    warning('di_cfg.analysis.plssvdParams provided but type/objective do not use PLS SVD. I will ignore di_cfg.analysis.plssvdParams')
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

% Default cluster forming p-threshold when clustering is used
if strcmp(di_cfg.analysis.type,'parametricFeature_inferenceCluster')
    if ~isfield(di_cfg.analysis,'clusterParams') || ~isfield(di_cfg.analysis.clusterParams,'clusterFormingPvalThreshold')
        if ~isfield(di_cfg.analysis,'clusterParams')
            di_cfg.analysis.clusterParams = struct();
        end
        di_cfg.analysis.clusterParams.clusterFormingPvalThreshold = 0.05;
        warning(['di_cfg.analysis.clusterParams.clusterFormingPvalThreshold not provided. I am using ' num2str(di_cfg.analysis.clusterParams.clusterFormingPvalThreshold) ' by default'])
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
