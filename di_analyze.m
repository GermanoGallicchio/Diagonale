function results = di_analyze(di_cfg, Y, X)
% Main statistical analysis wrapper
%
% Front-end function that orchestrates the entire analysis pipeline:
% design detection, validation, permutation/bootstrap resampling, and 
% statistical testing. Supports correlation, independent groups, repeated 
% measures, and mixed designs using regression or PLS-SVD.
%
% INPUT:
%   di_cfg      - struct with validated analysis configuration containing:
%                 .analysis.nIterations, .analysis.dataStruct, .dimensions
%                 .analysis.type ('empiricalFeature_inferenceFeature', 'parametricFeature_inferenceFeature', or 'parametricFeature_inferenceCluster')
%   X           - Primary data matrix (m x pX) where m = observations, pX = variables
%   Y           - Design/response matrix (m x pY). Content is flexible and analysis-dependent:
%                 * Continuous predictors for correlation/regression
%                 * Categorical codes for group comparisons [0 1]
%                 * Condition/contrast codes for repeated measures [1 0]
%                 * Multivariate data matrix for PLS-SVD
%                 * Can contain mixed types (e.g., continuous + categorical columns)
%
% OUTPUT:
%   results     - struct containing:
%                 .designCode: [repMeasures, indGroups] binary code
%                 .stat: observed test statistic matrix
%                 .pVal: p-value matrix
%                 .clusterMetrics: cluster-based inference results (if applicable)
%
% PREREQUISITES:
%   - di_cfg must be validated with di_validateAnalysis first
%   - di_cfg.analysis.dataStruct must contain: observationID, 
%     and optional indFactor# and repFactor# columns
%
% ANALYSIS TYPES SUPPORTED:
%   - Correlation (no factors)
%   - Independent samples t-tests / between-groups ANOVA
%   - Paired t-tests / repeated measures ANOVA
%   - Mixed designs (groups + repeated measures)
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% sanity checks

% make sure the analysis is validated prior to this function
if ~isfield(di_cfg,'validation')
    analysisUnvalidated = true;
else
    analysisUnvalidated = false;
    if ~isfield(di_cfg.validation,'analysis')
        analysisUnvalidated = true;
    else
        if ~di_cfg.validation.analysis
            analysisUnvalidated = true;
        end
    end
end
if analysisUnvalidated
    error('di_cfg.analysis not validated. use: di_cfg = di_validateAnalysis(di_cfg)')
end



% Y and X are matrices
if ~ismatrix(Y)
    error('Y must be a matrix')
end
if ~ismatrix(X)
    error('X must be a matrix')
end

% Y and X have the same num of rows
if size(Y,1)~=size(X,1)
    error('Y and X must have the same num of rows')
end

% X (and Y) must have the same num of rows as in di_cfg.analysis.dataStruct
if size(X,1)~=size(di_cfg.analysis.dataStruct,1)
    error('X and di_cfg.analysis.dataStruct must have the same num of rows')
end

% ignore_col if provided must be a row vector of same length as columns of X
if isfield(di_cfg.analysis, 'ignore_col')
    if ~isrow(di_cfg.analysis.ignore_col)
        error('di_cfg.analysis.ignore_col must be a row vector')
    else
        if size(di_cfg.analysis.ignore_col,2)~=size(X,2)
            error('di_cfg.analysis.ignore_col must have same length as columns of X')
        end
    end
end

% ignore_row if provided must be a column vector of same length as rows of X
if isfield(di_cfg.analysis, 'ignore_row')
    if ~iscolumn(di_cfg.analysis.ignore_row)
        error('di_cfg.analysis.ignore_row must be a column vector')
    else
        if size(di_cfg.analysis.ignore_row,1)~=size(X,1)
            error('di_cfg.analysis.ignore_row must have same length as rows of X')
        end
    end
end

% univariate analysis only does one comparison at the time
if ismember(di_cfg.analysis.type, ["empiricalFeature_inferenceFeature" "parametricFeature_inferenceFeature" "parametricFeature_inferenceCluster"])
    if size(Y,2)>1
        error(['more than one column in matrix Y (i.e., more than one comparison) not supported for ' di_cfg.analysis.type ])
    end
end

%% apply ignore_row, if provided

if isfield(di_cfg.analysis,'ignore_row')

    % take a subset of matrices Y and X and of the accompanying dataStruct
    Y = Y(~di_cfg.analysis.ignore_row,:);
    X = X(~di_cfg.analysis.ignore_row,:);
    di_cfg.analysis.dataStruct = di_cfg.analysis.dataStruct(~di_cfg.analysis.ignore_row,:);

end

%% keep a copy of the matrices Y and X

X_orig   = X;
Y_orig   = Y;

%% shortcuts

% num of iterations
nIterations = di_cfg.analysis.nIterations;

% size of matrices Y and X
[m, pY] = size(Y);
[~, pX] = size(X);

% num of dimensions and their numerosity 
% in d# order (fieldnames order)
dimKeys  = fieldnames(di_cfg.dimensions);
dimTypes = cellfun(@(k) di_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);
dimSizes = cellfun(@(k) length(di_cfg.dimensions.(k).vec), dimKeys);
nDims    = numel(dimKeys);
pX       = prod(dimSizes);  % total number of X features across all dimensions

contIdx = find(strcmp(dimTypes, 'continuous'));
sphIdx  = find(strcmp(dimTypes, 'spherical'));
catIdx  = find(strcmp(dimTypes, 'categorical'));

if numel(sphIdx) > 1
    error('Only one spherical dimension is supported');
end

% cluster analyses do not support categorical dimensions
if strcmp(di_cfg.analysis.type, "parametricFeature_inferenceCluster") && ~isempty(catIdx)
    error('categorical dimensions are not supported for cluster-based analyses');
end

%% parse design
% understand what analysis the user wants to do

designCode = di_parseDesign(di_cfg,Y);

% keep a copy of the designCode in the analysis
di_cfg.analysis.designCode = designCode;

%% perform the analysis

% get resampling indices
rowIdx = di_reorderRowsGenerate(di_cfg);

%% delegate to appropriate analysis function based on type and design code

% build key using analysis type and 2-digit design code
key = sprintf('%s & %d  %d', di_cfg.analysis.type, di_cfg.analysis.designCode(1), di_cfg.analysis.designCode(2));

switch key
    case 'empiricalFeature_inferenceFeature & 0  0'

        keyboard; % UNTIL HERE MAYBE OK

        % correlation
        % - empirical (via simulations) at feature level
        % - FDR correction
        % - cluster forming (descriptive)
        results = di_analysis_empiricalFeature_inferenceFeature_correlation(di_cfg, Y_orig, X_orig, rowIdx);

    case 'empiricalFeature_inferenceFeature & 0  1'
        
        % independent sample/groups t-test
        % - empirical (via simulations) at feature level
        % - FDR correction
        % - cluster forming (descriptive)
        results = di_analysis_empiricalFeature_inferenceFeature_ttestInd(di_cfg, Y_orig, X_orig, rowIdx);

    case 'empiricalFeature_inferenceFeature & 1  0'

        keyboard; % UNTIL HERE MAYBE OK

        % paired sample t-test
        % - empirical (via simulations) at feature level
        % - FDR correction
        % - cluster forming (descriptive)
        results = di_analysis_empiricalFeature_inferenceFeature_ttestPaired(di_cfg, Y_orig, X_orig, rowIdx);

    case 'parametricFeature_inferenceCluster & 0  0'
        
        keyboard; % UNTIL HERE MAYBE OK

        % correlation
        % - theoretical at feature level
        % - cluster forming (for inference)
        results = di_analysis_parametricFeature_inferenceCluster_correlation(di_cfg, Y_orig, X_orig, rowIdx);

    case 'parametricFeature_inferenceCluster & 0  1'
        % independent sample/groups t-test
        % - theoretical at feature level
        % - cluster forming (for inference)
        results = di_analysis_parametricFeature_inferenceCluster_ttestInd(di_cfg, Y_orig, X_orig, rowIdx);

    case 'parametricFeature_inferenceCluster & 1  0'

        keyboard; % UNTIL HERE MAYBE OK

        % paired sample t-test
        % - theoretical at feature level
        % - cluster forming (for inference)
        results = di_analysis_parametricFeature_inferenceCluster_ttestPaired(di_cfg, Y_orig, X_orig, rowIdx);

    case {'PLS_SVD & 0  0', 'PLS_SVD & 1  0', 'PLS_SVD & 0  1', 'PLS_SVD & 1  1'}

        keyboard; % UNTIL HERE MAYBE OK

        % PLS-SVD multivariate analysis (all design codes)
        results = di_analysis_plsSVD(di_cfg, Y_orig, X_orig, rowIdx);

    otherwise
        % any other analysis type/design code combination
        error(['not yet coded: ' di_cfg.analysis.type ' & ' num2str(di_cfg.analysis.designCode)])
end

% keep a copy of the analysis in the results
results.analysis = di_cfg.analysis;

%% compute inferential metrics
% compute p values for permutation testing (from H0 shuffling)
% compute BR and CI for bootstrap stability (from resampling)

results = di_inference(di_cfg,results);

%% cluster descriptive metrics
% this bit is only doing something if clusters are meaningful

results = di_describeClusters(di_cfg,results);


end
