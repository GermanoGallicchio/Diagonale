function results = di_analyze(di_cfg, Y, X)
% Front-end wrapper function that orchestrates the entire analysis pipeline:
% design detection, validation, permutation/bootstrap resampling, and 
% statistical testing. Supports correlation, independent groups, repeated 
% measures, and mixed designs using regression or PLS-SVD.
%
% INPUT:
%   di_cfg      - struct with validated analysis configuration containing:
%                 .analysis.nIterations, .analysis.dataStruct, .featureDims
%                 .analysis.inferenceLevel ('feature', 'cluster', or 'latent')
%   Y           - Primary data matrix (m x pData) where m = observations
%   X           - Design/response matrix (m x pDesign). Content is flexible and analysis-dependent:
%                 * Continuous predictors for correlation or PLS-SVD
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
%   - di_cfg.featureDims and di_cfg.analysis must be validated
%   - di_cfg.analysis.dataStruct must contain: observationID, 
%     indFactor# and repFactor# columns
%
% ANALYSIS TYPES SUPPORTED:
%   - Correlation (no factors)
%   - Independent samples t-tests / between-groups ANOVA
%   - Paired t-tests / repeated measures ANOVA
%   - Mixed designs (groups + repeated measures)
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% sanity checks and defaults

% sanity check: dimensions and analysis validated
dimensionValidation = false; % initialize it as false
analysisValidation  = false; % initialize it as false
if isfield(di_cfg,'validation')
    if isfield(di_cfg.validation,'featureDims')
        dimensionValidation = true;
    end
    if isfield(di_cfg.validation,'analysis')
        analysisValidation = true;
    end
end
if dimensionValidation==false
    error('di_cfg.featureDims not validated. use: di_cfg = di_validateFeatureDims(di_cfg)')
end
if analysisValidation==false
    error('di_cfg.analysis not validated. use: di_cfg = di_validateAnalysis(di_cfg)')
end

% sanity check: Y and X are matrices
if ~ismatrix(X)
    error('Y must be a matrix')
end
if ~ismatrix(Y)
    error('X must be a matrix')
end

% Y and X have the same num of rows
if size(X,1)~=size(Y,1)
    error('Y and X must have the same num of rows')
end

% X (and Y) must have the same num of rows as in di_cfg.analysis.dataStruct
if size(Y,1)~=size(di_cfg.analysis.dataStruct,1)
    error('X and di_cfg.analysis.dataStruct must have the same num of rows')
end

% ignore_col 
% if provided must be a row vector of same length as columns of Y (feature/data matrix)
% if not provided, nothing is ignored
if isfield(di_cfg.analysis, 'ignore_col')
    if ~isrow(di_cfg.analysis.ignore_col)
        error('\\ di_cfg.analysis.ignore_col must be a row vector')
    else
        if size(di_cfg.analysis.ignore_col,2)~=size(Y,2)
            error('\\ di_cfg.analysis.ignore_col must have same length as columns of Y')
        end
    end
end

% ignore_row 
% if provided must be a column vector of same length as rows of Y
% % if not provided, nothing is ignored
if isfield(di_cfg.analysis, 'ignore_row')
    if ~iscolumn(di_cfg.analysis.ignore_row)
        error('di_cfg.analysis.ignore_row must be a column vector')
    else
        if size(di_cfg.analysis.ignore_row,1)~=size(Y,1)
            error('di_cfg.analysis.ignore_row must have same length as rows of Y')
        end
    end
end

% univariate analysis only does one comparison at the time
if ismember(di_cfg.analysis.inferenceLevel, ["feature" "cluster"])
    if size(X,2)>1
        error(['more than one column in matrix X (i.e., more than one comparison) not supported for inferenceLevel = ' char(di_cfg.analysis.inferenceLevel) ])
    end
end

%% keep a copy of the matrices Y and X

Y_input = Y;
X_input = X;

%% apply ignore_row, if provided

if isfield(di_cfg.analysis,'ignore_row')

    % take a subset of matrices Y and X and of the accompanying dataStruct
    X = X(~di_cfg.analysis.ignore_row,:);
    Y = Y(~di_cfg.analysis.ignore_row,:);
    di_cfg.analysis.dataStruct = di_cfg.analysis.dataStruct(~di_cfg.analysis.ignore_row,:);

end

%% keep a copy of the matrices Y and X, after removal of of ignore_row

Y_input_pruned   = Y;
X_input_pruned   = X;

%% shortcuts

% matrix sizes (kept explicit for readability)
nRows_data   = size(Y, 1);
nCols_data   = size(Y, 2);
nCols_design = size(X, 2);

% num of dimensions and their numerosity 
% in d# order (fieldnames order)
dimKeys  = fieldnames(di_cfg.featureDims);
dimTypes = cellfun(@(k) di_cfg.featureDims.(k).type, dimKeys, 'UniformOutput', false);
dimSizes = cellfun(@(k) length(di_cfg.featureDims.(k).vec), dimKeys);
totalFeatureCount = prod(dimSizes);  % expected number of feature columns from dimensions

contIdx = find(strcmp(dimTypes, 'continuous'));
sphIdx  = find(strcmp(dimTypes, 'spherical'));
catIdx  = find(strcmp(dimTypes, 'categorical'));

% keep these explicit locals for lightweight runtime sanity readability
if nRows_data ~= size(X,1)
    error('Y and X must have the same num of rows')
end
if nCols_data ~= totalFeatureCount
    warning('Y columns (%d) differ from expected feature count from dimensions (%d)', nCols_data, totalFeatureCount);
end
if ismember(di_cfg.analysis.inferenceLevel, ["feature" "cluster"]) && nCols_design ~= 1
    error(['more than one column in matrix X (i.e., more than one comparison) not supported for inferenceLevel = ' char(di_cfg.analysis.inferenceLevel) ])
end

if numel(sphIdx) > 1
    error('\\ Only one spherical dimension is supported');
end

% cluster analyses do not support categorical dimensions
% because clusters are based on adjacency
if strcmp(di_cfg.analysis.inferenceLevel, "cluster") && ~isempty(catIdx)
    error('\\ categorical dimensions are not supported for cluster-based analyses');
end

%% parse design
% understand what analysis the user wants to do

designCode = di_parseDesign(di_cfg,X);

% keep a copy of the designCode in the analysis
di_cfg.analysis.designCode = designCode;

%% validate OLS parameters
% runs for both objectives (unlike di_validateH0, which is gated to
% permutationH0testing)

di_cfg = di_validateOLS(di_cfg);

%% validate null hypothesis
% if a null hypothesis is tested

if strcmp(di_cfg.analysis.objective, 'permutationH0testing') 
    di_cfg = di_validateH0(di_cfg);
end


%% perform the analysis

% get resampling indices
[rowIdx, di_cfg] = di_reorderRowsGenerate(di_cfg);

%% delegate to appropriate analysis function 
% based on type and design code

% build key using inference level and 2-digit design code
key = sprintf('%s & %d  %d', di_cfg.analysis.inferenceLevel, di_cfg.analysis.designCode(1), di_cfg.analysis.designCode(2));
switch key
    case {'feature & 0  0', ...
          'feature & 0  1', ...
          'feature & 1  0', ...
          'cluster & 0  0', ...
          'cluster & 0  1', ...
          'cluster & 1  0'}

        % unified OLS engine for all univariate designs and analysis types:
        %   design [0 0] : correlation (Pearson or Spearman via di_cfg.analysis.OLS.corrType)
        %   design [0 1] : independent-groups t-test (Welch or Student via di_cfg.analysis.OLS.varianceType)
        %   design [1 0] : paired t-test (subject fixed effects)
        %
        %   feature + empirical : beta statVal, empirical p-values via permutation
        %   feature + parametric: beta statVal, OLS/Welch SE -> parametric p-values
        %   cluster             : beta statVal, parametric p-values -> cluster forming
        % NOTE: convention here is Y=data matrix, X=design/contrast matrix.
        results = di_analysis_OLS(di_cfg, Y_input_pruned, X_input_pruned, rowIdx);


    case {'latent & 0  0', 'latent & 1  0', 'latent & 0  1', 'latent & 1  1'}
        % PLS-SVD multivariate analysis (all design codes)
        results = di_analysis_plsSVD(di_cfg, Y_input_pruned, X_input_pruned, rowIdx);


    otherwise
        % any other analysis type/design code combination
        error(['\\ not yet coded: ' char(di_cfg.analysis.inferenceLevel) ' & ' num2str(di_cfg.analysis.designCode)])

end

% keep a copy of the analysis in the results
results.analysis = di_cfg.analysis;

%% compute inferential metrics
% compute p values for permutation testing (from H0 shuffling)
% compute BR and CI for bootstrap stability (from resampling)

results = di_inference(di_cfg,results);


%% cluster descriptive metrics
% this bit is only doing something if clusters are meaningful
if isfield(di_cfg.analysis, 'clusterParams')
    results = di_clusterDescribe(di_cfg,results);
end

% keep full input matrices for downstream viewers after all processing
results.inputData.Y = Y_input;
results.inputData.X = X_input;

end
