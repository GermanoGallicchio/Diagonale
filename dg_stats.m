function results = dg_stats(dg_cfg, Y, X)
% Main statistical analysis wrapper
%
% Front-end function that orchestrates the entire analysis pipeline:
% design detection, validation, permutation/bootstrap resampling, and 
% statistical testing. Supports correlation, independent groups, repeated 
% measures, and mixed designs using regression or PLS-SVD.
%
% INPUT:
%   dg_cfg      - struct with validated analysis configuration containing:
%                 .analysis.nIterations, .analysis.dataStruct, .dimensions
%                 .analysis.type ('empirical_FDR' or 'theoretical_maxT')
%   X           - Primary data matrix (m x pX) where m = observations, pX = variables
%   Y           - Design/response matrix (m x pY). Role depends on analysis:
%                 * Correlation: Y is a secondary data matrix to correlate with X
%                 * Grouped comparisons [0 1]: Y contains group identifiers
%                 * Paired/RM comparisons [1 0]: Y contains condition/contrast codes
%                 * PLS-SVD: Y is a secondary data matrix for multivariate analysis
%
% OUTPUT:
%   results     - struct containing:
%                 .designCode: [repMeasures, indGroups] binary code
%                 .stat: observed test statistic matrix
%                 .pVal: p-value matrix
%                 .clusterMetrics: cluster-based inference results (if applicable)
%
% PREREQUISITES:
%   - dg_cfg must be validated with dg_validateAnalysis first
%   - dg_cfg.analysis.dataStruct must contain: observationID, 
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
if ~isfield(dg_cfg,'validation')
    analysisUnvalidated = true;
else
    analysisUnvalidated = false;
    if ~isfield(dg_cfg.validation,'analysis')
        analysisUnvalidated = true;
    else
        if ~dg_cfg.validation.analysis
            analysisUnvalidated = true;
        end
    end
end
if analysisUnvalidated
    error('dg_cfg.analysis not validated. use: dg_cfg = dg_validateAnalysis(dg_cfg)')
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

% X (and Y) must have the same num of rows as in dg_cfg.analysis.dataStruct
if size(X,1)~=size(dg_cfg.analysis.dataStruct,1)
    error('X and dg_cfg.analysis.dataStruct must have the same num of rows')
end

% ignore_col if provided must be a row vector of same length as columns of X
if isfield(dg_cfg.analysis, 'ignore_col')
    if ~isrow(dg_cfg.analysis.ignore_col)
        error('dg_cfg.analysis.ignore_col must be a row vector')
    else
        if size(dg_cfg.analysis.ignore_col,2)~=size(X,2)
            error('dg_cfg.analysis.ignore_col must have same length as columns of X')
        end
    end
end

% ignore_row if provided must be a column vector of same length as rows of X
if isfield(dg_cfg.analysis, 'ignore_row')
    if ~iscolumn(dg_cfg.analysis.ignore_row)
        error('dg_cfg.analysis.ignore_row must be a column vector')
    else
        if size(dg_cfg.analysis.ignore_row,1)~=size(X,1)
            error('dg_cfg.analysis.ignore_row must have same length as rows of X')
        end
    end
end

% univariate analysis only does one comparison at the time
if ismember(dg_cfg.analysis.type, ["empirical_FDR" "theoretical_maxT"])
    if size(Y,2)>1
        error(['more than one column in matrix Y (i.e., more than one comparison) not supported for ' dg_cfg.analysis.type ])
    end
end

%% apply ignore_row, if provided

if isfield(dg_cfg.analysis,'ignore_row')

    % take a subset of matrices Y and X and of the accompanying dataStruct
    Y = Y(~dg_cfg.analysis.ignore_row,:);
    X = X(~dg_cfg.analysis.ignore_row,:);
    dg_cfg.analysis.dataStruct = dg_cfg.analysis.dataStruct(~dg_cfg.analysis.ignore_row,:);

end

%% keep a copy of the matrices Y and X

X_orig   = X;
Y_orig   = Y;

%% shortcuts

% num of iterations
nIterations = dg_cfg.analysis.nIterations;

% size of matrices Y and X
[m, pY] = size(Y);
[~, pX] = size(X);

% num of dimensions and their numerosity 
% in d# order (fieldnames order)
dimKeys  = fieldnames(dg_cfg.dimensions);
dimTypes = cellfun(@(k) dg_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);
dimSizes = cellfun(@(k) length(dg_cfg.dimensions.(k).vec), dimKeys);
nDims    = numel(dimKeys);
Nall     = prod(dimSizes);

contIdx = find(strcmp(dimTypes, 'continuous'));
sphIdx  = find(strcmp(dimTypes, 'spherical'));
catIdx  = find(strcmp(dimTypes, 'categorical'));

if numel(sphIdx) > 1
    error('Only one spherical dimension is supported');
end

% cluster analyses do not support categorical dimensions
if strcmp(dg_cfg.analysis.type, "theoretical_maxT") && ~isempty(catIdx)
    error('categorical dimensions are not supported for cluster-based analyses');
end

%% parse design
% understand what analysis the user wants to do

designCode = dg_parseDesign(dg_cfg,Y);

% keep a copy of the designCode in the analysis and in the results
dg_cfg.analysis.designCode = designCode;
results.designCode = designCode;

keyboard; % UNTIL HERE MAYBE OK

%% perform the analysis

% get resampling indices
rowIdx = dg_reorderRowsGenerate(dg_cfg);

%% delegate to appropriate analysis function based on type and design code

% Build key using analysis type and 2-digit design code
key = sprintf('%s & %d  %d', dg_cfg.analysis.type, dg_cfg.designCode(1), dg_cfg.designCode(2));

switch key
    case 'empirical_FDR & 0  0'
        % Correlation analysis (no groups, no repeated measures)
        results = dg_analysis_empirical_correlation(dg_cfg, Y_orig, X_orig, rowIdx, nIterations, pX);

    case 'empirical_FDR & 0  1'
        % Independent groups t-test (independent groups, no repeated measures)
        results = dg_analysis_empirical_ttestInd(dg_cfg, Y_orig, X_orig, rowIdx, nIterations, pX);

    case 'empirical_FDR & 1  0'
        % Paired t-test (repeated measures, no independent groups)
        results = dg_analysis_empirical_ttestPaired(dg_cfg, Y_orig, X_orig, rowIdx, nIterations, pX);

    case 'theoretical_maxT & 0  0'
        % Cluster-based correlation analysis
        results = dg_analysis_theoretical_correlation(dg_cfg, Y_orig, X_orig, rowIdx, nIterations, Nall);

    case 'theoretical_maxT & 0  1'
        % Cluster-based independent groups t-test
        results = dg_analysis_theoretical_ttestInd(dg_cfg, Y_orig, X_orig, rowIdx, nIterations, Nall);

    case 'theoretical_maxT & 1  0'
        % Cluster-based paired t-test
        results = dg_analysis_theoretical_ttestPaired(dg_cfg, Y_orig, X_orig, rowIdx, nIterations, Nall);

    case {'PLS_SVD & 0  0', 'PLS_SVD & 1  0', 'PLS_SVD & 0  1', 'PLS_SVD & 1  1'}
        % PLS-SVD multivariate analysis (all design codes)
        results = dg_analysis_plsSVD(dg_cfg, Y_orig, X_orig, rowIdx, nIterations, Nall);

    otherwise
        % Any other analysis type/design code combination
        error(['not yet coded: ' dg_cfg.analysis.type ' & ' num2str(dg_cfg.designCode)])
end

%% compute inferential metrics
% compute p values for permutation testing
% compute BR and CI for bootstrap stability

results = dg_inference(dg_cfg,results);
       
%% cluster descriptive metrics

results = dg_describeClusters(dg_cfg,results);

end
