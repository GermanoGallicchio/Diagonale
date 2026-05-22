function designCode = di_parseDesign(di_cfg, Y)
% Deduces the experimental design from dataStruct configuration
%
% Determines whether the analysis involves repeated measures factors and/or 
% independent groups. Returns a 
% 2-digit binary code indicating the design structure.
%
% INPUT:
%   di_cfg      - struct with .analysis.dataStruct table containing:
%                 - observationID: unique identifier for each observation
%                 - indFactor#: independent group factors (optional)
%                 - repFactor#: repeated measure factors (optional)
%   Y           - data matrix (m x pY) used for validation only
%
% OUTPUT:
%   designCode  - 1x2 binary vector [repMeasures, indGroups]
%                 [0 0]: Correlation/association (no factors)
%                 [0 1]: At least one Independent observations factor
%                 [1 0]: At least one Repeated measures factor
%                 [1 1]: Mixed design (independent + repeated measures)
%
% watch out for:
%   - Unbalanced designs (subjects with different number of observations)
%     will trigger a warning but will use max count for structure detection
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% sanity checks

varLbl = di_cfg.analysis.dataStruct.Properties.VariableNames;

% if indFactor# not entered in dataStruct: error
if ~any(startsWith(varLbl,'indFactor'))
    error('\ "indFactor#" column(s) in di_cfg.analysis.dataStruct not entered')
end

% if repFactor# not entered in dataStruct: error
if ~any(startsWith(varLbl,'repFactor'))
    error('\\ "repFactor#" column(s) in di_cfg.analysis.dataStruct not entered')
end

% dataStruct must include one observationID column
if ~any(contains(varLbl,'observationID'))
    error('\\ "observationID" column in di_cfg.analysis.dataStruct not entered.')
end

%% pull information from the table

% num of rows
m = size(di_cfg.analysis.dataStruct,1); 

uniqueness = cellfun(@(k) length(unique(di_cfg.analysis.dataStruct.(k))), varLbl, 'UniformOutput', true);

% num of unique observations
uniqueObservationID_numerosity = uniqueness(strcmp(varLbl,'observationID'));

% num of indFactors
nUnique_indFactors = sum(startsWith(varLbl,'indFactor'));

% num of levels per each indFactor
nUnique_levelsPerIndFactors = uniqueness(startsWith(varLbl,'indFactor'));

% num of levels across all indFactors (pooled)
nUnique_indFactorLevels = prod(nUnique_levelsPerIndFactors);

% num of repFactors
nUnique_repFactors = sum(startsWith(varLbl,'repFactor'));

% num of levels per each repFactor
nUnique_levelsPerRepFactors = uniqueness(startsWith(varLbl,'repFactor'));

% num of levels across all repFactor (pooled)
nUnique_repFactorLevels = prod(nUnique_levelsPerRepFactors);

% check if there are repeated observations per observationID
obsID_counts = groupcounts(di_cfg.analysis.dataStruct.observationID);
repeatedPerObservationID_numerosity = max(obsID_counts);

% warning if unbalanced
if length(unique(obsID_counts)) > 1
    warning('unbalanced design detected: observations have different numbers of repeated measures (range: %d to %d).', ...
        min(obsID_counts), max(obsID_counts));
end

%% deduct intended design

% create a two digit binary number, indicating whether the analysis
% intends to compare across repeated measure factors and/or independent groups
% from left to right:
% 1st digit = repFactor (e.g., paired sample t-tests, repeated measures ANOVA)
% 2nd digit = indFactor (e.g., independent sample t-tests, between-groups ANOVA)
% it does not matter if there are multiple repFactor (or indFactor), the
% values are still 0 (no such comparison) or 1 (there are such comparisons)

designOptions = di_designOptions;

% Determine design code only based on data structure (dataStruct)
% (ie, not based on Y dimensions, since Y can have varying columns for PLS-SVD)
hasRepMeasures = (repeatedPerObservationID_numerosity > 1);
hasIndGroups = (nUnique_indFactorLevels > 1);

% initialize design code
designCode = [0 0];

if hasRepMeasures
    designCode(1) = 1;  % repeated measures present
end

if hasIndGroups
    designCode(2) = 1;  % independent groups of observations present
end


%% plot structure of Y

if di_cfg.analysis.figFlag
    figure()
    imagesc(normalize(Y,'range',[-1 1]))
    title('structure of Y')
    %axis equal
    box off
    axis off
    colormap(parula)
    colorbar;
end

%% 
if di_cfg.analysis.verbose
    rowIdx = (designOptions{:,"codeDec"}==bin2dec(num2str(designCode)))';
    fprintf(['\n' '\\ Design num is: ' num2str(designOptions{rowIdx,"codeDec"}) ', ' designOptions{rowIdx,"lbl"}{1} '\n'])
end

