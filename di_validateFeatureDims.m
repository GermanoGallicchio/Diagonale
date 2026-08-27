function di_cfg = di_validateFeatureDims(di_cfg)
% Validate di_cfg.featureDims
% 1. sanity checking required fields (and potentially add some as defaults)
% 2. add a "validated" flag

%% dimension types supported by Diagonale

dimensionTypes = {'continuous' 'spherical' 'categorical'};

%% validation checks

% sanity check: di_cfg.featureDims exists
if ~isfield(di_cfg,'dimensions')
    error('\\ di_cfg needs field: dimensions');
end

% DELETE THIS SECTION AS IS IT MADE OBSOLETE BY THE SECTION BELOW
% % sanity check: subfields of di_cfg.featureDims need to start with "d"
% if any(~startsWith(fieldnames(di_cfg.featureDims),'d'))
%     error('di_cfg.featureDims dimensions need to start with d');
% end

% sanity check: subfields of di_cfg.featureDims need to be of the form "d" followed by a number
pattern = "d" + digitsPattern();
if any(~contains(fieldnames(di_cfg.featureDims),pattern))
    error('\\ di_cfg.featureDims dimensions need to be labeled as d followed by a number (e.g., d1)');
end

dims = fieldnames(di_cfg.featureDims);
for dIdx = 1:length(dims)

    % --- type ---
    % sanity check: di_cfg.featureDims.d#.type needs to exist
    if ~isfield(di_cfg.featureDims.(dims{dIdx}),'type')
        error('\\ di_cfg.featureDims dimensions need to contain a type field');    
    end

    % sanity check: di_cfg.featureDims.d#.type needs to be one of the supported types
    if ~any(strcmp(di_cfg.featureDims.(dims{dIdx}).type,dimensionTypes))
        disp("allowed types: ")
        disp(dimensionTypes)
        error('\\ di_cfg.featureDims dimensions need to contain a type field of any allowed type');    
    end

    % --- vec ---
    % sanity check: di_cfg.featureDims.d#.vec needs to exist
    if ~isfield(di_cfg.featureDims.(dims{dIdx}),'vec')
        error('\\ di_cfg.featureDims dimensions need to contain a vec field');    
    end

    % sanity check: di_cfg.featureDims.d#.vec is an actual vector
    if sum(size(di_cfg.featureDims.(dims{dIdx}).vec)~=1)~=1
        error('\\ di_cfg.featureDims.d*.vec must be a vector of non unit length')
    end

    % --- coord --- only for spherical dimensions
    if strcmp(di_cfg.featureDims.(dims{dIdx}).type,'spherical')
        if ~isfield(di_cfg.featureDims.(dims{dIdx}),'coord')
            error('di_cfg.featureDims.d* of spherical type must have a coord field')
        else
            if ~isfield(di_cfg.featureDims.(dims{dIdx}).coord,'sphTheta')
                error('\\ di_cfg.featureDims.d* of spherical type must have a coord field with sphTheta subfield')
            end
            if ~isfield(di_cfg.featureDims.(dims{dIdx}).coord,'sphPhi')
                error('\\ di_cfg.featureDims.d* of spherical type must have a coord field with sphPhi subfield')
            end
        end
        n_sphTheta = length(di_cfg.featureDims.(dims{dIdx}).coord.sphTheta);
        n_sphPhi = length(di_cfg.featureDims.(dims{dIdx}).coord.sphPhi);
        if ~isequal(n_sphTheta, n_sphPhi, length(di_cfg.featureDims.(dims{dIdx}).vec))
            error('\\ di_cfg.featureDims.d* of spherical type must have a coord field with sphTheta and sphPhi subfield with same numerosity as di_cfg.featureDims.d*.vec')
        end
        if isfield(di_cfg.featureDims.(dims{dIdx}).coord,'sphRadius')
            n_sphRadius = length(di_cfg.featureDims.(dims{dIdx}).coord.sphRadius);
            if ~isequal(n_sphRadius, length(di_cfg.featureDims.(dims{dIdx}).vec))
                error('\\ di_cfg.featureDims.d* of spherical type must have a coord field, and it can have a sphRadius field--if it does, it must be same numerosity as di_cfg.featureDims.d*.vec')
            end
        end
    end
end

% spherical dimension must be declared last
% keeping the exception order explicit: "at most one spherical" is checked first
% and the spatial-order error is raised only when the declaration order is wrong

dimKeys = fieldnames(di_cfg.featureDims);
dimTypes = cellfun(@(k) di_cfg.featureDims.(k).type, dimKeys, 'UniformOutput', false);
sphIdx = find(strcmp(dimTypes,'spherical'));
if numel(sphIdx) > 1
    error('\ Only one spherical dimension is supported');
end
if ~isempty(sphIdx) && sphIdx ~= numel(dimKeys)
    error(['Spherical dimension must be declared last. Found ''%s'' (type ' ...
        '''spherical'') at position %d of %d.\n' ...
        'Reorder the d# fields of di_cfg.dimensions so the spherical ' ...
        'dimension comes last (e.g. d1 = time, d2 = frequency, ' ...
        'd3 = channel).\n' ...
        'This convention keeps figure data and figure axis labels derived ' ...
        'from the same declaration order.'], ...
        dimKeys{sphIdx}, sphIdx, numel(dimKeys));
end

%% defaults

for dIdx = 1:length(dims)

    % di_cfg.featureDims.d#.lbl = d# + type
    if ~isfield(di_cfg.featureDims.(dims{dIdx}),'lbl')
        dimLbl = [dims{dIdx}  di_cfg.featureDims.(dims{dIdx}).type];
        warning(['\\ label not provided by user for dimension ' dims{dIdx} '. Using this label: ' dimLbl])
        di_cfg.featureDims.(dims{dIdx}).lbl = dimLbl;
    end

    % di_cfg.featureDims.d#.units
    if ~isfield(di_cfg.featureDims.(dims{dIdx}),'units')
        dimUnits = 'undefinedunits';
        warning(['\\ units not provided by user for dimension ' dims{dIdx} ' (' di_cfg.featureDims.(dims{dIdx}).lbl ') . Using this label: ' dimUnits])
        di_cfg.featureDims.(dims{dIdx}).units = dimUnits;
    end

    % di_cfg.featureDims.d#.coord.sphRadius = 1
    if strcmp(di_cfg.featureDims.(dims{dIdx}).type,'spherical')
        if ~isfield(di_cfg.featureDims.(dims{dIdx}).coord,'sphRadius')
            di_cfg.featureDims.(dims{dIdx}).coord.sphRadius = ones(size(di_cfg.featureDims.(dims{dIdx}).vec));
            warning('\\ sphRadius not provided. using 1 by default')
        end
    end
end

%% add 'num' subfield of each dimension for numerosity

for dIdx = 1:length(dims)
    nVec = length(di_cfg.featureDims.(dims{dIdx}).vec);
    di_cfg.featureDims.(dims{dIdx}).num = nVec;
end

%% validation done

% add a validated flag
di_cfg.validation.featureDims = true;

disp('\\ dimensions validated')
end
