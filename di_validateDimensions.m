function di_cfg = di_validateDimensions(di_cfg)
% Validate di_cfg.dimensions
% 1. sanity checking required fields (and potentially add some as defaults)
% 2. add a "validated" flag

%% dimension types supported by Diagonale

dimensionTypes = {'continuous' 'spherical' 'categorical'};

%% validation checks

% sanity check: di_cfg.dimensions exists
if ~isfield(di_cfg,'dimensions')
    error('\\ di_cfg needs field: dimensions');
end

% DELETE THIS SECTION AS IS IT MADE OBSOLETE BY THE SECTION BELOW
% % sanity check: subfields of di_cfg.dimensions need to start with "d"
% if any(~startsWith(fieldnames(di_cfg.dimensions),'d'))
%     error('di_cfg.dimensions dimensions need to start with d');
% end

% sanity check: subfields of di_cfg.dimensions need to be of the form "d" followed by a number
pattern = "d" + digitsPattern();
if any(~contains(fieldnames(di_cfg.dimensions),pattern))
    error('\\ di_cfg.dimensions dimensions need to be labeled as d followed by a number (e.g., d1)');
end

dims = fieldnames(di_cfg.dimensions);
for dIdx = 1:length(dims)

    % --- type ---
    % sanity check: di_cfg.dimensions.d#.type needs to exist
    if ~isfield(di_cfg.dimensions.(dims{dIdx}),'type')
        error('\\ di_cfg.dimensions dimensions need to contain a type field');    
    end

    % sanity check: di_cfg.dimensions.d#.type needs to be one of the supported types
    if ~any(strcmp(di_cfg.dimensions.(dims{dIdx}).type,dimensionTypes))
        disp("allowed types: ")
        disp(dimensionTypes)
        error('\\ di_cfg.dimensions dimensions need to contain a type field of any allowed type');    
    end

    % --- vec ---
    % sanity check: di_cfg.dimensions.d#.vec needs to exist
    if ~isfield(di_cfg.dimensions.(dims{dIdx}),'vec')
        error('\\ di_cfg.dimensions dimensions need to contain a vec field');    
    end

    % sanity check: di_cfg.dimensions.d#.vec is an actual vector
    if sum(size(di_cfg.dimensions.(dims{dIdx}).vec)~=1)~=1
        error('\\ di_cfg.dimensions.d*.vec must be a vector of non unit length')
    end

    % --- coord --- only for spherical dimensions
    if strcmp(di_cfg.dimensions.(dims{dIdx}).type,'spherical')
        if ~isfield(di_cfg.dimensions.(dims{dIdx}),'coord')
            error('di_cfg.dimensions.d* of spherical type must have a coord field')
        else
            if ~isfield(di_cfg.dimensions.(dims{dIdx}).coord,'sphTheta')
                error('\\ di_cfg.dimensions.d* of spherical type must have a coord field with sphTheta subfield')
            end
            if ~isfield(di_cfg.dimensions.(dims{dIdx}).coord,'sphPhi')
                error('\\ di_cfg.dimensions.d* of spherical type must have a coord field with sphPhi subfield')
            end
        end
        n_sphTheta = length(di_cfg.dimensions.(dims{dIdx}).coord.sphTheta);
        n_sphPhi = length(di_cfg.dimensions.(dims{dIdx}).coord.sphPhi);
        if ~isequal(n_sphTheta, n_sphPhi, length(di_cfg.dimensions.(dims{dIdx}).vec))
            error('\\ di_cfg.dimensions.d* of spherical type must have a coord field with sphTheta and sphPhi subfield with same numerosity as di_cfg.dimensions.d*.vec')
        end
        if isfield(di_cfg.dimensions.(dims{dIdx}).coord,'sphRadius')
            n_sphRadius = length(di_cfg.dimensions.(dims{dIdx}).coord.sphRadius);
            if ~isequal(n_sphRadius, length(di_cfg.dimensions.(dims{dIdx}).vec))
                error('\\ di_cfg.dimensions.d* of spherical type must have a coord field, and it can have a sphRadius field--if it does, it must be same numerosity as di_cfg.dimensions.d*.vec')
            end
        end
    end
end

%% defaults

for dIdx = 1:length(dims)

    % di_cfg.dimensions.d#.lbl = d# + type
    if ~isfield(di_cfg.dimensions.(dims{dIdx}),'lbl')
        dimLbl = [dims{dIdx}  di_cfg.dimensions.(dims{dIdx}).type];
        warning(['\\ label not provided by user for dimension ' dims{dIdx} '. Using this label: ' dimLbl])
        di_cfg.dimensions.(dims{dIdx}).lbl = dimLbl;
    end

    % di_cfg.dimensions.d#.units
    if ~isfield(di_cfg.dimensions.(dims{dIdx}),'units')
        dimUnits = 'undefinedunits';
        warning(['\\ units not provided by user for dimension ' dims{dIdx} ' (' di_cfg.dimensions.(dims{dIdx}).lbl ') . Using this label: ' dimUnits])
        di_cfg.dimensions.(dims{dIdx}).units = dimUnits;
    end

    % di_cfg.dimensions.d#.coord.sphRadius = 1
    if strcmp(di_cfg.dimensions.(dims{dIdx}).type,'spherical')
        if ~isfield(di_cfg.dimensions.(dims{dIdx}).coord,'sphRadius')
            di_cfg.dimensions.(dims{dIdx}).coord.sphRadius = ones(size(di_cfg.dimensions.(dims{dIdx}).vec));
            warning('\\ sphRadius not provided. using 1 by default')
        end
    end
end

%% add 'num' subfield of each dimension for numerosity

for dIdx = 1:length(dims)
    nVec = length(di_cfg.dimensions.(dims{dIdx}).vec);
    di_cfg.dimensions.(dims{dIdx}).num = nVec;
end

%% validation done

% add a validated flag
di_cfg.validation.dimensions = true;

disp('\\ dimensions validated')
end
