function pe_cfg = pe_validateConfig(pe_cfg)
% Validate presence of required fields in pe_cfg.dimensions and include basic computations for shortcuts later on.

% Usage:
%   pe_cfg = pe_validateConfig(pe_cfg)
%
% This function throws informative errors if fields are missing or invalid.
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%%

dimensionTypes = {'continuous' 'spherical' 'categorical'};

%% input validation



% pe_cfg.dimensions exists
if ~isfield(pe_cfg,'dimensions')
    error('pe_cfg needs field: dimensions');
end

% pe_cfg.dimensions fields need to start with d
if any(~startsWith(fieldnames(pe_cfg.dimensions),'d'))
    error('pe_cfg.dimensions dimensions need to start with d');
end

% pe_cfg.dimensions fields need to be of pattern 'd' followed by a number
pattern = "d" + digitsPattern();
if any(~contains(fieldnames(pe_cfg.dimensions),pattern))
    error('pe_cfg.dimensions dimensions need to be labeled as d followed by a number');
end



dims = fieldnames(pe_cfg.dimensions);
for dIdx = 1:length(dims)

    % pe_cfg.dimensions.d* must have a type field
    if ~any(strcmp(fieldnames(pe_cfg.dimensions.(dims{dIdx})),'type'))
        error('pe_cfg.dimensions dimensions need to contain a type field');    
    end

    % pe_cfg.dimensions need to belong to one of the allowed types
    if ~any(strcmp(pe_cfg.dimensions.(dims{dIdx}).type,dimensionTypes))
        disp("allowed types: ")
        disp(dimensionTypes)
        error('pe_cfg.dimensions dimensions need to contain a type field of any allowed type');    
    end

    % there can only be one spherical type
    % TO DO

    
    % pe_cfg.dimensions.d* must have a vec field
    if ~any(strcmp(fieldnames(pe_cfg.dimensions.(dims{dIdx})),'vec'))
        error('pe_cfg.dimensions dimensions need to contain a vec field');    
    end

    % pe_cfg.dimensions.d*.vec must be a vector 
    % (i.e., one and only one dimension different from 1)
    if sum(size(pe_cfg.dimensions.(dims{dIdx}).vec)~=1)~=1
        error('pe_cfg.dimensions.d*.vec must be a vector of non unit length')
    end

    % pe_cfg.dimensions.d* of spherical type must have a coord field with
    % subfields sphTheta and sphPhi of same size as vec
    if strcmp(pe_cfg.dimensions.(dims{dIdx}).type,'spherical')

        if ~isfield(pe_cfg.dimensions.(dims{dIdx}),'coord')
            error('pe_cfg.dimensions.d* of spherical type must have a coord field')
        else
            if ~isfield(pe_cfg.dimensions.(dims{dIdx}).coord,'sphTheta')
                error('pe_cfg.dimensions.d* of spherical type must have a coord field with sphTheta subfield')
            end
            if ~isfield(pe_cfg.dimensions.(dims{dIdx}).coord,'sphPhi')
                error('pe_cfg.dimensions.d* of spherical type must have a coord field with sphPhi subfield')
            end
        end

        % sphTheta and sphPhi must be same length as vec
        n_sphTheta = length(pe_cfg.dimensions.(dims{dIdx}).coord.sphTheta);
        n_sphPhi = length(pe_cfg.dimensions.(dims{dIdx}).coord.sphPhi);
        if ~isequal(n_sphTheta, n_sphPhi, length(pe_cfg.dimensions.(dims{dIdx}).vec))
            error('pe_cfg.dimensions.d* of spherical type must have a coord field with sphTheta and sphPhi subfield with same numerosity as pe_cfg.dimensions.d*.vec')
        end

        % if sphRadius is provded (not compulsory), it also must match the
        % same length as vec
        if isfield(pe_cfg.dimensions.(dims{dIdx}).coord,'sphRadius')
            n_sphRadius = length(pe_cfg.dimensions.(dims{dIdx}).coord.sphRadius);
            if ~isequal(n_sphRadius, length(pe_cfg.dimensions.(dims{dIdx}).vec))
                error('pe_cfg.dimensions.d* of spherical type must have a coord field, and it can have a sphRadius field--if it does, it must be same numerosity as pe_cfg.dimensions.d*.vec')
            end
        end
        

    end
end

%% defaults

for dIdx = 1:length(dims)

    % dimension label
    if ~isfield(pe_cfg.dimensions.(dims{dIdx}),'lbl')
        dimLbl = [dims{dIdx}  pe_cfg.dimensions.(dims{dIdx}).type];
        warning(['dimension label not provided by user for ' dims{dIdx} '. Using this label: ' dimLbl])
        pe_cfg.dimensions.(dims{dIdx}).lbl = dimLbl;
    end

    % units
    if ~isfield(pe_cfg.dimensions.(dims{dIdx}),'units')
        dimUnits = 'undefinedunits';
        warning(['dimension units not provided by user for ' dims{dIdx} '. Using this label: ' dimUnits])
        pe_cfg.dimensions.(dims{dIdx}).units = dimUnits;
    end


    % if spherical, is sphRadius not provided, use 1
    if strcmp(pe_cfg.dimensions.(dims{dIdx}).type,'spherical')
        if ~isfield(pe_cfg.dimensions.(dims{dIdx}).coord,'sphRadius')
            pe_cfg.dimensions.(dims{dIdx}).coord.sphRadius = ones(size(pe_cfg.dimensions.(dims{dIdx}).vec));
            warning('sphRadius not provided. using 1 by default')
                    
        end
    end

end



%% add fields

for dIdx = 1:length(dims)

    % include a num field identifying the numerosity of this dimension
    nVec = length(pe_cfg.dimensions.(dims{dIdx}).vec);
    pe_cfg.dimensions.(dims{dIdx}).num = nVec;

end


disp('Diagonale: dimensions validated')
end
