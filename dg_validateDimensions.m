function dg_cfg = dg_validateDimensions(dg_cfg)
% Validate presence of required fields in dg_cfg.dimensions and include basic computations for shortcuts later on.

dimensionTypes = {'continuous' 'spherical' 'categorical'};

if ~isfield(dg_cfg,'dimensions')
    error('dg_cfg needs field: dimensions');
end

if any(~startsWith(fieldnames(dg_cfg.dimensions),'d'))
    error('dg_cfg.dimensions dimensions need to start with d');
end

pattern = "d" + digitsPattern();
if any(~contains(fieldnames(dg_cfg.dimensions),pattern))
    error('dg_cfg.dimensions dimensions need to be labeled as d followed by a number');
end

dims = fieldnames(dg_cfg.dimensions);
for dIdx = 1:length(dims)
    if ~any(strcmp(fieldnames(dg_cfg.dimensions.(dims{dIdx})),'type'))
        error('dg_cfg.dimensions dimensions need to contain a type field');    
    end
    if ~any(strcmp(dg_cfg.dimensions.(dims{dIdx}).type,dimensionTypes))
        disp("allowed types: ")
        disp(dimensionTypes)
        error('dg_cfg.dimensions dimensions need to contain a type field of any allowed type');    
    end
    if ~any(strcmp(fieldnames(dg_cfg.dimensions.(dims{dIdx})),'vec'))
        error('dg_cfg.dimensions dimensions need to contain a vec field');    
    end
    if sum(size(dg_cfg.dimensions.(dims{dIdx}).vec)~=1)~=1
        error('dg_cfg.dimensions.d*.vec must be a vector of non unit length')
    end
    if strcmp(dg_cfg.dimensions.(dims{dIdx}).type,'spherical')
        if ~isfield(dg_cfg.dimensions.(dims{dIdx}),'coord')
            error('dg_cfg.dimensions.d* of spherical type must have a coord field')
        else
            if ~isfield(dg_cfg.dimensions.(dims{dIdx}).coord,'sphTheta')
                error('dg_cfg.dimensions.d* of spherical type must have a coord field with sphTheta subfield')
            end
            if ~isfield(dg_cfg.dimensions.(dims{dIdx}).coord,'sphPhi')
                error('dg_cfg.dimensions.d* of spherical type must have a coord field with sphPhi subfield')
            end
        end
        n_sphTheta = length(dg_cfg.dimensions.(dims{dIdx}).coord.sphTheta);
        n_sphPhi = length(dg_cfg.dimensions.(dims{dIdx}).coord.sphPhi);
        if ~isequal(n_sphTheta, n_sphPhi, length(dg_cfg.dimensions.(dims{dIdx}).vec))
            error('dg_cfg.dimensions.d* of spherical type must have a coord field with sphTheta and sphPhi subfield with same numerosity as dg_cfg.dimensions.d*.vec')
        end
        if isfield(dg_cfg.dimensions.(dims{dIdx}).coord,'sphRadius')
            n_sphRadius = length(dg_cfg.dimensions.(dims{dIdx}).coord.sphRadius);
            if ~isequal(n_sphRadius, length(dg_cfg.dimensions.(dims{dIdx}).vec))
                error('dg_cfg.dimensions.d* of spherical type must have a coord field, and it can have a sphRadius field--if it does, it must be same numerosity as dg_cfg.dimensions.d*.vec')
            end
        end
    end
end

for dIdx = 1:length(dims)
    if ~isfield(dg_cfg.dimensions.(dims{dIdx}),'lbl')
        dimLbl = [dims{dIdx}  dg_cfg.dimensions.(dims{dIdx}).type];
        warning(['Diagonale: label not provided by user for dimension ' dims{dIdx} '. Using this label: ' dimLbl])
        dg_cfg.dimensions.(dims{dIdx}).lbl = dimLbl;
    end
    if ~isfield(dg_cfg.dimensions.(dims{dIdx}),'units')
        dimUnits = 'undefinedunits';
        warning(['Diagonale: units not provided by user for dimension ' dims{dIdx} ' (' dg_cfg.dimensions.(dims{dIdx}).lbl ') . Using this label: ' dimUnits])
        dg_cfg.dimensions.(dims{dIdx}).units = dimUnits;
    end
    if strcmp(dg_cfg.dimensions.(dims{dIdx}).type,'spherical')
        if ~isfield(dg_cfg.dimensions.(dims{dIdx}).coord,'sphRadius')
            dg_cfg.dimensions.(dims{dIdx}).coord.sphRadius = ones(size(dg_cfg.dimensions.(dims{dIdx}).vec));
            warning('sphRadius not provided. using 1 by default')
        end
    end
end

for dIdx = 1:length(dims)
    nVec = length(dg_cfg.dimensions.(dims{dIdx}).vec);
    dg_cfg.dimensions.(dims{dIdx}).num = nVec;
end


% add a validated flag
dg_cfg.validation.dimensions = true;

disp('Diagonale: dimensions validated')
end
