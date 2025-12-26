function di_cfg = di_validateDimensions(di_cfg)
% Validate presence of required fields in di_cfg.dimensions and include basic computations for shortcuts later on.

dimensionTypes = {'continuous' 'spherical' 'categorical'};

if ~isfield(di_cfg,'dimensions')
    error('di_cfg needs field: dimensions');
end

if any(~startsWith(fieldnames(di_cfg.dimensions),'d'))
    error('di_cfg.dimensions dimensions need to start with d');
end

pattern = "d" + digitsPattern();
if any(~contains(fieldnames(di_cfg.dimensions),pattern))
    error('di_cfg.dimensions dimensions need to be labeled as d followed by a number');
end

dims = fieldnames(di_cfg.dimensions);
for dIdx = 1:length(dims)
    if ~any(strcmp(fieldnames(di_cfg.dimensions.(dims{dIdx})),'type'))
        error('di_cfg.dimensions dimensions need to contain a type field');    
    end
    if ~any(strcmp(di_cfg.dimensions.(dims{dIdx}).type,dimensionTypes))
        disp("allowed types: ")
        disp(dimensionTypes)
        error('di_cfg.dimensions dimensions need to contain a type field of any allowed type');    
    end
    if ~any(strcmp(fieldnames(di_cfg.dimensions.(dims{dIdx})),'vec'))
        error('di_cfg.dimensions dimensions need to contain a vec field');    
    end
    if sum(size(di_cfg.dimensions.(dims{dIdx}).vec)~=1)~=1
        error('di_cfg.dimensions.d*.vec must be a vector of non unit length')
    end
    if strcmp(di_cfg.dimensions.(dims{dIdx}).type,'spherical')
        if ~isfield(di_cfg.dimensions.(dims{dIdx}),'coord')
            error('di_cfg.dimensions.d* of spherical type must have a coord field')
        else
            if ~isfield(di_cfg.dimensions.(dims{dIdx}).coord,'sphTheta')
                error('di_cfg.dimensions.d* of spherical type must have a coord field with sphTheta subfield')
            end
            if ~isfield(di_cfg.dimensions.(dims{dIdx}).coord,'sphPhi')
                error('di_cfg.dimensions.d* of spherical type must have a coord field with sphPhi subfield')
            end
        end
        n_sphTheta = length(di_cfg.dimensions.(dims{dIdx}).coord.sphTheta);
        n_sphPhi = length(di_cfg.dimensions.(dims{dIdx}).coord.sphPhi);
        if ~isequal(n_sphTheta, n_sphPhi, length(di_cfg.dimensions.(dims{dIdx}).vec))
            error('di_cfg.dimensions.d* of spherical type must have a coord field with sphTheta and sphPhi subfield with same numerosity as di_cfg.dimensions.d*.vec')
        end
        if isfield(di_cfg.dimensions.(dims{dIdx}).coord,'sphRadius')
            n_sphRadius = length(di_cfg.dimensions.(dims{dIdx}).coord.sphRadius);
            if ~isequal(n_sphRadius, length(di_cfg.dimensions.(dims{dIdx}).vec))
                error('di_cfg.dimensions.d* of spherical type must have a coord field, and it can have a sphRadius field--if it does, it must be same numerosity as di_cfg.dimensions.d*.vec')
            end
        end
    end
end

for dIdx = 1:length(dims)
    if ~isfield(di_cfg.dimensions.(dims{dIdx}),'lbl')
        dimLbl = [dims{dIdx}  di_cfg.dimensions.(dims{dIdx}).type];
        warning(['Diagonale: label not provided by user for dimension ' dims{dIdx} '. Using this label: ' dimLbl])
        di_cfg.dimensions.(dims{dIdx}).lbl = dimLbl;
    end
    if ~isfield(di_cfg.dimensions.(dims{dIdx}),'units')
        dimUnits = 'undefinedunits';
        warning(['Diagonale: units not provided by user for dimension ' dims{dIdx} ' (' di_cfg.dimensions.(dims{dIdx}).lbl ') . Using this label: ' dimUnits])
        di_cfg.dimensions.(dims{dIdx}).units = dimUnits;
    end
    if strcmp(di_cfg.dimensions.(dims{dIdx}).type,'spherical')
        if ~isfield(di_cfg.dimensions.(dims{dIdx}).coord,'sphRadius')
            di_cfg.dimensions.(dims{dIdx}).coord.sphRadius = ones(size(di_cfg.dimensions.(dims{dIdx}).vec));
            warning('sphRadius not provided. using 1 by default')
        end
    end
end

for dIdx = 1:length(dims)
    nVec = length(di_cfg.dimensions.(dims{dIdx}).vec);
    di_cfg.dimensions.(dims{dIdx}).num = nVec;
end


% add a validated flag
di_cfg.validation.dimensions = true;

disp('Diagonale: dimensions validated')
end
