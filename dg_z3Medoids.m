function [z3_MedIdx, z3_wMedIdx] = dg_z3Medoids(dg_parameters, statMat, clustIdxMat,angularDistMat)
% Identifies medoids in sensor space for EEG/MEG data analysis

% (Body identical to original PE_z3Medoids_y1x2z3 with updated names.)

% input validation
narginchk(4,4);
if isempty(statMat) || isempty(clustIdxMat) || isempty(angularDistMat)
    error('input variables mst not be empty');
end
if ~isnumeric(statMat)
    error('statMat must be numeric');
end
if ~(isnumeric(clustIdxMat) || islogical(clustIdxMat))
    error('clustIdxMat must be numeric or logical');
end
if ~isnumeric(angularDistMat)
    error('angularDistMat must be numeric');
end
if any(isnan(statMat(:))) || any(isinf(statMat(:)))
    error('statMat contains NaN or Inf values');
end
if any(isnan(angularDistMat(:))) || any(isinf(angularDistMat(:)))
    error('angularDistMat contains NaN or Inf values');
end
if ~all(ismember(clustIdxMat(:), [0 1]))
    error('clustIdxMat must contain only 0s and 1s');
end
if ~isequal(size(clustIdxMat),size(statMat))
    error('clustIdxMat and statMat must have the same dimensionality')
end
if ~isequal(size(statMat,3),size(angularDistMat,1))
    error('third dimension of both statMat and clustIdxMat should be z3')
end

[ny1, nx2, nz3] = size(statMat);
[y1Matrix, x2Matrix, z3Matrix] = ndgrid(1:ny1, 1:nx2, 1:nz3);

chansInCluster = unique(z3Matrix(clustIdxMat)); 
nChansInCluster = length(chansInCluster);

if dg_parameters.verbose
    if nChansInCluster == 0
        warning('No channels found in cluster');
        z3_MedIdx = [];
        z3_wMedIdx = [];
        return;
    end
end

if nChansInCluster == 1
    z3_MedIdx = chansInCluster;
    z3_wMedIdx = chansInCluster;
    return;
end

angDistSubsetMat = angularDistMat(chansInCluster,chansInCluster);
angDistSubsetVec = sum(angDistSubsetMat,1);

clustStat = statMat(clustIdxMat);
clustChan = z3Matrix(clustIdxMat);
chanRegularMass = nan(1,nChansInCluster);
for chanIdx = 1:nChansInCluster
    chanStat = clustStat(clustChan==chansInCluster(chanIdx));
    chanRegularMass(chanIdx) = sum(chanStat);
end
chanMass = chanRegularMass;
if ~isequal(abs(sum(sign(chanMass))),nChansInCluster)
    warning('mass sign check mismatch; OK if using geometric approach')
end

H = normalize(-abs(chanMass), 'range', [min(angDistSubsetVec) max(angDistSubsetVec)]);

% medoid index
[~, minIdx] = min(angDistSubsetVec);
z3_MedIdx = chansInCluster(minIdx);

% weighted medoid index
[~, minWIdx] = min(angDistSubsetVec .* H);
z3_wMedIdx = chansInCluster(minWIdx);

end
