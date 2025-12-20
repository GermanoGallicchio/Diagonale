function rowIdx = dg_resample(dg_cfg)
% Build resampling/permutation indices for rows in designTbl based on the
% objective and design code.

if ~isfield(dg_cfg, 'designTbl') || ~istable(dg_cfg.designTbl)
    error('dg_cfg.designTbl must be a table');
end
designTbl = dg_cfg.designTbl;

if ~all(ismember({'subjID','groupID'}, designTbl.Properties.VariableNames))
    error('designTbl must contain subjID and groupID columns');
end

if ~isfield(dg_cfg,'nIterations') || ~isnumeric(dg_cfg.nIterations) || ~isscalar(dg_cfg.nIterations) || dg_cfg.nIterations<1
    error('dg_cfg.nIterations must be a positive scalar');
end

if ~isfield(dg_cfg,'objective') || (~ischar(dg_cfg.objective) && ~isstring(dg_cfg.objective))
    error('dg_cfg.objective must be a string/char');
end

if ~isfield(dg_cfg,'designCode') || ~isnumeric(dg_cfg.designCode) || numel(dg_cfg.designCode)~=3
    error('dg_cfg.designCode must be a 1x3 numeric vector');
end

nIterations = double(dg_cfg.nIterations);
nRow = height(designTbl);

if isfield(dg_cfg,'randomSeed') && ~isempty(dg_cfg.randomSeed)
    rng(dg_cfg.randomSeed, 'twister');
end

[unique_subjID,~,subjIdxPerRow] = unique(designTbl.subjID,'stable');
nUnique_subjID = numel(unique_subjID);

unique_groupID = string(unique(designTbl.groupID,'stable'));
nUnique_groupID = numel(unique_groupID);

rmF_col_idx = contains(designTbl.Properties.VariableNames,'rmFactor');
nUnique_RMFactors = nnz(rmF_col_idx);
nUnique_RMLevelsPerFactors = nan(1,nUnique_RMFactors);
colCounter = 0;
for rmIdx = find(rmF_col_idx)
    colCounter = colCounter + 1;
    nUnique_RMLevelsPerFactors(1,colCounter) = numel(unique(designTbl{:,rmIdx}));
end
if isempty(nUnique_RMLevelsPerFactors)
    nUnique_RMLevels = 1;
else
    nUnique_RMLevels = prod(nUnique_RMLevelsPerFactors);
end

groupBelonging = strings(nUnique_subjID,1);
for subjIdx = 1:nUnique_subjID
    mask       = designTbl.subjID == unique_subjID(subjIdx);
    groupBelonging(subjIdx) = string(designTbl.groupID(find(mask,1)));
end

rowIdxBySubj = accumarray(subjIdxPerRow, (1:nRow)', [], @(x){x});

groupRowsByG    = cell(nUnique_groupID,1);
origGroupRowCnt = zeros(nUnique_groupID,1);
for gIdx = 1:nUnique_groupID
    subs = find(groupBelonging==unique_groupID(gIdx));
    groupRowsByG{gIdx}    = vertcat(rowIdxBySubj{subs});
    origGroupRowCnt(gIdx) = numel(groupRowsByG{gIdx});
end

key = [char(dg_cfg.objective) ' & ' num2str(dg_cfg.designCode)];

switch key
    case {'permutationH0testing & 1  0  0', ...
          'permutationH0testing & 0  1  0', ...
          'permutationH0testing & 1  1  0'}
        rowIdx = zeros(nRow, nIterations);
        rowIdx(:,1) = (1:nRow)';
        for itIdx = 2:nIterations
            rowIdx(:,itIdx) = randperm(nRow)';
        end

    case {'permutationH0testing & 0  0  1'}
        rowIdx = zeros(nRow, nIterations);
        rowIdx(:,1) = (1:nRow)';
        for itIdx = 2:nIterations
            rowIdx(:,itIdx) = (1:nRow)';
            for subjIdx = 1:nUnique_subjID
                idx = rowIdxBySubj{subjIdx};
                half = numel(idx)/2;
                idx1 = idx(1:half);
                idx2 = idx(half+1:end);
                if rand>0.5
                    rowIdx(idx1,itIdx) = idx2;
                    rowIdx(idx2,itIdx) = idx1;
                end
            end
        end

    case {'bootstrapStability & 1  0  0', ...
          'bootstrapStability & 0  1  0', ...
          'bootstrapStability & 1  1  0', ...
          'bootstrapStability & 0  0  1', ...
          'bootstrapStability & 0  1  1', ...
          'bootstrapStability & 1  0  1', ...
          'bootstrapStability & 1  1  1'}
        rowIdx = zeros(nRow, nIterations);
        rowIdx(:,1) = (1:nRow)';
        for itIdx = 2:nIterations
            rows = [];
            for gIdx = 1:nUnique_groupID
                subs = find(groupBelonging==unique_groupID(gIdx));
                resampSubs = subs(randi(numel(subs), [1 numel(subs)]));
                resampRows = vertcat(rowIdxBySubj{resampSubs});
                rows = [rows; resampRows];
            end
            rowIdx(:,itIdx) = rows;
        end
end

end
