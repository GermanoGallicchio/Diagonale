function [rowIdx, di_cfg] = di_reorderRowsGenerate(di_cfg)
% Generate row reordering indices for 
% - permutation for inferential testing 
% - bootstrap for stability / precision estimation 
%
% TO DO (not yet available):
% - cross-validation for generalizability estimation // do this in another
% function. for example, the fact that iteration 1 is identical to the
% origianl won't work for cross validation
%
% row indices are based on dataStruct, analysis objective, and (for
% permutation) an explicit H0 hypothesis.
%
% the rowIdx output is then verified, and its verification checks can be seen when di_cfg.analysis.verbose = true
%
% INPUT:
%   di_cfg              - struct with fields:
%     .analysis.dataStruct - table with columns: observationID, indFactor#, repFactor#
%     .analysis.nIterations - positive integer, number of resampling iterations
%     .analysis.objective   - 'permutationH0testing' or 'bootstrapStability'
%     .analysis.randomSeed  - (optional) integer for reproducibility
%     .analysis.H0hypothesis - required when objective='permutationH0testing'
%                              allowed values:
%                              'between', 'within', 'within-by-group',
%                              'between-by-condition'
%
% OUTPUT:
%   rowIdx              - (nRows x nIterations) matrix of row indices
%                         Column 1 is always unchanged (original data)
%                         Columns 2-nIterations contain permuted/resampled indices
%
% Notes:
% 1. The first iteration is identical to the original data (no permutation
% or resampling)
% 2. PERMUTATION is controlled by di_cfg.analysis.H0hypothesis
%    'between'               : shuffle entire observations across positions
%                              (preserve repeated measures, if there are any, within each observation)
%    'within'                : shuffle repeated measures within each observation
%                              (preserve observation identities)
%    'within-by-group'       : same permutation mechanics as 'between'
%    'between-by-condition'  : same permutation mechanics as 'within'
%
% 3. BOOTSTRAP is always done by resamples with replacement full observation
% units (eg, a subject with all its conditions)
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% sanity checks

% there is a random seed
if ~isfield(di_cfg.analysis, 'randomSeed')
    error('di_cfg.analysis.randomSeed is missing. Run di_validateAnalysis first.');
end

% permutation requires explicit H0 hypothesis
if ~isfield(di_cfg.validation,'H0hypothesis') || di_cfg.validation.H0hypothesis~=1
    error('null hypothesis not validated. Likely a coding bug due to misplaced di_validateH0 in the pipeline');
end
%% shortcuts

dataStruct = di_cfg.analysis.dataStruct;
nIterations = double(di_cfg.analysis.nIterations);
nRow = height(dataStruct);
H0hypothesis = di_cfg.analysis.H0hypothesis;

%% set random seed
rng(double(di_cfg.analysis.randomSeed), 'twister');

if isfield(di_cfg.analysis, 'verbose') && di_cfg.analysis.verbose
    fprintf('random seed used for row reordering: %d\n', double(di_cfg.analysis.randomSeed));
end

%% independent observation factors

% independent observations
[unique_obsID,~,obsIdxPerRow] = unique(di_cfg.analysis.dataStruct.observationID,'stable');
nUnique_obsID = numel(unique_obsID);

% independent group factor columns (if any)
indF_col_idx = startsWith(di_cfg.analysis.dataStruct.Properties.VariableNames, 'indFactor');
nUnique_indFactors = nnz(indF_col_idx);

if nUnique_indFactors > 0
    % create group identifier from indFactor columns
    indFactorCols = di_cfg.analysis.dataStruct(:, indF_col_idx);
    unique_groups = unique(indFactorCols, 'rows', 'stable');
    nUnique_groups = height(unique_groups);
else
    % no groups of independent observatioons. treat all observations as if they belong to the same group
    nUnique_groups = 1;
end

%% repeated measure factors

rmF_col_idx = startsWith(di_cfg.analysis.dataStruct.Properties.VariableNames, 'repFactor');
nUnique_repFactors = nnz(rmF_col_idx);
nUnique_RMLevelsPerFactors = nan(1, nUnique_repFactors);
colCounter = 0;
for rmIdx = find(rmF_col_idx)
    colCounter = colCounter + 1;
    nUnique_RMLevelsPerFactors(1, colCounter) = numel(unique(di_cfg.analysis.dataStruct{:, rmIdx}));
end
if isempty(nUnique_RMLevelsPerFactors)
    nUnique_RMLevels = 1;
else
    nUnique_RMLevels = prod(nUnique_RMLevelsPerFactors);
end

%% belonging of each observation unit (eg, subject) to its respective group
% only for bootstrap (not for permutation)
% this cell enables stratified bootstrap by identifying which group each
% observation belongs to, so each observation (eg, subject) is resampled within their respective group

% map each observation (e.g., subject) to its group membership based on 
% indFactor columns (e.g., group).
% each cell contains the indFactor values (group identifier)
% for one observation, or [] if no independent grouping factors exist.
groupBelonging = cell(nUnique_obsID, 1);
for obsIdx = 1:nUnique_obsID
    obsIS_rowIdx = dataStruct.observationID == unique_obsID(obsIdx);
    if nUnique_indFactors > 0
        % store the group identifier (indFactor values) for this observation
        % take first row since all rows for an observation have identical group membership
        groupBelonging{obsIdx} = dataStruct{find(obsIS_rowIdx, 1), indF_col_idx};
    else
        % No groups
        groupBelonging{obsIdx} = [];
    end
end

%% belonging of each repeated measure row to its observation unit (eg, subject)
% create a lookup cell mapping each observation ID to the row indices in
% dataStruct associated with it.
% useful for:
% 1- permutation with within-observation H0 hypotheses 
%   to shuffle conditions within each subject
% 2- bootstrap resampling (all design codes)
%   to (re)sample entire observations and letting each observation come
%   along all their repeated measures
rowIdxByObs = accumarray(obsIdxPerRow, (1:nRow)', [], @(x){x});
% obsIdxPerRow maps each row of dataStruct to the ordinal index of its unique observationID (1..nUnique_obsID),
% rowIdxByObs{idx} contains a vector of row indices belonging to each unique observation ID.
% eg: if subject 1 has rows [1,2,3] and subject 2 has rows [4,5,6], then rowIdxByObs{1} = [1;2;3] and rowIdxByObs{2} = [4;5;6]



%% delete this entire cell // obsolete/superseded by previous cell
% organize row indices by group membership for stratified resampling.
% legacy code from a previous implementation
% groupRowsByG    = cell(nUnique_groups, 1);
% origGroupRowCnt = zeros(nUnique_groups, 1);
% if nUnique_indFactors > 0
%     for gIdx = 1:nUnique_groups
%         % Find observations belonging to this group by comparing their 
%         % indFactor values (stored in groupBelonging) with the current group's identifier
%         group_mask = true(nUnique_obsID, 1);
%         for obsIdx = 1:nUnique_obsID
%             if isempty(groupBelonging{obsIdx})
%                 group_mask(obsIdx) = false;
%             else
%                 % Check if this observation's group matches current group
%                 match = all(groupBelonging{obsIdx} == unique_groups{gIdx,:}, 2);
%                 group_mask(obsIdx) = any(match);
%             end
%         end
%         % Collect all row indices from observations in this group
%         obs_in_group = find(group_mask);
%         groupRowsByG{gIdx} = vertcat(rowIdxByObs{obs_in_group});
%         origGroupRowCnt(gIdx) = numel(groupRowsByG{gIdx});
%     end
% else
%     % Single group: all observations belong to one group
%     groupRowsByG{1} = (1:nRow)';
%     origGroupRowCnt(1) = nRow;
% end

%% generate indices

% note: 
% permutation adn bootstrap leave the first iteration as the original 
% this is intentional, to make the code work also to get results on the
% original data (itIdx==1) to compare against all iterations
% (itIdx==1:nIterations)

switch char(di_cfg.analysis.objective)
    case 'permutationH0testing'
        switch H0hypothesis
            case {'between', 'within-by-group'}
                % Shuffle entire observations across positions while preserving
                % each observation's full repeated-measures structure.
                rowIdx = zeros(nRow, nIterations);
                rowIdx(:,1) = (1:nRow)';
                for itIdx = 2:nIterations
                    permuted_obsIdx = randperm(nUnique_obsID);
                    rowIdx(:,itIdx) = vertcat(rowIdxByObs{permuted_obsIdx});
                end

            case {'within', 'between-by-condition'}
                % Shuffle repeated measures within each observation; observation
                % identities are preserved.
                rowIdx = zeros(nRow, nIterations);
                rowIdx(:,1) = (1:nRow)';
                for itIdx = 2:nIterations
                    rowIdx(:,itIdx) = (1:nRow)';
                    for obsIdx = 1:nUnique_obsID
                        idx = rowIdxByObs{obsIdx};
                        idx_shuffled = idx(randperm(length(idx)));
                        rowIdx(idx,itIdx) = idx_shuffled;
                    end
                end
        end

    case 'bootstrapStability'
        rowIdx = zeros(nRow, nIterations); % initialize
        rowIdx(:,1) = (1:nRow)'; % original
        for itIdx = 2:nIterations
            rows = []; % initialize
            for gIdx = 1:nUnique_groups % loop through groups for stratified bootstrap
                if nUnique_indFactors > 0
                    % filter to observations that belong to this specific group
                    targetGroup = unique_groups{gIdx,:};
                    obs_in_group = find(cellfun(@(x) ~isempty(x) && isequal(x, targetGroup), groupBelonging));
                else
                    % no independent-group factor: all observations form one pool
                    obs_in_group = (1:nUnique_obsID)';
                end
                resampObs = obs_in_group(randi(numel(obs_in_group), [1 numel(obs_in_group)])); % sample with replacement among obsIDs (eg, subjects) belonging to group gIdx
                resampRows = vertcat(rowIdxByObs{resampObs}); % make sure each obsID carries their repeated measures
                rows = [rows; resampRows]; % concatenate rows
            end
            rowIdx(:,itIdx) = rows;
        end

    otherwise
        error('Unsupported analysis objective: %s', char(di_cfg.analysis.objective));
end

%% sanity checks: Verify resampling respects data structure

% sanity check: first iteration must be original data (unchanged)
if ~isequal(rowIdx(:,1), (1:nRow)')
    error('First iteration (rowIdx(:,1)) is not the original data sequence (1:nRow)');
end

% sanity check: for permutation, all rows are used exactly once per iteration
if strcmp(di_cfg.analysis.objective, 'permutationH0testing')
    for itIdx = 1:nIterations
        used_rows = rowIdx(:,itIdx);
        if length(unique(used_rows)) ~= nRow || length(used_rows) ~= nRow
            error('Iteration %d: Permutation did not use all rows exactly once. Duplicate or missing rows detected.', itIdx);
        end
    end
    if di_cfg.analysis.verbose
    fprintf('permutation verification: all rows used exactly once per iteration\n');
    end
end

% sanity check: for permutation with repeated measures, repeated measures still belong to the original observation
% this verifies that within-subject shuffling preserves observation identity
% (conditions are shuffled within each subject, but no condition moves to a different subject)
if strcmp(di_cfg.analysis.objective, 'permutationH0testing')
    is_within_shuffle = any(strcmp(di_cfg.analysis.H0hypothesis, {'within', 'between-by-condition'}));
    if is_within_shuffle % H0 implies within-observation shuffling
    for itIdx = 2:nIterations % start from 2 since itIdx 1 is original data
        for obsIdx = 1:nUnique_obsID
            % get the row indices that originally belong to this observation
            original_obs_rows = rowIdxByObs{obsIdx};
            
            % get the permuted row indices for this observation in this iteration
            % rowIdx(original_obs_rows, itIdx) tells us which rows were shuffled into these positions
            permuted_obs_rows = rowIdx(original_obs_rows, itIdx);
            
            % check 1: all rows for this observation should still belong to it (same count)
            % if permutation is correct, the number of rows should not change
            if length(permuted_obs_rows) ~= length(original_obs_rows)
                error('Iteration %d, Observation %d: Permutation broke observation structure. Expected %d rows, got %d.', ...
                    itIdx, obsIdx, length(original_obs_rows), length(permuted_obs_rows));
            end
            
            % check 2: verify permuted rows are all from the same observation in dataStruct
            % extract the observationIDs from the permuted row positions
            obs_ids_in_permuted = dataStruct.observationID(permuted_obs_rows);
            % all should match this observation's ID (within-subject shuffle only)
            if ~all(obs_ids_in_permuted == unique_obsID(obsIdx))
                error('Iteration %d, Observation %d: Permutation violated within-observation structure.', itIdx, obsIdx);
            end
        end
    end
    if di_cfg.analysis.verbose
    fprintf('permutation with repeated measures verification: observation structure preserved\n');
    end
    end
end

% sanity check: for between-observation permutation hypotheses
% verify that each observation appears exactly once and its full set of repeated measures stays intact
if strcmp(di_cfg.analysis.objective, 'permutationH0testing')
    if any(strcmp(di_cfg.analysis.H0hypothesis, {'between', 'within-by-group'}))
        for itIdx = 2:nIterations % start from 2 since itIdx 1 is original data
            permuted_obsIDs = dataStruct.observationID(rowIdx(:, itIdx));
            for obsIdx = 1:nUnique_obsID
                original_count = numel(rowIdxByObs{obsIdx});
                current_count = sum(permuted_obsIDs == unique_obsID(obsIdx));
                if current_count ~= original_count
                    error('Iteration %d, Observation %d: H0_between permutation broke observation integrity. Expected %d rows, got %d.', ...
                        itIdx, obsIdx, original_count, current_count);
                end
            end
        end
        if di_cfg.analysis.verbose
            fprintf('permutation mixed design (H0_between) verification: observations preserved with full repeated measures\n');
        end
    end
end

% sanity check: for bootstrap, verify entire observations are (re)sampled together (with replacement)
% this ensures that when an observation (eg, subject) is resampled, ALL its conditions come along.
% observations(eg, subjects) can appear 0, 1, or multiple times (with replacement), but must always be complete
if strcmp(di_cfg.analysis.objective, 'bootstrapStability')
    for itIdx = 2:nIterations % start from 2 since itIdx 1 is original data
        % get the resampled row indices and their corresponding observation IDs
        resampled_rows = rowIdx(:, itIdx);
        resampled_obsIDs = dataStruct.observationID(resampled_rows);
        
        % check that observations appear in contiguous blocks (all conditions together)
        % when an observation (eg, subject) is sampled, all its repeated measure rows (eg, conditions) should appear consecutively,
        % just because they were done one obsID at the time and concatenating their row indices.
        current_obsID = resampled_obsIDs(1); % start with first observation in resampled data
        for row_idx = 2:nRow
            obs_in_current_row = resampled_obsIDs(row_idx);
            
            % if observation changes, update the tracker
            % (we're moving to a different observation in the resampled sequence)
            if obs_in_current_row ~= current_obsID
                current_obsID = obs_in_current_row;
            end
            
            % Verify this observation's rows appear contiguously.
            % find all positions where this observation appears in the resampled data
            obs_row_locs = find(resampled_obsIDs == current_obsID);
            
            % criterion 1: check if rows are not all contiguous (diff not all equal to 1)
            % this detects when an observation appears in non-consecutive row positions
            criterion1 = ~all(diff(obs_row_locs) == 1);
            
            % criterion 2: check if observation appears more than once in resampled data
            % this detects when an observation was sampled multiple times (with replacement)
            criterion2 = length(obs_row_locs) > 1;
            
            if criterion1 && criterion2
                % both criteria true: observation appears non-contiguously and multiple times
                % this is allowed because with replacement resampling, an observation can appear 0, 1, or k times
                % (eg, subject 1 with 3 conditions could appear as rows [1,2,3] and [7,8,9])
                % however, we must validate criteria 3 and 4 to ensure resampling integrity
                
                % get the original number of rows (conditions) for this observation
                original_count = length(rowIdxByObs{current_obsID == unique_obsID});
                % get how many times this observation appears in total in resampled data
                current_count = length(obs_row_locs);
                
                % criterion 3: current_count must be a multiple of original_count (violation if NOT)
                % with replacement resampling: if observation sampled k times, expect k x original_count rows total
                % example: if original_count=3 (3 conditions), current_count should be 3, 6, 9, 12, ... (not 4, 5, 7, ...)
                % a non-multiple value indicates rows were lost, duplicated incorrectly, or misaligned
                criterion3 = mod(current_count, original_count) ~= 0;
                
                if criterion3
                    error('Iteration %d: Bootstrap violated observation integrity for observationID %d.', ...
                        itIdx, current_obsID);
                end
                
                % criterion 4: all original conditions must be present in resampled data (violation if missing any)
                % checks that no condition was omitted or oversampled exclusively
                % rationale: even if counts are correct (criterion3), a condition could be completely missing
                % example: 3 conditions {A,B,C} resampled could mistakenly return {A,B,B,C,C,C} (C oversampled and A undersampled)
                if nUnique_repFactors > 0
                    % extract all conditions present in original observation by finding its rows and getting unique repFactor values
                    original_repFactorRows = find(dataStruct.observationID == current_obsID);
                    orig_repFactorVals = dataStruct(original_repFactorRows, rmF_col_idx);
                    orig_repFactorUnique = unique(orig_repFactorVals, 'rows', 'stable');
                    
                    % extract all conditions present in resampled observation using the row positions found earlier
                    % resampled_rows(obs_row_locs) gives the actual row indices in resampled data for this observation
                    resampled_repFactorRows = resampled_rows(obs_row_locs);
                    resampled_repFactorVals = dataStruct(resampled_repFactorRows, rmF_col_idx);
                    resampled_repFactorUnique = unique(resampled_repFactorVals, 'rows', 'stable');
                    
                    % compare as tables to support numeric and non-numeric
                    % repeated-factor values (e.g., strings/cellstr) without
                    % triggering ismember(...,'rows') warnings on cell arrays.
                    nMismatch = height(orig_repFactorUnique) ~= height(resampled_repFactorUnique);
                    missingInResampled = ~isempty(setdiff(orig_repFactorUnique, resampled_repFactorUnique, 'rows'));
                    % criterion4 is a violation if:
                    % - number of unique conditions differs: original had N, resampled has M ≠ N (missing or extra conditions)
                    % - OR: at least one original condition is absent from resampled data
                    criterion4 = nMismatch || missingInResampled;
                else
                    % no repeated measures (only 1 condition per observation), so all conditions trivially present
                    criterion4 = false;
                end
                
                if criterion4
                    error('Iteration %d: Bootstrap missing conditions for observationID %d.', ...
                        itIdx, current_obsID);
                end
            end
        end
    end
    if di_cfg.analysis.verbose
    fprintf('bootstrap verification: entire observations sampled with replacement\n');
    end
end

end
