function results = di_analysis_plsSVD(di_cfg, Y_orig, X_orig, rowIdx)
% Partial Least Squares via Singular Value Decomposition with
% permutation testing or bootstrap stability estimation.
% [U S V] = svd(X' * Y), where X is design/contrast/data and Y is main data
% note: I reversed what McIntosh does to keep a similar notation with OLS
%
% Cross-covariance matrix between X and Y is first computed. both Y and X
% contribute to latent structure. That matrix is decomposed via SVD into
% multiple modes. Each mode has a singular value (its magnitude) and two
% singular vectors (one for X and one for Y mapping the contribution of
% those variables to the latent structure). Permutation testing is used to
% test the strength of the modes' singular value. Bootstrap estimation is
% used, only for significant modes, to identify original features of Y
% (potentially of X too) that are stable under sampling variability.
%
% notes:
% 1. Y and X are always column-mean-centered at each iteration. Z-scoring
% is optional and controlled by di_cfg.analysis.plssvdParams.zscoringVec.
% 2. If X is a matrix of contrasts, it's useful to use orthogonal (full rank) and possibly
% orthonormal (each column has same Euclidean length) so that the
% cross-covariance matrix with Y represents various directions with no bias
% for one or another.
% 2b. If the directions of X were not really the focus on the study, the
% results (U loadings) can be projected to a matrix with contrasts that
% were more interesting to the user. This currently has to be done outside
% Diagonale... not a big job in terms of computation and code. I might
% embed this as a feature in a future iteration
%
% INPUT:
%   di_cfg        - validated analysis configuration structure
%   Y_orig        - data matrix (m x pY) where pY = product of all dimensions
%   X_orig        - codes / design matrix
%   rowIdx        - resampling row indices from di_reorderRowsGenerate
%
% OUTPUT:
%   results       - unified results structure:
%
%     .observed.statVal       (nModes x 1) observed singular values
%     .observed.clusters      empty struct (PLS operates on modes, not spatial features)
% TO DO: I might want to create clusters from reliable features (e.g.,|BR|>2)
%
%     .PLS_SVD      loading level (feature level)
%       .loadings.U_obs       (pX x nModes) X-side loadings in X-variable space
%       .loadings.V_obs       (pY x nModes) Y-side loadings in Y-variable space
%       .loadings.XU_obs      (m x nModes) X latent scores
%       .loadings.YV_obs      (m x nModes) Y latent scores
%
%     .PLS_SVD      mode level
%       .modes.nModes         number of extracted latent variables
%       .modes.s              (1 x nModes) singular values
%       .modes.p              (1 x nModes) proportion of covariance per mode
%       .modes.wilk           (1 x nModes) cumulative variance from end
%       .modes.inertia        (1 x 1) total covariance/inertia
%       .modes.r              (1 x nModes) correlation between Y and X scores
%
%     .resampling   mode level inference (permutation)
%         .permutationH0.modes.s          (nIterations x nModes) singular values
%         .permutationH0.modes.wilk       (nIterations x nModes) Wilks lambda
%
%     .resampling   loading level (feature level) inference (bootstrap)
%         .bootstrapStability.loadings.U_boot   (pX x nModes x nIterations)
%         .bootstrapStability.loadings.V_boot   (pY x nModes x nIterations)
%         .bootstrapStability.loadings.XU_boot  (m x nModes x nIterations)
%         .bootstrapStability.loadings.YV_boot  (m x nModes x nIterations)
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% shortcuts

% number of iterations
nIterations = di_cfg.analysis.nIterations;

% total number of features (product of all dimension sizes)
dimKeys  = fieldnames(di_cfg.dimensions);
dimSizes = cellfun(@(k) length(di_cfg.dimensions.(k).vec), dimKeys);
pY       = prod(dimSizes);  % total number of Y features

% analysis objective and type
analysisObjective = di_cfg.analysis.objective;
inferenceLevel    = di_cfg.analysis.inferenceLevel;
designCode        = di_cfg.analysis.designCode;

% number of observations (rows in Y and X)
m = size(X_orig, 1);

% subject-centering mode shortcut
if isfield(di_cfg.analysis, 'subjectCentering')
    subjectCentering = di_cfg.analysis.subjectCentering;
else
    subjectCentering = '';
end

%% sanity checks
% validate input dimensions

% --- general parameters ---

% num of dimensions must be at least 1
if pY < 1
    error('pY must be at least 1 (total number of features)');
end

% num of dimensions must be same as num of features in Y matrix
if size(Y_orig, 2) ~= pY
    error('Y_orig must have pY columns, same number as features');
end

% num of iterations must be at least 1
if nIterations < 1
    error('nIterations must be at least 1');
end

% num of iterations must be same as num of columns of rowIdx (resampling indices)
if size(rowIdx, 2) ~= nIterations
    error('rowIdx must have nIterations columns');
end

% num of rows in Y must be same as num of rows of rowIdx (resampling indices)
if size(rowIdx, 1) ~= size(Y_orig, 1)
    error('rowIdx must have same number of rows as Y_orig');
end

% X and Y same number of rows
if size(X_orig, 1) ~= size(Y_orig, 1)
    error('Y_orig and X_orig must have same number of rows');
end

% a shortcut
obsID_base = di_cfg.analysis.dataStruct.observationID;

% --- specific to this analysis ---

% shortcuts
isSubjectCentering = strcmp(subjectCentering, 'center');
isBootstrapStability = strcmp(analysisObjective, 'bootstrapStability');

% make sure X_orig is numeric and not logic (it could happen for certain
% contrasts)
X_orig = double(X_orig);

% check that PLS SVD parameters are configured
if ~isfield(di_cfg.analysis, 'plssvdParams')
    error('di_cfg.analysis.plssvdParams must be defined for %s + %s', analysisObjective, inferenceLevel);
end

% z-scoring options used throughout
zscoringVec = di_cfg.analysis.plssvdParams.zscoringVec;

% User-controlled subject-centering mode (for Y side).
% For repeated-measures designs, this must be set explicitly.
hasRepeatedComponent = (designCode(1) == 1);
if isempty(subjectCentering)
    if hasRepeatedComponent
        error(['PLS-SVD with repeated-measures designs requires an explicit subject-centering choice. ' ...
            'Set di_cfg.analysis.subjectCentering to ''center'' or ''noCenter''.']);
    else
        subjectCentering = 'noCenter';
    end
end

if ~any(strcmp(subjectCentering, {'center', 'noCenter'}))
    error('\\ Invalid di_cfg.analysis.subjectCentering in di_analysis_plsSVD. Use ''center'' or ''noCenter''.');
end

if ~hasRepeatedComponent && strcmp(subjectCentering, 'center')
    error('\\ di_cfg.analysis.subjectCentering=''center'' requires a repeated-measures component in the design.');
end

%% computed expected number of modes

% TO DO: improve this section (low priority) as a useful sanity check.
% currently commented out because I might need to add some tolerance to
% match size(S,1) below

% % compute mode count on preprocessed matrices (same transform logic used in-loop)
% % Y side = data matrix (features), X side = design/contrast matrix
% Y_rank_basis = normalize(Y_orig, 'center');
% if zscoringVec(1)
%     Y_rank_basis = zscore(Y_orig);
% end
% 
% X_rank_basis = normalize(X_orig, 'center');
% if zscoringVec(2)
%     X_rank_basis = zscore(X_orig);
% end
% rankX = rank(X_rank_basis);
% 
% rankY = rank(Y_rank_basis);
% 
% % compute num of modes from both sides
% nModes = min(rankX, rankY);

%% main analysis

% For subject-centering in permutation mode, observation mapping is fixed.
if isSubjectCentering && ~isBootstrapStability
    obsID_for_Y_const = obsID_base;
end

for itIdx = 1:nIterations

    % apply row reordering based on indices in rowIdx
    [Y, X] = di_reorderRowsApply(di_cfg, Y_orig, X_orig, rowIdx, itIdx);

    % optional subject-centering of Y (user-controlled via subjectCentering flag).
    if isSubjectCentering

        if isBootstrapStability
            obsID_for_Y = obsID_base(rowIdx(:, itIdx));
        else
            % For permutationH0testing X is not reordered by design; keep original mapping.
            obsID_for_Y = obsID_for_Y_const;
        end

        uniqueObs = unique(obsID_for_Y, 'stable');
        for obsIdx = 1:numel(uniqueObs)
            maskObs = (obsID_for_Y == uniqueObs(obsIdx));
            Y(maskObs, :) = Y(maskObs, :) - mean(Y(maskObs, :), 1);
        end
    end

    % mean center columns of Y and X
    Yz = normalize(Y,'center');
    Xz = normalize(X,'center');

    % zscore columns of Y, depending on zscore choice
    if zscoringVec(1)
        Yz = zscore(Y,0,1);
    end
    % zscore columns of X, depending on zscore choice
    if zscoringVec(2)
        Xz = zscore(X,0,1);
    end

    % apply ignore mask (features / columns of Y)
    Yz(:,logical(di_cfg.analysis.ignore_col)) = 0;

    % temporary section for debugging // can be kept of removed...
    if di_cfg.analysis.figFlag  &&  itIdx==1
        figure()
        imagesc(Yz)
        title('structure of X at itIdx==1 after preprocessing')
        %axis equal
        box off
        axis off
        colormap(parula)
        colorbar
        drawnow
    end

    % cross-product (product of design/data transpose and data)
    C = (Xz' * Yz)/(m-1);

    % SVD
    % "econ" for reduced SVD (only non-zero singular values) 
    [U,S,V] = svd(C,"econ");

    % number of modes as computed in SVD
    nModes = size(S,1);

    % loadings sign and direction alignment (only for bootstrap)
    if strcmp(analysisObjective, 'bootstrapStability') && nModes > 1
        if itIdx == 1 % reference loadings defined locally in this section.
            U_ref = U;
            V_ref = V;
        else

            % SVD has sign ambiguity in the singular vectors (loadings):
            % [U,S,V] and [-U,S,-V] are equivalent solutions. The choice
            % can be taken intentionally (by convention) based on the observed data
            % (1st iteration).
            % for permutation it's not necessary to align the sign
            % ebcaseu we only care about the singular values, which are unsigned
            % and untouched by the sign alignment anyway. Also, due to the
            % permutation before the computation of the cross-covariance, we don't
            % necessarily expect similar singular vectors across iterations.
            % For bootstrap it is critical as this analysis only looks at the
            % singular vectors (loadings) and looks at the consistency in the
            % values, therefore an unaccounted sign change from x to -x would be a
            % big problem. Also, we can expect similar singular vectors across
            % multiple bootstrap iterations as we are not generating a null
            % hypothesis distribution and instead maintaining the mapping between X
            % and Y.
            % implementation: use the 1st iteration loadings as reference, and for
            % subsequent iterations, match the signs to that of the reference
            % loading. SVD sign ambiguity is coupled: to maintain the same dot-product,
            % if U and V need flipping, both need to flip. under no case, only one is flipped
            [U, V] = di_signAlign(U, V, U_ref, V_ref);

            % for bootstrap with nModes>1, we also need to match direction across iterations.
            % This can easily be done through another SVD ("procrustes" rotation... terrible name)
            [U, V] = di_procrustesRotate(U_ref, V_ref, U, V);

        end
    end

    % metrics: projection of X and Y on each mode (i.e., latent scores)
    XU = Xz * U;    % X latent scores (one per mode)
    YV = Yz * V;    % Y latent scores (one per mode)

    % metrics: correlation between the two latent scores, per mode
    r = diag(corr(XU, YV)); % Pearson correlation between X and Y latent variables (i.e., scores = data * singular vectors) (one per mode)
    r = r(:); % enforce column vector

    % metrics: singular value based
    s = diag(S); % singular value (one per mode)
    s = s(:); % enforce column vector
    p = s.^2 / sum(s.^2); % proportion of variance explained by each mode (only meaningful for multiple modes)
    p = p(:); % enforce column vector

    % store original-data values (observed)
    if itIdx==1
        U_obs  = U;   % X-side loadings (U singular vectors)
        V_obs  = V;   % Y-side loadings (V singular vectors)
        XU_obs = XU;  % latent scores: X * singular vectors for X
        YV_obs = YV;  % latent scores: Y * singular vectors for Y
        s_obs = s';   % as row (one per mode)
        p_obs = p';   % as row (one per mode)
        r_obs = r';   % as row (one per mode)
    end


    % store output of each iteration depending on objective
    switch analysisObjective
        case 'permutationH0testing'

            % initialize mode-level metrics
            if itIdx==1
                resampling_s          = zeros(nIterations, nModes);  % singular values per mode
                resampling_wilk       = zeros(nIterations, nModes);  % Wilks lambda
            end

            % fill each iteration
            resampling_s(itIdx, :)          = s';
            resampling_wilk(itIdx, :)       = flipud(cumsum(flipud(s(:).^2)))';

        case 'bootstrapStability'

            % initialize loading wise metrics
            if itIdx==1
                U_boot  = zeros(size(X_orig,2),nModes,nIterations);
                V_boot  = zeros(size(Y_orig,2),nModes,nIterations);
                XU_boot = zeros(size(X_orig,1),nModes,nIterations);
                YV_boot = zeros(size(Y_orig,1),nModes,nIterations);
            end

            % fill each iteration
            U_boot(:,:,itIdx) = U;
            V_boot(:,:,itIdx) = V;
            XU_boot(:,:,itIdx) = XU;
            YV_boot(:,:,itIdx) = YV;
    end

    di_counter(itIdx,nIterations)  % iteration counter
end


%% collate results into unified structure

% --- observed data ---

results.observed.statVal  = [];             % empty: PLS doesn't produce univariate statistics
results.observed.clusters = struct();       % empty for now...
% TO DO: I might want to embed some post-hoc clustering based on |BR| > 2 purely for descriptive reasons

% feature / loading level (contributions to latent variables)
results.PLS_SVD.loadings.U_obs    = U_obs;          % (pX x nModes) X-side loadings in X-variable space
results.PLS_SVD.loadings.V_obs    = V_obs;          % (pY x nModes) Y-side loadings in Y-variable space
results.PLS_SVD.loadings.XU_obs   = XU_obs;         % (m x nModes) X latent scores
results.PLS_SVD.loadings.YV_obs   = YV_obs;         % (m x nModes) Y latent scores

% singular value level (overall latent variable metrics)
results.PLS_SVD.modes.nModes   = nModes;         % number of latent variables
results.PLS_SVD.modes.s       = s_obs;          % (1 x nModes) singular values
results.PLS_SVD.modes.p       = p_obs;          % (1 x nModes) proportion of covariance explained
results.PLS_SVD.modes.wilk    = flipud(cumsum(flipud(s_obs(:).^2)))'; % (1 x nModes) cumulative variance from end
results.PLS_SVD.modes.inertia = sum(s_obs.^2); % (1 x 1) total covariance/inertia in observed data
results.PLS_SVD.modes.r       = r_obs;          % (1 x nModes) correlation between Y and X scores


% --- simulated data ---

switch analysisObjective
    case 'permutationH0testing'
        % permutation testing: SV-level inference
        % store SV-level metrics from each null iteration (for downstream
        % maxT correction and empirical p value computation).

        results.simulated.permutationH0.modes.s          = resampling_s;           % (nIterations x nModes)
        results.simulated.permutationH0.modes.wilk       = resampling_wilk;        % (nIterations x nModes)

    case 'bootstrapStability'
        % bootstrap stability: singular vector-level 
        % store bootstrap resamples of singular vectors and data projection to them 
        % to later (downstream) assess their stability under sampling
        % variaiblity. Under .singularVectors to reflect 
        results.simulated.bootstrapStability.loadings.U_boot   = U_boot(:,:,:);   % (pX x nModes x nIterations)
        results.simulated.bootstrapStability.loadings.V_boot   = V_boot(:,:,:);   % (pY x nModes x nIterations)
        results.simulated.bootstrapStability.loadings.XU_boot  = XU_boot(:,:,:);  % (m x nModes x nIterations)
        results.simulated.bootstrapStability.loadings.YV_boot  = YV_boot(:,:,:);  % (m x nModes x nIterations)
end

% TO DO: the fields could reflect a more linear algebra language. "modes"
% is not good because everything here is mode-wise. what currently is .modes
% could be .sigma or singularValuesMat to indicate it is based on singular value (U Sigma V).
% what currently is .loadings could be relabelled singularVectorsMat to
% clarify it is based on those

%% sanity checks: output

% PLS_SVD structure consistency
if results.PLS_SVD.modes.nModes < 1
    error('nModes must be at least 1');
end
if size(results.PLS_SVD.modes.s, 2) ~= results.PLS_SVD.modes.nModes
    error('singular values .modes.s does not match nModes');
end
if size(results.PLS_SVD.loadings.U_obs, 2) ~= results.PLS_SVD.modes.nModes
    error('X-side loadings .loadings.U_obs does not match nModes');
end
if size(results.PLS_SVD.loadings.V_obs, 2) ~= results.PLS_SVD.modes.nModes
    error('Y singular vectors .loadings.V_obs does not match nModes');
end
if size(results.PLS_SVD.loadings.V_obs, 1) ~= pY
    error('Y singular vectors .loadings.V_obs does not have pY rows');
end

end
