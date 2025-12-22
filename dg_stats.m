function results = dg_stats(dg_cfg, Y, X)

% This is the front-end wrapper function that reads the input, calls the
% "back-end" functions to do the jobs (e.g., understand what design the
% user wants, do sanity checks, do the analysis) and presents the output to
% the user. 
%
% [results, clusterMetrics, clustThreshold] = PE_Stats_y1x2z3(dataArray, PE_parameters, DimStruct)
%
% INPUT:
%   Y           -
%
%   X           -
%
%   dg_cfg
%
% OUTPUT:
%
%
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% sanity checks

% make sure the analysis is validated prior to this function
if ~isfield(dg_cfg,'validation')
    analysisUnvalidated = true;
else
    analysisUnvalidated = false;
    if ~isfield(dg_cfg.validation,'analysis')
        analysisUnvalidated = true;
    else
        if ~dg_cfg.validation.analysis
            analysisUnvalidated = true;
        end
    end
end
if analysisUnvalidated
    error('dg_cfg.analysis not validated. use: dg_cfg = dg_validateAnalysis(dg_cfg)')
end



% Y and X are matrices
if ~ismatrix(Y)
    error('Y must be a matrix')
end
if ~ismatrix(X)
    error('X must be a matrix')
end

% Y and X have the same num of rows
if size(Y,1)~=size(X,1)
    error('Y and X must have the same num of rows')
end

% X (and Y) must have the same num of rows as in dg_cfg.analysis.dataStruct
if size(X,1)~=size(dg_cfg.analysis.dataStruct,1)
    error('X and dg_cfg.analysis.dataStruct must have the same num of rows')
end

% ignore_col if provided must be a row vector of same length as columns of X
if isfield(dg_cfg.analysis, 'ignore_col')
    if ~isrow(dg_cfg.analysis.ignore_col)
        error('dg_cfg.analysis.ignore_col must be a row vector')
    else
        if size(dg_cfg.analysis.ignore_col,2)~=size(X,2)
            error('dg_cfg.analysis.ignore_col must have same length as columns of X')
        end
    end
end

% ignore_row if provided must be a column vector of same length as rows of X
if isfield(dg_cfg.analysis, 'ignore_row')
    if ~iscolumn(dg_cfg.analysis.ignore_row)
        error('dg_cfg.analysis.ignore_row must be a column vector')
    else
        if size(dg_cfg.analysis.ignore_row,1)~=size(X,1)
            error('dg_cfg.analysis.ignore_row must have same length as rows of X')
        end
    end
end

% univariate analysis only does one comparison at the time
if ismember(dg_cfg.analysis.type, ["empiricalL1_FDR" "theoreticalL1_clusterMaxT"])
    if size(Y,2)>1
        error(['more than one column in matrix Y (i.e., more than one comparison) not supported for ' dg_cfg.analysis.type ])
    end
end

%% apply ignore_row, if provided

if isfield(dg_cfg.analysis,'ignore_row')

    % take a subset of matrices Y and X and of the accompanying dataStruct
    Y = Y(~dg_cfg.analysis.ignore_row,:);
    X = X(~dg_cfg.analysis.ignore_row,:);
    dg_cfg.analysis.dataStruct = dg_cfg.analysis.dataStruct(~dg_cfg.analysis.ignore_row,:);

end

%% keep a copy of the matrices Y and X

X_orig   = X;
Y_orig   = Y;

%% shortcuts

% num of iterations
nIterations = dg_cfg.analysis.nIterations;

% size of matrices Y and X
[m, pY] = size(Y);
[~, pX] = size(X);

% num of dimensions and their numerosity 
% in d# order (fieldnames order)
dimKeys  = fieldnames(dg_cfg.dimensions);
dimTypes = cellfun(@(k) dg_cfg.dimensions.(k).type, dimKeys, 'UniformOutput', false);
dimSizes = cellfun(@(k) length(dg_cfg.dimensions.(k).vec), dimKeys);
nDims    = numel(dimKeys);
Nall     = prod(dimSizes);

contIdx = find(strcmp(dimTypes, 'continuous'));
sphIdx  = find(strcmp(dimTypes, 'spherical'));
catIdx  = find(strcmp(dimTypes, 'categorical'));

if numel(sphIdx) > 1
    error('Only one spherical dimension is supported');
end

% cluster analyses do not support categorical dimensions
if strcmp(dg_cfg.analysis.type, "theoreticalL1_clusterMaxT") && ~isempty(catIdx)
    error('categorical dimensions are not supported for cluster-based analyses');
end

%% parse design
% understand what analysis the user wants to do

keyboard; % UNTIL HERE OK

designCode = dg_parseDesign(dg_cfg,Y);
dg_cfg.designCode = designCode;
results.designCode = designCode;



%% perform the analysis

% get resampling indices
rowIdx = dg_resample(dg_cfg);

switch [dg_cfg.analysis ' & ' num2str(dg_cfg.designCode)]
    case 'empiricalL1_FDR & 1  0  0' % -- correlation --

        % choose analysis subtype
        list = ["Pearson" "Spearman" "Kendall" "cylindrical"];
        [idx,tf] = listdlg('ListString',list,'SelectionMode','single','ListSize',[160 100],'PromptString','choose correlation type');
        if tf==0; idx=1; warning('You did not choose correlation type. I chose for you: Pearson'); end
        corrType = list(idx);
        if strcmp(corrType,'cylindrical')
            error('not coded yet') % it will require its own function to keep things tidy
        end

        for itIdx = 1:nIterations

            % sort rows as appropriate
            [Y,X] = dg_sortRows(dg_cfg,Y_orig,X_orig,rowIdx,itIdx);

            % perform test
            [statVal, pVal] = corr(Y, X, 'type', corrType); % test
  
            if itIdx==1
                statVal_obs           = statVal;
                pVal_obs              = pVal;
            end

            if itIdx==1
                statVal_resamp = zeros(nIterations,pX);
            end
            statVal_resamp(itIdx,:) = statVal;

            dg_counter(itIdx,nIterations)  % iteration counter
        end

        % collate results
        results.statVal_obs        = statVal_obs;
        results.pVal_obs           = pVal_obs;
        results.resampling.statVal_resamp = statVal_resamp;

    case 'empiricalL1_FDR & 0  1  0' % -- independent sample t-test --
        for itIdx = 1:nIterations
            % sort rows as appropriate
            [Y,X] = dg_sortRows(dg_cfg,Y_orig,X_orig,rowIdx,itIdx);

            % perform test
            varType = 'unequal';  % equal | unequal (for info see doc ttest2)
            [~,p,~,stats] = ttest2(X(Y==max(Y),:),X(Y==min(Y),:),'Vartype',varType);
            statVal = stats.tstat;
            pVal = p;

            if itIdx==1
                statVal_obs           = statVal;
                pVal_obs              = pVal;
            end

            if itIdx==1
                statVal_resamp = zeros(nIterations,pX);
            end
            statVal_resamp(itIdx,:) = statVal;

            dg_counter(itIdx,nIterations)  % iteration counter
        end

        % collate results
        results.statVal_obs        = statVal_obs;
        results.pVal_obs           = pVal_obs;
        results.resampling.statVal_resamp = statVal_resamp;

    case 'empiricalL1_FDR & 0  0  1' % -- paired sample t-test --
        for itIdx = 1:nIterations

            % sort rows as appropriate
            [Y,X] = dg_sortRows(dg_cfg,Y_orig,X_orig,rowIdx,itIdx);

            % rows of conditions belonging to first and second halves
            cond_firstHalf = 1:length(Y)/2;
            cond_secondHalf = 1+length(Y)/2:length(Y);
            % rows of conditions with the largest and lowest constrast code
            % ie, find which hald corresponds with the largest constrast code
            contCodes = [unique(Y(cond_firstHalf)) unique(Y(cond_secondHalf))];
            [~, sortIdx] = sort(contCodes,'descend');
            if diff(sortIdx)<0
                cond_maxHalf = cond_secondHalf;
                cond_minHalf = cond_firstHalf;
            else
                cond_maxHalf = cond_firstHalf;
                cond_minHalf = cond_secondHalf;
            end

            % perform test
            [~,p,~,stats] = ttest(X(cond_maxHalf,:),X(cond_minHalf,:));
            statVal = stats.tstat;
            pVal = p;

            if itIdx==1
                statVal_obs           = statVal;
                pVal_obs              = pVal;
            end

            if itIdx==1
                statVal_resamp = zeros(nIterations,pX);
            end
            statVal_resamp(itIdx,:) = statVal;

            dg_counter(itIdx,nIterations)  % iteration counter
        end

        % collate results
        results.statVal_obs        = statVal_obs;
        results.pVal_obs           = pVal_obs;
        results.resampling.statVal_resamp = statVal_resamp;

    case 'theoreticalL1_clusterMaxT & 1  0  0' % -- correlation --

        % initialize lvl2 (cluster) metrics
        resampling  = repmat(struct('id', [], 'size', [], 'mass', []), 1, nIterations);

        % choose analysis subtype
        list = ["Pearson" "Spearman" "Kendall" "cylindrical"];
        [idx,tf] = listdlg('ListString',list,'SelectionMode','single','ListSize',[160 100],'PromptString','choose correlation type');
        if tf==0; idx=1; warning('You did not choose correlation type. I chose for you: Pearson'); end
        corrType = list(idx);
        if strcmp(corrType,'cylindrical')
            error('not coded yet') % it will require its own function to keep things tidy
        end

        for itIdx = 1:nIterations

            % sort rows as appropriate
            [Y,X] = dg_sortRows(dg_cfg,Y_orig,X_orig,rowIdx,itIdx);

            % perform test
            [statVal, pVal] = corr(Y, X, 'type', corrType); % test

            % form clusters
            [clusterMembership, clustIDList, metrics] = dg_clusterForming(dg_cfg, statVal ,pVal);

            if itIdx==1
                statVal_obs           = statVal;
                pVal_obs              = pVal;
                clusterMembership_obs = clusterMembership ;
                clustIDList_obs       = clustIDList ;
                metrics_obs           = metrics;
            end

            switch dg_cfg.objective
                case 'permutationH0testing'
                    resampling(1,itIdx) = metrics;

                case 'bootstrapStability'
                    if itIdx==1
                        statVal_boot = zeros(nIterations,Nall);
                    end
                    statVal_boot(itIdx,:) = statVal;
            end

            dg_counter(itIdx,nIterations)  % iteration counter
        end

        % collate results
        results.statVal_obs                    = statVal_obs;
        results.pVal_obs                       = pVal_obs;
        results.clusters.clusterMembership_obs = clusterMembership_obs;
        results.clusters.clustIDList_obs       = clustIDList_obs;
        results.clusters.metrics_obs           = metrics_obs;
        
        switch dg_cfg.objective
            case 'permutationH0testing'
                results.resampling.metrics           = resampling;
            case 'bootstrapStability'
                results.resampling.statVal_boot      = statVal_boot(:,:);
        end

    case 'theoreticalL1_clusterMaxT & 0  1  0' % -- independent sample t-test --

        % initialize lvl2 (cluster) metrics
        resampling  = repmat(struct('id', [], 'size', [], 'mass', []), 1, nIterations);

        for itIdx = 1:nIterations

            % sort rows as appropriate
            [Y,X] = dg_sortRows(dg_cfg,Y_orig,X_orig,rowIdx,itIdx);

            % perform test
            varType = 'unequal';  % equal | unequal (for info see doc ttest2)
            [~,p,~,stats] = ttest2(X(Y==max(Y),:),X(Y==min(Y),:),'Vartype',varType);
            statVal = stats.tstat;
            pVal = p;

            % form clusters
            [clusterMembership, clustIDList, metrics] = dg_clusterForming(dg_cfg, statVal ,pVal);

            if itIdx==1
                statVal_obs           = statVal;
                pVal_obs              = pVal;
                clusterMembership_obs = clusterMembership ;
                clustIDList_obs       = clustIDList ;
                metrics_obs           = metrics;
            end

            switch dg_cfg.objective
                case 'permutationH0testing'
                    resampling(1,itIdx) = metrics;

                case 'bootstrapStability'
                    if itIdx==1
                        statVal_boot = zeros(nIterations,Nall);
                    end
                    statVal_boot(itIdx,:) = statVal;
            end

            dg_counter(itIdx,nIterations)  % iteration counter
        end
                
        % collate results
        results.statVal_obs                    = statVal_obs;
        results.pVal_obs                       = pVal_obs;
        results.clusters.clusterMembership_obs = clusterMembership_obs;
        results.clusters.clustIDList_obs       = clustIDList_obs;
        results.clusters.metrics_obs           = metrics_obs;

        switch dg_cfg.objective
            case 'permutationH0testing'
                results.resampling.metrics           = resampling;
            case 'bootstrapStability'
                results.resampling.statVal_boot      = statVal_boot(:,:);
        end

    case 'theoreticalL1_clusterMaxT & 0  0  1' % -- paired sample t-test --

        % initialize lvl2 (cluster) metrics
        resampling  = repmat(struct('id', [], 'size', [], 'mass', []), 1, nIterations);

        for itIdx = 1:nIterations

            % sort rows as appropriate
            [Y,X] = dg_sortRows(dg_cfg,Y_orig,X_orig,rowIdx,itIdx);

             % rows of conditions belonging to first and second halves
            cond_firstHalf = 1:length(Y)/2;
            cond_secondHalf = 1+length(Y)/2:length(Y);
            % rows of conditions with the largest and lowest constrast code
            % ie, find which hald corresponds with the largest constrast code
            contCodes = [unique(Y(cond_firstHalf)) unique(Y(cond_secondHalf))];
            [~, sortIdx] = sort(contCodes,'descend');
            if diff(sortIdx)<0
                cond_maxHalf = cond_secondHalf;
                cond_minHalf = cond_firstHalf;
            else
                cond_maxHalf = cond_firstHalf;
                cond_minHalf = cond_secondHalf;
            end

            % perform test
            [~,p,~,stats] = ttest(X(cond_maxHalf,:),X(cond_minHalf,:));
            statVal = stats.tstat;
            pVal = p;

            % form clusters
            [clusterMembership, clustIDList, metrics] = dg_clusterForming(dg_cfg, statVal ,pVal);

            if itIdx==1
                statVal_obs           = statVal;
                pVal_obs              = pVal;
                clusterMembership_obs = clusterMembership ;
                clustIDList_obs       = clustIDList ;
                metrics_obs           = metrics;
            end

            switch dg_cfg.objective
                case 'permutationH0testing'
                    resampling(1,itIdx) = metrics;

                case 'bootstrapStability'
                    if itIdx==1
                        statVal_boot = zeros(nIterations,Nall);
                    end
                    statVal_boot(itIdx,:) = statVal;
            end

            dg_counter(itIdx,nIterations)  % iteration counter
        end
                
        % collate results
        results.statVal_obs                    = statVal_obs;
        results.pVal_obs                       = pVal_obs;
        results.clusters.clusterMembership_obs = clusterMembership_obs;
        results.clusters.clustIDList_obs       = clustIDList_obs;
        results.clusters.metrics_obs           = metrics_obs;

        switch dg_cfg.objective
            case 'permutationH0testing'
                results.resampling.metrics           = resampling;
            case 'bootstrapStability'
                results.resampling.statVal_boot      = statVal_boot(:,:);
        end

    case {'PLS_SVD & 1  0  0'   'PLS_SVD & 0  1  0'   'PLS_SVD & 1  1  0'}
        nModes = min(rank(Y),rank(X));

        for itIdx = 1:nIterations

            % sort rows as appropriate
            [Y,X] = dg_sortRows(dg_cfg,Y_orig,X_orig,rowIdx,itIdx);

            % mean center columns of Y and X
            Yz = normalize(Y,'center');
            Xz = normalize(X,'center');

            % zscore columns of Y (optional)
            if dg_cfg.pls_svdParams.zscoringVec(1)
                Yz = zscore(Y);
            end
            % zscore columns of X (optional)
            if dg_cfg.pls_svdParams.zscoringVec(2)
                Xz = zscore(X);
            end

            % apply ignore mask (features in X)
            Xz(:,logical(dg_cfg.R_ignore)) = 0;

            % scale by frobenius norm (optional)
            % DELETE THIS SECTION UNLESS JUSTIFIED
            %                     if logical(dg_cfg.pls_svdParams.froFlag)
            %                         Lz = Lz / norm(Lz,'fro');
            %                         Rz = Rz / norm(Rz,'fro');
            %                     end

            % cross-product
            %                     C = (Yz'*Xz_pca)/(m-1);  % for PCA
            C = (Yz'*Xz)/(m-1);

            % SVD
            [U,S,V] = svd(C,"econ");

            % sanity check: nModes as expected
            if nModes~=size(S,1)
                error('unexpected number of SVD modes')
            end

            % resolve sign uncertainty
            % match SVD outcome with observed data
            % only for last iteration
            if itIdx==1
                [U, V] = dg_signConvention(Yz,Xz,U,V);
            end

            % metrics: level-1 (feature based)
            YU = Yz*U;
            XV = Xz*V;

            % metrics: level-2 (SV based)
            s = diag(S); % singular value of each mode
            p = s.^2 / sum(s.^2); % proportion of covariance explained by each mode
            r = diag(corr(YU, XV)); % Pearson correlation between left and right latent variables (i.e., scores = data * singular vectors) per mode
            % store the above as col vectors
            s = s(:);
            p = p(:);
            r = r(:);


            % store original-data values (observed)
            if itIdx==1
                U_obs = U;
                V_obs = V;
                YU_obs = YU;
                XV_obs = XV;
                s_obs = s';
                p_obs = p';
                r_obs = r';
            end

            % store output of each iteration as appropriate
            switch dg_cfg.objective
                case 'permutationH0testing'
                    % initialize lvl2 (mode wise) metrics
                    if itIdx==1
                        resampling  = repmat(struct('s', [], 'inertia', [], 'wilk', [], 'sequential', []), 1, nIterations);
                    end
                    % fill each field per iteration
                    resampling(:,itIdx).s = s;
                    resampling(:,itIdx).inertia = sum(s);
                    resampling(:,itIdx).wilk = flipud(cumsum(flipud(s(:).^2)));
                    resampling(:,itIdx).sequential = cumsum(s(:).^2);

                case 'bootstrapStability'
                    % initialize lvl1 (loading wise) metrics
                    % store singular and latent vectors (all iterations)
                    if itIdx==1
                        U_boot = zeros(size(Y,2),nModes,nIterations);
                        V_boot = zeros(size(X,2),nModes,nIterations);
                        YU_boot = zeros(size(Y,1),nModes,nIterations);
                        XV_boot = zeros(size(X,1),nModes,nIterations);
                    end
                    U_boot(:,:,itIdx) = U;
                    V_boot(:,:,itIdx) = V;
                    YU_boot(:,:,itIdx) = YU;
                    XV_boot(:,:,itIdx) = XV;
            end

            dg_counter(itIdx,nIterations)  % iteration counter
        end

        % collate results
        % observed data
        results.PLS_SVD.nModes = nModes;
        results.PLS_SVD.s_obs = s_obs;
        results.PLS_SVD.p_obs = p_obs;
        results.PLS_SVD.r_obs = r_obs;
        results.PLS_SVD.U_obs = U_obs;
        results.PLS_SVD.V_obs = V_obs;
        results.PLS_SVD.YU_obs = YU_obs;
        results.PLS_SVD.XV_obs = XV_obs;

        switch dg_cfg.objective
            case 'permutationH0testing'
                results.resampling     = resampling;
            case 'bootstrapStability'
                % all iterations at lvl1
                results.resampling.U_boot   = U_boot(:,:,:);
                results.resampling.V_boot   = V_boot(:,:,:);
                results.resampling.YU_boot  = YU_boot(:,:,:);
                results.resampling.XV_boot  = XV_boot(:,:,:);
        end

    otherwise % -- any other design --
        error(['not yet coded: ' dg_cfg.objective ' ' num2str(dg_cfg.designCode)])
end

%% compute inferential metrics
% compute p values for permutation testing
% compute BR and CI for bootstrap stability

results = dg_inference(dg_cfg,results);
       
%% cluster descriptive metrics

results = dg_describeClusters(dg_cfg,results);

end
