function P = F_DiscretePoissonMultinomialMultiAny_V2(p, mode,base,opts)
% Poisson-Multinomial distribution via dynamic programming
% Supports log-space (L,S) with auto-switch from P-space.
%
% INPUT:
%   p    : (n x m x latent...) array of outcome probabilities
%   mode : 'P' (default) probability-space
%          'LS' log-space with support mask
%
% OUTPUT:
%   P    : distribution over counts
if nargin<4
    opts.PrecisionCheckUpdate = 0;
end
if nargin<3
    base = 'e';
end
if nargin<2
    mode = 'P';
end


numberOfDimensions = ndims(p);
latentDimensions= cell(1,numberOfDimensions-2);
[num_sites, num_outcomes, latentDimensions{:}] = size(p);
latentDims = [latentDimensions{:}];

% Collapse latent dims to check support
if numberOfDimensions > 2
    pmax = max(reshape(p,[num_sites,num_outcomes,numel(p)/(num_sites*num_outcomes)]),[],3);
else
    pmax = p;
end
K = sum(pmax>0,1); % max counts

sz = [K+1, latentDims];

% --- Initialization
if strcmp(mode,'P')
    P = ndSparse.build(sz,0);
    idx0 = num2cell(ones(1,num_outcomes));%incorrect
    P(idx0{:},:) = 1;
else
    % log-space
    L = ndSparse.build(sz,0); % 1 at zero state
    S = ndSparse.build(sz,0);
    idx0 = num2cell(ones(1,num_outcomes));%incorrect
    L(idx0{:},:) = 0;  % log(1)=0
    S(idx0{:},:) = 1;  % true one
    mode = 'LS';
end

switchedToLS = false;
updateBaseNow = 0;
% --- Loop over trials
for i = 1:num_sites
    % Extract probabilities
    idxSite = [num2cell(i), repmat({':'},1,length(latentDims)+1)];
    pSite = reshape(p(idxSite{:}),[num_outcomes latentDims]);
    if strcmp(mode,'P') && ~switchedToLS
        % ----- Probability-space update -----
        [Pnew, switchNow] = updateTrialP(P,pSite);
        if switchNow
            % convert into (L,S) representation
            [L,S] = convertPtoLS(P,base);
            [L,S,updateBaseNow] = updateTrialLS(L,S,pSite,opts,base);
            %[L,S] = normalizeLS(L,S,K,base);
            mode = 'LS'; switchedToLS = true;
        else
            [Z,~] = rs_sum(Pnew,1:length(K));
            Pnew = Pnew ./ Z;
            P = Pnew;
        end
    else
        % ----- Log-space update -----
        [L,S,updateBaseNow] = updateTrialLS(L,S,pSite,opts,base);
        %[L,S] = normalizeLS(L,S,K,base);
    end
    if updateBaseNow
        error('base');
    end
end

% Back convert if log-space
if strcmp(mode,'LS')
    [L,S] = normalizeLS(L,S,K,base);
    P = backConvertLS(L,S,K,base);
end
end


function F = Fmask(L,S)
% Mask: 1 if not a true zero, 0 if (L=0,S=0)
F = ~(L==0 & S==0);
end

function P = backConvertLS(L,S,outcomeDims,base)
% Convert (L,S) -> probability
% backConvertLS : convert (L,S) representation back to probabilities P.
%
% Inputs:
%   L    : array of log-values (same shape as probability array).
%   S    : array of flags (same shape) — disambiguates L==0.
%   base : (optional) log base, default = exp(1).
%
% Output:
%   P    : probabilities in [0,1].

% Define F(L,S) = 0 if (L==0 & S==0), else 1
isSparse = issparse(L) || isa(L,'ndSparse');
if isSparse
    if nargin < 4 || strcmp(base,'e') % natural log
        P = Fmask(L,S) .* spfun(@exp,L);
    else
        P = Fmask(L,S) .*  spfun(@(x) base.^x, L);
    end
    P((L==0).*(S==1)==1) = 1;
else
    if nargin < 4 || strcmp(base,'e') % natural log
        P = double(simplify(Fmask(L,S))) .* exp(L);
    else
        P = double(simplify(Fmask(L,S))) .* base.^L;
    end
end
  [Z,~] = rs_sum(P,1:length(outcomeDims));
if isSparse
   Pind = num2cell(P, sort(full(outcomeDims)));   
   P = P./ repmat(Z,size(P)-size(Pind)+1);
else
    P = P ./ Z;
end

end

function [Lnew,Snew] = multiplyLS(Li,Si,Lj,Sj)
% Multiply two probabilities in log domain
Lnew = Fmask(Li,Si).*Fmask(Lj,Sj).*(Li+Lj);
Snew = (Si==1 & Sj==1); % only true 1*1 stays true 1
end

function [Lmerge,Smerge] = mergeLS(Lg,Sg,base)
% Merge by log-sum-exp, ignoring true zeros
if nargin < 3 || isempty(base)
    base = 'e'; % default = natural log
end

% Detect sparse vs full
% Handle empty input
if isempty(Lg)
    Lmerge = 0;
    Smerge = 0;
    return;
end

isSparse = issparse(Lg) || isa(Lg,'ndSparse');
if isSparse
    % Build mask = not true zeros
    mask = Fmask(Lg,Sg);
    if ~nnz(mask)
        Lmerge = 0;
        Smerge = 0;
        return;
    end
    % take only valid entries
    Lg = Lg(mask);
    Lmax = max(Lg(:));
    if nargin < 3 || strcmp(base,'e')
        contrib = spfun(@(x) exp(x - Lmax), Lg);
        sumContrib = sum(contrib(:));
        Lmerge = Lmax + spfun(@log,sumContrib);
    else
        contrib = spfun(@(x) base.^(x - Lmax), Lg);
        sumContrib = sum(contrib(:));
        Lmerge = Lmax + spfun(@log,sumContrib) / log(base);
    end
    Smerge = (Lmerge==0);
else
    mask = Fmask(Lg,Sg);
    if ~any(mask)
        Lmerge = 0;
        Smerge = 0;
        if ~(ismcc || isdeployed)
            %#exclude sym
            %#exclude eval
            if (isa(mask,'sym'))
                Lmerge = sym(Lmerge);
                Smerge = sym(Smerge);
            end
        else
            if (isa(mask,'hpf'))
                Lmerge = hpf(Lmerge);
                Smerge = hpf(Smerge);
            end
        end
        return;
    end
    Lvalid = Lg(mask);
    Lmax = max(Lvalid);
    if isequal(base,'e')
        contrib = sum(exp(Lvalid - Lmax));
        Lmerge = Lmax + log(contrib);
    else
        contrib = sum(base.^(Lvalid - Lmax));
        Lmerge = Lmax + log(contrib)/log(base);
    end
    % S=1 only if exactly equal to true 1
    Smerge = (Lmerge==0);
    % Smerge = 1 if exactly one contributor is a true 1
    %Smerge = (numel(Lvalid)==1 && all(Lvalid==0) && any(Svalid==1));
end
% % Decide Smerge:
% % Smerge should be true only if merged probability equals exactly 1.
% % That happens only if:
% %   - there is at least one contributor with (L==0 & S==1)
% %   - and all other contributors are true zeros (i.e. excluded)
% % or all contributors
% % Check for the unique exact-one contributor case:
% isOneContrib = (L == 0) & (S == 1);
% % compute number of valid contributions (after K filter)
% numValid = numel(Lvalid);
% if any(isOneContrib) && (sum(isOneContrib) == 1) && (numValid == 1)
%     % exactly one contributor with value 1 and everyone else is excluded,
%     % merged sum is exactly 1
%     Smerge = true;
% else
%     % otherwise Smerge true only if Lmerge numerically equals zero:
%     % (this happens rarely due to numeric rounding; use tolerance)
%         Smerge = abs(Lmerge) <= tol;
% end
end
function [S,D] = rs_sum(S,D)
if ~isempty(D)
    S = sum(S,D(1));
    D(1) = [];
    [S,D] = rs_sum(S,D);
end
end

function [S,D] = rs_max(S,D)
if ~isempty(D)
    S = max(S,[],D(1));
    D(1) = [];
    [S,D] = rs_max(S,D);
end
end

function [Lnorm,Snorm] = normalizeLS(L,S,outcomeDims,base)
if nargin < 4
    base = 'e';
end

% Normalize probabilities across outcome counts per latent slice
% Pnorm2 = backConvertLS(L,S);
% [Lnorm2,Snorm2] = convertPtoLS(Pnorm2,base);
NumberOfLatentDims = ndims(L)-length(outcomeDims);

F = Fmask(L,S);
Lind = num2cell(L, 1:length(outcomeDims));      
Find = num2cell(F, 1:length(outcomeDims));   
[Zind,~] = cellfun(@(x,y) Lmerge_Func(x,y,base),Lind,Find,'Un',0);
%Znorm = cell2mat(Zind);
Zind = cellfun(@full,Zind,'Un',0);
Zist = squeeze(Zind);
Zflat = [Zist{:}];
Zcoord = cell(1,NumberOfLatentDims);
[Zcoord{:}] = ind2sub(size(Zist),1:numel(Zist));
Zcoord_transpose = cellfun(@transpose,Zcoord,'Un',0);
Zcoord_table = [Zcoord_transpose{:}];
Znorm = ndSparse.build([ones(numel(Zist),length(outcomeDims)) Zcoord_table],Zflat,[size(Zind)]);
Lnorm = L - repmat(Znorm,size(L)-size(Lind)+1).*F;
Snorm = S;
%Pnorm = backConvertLS(Lnorm,Snorm);  
% Lm = rs_max(L,1:length(outcomeDims));
% 
% Lmax = max(L(L~=0));
% if  strcmp(base,'e')
%     F_Latent = rs_sum(F,1:length(outcomeDims));
%     Lmax = max(L(F));
%     contrib = spfun(@(x) exp(x - Lmax), L);
%     contrib((L==0).*(S==1)==1) = exp(-full(Lmax));
%     sumContrib = rs_sum(contrib,1:length(outcomeDims));
%     Z_latent = Lmax + spfun(@log,sumContrib);
%     Z_latent(F_Latent==0) = 0;
%     Lmax2 = 0.5*max(L(L~=0));
%     contrib2 = spfun(@(x) exp(x - Lmax2), L);
%     contrib2((L==0).*(S==1)==1) = exp(-full(Lmax2));
%     sumContrib2 = rs_sum(contrib2,1:length(outcomeDims));
%     Z_latent2 = Lmax2 + spfun(@log,sumContrib2);
% 
% else
%     contrib = spfun(@(x) base.^(x - Lmax), L);
%     contrib(L==0) = 1;
%     sumContrib = sum(F.*contrib);
%     Z_latent = Lmax + spfun(@log,sumContrib) / log(base);
% end
% 
% if strcmp(base,'e')
%     Lmerged_G= spfun(@(x) log(x)+max(Lg(mask)),sum(spfun(@(x) exp(x - max(Lg(mask))), Lg(mask))));
% else
%     Lmerged_G = spfun(@(x) log(x)/log(base) +max(Lg(mask)),sum(spfun(@(x) base.^(x - max(Lg(mask))), Lg(mask))));
% end
% Re-encode into (L,S)
% [Lnorm2sym,Snorm2sym] = convertPtoLS(Pnorm2sym,base);
% Pnormsym = exp(Lnorm2sym);
% 
% [Z2sym,~] = rs_sum(Pnormsym,1:length(outcomeDims));
% [Z3sym,~] = rs_sum(Pnorm2sym,1:length(outcomeDims));
% 
% [Zind_sym,~] = cellfun(@(x,y) Lmerge_Sym_Func(sym(sparse(x)),sym(sparse(y)),base),Lind,Find,'Un',0);
% Znorm_sym = cell2mat(Zind_sym);
% Lnorm_sym = reshape(sym(sparse(L)),size(L)) - repmat(Znorm_sym,size(L)-size(Lind)+1).*reshape(sym(sparse(F)),size(F));
% Snorm_sym = reshape(sym(sparse(Snorm)),size(L));
% Pnorm_sym = backConvertLS(Lnorm_sym,Snorm_sym);
% [Z4sym,~] = rs_sum(Pnorm_sym,1:length(outcomeDims));
end
% function [Lnorm, Snorm] = normalizeLS(L, S)
% % normalizeLS : normalize (L,S) representation so that P sums to 1
% %
% % Inputs:
% %   L : log-values (array)
% %   S : flag array, distinguishes true zero vs true one
% %
% % Outputs:
% %   Lnorm, Snorm : normalized version
%
%     % Case A: If there are true ones
%     if any(S(:) == 1)
%         % Only those entries survive
%         Snorm = (S == 1);
%         Lnorm = zeros(size(L));   % log(1) = 0
%         return;
%     end
%
%     % Case B: Otherwise, normalize the "soft" distribution in log-space
%     % Mask out true zeros
%     mask = ~( (L==0) & (S==0) );
%     active = find(mask);
%
%     if isempty(active)
%         % Degenerate case: all zeros -> return same
%         Lnorm = L;
%         Snorm = S;
%         return;
%     end
%
%     % Use log-sum-exp trick over active entries
%     Lmax = max(L(active));
%     Z = Lmax + log(sum(exp(L(active) - Lmax)));
%
%     % Normalize
%     Lnorm = L - Z;
%     Snorm = S; % unchanged
% end

function [L,S] = convertPtoLS(P,base)
% Convert probability distribution into (L,S)
if nargin <2
    base = 'e';
end
isSparse = issparse(P) || isa(P,'ndSparse');
if isSparse
    S = ndSparse.build(size(P),0);
    if nargin < 2 || strcmp(base,'e')
        L = spfun(@log,P);
    else
        L = spfun(@log,P)/log(base);
    end
    maskOnes = (P==1);
    L(maskOnes) = 0;
    S(maskOnes) = 1;
else
    S = ndSparse.build(size(P),0);
    if nargin < 2 || strcmp(base,'e')
        L = log(P);
    else
        L = log(P)/log(base);
    end
    maskOnes = logical(P==1);
    L(find(maskOnes)) = 0;
    S(find(maskOnes)) = 1;
end
end

function [Pnew, switchNow] = updateTrialP(P,pSite)
% Update in probability-space, check for underflow
switchNow = false;
num_outcomes = size(pSite,1);
latentDims = size(pSite);
latentDims = latentDims(2:end);
sz = size(P);% full size of the distribution
Pnew = ndSparse.build(sz,0);
% --- Compute "none" probability ---
piNone = ndSparse(1 - sum(pSite,1));% [latentDims...]

% --- Find all active states in current P ---
subsCell = cell(1,ndims(P));% turn into cell array of coordinates if needed
linear_indx = find(P);
if isempty(linear_indx)
    % only zero state exists
    linear_indx = sub2ind(sz, ones(1,num_outcomes+numel(latentDims)));
end
[subsCell{:}] = ind2sub(size(P), linear_indx);  % [numActive, m+numLatent]

%count outcomes
counts = cell2mat(subsCell(:,1:num_outcomes)) - 1;   % zero-based counts
current_Probabilities = P(linear_indx);
current_Probabilities = current_Probabilities(:);              % [numActive,1]

probsNo_BindingEvents = reshape(piNone,[],1);
if ~isempty(latentDims)
    latentIdx = cell2mat(subsCell(:,num_outcomes+1:end));   % latent subs
    latentSubs = num2cell(latentIdx,1);
    linLatent = sub2ind(latentDims,latentSubs{:});
    probsNo_BindingEvents_Latent = probsNo_BindingEvents(linLatent);
else
    latentIdx = zeros(numel(linear_indx),0);
    latentSubs = {};
    probsNo_BindingEvents_Latent = probsNo_BindingEvents;
end
wNone = current_Probabilities.*probsNo_BindingEvents_Latent;
% none outcome
argsNone = [num2cell(counts + 1,1),latentSubs];
linIdxNone = sub2ind(sz,argsNone{:});
% Detect underflow: nonzero*nonzero→0
if ~(ismcc || isdeployed)
    %#exclude sym
    diff_None_Error = sym(sparse(wNone)) - sym(sparse(current_Probabilities)).*sym(sparse(probsNo_BindingEvents_Latent));
    diff_Reference = sym(zeros(size(wNone)));
else
    diff_None_Error = hpf(full(wNone)) - hpf(full(current_Probabilities)).*hpf(full(probsNo_BindingEvents_Latent));
    diff_Reference = hpf(zeros(size(wNone)));
end

if sum((current_Probabilities>0).*(probsNo_BindingEvents_Latent>0).*(wNone==0))>0 || ~isequal(diff_None_Error,diff_Reference)
    %any(current_Probabilities>0 & probsNo_BindingEvents_Latent>0 & wNone==0)
    switchNow = true;
    return;
end

tmpNone = accumarray(linIdxNone, full(wNone), [prod(sz),1]);
Pnew = Pnew + ndSparse(reshape(tmpNone,sz));

for j = 1:num_outcomes
    mask = counts(:,j) < sz(j)-1;
    if ~any(mask)
        continue;
    end
    countsInc = counts(mask,:);
    countsInc(:,j) = countsInc(:,j) + 1; % Increment outcome j
    valuesToIncrement_Probabilities = current_Probabilities(mask);

    % Probability for this outcome
    idxOutcome = [num2cell(j), repmat({':'},1,length(latentDims))];
    probs_BindingEvent_J = reshape(pSite(idxOutcome{:}),[],1);

    % Repeat along latent dims
    if ~isempty(latentDims)
        latentIdxInc = latentIdx(mask,:);
        latentSubsInc = num2cell(latentIdxInc,1);
        latentLinInc = sub2ind(latentDims,latentSubsInc{:});
        probs_BindingEvent_J_Latent  = probs_BindingEvent_J(latentLinInc);
    else
        latentSubsInc = {};
        probs_BindingEvent_J_Latent  = probs_BindingEvent_J;
    end
    wJ = valuesToIncrement_Probabilities .*probs_BindingEvent_J_Latent;

    % Detect underflow: nonzero*nonzero→0
    if ~(ismcc || isdeployed)
        %#exclude sym
        diff_BindingEvent_J_Error = sym(sparse(wJ)) - sym(sparse(valuesToIncrement_Probabilities)).*sym(sparse(probs_BindingEvent_J_Latent));
        diff_Reference = sym(zeros(size(wJ)));
    else
        diff_BindingEvent_J_Error = hpf(full(wJ)) - hpf(full(valuesToIncrement_Probabilities)).*hpf(full(probs_BindingEvent_J_Latent));
        diff_Reference = hpf(zeros(size(wJ)));
    end

    if sum((valuesToIncrement_Probabilities>0).*(probs_BindingEvent_J_Latent>0).*(wJ==0))>0  || ~isequal(diff_BindingEvent_J_Error,diff_Reference)
        switchNow = true;
        return;
    end

    % Linear indices in Pnew
    countsIncPlus = countsInc + 1;
    subsColsInc = num2cell(countsIncPlus,1);
    argsInc = [subsColsInc, latentSubsInc];
    linIdxInc = sub2ind(sz,argsInc{:});

    % Accumulate probabilities
    tmpJ = accumarray(linIdxInc,full(wJ),[prod(sz),1]);
    Pnew = Pnew + ndSparse(reshape(tmpJ,sz));
end
end

function [Lmerged_G,Smerged_G] = Lmerge_Func(Lg,mask,base)
if (nnz(mask)==0)
    Lmerged_G = 0;
    Smerged_G = 0;
    return;
end
if strcmp(base,'e')
    Lmerged_G= spfun(@(x) log(x)+max(Lg(mask)),sum(spfun(@(x) exp(x - max(Lg(mask))), Lg(mask))));
else
    Lmerged_G = spfun(@(x) log(x)/log(base) +max(Lg(mask)),sum(spfun(@(x) base.^(x - max(Lg(mask))), Lg(mask))));
end
Smerged_G = (Lmerged_G==0);
end
function [Lmerged_G_Sym,Smerged_G_Sym] = Lmerge_Sym_Func(Lg,mask,base)
if (sum(double(simplify(mask)))==0)
    if (isa(Lg,'hpf'))
    Lmerged_G_Sym = hpf(0);
    Smerged_G_Sym = 0;
    end
    if (isa(Lg,'sym'))
    Lmerged_G_Sym = sym(0);
    Smerged_G_Sym = 0;
    end
    return;
end
mask = find(double(simplify(mask)));
if strcmp(base,'e')
    Lmerged_G_Sym = max(Lg(mask)) + log(full(sum(exp(Lg(mask) - max(Lg(mask)))))); 
else
    Lmerged_G_Sym = max(Lg(mask)) + log(full(sum(base.^(Lg(mask)  - max(Lg(mask))))))/log(base);
end
if (isa(Lg,'hpf'))
    Lmerged_G_Sym = Lmerged_G_Sym-log(hpf(1));
end
if (isa(Lg,'sym'))
    Smerged_G_Sym = double(Lmerged_G_Sym==0);
end
if (isa(Lg,'hpf'))
    Smerged_G_Sym = (Lmerged_G_Sym==0);
end
end


function [Lnew, Snew,updateBaseNow] = updateTrialLS(L, S, pSite, opts, base)
% Update a single trial in log-space using L & S ndSparse arrays
% L : ndSparse log-probabilities
% S : ndSparse flags for true 1s
% pSite : [m x latentDims...] probability for current trial
% base : optional log base (default: e)
updateBaseNow = false;
% detect inaccuracy in log status and increment base to larger value
if nargin < 5
    base = 'e';
end
if nargin < 4
    opts.PrecisionCheckUpdate = 0;
end
num_outcomes = size(pSite,1);
latentDims = size(pSite);
latentDims = latentDims(2:end);
% Prepare output
sz = size(L);
Lnew = ndSparse.build(sz, 0);
Snew = ndSparse.build(sz, 0);
% Compute "none" probability
piNone = ndSparse(1 - sum(pSite,1));% [latentDims...] treat as sparse
% Get all active indices in current L and S
linear_idx = find(Fmask(L,S)); % only consider non-zero / non-true-zero states
if isempty(linear_idx)
    % Only zero states exist, start from zero counts
    linear_idx = sub2ind(sz, ones(1,num_outcomes+numel(latentDims)));
end
subsCell = cell(1, ndims(L));
[subsCell{:}] = ind2sub(sz, linear_idx);
counts = cell2mat(subsCell(:,1:num_outcomes)) - 1;
current_LogProbabilities_L = L(linear_idx);
current_LogProbabilities_S = S(linear_idx);
current_LogProbabilities_L  = current_LogProbabilities_L(:);              % [numActive,1]
current_LogProbabilities_S = current_LogProbabilities_S(:);              % [numActive,1]
probsNo_BindingEvents = reshape(piNone,[],1); % flatten latent dims
if ~isempty(latentDims)
    latentIdx = cell2mat(subsCell(:,num_outcomes+1:end));   % latent subs
    latentSubs = num2cell(latentIdx,1);
    linLatent = sub2ind(latentDims,latentSubs{:});
    [L_No_BindingEvents_Latent,S_No_BindingEvents_Latent] = convertPtoLS(probsNo_BindingEvents(linLatent),base);
else
    latentIdx = zeros(numel(linear_idx),0);
    latentSubs = {};
    [L_No_BindingEvents_Latent,S_No_BindingEvents_Latent] = convertPtoLS(probsNo_BindingEvents,base);
end

%% --- Update "none" outcome ---
% Compute wNone = multiplyLS(valsL, valsS, log(piNone))
[Lnone, Snone] = multiplyLS(current_LogProbabilities_L, current_LogProbabilities_S, L_No_BindingEvents_Latent, S_No_BindingEvents_Latent);
if opts.PrecisionCheckUpdate
    if ~(ismcc || isdeployed)
        %#exclude sym
        %#exclude eval
        [Lnone_Sym, Snone_Sym] = multiplyLS(sym(sparse(current_LogProbabilities_L)), sym(sparse(current_LogProbabilities_S)),...
            sym(sparse(L_No_BindingEvents_Latent)), sym(sparse(S_No_BindingEvents_Latent)));
        diff_None_Error = (sym(sparse(Lnone))-Lnone_Sym).*(sym(sparse(Snone)) - sym(double(simplify(Snone_Sym))));
        diff_Reference = sym(zeros(size(Lnone)));
    else
        [Lnone_Sym, Snone_Sym] = multiplyLS(hpf(full(current_LogProbabilities_L)), hpf(full(current_LogProbabilities_S)),...
            hpf(full(L_No_BindingEvents_Latent)), hpf(full(S_No_BindingEvents_Latent)));
        diff_None_Error = (hpf(full(Lnone))-Lnone_Sym).*hpf(full(Snone-Snone_Sym));
        diff_Reference = hpf(zeros(size(Lnone)));
    end
    if sum((current_LogProbabilities_S==0).*(S_No_BindingEvents_Latent==0).*(Snone==1))>0 || ~isequal(diff_None_Error,diff_Reference)
        updateBaseNow = true;isAlways(diff_None_Error== diff_Reference)
        return;
    end
end
% Linear indices in Lnew
argsNone = [num2cell(counts + 1,1),latentSubs];
linIdxNone = sub2ind(sz,argsNone{:});
uIdxNone  = unique(linIdxNone);
linIdxNone_All = [uIdxNone ;linIdxNone];
Lnone_All = [Lnew(uIdxNone); Lnone];
Snone_All = [Snew(uIdxNone); Snone];
[None_Groups,None_Group_ID] = findgroups(linIdxNone_All);
mask_None_All = Fmask(Lnone_All,Snone_All);
num_in_mask = splitapply(@sum,full(mask_None_All),None_Groups);
if opts.PrecisionCheckUpdate
    if ~(ismcc || isdeployed)
        %#exclude sym
        %#exclude eval
        mask_None_All_Sym = Fmask(sym(sparse(Lnone_All)),sym(sparse(Snone_All)));
        num_in_mask_From_Sym = splitapply(@sum,double(simplify(mask_None_All_Sym)),None_Groups);
    else
        mask_None_All_Sym = Fmask(hpf(full(Lnone_All)),hpf(full(Snone_All)));
        num_in_mask_From_Sym = splitapply(@sum,mask_None_All_Sym,None_Groups);
    end
    if ~isequal(num_in_mask,num_in_mask_From_Sym)
        updateBaseNow = true;
        return;
    end
end
[Lmerged_G, Smerged_G] = splitapply(@(Lg,Mask) Lmerge_Func(Lg,Mask,base),full(Lnone_All),full(mask_None_All),None_Groups);
if opts.PrecisionCheckUpdate
    if ~(ismcc || isdeployed)
        %#exclude sym
        %#exclude eval
        [Lmerged_G_Sym, Smerged_G_Sym] = splitapply(@(Lg,Mask) Lmerge_Sym_Func(Lg,Mask,base),sym(sparse(Lnone_All)),mask_None_All_Sym,None_Groups);
        diff_None_Error = (sym(sparse(Lmerged_G))-Lmerged_G_Sym).*sym(Smerged_G - Smerged_G_Sym);
        diff_Reference = sym(zeros(size(Lmerged_G)));
    else
        [Lmerged_G_Sym, Smerged_G_Sym] = splitapply(@(Lg,Mask) Lmerge_Sym_Func(Lg,Mask,base),hpf(full(Lnone_All)),mask_None_All_Sym,None_Groups);
        diff_None_Error = (hpf(full(Lmerged_G))-Lmerged_G_Sym).*hpf(full(Smerged_G-Smerged_G_Sym));
        diff_Reference = hpf(zeros(size(Lmerged_G)));
    end
    if ~isequal(diff_None_Error,diff_Reference)
        updateBaseNow = true;
        return;
    end
end
Lnew(None_Group_ID(num_in_mask>0)) = Lmerged_G(num_in_mask>0);
Snew(None_Group_ID(num_in_mask>0)) = Smerged_G(num_in_mask>0);
if sum((num_in_mask==0))>0
    Lnew(None_Group_ID(num_in_mask==0)) = 0;
    Snew(None_Group_ID(num_in_mask==0)) = 0;
end
%% --- Update each outcome ---
for j = 1:num_outcomes
    mask = counts(:,j) < sz(j)-1; % max counts
    if ~any(mask)
        continue;
    end

    countsInc = counts(mask,:);
    countsInc(:,j) = countsInc(:,j)+1; % increment outcome j
    valuesToIncrement_LogProbabilities_L = current_LogProbabilities_L(mask);
    valuesToIncrement_LogProbabilities_S = current_LogProbabilities_S(mask);

    % Probability for this outcome
    idxOutcome = [num2cell(j), repmat({':'},1,length(latentDims))];
    probs_BindingEvent_J = reshape(pSite(idxOutcome{:}),[],1);

    % Repeat along latent dims
    if ~isempty(latentDims)
        latentIdxInc = latentIdx(mask,:);
        latentSubsInc = num2cell(latentIdxInc,1);
        latentLinInc = sub2ind(latentDims,latentSubsInc{:});
        [L_BindingEvent_J_Latent,S_BindingEvent_J_Latent] = convertPtoLS(probs_BindingEvent_J(latentLinInc),base);
    else
        latentSubsInc = {};
        [L_BindingEvent_J_Latent,S_BindingEvent_J_Latent] = convertPtoLS(probs_BindingEvent_J,base);
    end

    [Lprod, Sprod] = multiplyLS(valuesToIncrement_LogProbabilities_L, valuesToIncrement_LogProbabilities_S,...
        L_BindingEvent_J_Latent, S_BindingEvent_J_Latent);
    if opts.PrecisionCheckUpdate
        if ~(ismcc || isdeployed)
            %#exclude sym
            %#exclude eval
            [Lprod_Sym, Sprod_Sym] = multiplyLS(sym(sparse(valuesToIncrement_LogProbabilities_L)), sym(sparse(valuesToIncrement_LogProbabilities_S)),...
                sym(sparse(L_BindingEvent_J_Latent)), sym(sparse(S_BindingEvent_J_Latent)));
            diff_IncJ_Error = (sym(sparse(Lprod))-Lprod_Sym).*sym(Sprod - double(simplify(Sprod_Sym)));
            diff_Reference = sym(zeros(size(Lprod)));
        else
            [Lprod_Sym, Sprod_Sym] = multiplyLS(hpf(full(valuesToIncrement_LogProbabilities_L)), hpf(full(valuesToIncrement_LogProbabilities_S)),...
                hpf(full(L_BindingEvent_J_Latent)), hpf(full(S_BindingEvent_J_Latent)));
            diff_IncJ_Error = (hpf(full(Lprod))-Lprod_Sym).*hpf(full(Sprod)-Sprod_Sym);
            diff_Reference = hpf(zeros(size(Lprod)));
        end

        if sum((valuesToIncrement_LogProbabilities_S==0).*(S_BindingEvent_J_Latent==0).*(Sprod==1))>0 || ~isequal(diff_IncJ_Error,diff_Reference)
            updateBaseNow = true;
            return;
        end
    end
    % Linear indices in Lnew
    countsIncPlus = countsInc + 1;
    subsColsInc = num2cell(countsIncPlus,1);
    argsInc = [subsColsInc, latentSubsInc];
    linIdxInc = sub2ind(sz,argsInc{:});
    uIdxInc  = unique(linIdxInc);
    linIdxInc_All = [uIdxInc;linIdxInc];
    Lprod_All = [Lnew(uIdxInc); Lprod];
    Sprod_All = [Snew(uIdxInc); Sprod];
    [IncJ_Groups,IncJ_Group_ID] = findgroups(linIdxInc_All);
    mask_Prod_All = Fmask(Lprod_All,Sprod_All);
    num_in_mask = splitapply(@sum,full(mask_Prod_All),IncJ_Groups);
    if opts.PrecisionCheckUpdate
        if ~(ismcc || isdeployed)
            %#exclude sym
            %#exclude eval
            mask_Prod_All_Sym = Fmask(sym(sparse(Lprod_All)),sym(sparse(Sprod_All)));
            num_in_mask_From_Sym = splitapply(@sum,double(simplify(mask_Prod_All_Sym)),IncJ_Groups);
        else
            mask_Prod_All_Sym = Fmask(hpf(full(Lprod_All)),hpf(full(Sprod_All)));
            num_in_mask_From_Sym = splitapply(@sum,mask_Prod_All_Sym,IncJ_Groups);
        end
        if ~isequal(num_in_mask,num_in_mask_From_Sym)
            updateBaseNow = true;
            return;
        end
    end
    [Lmerged_G, Smerged_G] = splitapply(@(Lg,Mask) Lmerge_Func(Lg,Mask,base),full(Lprod_All),full(mask_Prod_All),IncJ_Groups);
    if opts.PrecisionCheckUpdate
        if ~(ismcc || isdeployed)
            %#exclude sym
            %#exclude eval
            [Lmerged_G_Sym,Smerged_G_Sym] = splitapply(@(Lg,Mask) Lmerge_Sym_Func(Lg,Mask,base),sym(sparse(Lprod_All)),simplify(mask_Prod_All_Sym),IncJ_Groups);
            diff_IncJ_Error = (sym(sparse(Lmerged_G))-Lmerged_G_Sym).*sym(Smerged_G - Smerged_G_Sym);
            diff_Reference = sym(zeros(size(Lmerged_G)));
        else
            [Lmerged_G_Sym,Smerged_G_Sym] = splitapply(@(Lg,Mask) Lmerge_Sym_Func(Lg,Mask,base),hpf(full(Lprod_All)),mask_Prod_All_Sym,IncJ_Groups);
            diff_IncJ_Error = (hpf(full(Lmerged_G))-Lmerged_G_Sym).*hpf(full(Smerged_G-Smerged_G_Sym));
            diff_Reference = hpf(zeros(size(Lmerged_G)));
        end
        if ~isequal(diff_IncJ_Error,diff_Reference)
            updateBaseNow = true;
            return;
        end
    end
    Lnew(IncJ_Group_ID(num_in_mask>0)) = Lmerged_G(num_in_mask>0);
    Snew(IncJ_Group_ID(num_in_mask>0)) = Smerged_G(num_in_mask>0);
    if sum((num_in_mask==0))>0
        Lnew(IncJ_Group_ID(num_in_mask==0)) = 0;
        Snew(IncJ_Group_ID(num_in_mask==0)) = 0;
    end
end
end
