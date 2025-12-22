function metrics = F_PMF_generalMetrics(P, N, M, L, R, K, S, Q,metrics_chosen)
% COMPUTEPMFMETRICS_GENERAL General metric computation from joint PMF P
%
% Inputs:
%   P - ndSparse or dense tensor representing the joint PMF
%   N - indices of primary outcomes (metrics computed over these)
%   M - indices of conditional outcomes
%   L - indices of conditional latent dimensions
%   R - indices of marginal outcomes (sum over before metric)
%   K - indices of marginal latent dims (sum over before metric)
%   S - indices of sliced outcomes (compute metric per combination)
%   Q - indices of sliced latent dims (compute metric per combination)
%
% Outputs:
%   metrics - struct with fields:
%       PCC, RCC, MI, NMI, PMI, TC, DTC, InteractionInfo
%
% Notes:
%   - Fully general: supports any number of N, M, L, R, K, S, Q
%   - Conditional and sliced dimensions handled separately
%   - Marginalized dimensions summed over
%   - Sparse tensors supported via ndSparse

% Step 0: defaults
if nargin < 9, Q = []; end
if nargin < 8, S = []; end
if nargin < 7, K = []; end
if nargin < 6, R = []; end
if nargin < 5, L = []; end
if nargin < 4, M = []; end
if nargin < 3, N = []; end

%% Step 0: Validate and permute dimensions
dimsAll = [N(:)', M(:)', L(:)', R(:)', K(:)', S(:)', Q(:)'];
dimsAll = unique(dimsAll,'stable'); % ensure no duplicates
P = permute(P, dimsAll);
targetDims = N;
sliceDims = [S(:)', Q(:)'];
condDims = [M(:)', L(:)'];
margDims = [R(:)', K(:)'];
%% Step 1: Marginalize over R and K
permOrder = [targetDims, setdiff(dimsAll, targetDims)];
P = permute(P, permOrder);
D = numel(targetDims);
%% Step 2: Generate slice combinations over S and Q, Enumerate slices (S,Q)
if isempty(sliceDims)
    sliceGroups = {[]};
else
    sz = size(P);
    szSlice = sz(sliceDims);
    sliceSubs = cell(1,numel(sliceDims));
    if length(sliceDims) >1
        [sliceSubs{:}] = ind2sub(szSlice, 1:prod(szSlice));
    else
        [sliceSubs{:}] = ind2sub([szSlice 1], 1:prod(szSlice));
    end
    sliceGroups = mat2cell(cell2mat(sliceSubs'), numel(sliceDims), ones(1,prod(szSlice)));
end
%% Step 3: Identify unique conditional slices for M+L
if isempty(condDims)
    condGroups = {[]};
else
    sz = size(P);
    szCond = sz(condDims);
    condSubs = cell(1,numel(condDims));
    if length(condDims) >1
    [condSubs{:}] = ind2sub(szCond, 1:prod(szCond));
    else
    [condSubs{:}] = ind2sub([szCond 1], 1:prod(szCond));
    end
    condGroups = mat2cell(cell2mat(condSubs'), numel(condDims), ones(1,prod(szCond)));
end
%% Step 4: Preallocate metric storage
metrics = struct();
metrics.PCC = nan([numel(N), numel(N), max(1,length(condGroups)), max(1,length(sliceGroups))]);
metrics.RCC = nan([numel(N), numel(N), max(1,length(condGroups)), max(1,length(sliceGroups))]);
metrics.SCC = nan([numel(N), numel(N), max(1,length(condGroups)), max(1,length(sliceGroups))]);
metrics.MI  = nan(max(1,length(condGroups)), max(1,length(sliceGroups)));
metrics.NMI = nan(max(1,length(condGroups)), max(1,length(sliceGroups)));
metrics.PMI = cell(max(1,length(condGroups)), max(1,length(sliceGroups)));
metrics.TC  = nan(max(1,length(condGroups)), max(1,length(sliceGroups)));
metrics.DTC = nan(max(1,length(condGroups)), max(1,length(sliceGroups)));
metrics.InteractionInfo = nan(max(1,length(condGroups)), max(1,length(sliceGroups)));
metrics.SI = cell(max(1,length(condGroups)), max(1,length(sliceGroups)));
metrics.scores = cell(max(1,length(condGroups)), max(1,length(sliceGroups)));
%% Step 5: Loop over slices and conditional groups
for cg = 1:numel(condGroups)
    condVals = condGroups{cg};
    % extract conditional slice
    if isempty(condVals)
        Pcond = P;
    else
        idx = repmat({':'}, 1, ndims(P));
        for d = 1:numel(condDims)
            idx{condDims(d)} = condVals(d);
        end
        Pcond = P(idx{:});
    end

    % normalize conditional slice
    if sum(Pcond(:)) > 0
        Pcond = Pcond ./ sum(Pcond(:));
    else
        continue;
    end

    % loop over slice groups
    for sg = 1:numel(sliceGroups)
        sliceVals = sliceGroups{sg};

        if isempty(sliceVals)
            Pslice = Pcond;
        else
            idx = repmat({':'}, 1, ndims(Pcond));
            for d = 1:numel(sliceDims)
                idx{sliceDims(d)} = sliceVals(d);
            end
            Pslice = Pcond(idx{:});
        end

        % marginalize if needed
        if ~isempty(margDims)
            Pslice = sumOverOtherDimsSparse(Pslice, setdiff(1:ndims(Pslice), margDims));
        end

        if sum(Pslice(:)) > 0
            Pslice = Pslice ./ sum(Pslice(:));
        else
            continue;
        end

        %% Step 5a: Linear metrics
        if (ismatrix(Pslice))
            opts.alpha = 2; opts.beta = 2; opts.rho = 1; opts.gamma = 2;
            opts.kappa = 4; opts.w00 = 0; opts.kernel='exp';
            opts.eps = 1e-6;opts.nbins=50;
            scores = compute_colabeling_scores_xy_full(reshape(1:length(find(sum(Pslice,2)>0)),[],1), reshape(1:length(find(sum(Pslice,1)>0)),1,[]), Pslice(1:length(find(sum(Pslice,2)>0)),1:length(find(sum(Pslice,1)>0))), opts);       
            metrics.scores{cg,sg} = scores;
        end
        [PCC,RCC,SCC] = computePCC_RCC_SCC(Pslice,D);
        metrics.PCC(:,:,cg,sg) = PCC;
        metrics.RCC(:,:,cg,sg) = RCC;
        metrics.SCC(:,:,cg,sg) = SCC;
        %% Step 5b: Information metrics
        metrics.MI(cg,sg)  = computeMutualInformationSparse(Pslice);
        metrics.NMI(cg,sg) = computeNormalizedMutualInformationSparse(Pslice);
        metrics.PMI{cg,sg} = computePointwiseMutualInformationSparse(Pslice);
        metrics.TC(cg,sg)  = computeTotalCorrelationSparse(Pslice);
        metrics.DTC(cg,sg) = computeDualTotalCorrelationSparse(Pslice);
        metrics.NTC(cg,sg) = computeNormalizedTotalCorrelationSparse(Pslice);
        metrics.InteractionInfo(cg,sg) = computeInteractionInformationSparse(Pslice);
        metrics.SI{cg,sg} = computeSpecificCorrelationSparse(Pslice);
    end
end
end




function H = computeEntropySparse(P)
% H = -sum p log2 p  (bits). P is normalized (sum(P)=1).
P = P/sum(P(:));
R = spfun(@log2,P);
Q = P.*R;
H = -sum(Q(:));
end
function Psub = sumOverOtherDimsSparse(P, dimsToKeep)
if isempty(dimsToKeep)
    Psub = sum(P(:));
    return
end
allDims = 1:ndims(P);
dimsToSum = setdiff(allDims, dimsToKeep);
Psub = P;
for d = sort(dimsToSum,'descend')
    Psub = sum(Psub, d);
end
end
function [E, EX2, covMat] = computeLinearMomentsSparse(P, D)
subsCell = cell(1,ndims(P)); 
linIndx = find(P);
[subsCell{:}] = ind2sub(size(P), linIndx);
vals = P(linIndx);
if isempty(vals)
    E = zeros(1,D); EX2 = zeros(1,D); covMat = zeros(D);
    return
end
coords = cell2mat(subsCell);
X = double(coords(:,1:D)) - 1;
vals = double(vals(:));
%sum less than equal to X, for marginals?
%average rank get ties
% example rank of P(x=xi,Y=yj) is rank (Ri+Rj)/2;
%0 = R1
E = sum(bsxfun(@times, vals, X), 1);
EX2 = sum(bsxfun(@times, vals, X.^2), 1);
EXY = sum(bsxfun(@times, vals, prod(X,2)), 1);
covMat = zeros(D);
for i = 1:D
    for j = i:D
        m = sum(vals .* (X(:,i) .* X(:,j)));
        covMat(i,j) = m - E(i)*E(j);
        covMat(j,i) = covMat(i,j);
    end
end


end
function [PCC, RCC,SCC] = computePCC_RCC_SCC(P, D)
[E, EX2, covMat] = computeLinearMomentsSparse(P, D);
varVec = EX2 - E.^2;
Dmat = sqrt(varVec(:) * varVec(:)');
PCC = covMat ./ (Dmat + eps);
Ecross = getCrossMoment(P,D);
EX2mat = sqrt(EX2(:) * EX2(:)');
RCC = Ecross ./ (EX2mat + eps);
SCC = NaN*ones(size(PCC));

Px = sum(P,2);
Py = sum(P,1);
Rx_equal = spfun(@cumsum,Px(find(Px>0)));%P(X<=x)
Ry_equal = spfun(@cumsum,Py(find(Py>0)));%P(Y<=y)
Rx_minus = [0; Rx_equal(1:end-1)];
Ry_minus = [0 Ry_equal(1:end-1)];

Ry = (Ry_minus+Ry_equal)/2;
Rx = (Rx_minus+Rx_equal)/2;
Ry = sym(full(Ry_equal - Py(find(Py>0))))+1/2;
Rx = sym(full(Rx_equal - Px(find(Px>0))))+1/2;

Txy = sum((Rx.*Ry).*sym(full(P(find(Px>0),find(Py>0)))),'all')-sum(Rx.*sym(full(Px(find(Px>0)))),'all')*sum(Ry.*sym(full(Py(find(Py>0)))),'all');

Bxy = sqrt((sum(Rx.^2.*sym(full(Px(find(Px>0)))),'all')-sum(Rx.*sym(full(Px(find(Px>0)))),'all')^2).*(sum(Ry.^2.*sym(full(Py(find(Py>0)))),'all')-sum(Ry.*sym(full(Py(find(Py>0)))),'all')^2));

SCC(1,2) = double(Txy/Bxy);
SCC(2,1) = double(Txy/Bxy);
% F = P.*((sumOverOtherDimsSparse(P, 1)-sumOverOtherDimsSparse(P, 2)).^2);
% x = 0:size(P,1)-1; 
% y = 0:size(P,2)-1; 
% 
% CMF1 = cumsum(sumOverOtherDimsSparse(P, 1));
% CMF2= cumsum(sumOverOtherDimsSparse(P, 2));
% CMF12 = cumsum(cumsum(P,1),2);
% SCC  = 1-6*trapz(y,trapz(x, full(F), 1)); 
% 
% PSRC = 1 %12*integral(copula) -3
% %\(\tau =4\iint F(x,y)dF(x,y)-1\)



end
function Ecross = getCrossMoment(P, D)
subsCell = cell(1,ndims(P)); 
linIndx = find(P);
[subsCell{:}] = ind2sub(size(P), linIndx);
vals = P(linIndx);
if isempty(vals)
    Ecross = zeros(D);
    return
end
coords = cell2mat(subsCell);
X = double(coords(:,1:D)) - 1;
vals = double(vals(:));
Ecross = zeros(D);
for i = 1:D
    for j = 1:D
        Ecross(i,j) = sum(vals .* (X(:,i) .* X(:,j)));
    end
end
end
function MI = computeMutualInformationSparse(P)
% MI(X1;X2;...;XN) = sum_i H(Xi) - H(X1..XN)
D = ndims(P);
Hjoint = computeEntropySparse(P);
Hsum = 0;
for d = 1:D
    % marginal over all dims except d
    Pm = sumOverOtherDimsSparse(P, d);
    Hm = computeEntropySparse(Pm);
    Hsum = Hsum + Hm;
end
MI = Hsum - Hjoint;
if (MI<0)
    MI1 = MI;
    PmargProd = sumOverOtherDimsSparse(P, 1).*sumOverOtherDimsSparse(P, 2);
    R = spfun(@log2,P(~iszero(PmargProd))./PmargProd(~iszero(PmargProd)));
    Q = P(~iszero(PmargProd)).*R;
    MI = -sum(Q(:));
end
end
function NMI = computeNormalizedMutualInformationSparse(P)
MI = computeMutualInformationSparse(P);
H = computeEntropySparse(P);
D = ndims(P);
Hsum = 0;
Hind = zeros(1,D);
for d = 1:D
    Pm = sumOverOtherDimsSparse(P, d);
    Hsum = Hsum + computeEntropySparse(Pm);
    Hind(d) = computeEntropySparse(Pm);
end
NMI = MI ./ (Hsum + eps);
DH = 1 - MI./H;
Dmin = MI./min(Hind);
IQR = MI./H - 1;
distMI = D*H - Hsum;

Q = sum(P,1).*sum(P,2);
R = spfun(@log2,P./Q);
R(isnan(R)) = 0;
I = sum(P.*R);
try
NMI_corr = MI./nthroot(prod(Hind),D);
catch
eee = 1;
end
end
function PMI = computePointwiseMutualInformationSparse(P)
subsCell = cell(1,ndims(P)); 
linIndx = find(P);
[subsCell{:}] = ind2sub(size(P), linIndx);
vals = P(linIndx);
if isempty(vals)
    PMI.coords = []; PMI.vals = [];
    return
end
coords = cell2mat(subsCell);
D = size(coords,2);
margMaps = cell(1,D);
for d = 1:D
    [~, ~, ic] = unique(coords(:,d));
    marg = accumarray(ic, full(vals));
    margMaps{d} = marg(ic);
end
prodMarg = ones(size(vals));
for d = 1:D
    prodMarg = prodMarg .* margMaps{d};
end
PMIvals = log2( (vals + eps) ./ (prodMarg + eps) );
PMI.coords = coords;
PMI.vals = PMIvals;
PMI.matrix = ndSparse.build(PMI.coords,PMI.vals,size(P));
end
function TC = computeTotalCorrelationSparse(P)
Hjoint = computeEntropySparse(P);
D = ndims(P);
Hsum = 0;
for d = 1:D
    Pm = sumOverOtherDimsSparse(P, d);
    Hsum = Hsum + computeEntropySparse(Pm);
end
TC = Hsum - Hjoint;
end
function NTC = computeNormalizedTotalCorrelationSparse(P)
TC = computeTotalCorrelationSparse(P);
H = computeEntropySparse(P);
NTC = TC./H;
end

function DTC = computeDualTotalCorrelationSparse(P)
D = ndims(P);
Hjoint = computeEntropySparse(P);
HothersSum = 0;
for d = 1:D
    others = setdiff(1:D, d);
    Pothers = sumOverOtherDimsSparse(P, others);
    HothersSum = HothersSum + computeEntropySparse(Pothers);
end
DTC = (1 - D) * Hjoint + HothersSum;
end
function II = computeInteractionInformationSparse(P)
D = ndims(P);
II = 0;
for k = 1:D
    subsets = nchoosek(1:D, k);
    sign = (-1)^(k+1);
    for r = 1:size(subsets,1)
        sel = subsets(r,:);
        Psel = sumOverOtherDimsSparse(P, sel);
        Hsel = computeEntropySparse(Psel);
        II = II + sign * Hsel;
    end
end
end
function SI = computeSpecificCorrelationSparse(P)
% Computes Specific Correlation (multi-variable PMI generalization)
% Returns struct with coords and SI values
subsCell = cell(1,ndims(P)); 
linIndx = find(P);
[subsCell{:}] = ind2sub(size(P), linIndx);
vals = P(linIndx);
if isempty(vals)
    SI.coords = []; SI.vals = [];
    return
end
coords = cell2mat(subsCell);
D = size(coords,2);

% compute marginals for each dimension
margMaps = cell(1,D);
for d = 1:D
    [~, ~, ic] = unique(coords(:,d));
    marg = accumarray(ic, full(vals));
    margMaps{d} = marg(ic);
end

% product of marginals
prodMarg = ones(size(vals));
for d = 1:D
    prodMarg = prodMarg .* margMaps{d};
end

% specific correlation
SIvals = log2( (vals + eps) ./ (prodMarg + eps) );

SI.coords = coords;
SI.vals   = SIvals;
SI.matrix = ndSparse.build(SI.coords,SI.vals,size(P));
end
