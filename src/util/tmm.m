<<<<<<< HEAD
<<<<<<< HEAD
function [norm_factors, norm_matrix] = tmm(tpm_matrix, logratioTrim, sumTrim,Acutoff,doWeighting,ref_idx)
% TMM normalization inspired by NOISeq::tmm with doWeighting = TRUE
% tpm_matrix: G x S matrix (genes x samples) of TPM values
% logratioTrim: proportion to trim from each tail of M-values (default 0.3)
if nargin < 2
    logratioTrim = 0.3;
end
if nargin < 3
    sumTrim = 0.05;
end
if nargin < 4
    Acutoff = -1e10;
end
if nargin <5
    doWeighting = 1;
end
if nargin <6
    ref_idx = 1;
end
S = size(tpm_matrix,2);
lib_sizes = sum(tpm_matrix, 1);

%remove all zero rows
non_zero_tpm_matrix = tpm_matrix;
non_zero_tpm_matrix(sum(tpm_matrix,2)==0,:) = []; 

% Choose reference as sample with median library size
f75 = ones(1,S);
f75 = quantile(non_zero_tpm_matrix./lib_sizes,0.75,1);
if(median(f75) < 1e-20)
    [~,ref_idx] = max(sum(sqrt(non_zero_tpm_matrix),1));
else
    [~,ref_idx] = min(abs(f75-mean(f75)));
end

ref = non_zero_tpm_matrix(:, ref_idx);
nR = sum(ref);
% Initialize output
norm_factors = ones(1, S);
for j = 1:S
    if j == ref_idx
        continue;  % Reference gets a scaling factor of 1
    end
    obs = non_zero_tpm_matrix(:, j);
    nO = sum(obs);
    logR = log2((obs/nO)./(ref/nR));
    absE = (log2(obs/nO) + log2(ref/nR))/2;
    asymptotic_v = (nO-obs)/nO./obs + (nR-ref)./ref;
    keep = isfinite(logR).*isfinite(absE).*double(absE>Acutoff);
    logR = logR(keep==1);
    absE = absE(keep==1);
    asymptotic_v = asymptotic_v(keep==1);
    n = sum(keep);
    loL = floor(n*logratioTrim)+1;
    hiL = n + 1 - loL;
    loS = floor(n*sumTrim)+1;
    hiS = n + 1 - loS;
    keep = double(tiedrank(logR)>=loL).*double(tiedrank(logR)<=hiL).*double(tiedrank(absE)>=loS).*double(tiedrank(absE)<hiS);
    if (doWeighting)
    norm_factors(j) = 2^(sum(logR(keep==1)./asymptotic_v(keep==1),'omitnan')/sum(1./asymptotic_v(keep==1),'omitnan'));
    else
    norm_factors(j) = 2^mean(logR(keep==1),'omitnan');    
    end
end
norm_factors = norm_factors./exp(mean(log(norm_factors)));

% Normalize TPM matrix
norm_matrix = tpm_matrix ./ norm_factors;

=======
function [norm_factors, norm_matrix] = tmm(tpm_matrix, logratioTrim, sumTrim,Acutoff,doWeighting,ref_idx)
% TMM normalization inspired by NOISeq::tmm with doWeighting = TRUE
% tpm_matrix: G x S matrix (genes x samples) of TPM values
% logratioTrim: proportion to trim from each tail of M-values (default 0.3)
if nargin < 2
    logratioTrim = 0.3;
end
if nargin < 3
    sumTrim = 0.05;
end
if nargin < 4
    Acutoff = -1e10;
end
if nargin <5
    doWeighting = 1;
end
if nargin <6
    ref_idx = 1;
end
S = size(tpm_matrix,2);
lib_sizes = sum(tpm_matrix, 1);

%remove all zero rows
non_zero_tpm_matrix = tpm_matrix;
non_zero_tpm_matrix(sum(tpm_matrix,2)==0,:) = []; 

% Choose reference as sample with median library size
f75 = ones(1,S);
f75 = quantile(non_zero_tpm_matrix./lib_sizes,0.75,1);
if(median(f75) < 1e-20)
    [~,ref_idx] = max(sum(sqrt(non_zero_tpm_matrix),1));
else
    [~,ref_idx] = min(abs(f75-mean(f75)));
end

ref = non_zero_tpm_matrix(:, ref_idx);
nR = sum(ref);
% Initialize output
norm_factors = ones(1, S);
for j = 1:S
    if j == ref_idx
        continue;  % Reference gets a scaling factor of 1
    end
    obs = non_zero_tpm_matrix(:, j);
    nO = sum(obs);
    logR = log2((obs/nO)./(ref/nR));
    absE = (log2(obs/nO) + log2(ref/nR))/2;
    asymptotic_v = (nO-obs)/nO./obs + (nR-ref)./ref;
    keep = isfinite(logR).*isfinite(absE).*double(absE>Acutoff);
    logR = logR(keep==1);
    absE = absE(keep==1);
    asymptotic_v = asymptotic_v(keep==1);
    n = sum(keep);
    loL = floor(n*logratioTrim)+1;
    hiL = n + 1 - loL;
    loS = floor(n*sumTrim)+1;
    hiS = n + 1 - loS;
    keep = double(tiedrank(logR)>=loL).*double(tiedrank(logR)<=hiL).*double(tiedrank(absE)>=loS).*double(tiedrank(absE)<hiS);
    if (doWeighting)
    norm_factors(j) = 2^(sum(logR(keep==1)./asymptotic_v(keep==1),'omitnan')/sum(1./asymptotic_v(keep==1),'omitnan'));
    else
    norm_factors(j) = 2^mean(logR(keep==1),'omitnan');    
    end
end
norm_factors = norm_factors./exp(mean(log(norm_factors)));

% Normalize TPM matrix
norm_matrix = tpm_matrix ./ norm_factors;

>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
=======
function [norm_factors, norm_matrix] = tmm(tpm_matrix, logratioTrim, sumTrim,Acutoff,doWeighting,ref_idx)
% TMM normalization inspired by NOISeq::tmm with doWeighting = TRUE
% tpm_matrix: G x S matrix (genes x samples) of TPM values
% logratioTrim: proportion to trim from each tail of M-values (default 0.3)
if nargin < 2
    logratioTrim = 0.3;
end
if nargin < 3
    sumTrim = 0.05;
end
if nargin < 4
    Acutoff = -1e10;
end
if nargin <5
    doWeighting = 1;
end
if nargin <6
    ref_idx = 1;
end
S = size(tpm_matrix,2);
lib_sizes = sum(tpm_matrix, 1);

%remove all zero rows
non_zero_tpm_matrix = tpm_matrix;
non_zero_tpm_matrix(sum(tpm_matrix,2)==0,:) = []; 

% Choose reference as sample with median library size
f75 = ones(1,S);
f75 = quantile(non_zero_tpm_matrix./lib_sizes,0.75,1);
if(median(f75) < 1e-20)
    [~,ref_idx] = max(sum(sqrt(non_zero_tpm_matrix),1));
else
    [~,ref_idx] = min(abs(f75-mean(f75)));
end

ref = non_zero_tpm_matrix(:, ref_idx);
nR = sum(ref);
% Initialize output
norm_factors = ones(1, S);
for j = 1:S
    if j == ref_idx
        continue;  % Reference gets a scaling factor of 1
    end
    obs = non_zero_tpm_matrix(:, j);
    nO = sum(obs);
    logR = log2((obs/nO)./(ref/nR));
    absE = (log2(obs/nO) + log2(ref/nR))/2;
    asymptotic_v = (nO-obs)/nO./obs + (nR-ref)./ref;
    keep = isfinite(logR).*isfinite(absE).*double(absE>Acutoff);
    logR = logR(keep==1);
    absE = absE(keep==1);
    asymptotic_v = asymptotic_v(keep==1);
    n = sum(keep);
    loL = floor(n*logratioTrim)+1;
    hiL = n + 1 - loL;
    loS = floor(n*sumTrim)+1;
    hiS = n + 1 - loS;
    keep = double(tiedrank(logR)>=loL).*double(tiedrank(logR)<=hiL).*double(tiedrank(absE)>=loS).*double(tiedrank(absE)<hiS);
    if (doWeighting)
    norm_factors(j) = 2^(sum(logR(keep==1)./asymptotic_v(keep==1),'omitnan')/sum(1./asymptotic_v(keep==1),'omitnan'));
    else
    norm_factors(j) = 2^mean(logR(keep==1),'omitnan');    
    end
end
norm_factors = norm_factors./exp(mean(log(norm_factors)));

% Normalize TPM matrix
norm_matrix = tpm_matrix ./ norm_factors;

>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
end