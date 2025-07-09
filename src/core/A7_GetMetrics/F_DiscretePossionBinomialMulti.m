function [Px1] = F_DiscretePossionBinomialMulti(p)
%% This function generates the poission binomial distribution for the odds
% that any number of targets binding sites are bound by probes given the
% probability that each binding site is bound by a probe. This distribution
% is computed parallelized for every on/off-target using the possion binomials
% probability generating function. The Poisson binomial is a distribution
% of the sum of independent Bernoulli trials that can have distinct
% success probability. This is for Probability Targets x Sites.
n_genes = size(p,1);n_sites = size(p,2);
p1 = p;
p1(p1==0) = 1;
alphas = repmat(prod(full(p1),2),[1 n_sites+1]);%is alpha needed if i just normalize it in the end
s = -(1-p)./p;
e1 = permute(s,[2 1]);
e = permute(s,[2 1]);
% Strip out infinities
e(isnan(e)) = 0;
e(isinf(e)) = 0;
%diff between being empty no finite.
% Expand recursion formula
q0 = CATnWrapper(arrayfun(@(C) [1 zeros(1,size(e,1))]',1:size(s,1),'Un',0),2);
for j=1:n_sites
    q0(2:(j+1),:) = q0(2:(j+1),:) - e(j,:).*q0(1:j,:);
end
% if only one site is bound??
if isequal(sort(e(imag(e)>0)),sort(conj(e(imag(e)<0))))
    q0 = real(q0);
end
flip_function = @(Row,Status) flip([flip(Row(1:1+sum(~isinf(Status)))); Row(2+sum(~isinf(Status)):end)]);
if (sum(isinf(q0)+isnan(q0),'all')+sum(alphas==0,'all')==0)
    q = permute(CATnWrapper(arrayfun(@(target) flip_function(q0(:,target),e1(:,target)),1:n_genes,'Un',0),2),[2 1]);
    Px1 = fliplr(alphas.*q);%Probability distribution added in reverse order of coefficients
    Px1 = Px1./repmat(sum(Px1,2),[1 size(Px1,2)]);
else
    if ~(ismcc || isdeployed)
        %#exclude sym
        alphas_sym = repmat(prod(sym(full(p1)),2),[1 n_sites+1]);%is alpha needed if i just normalize it in the end
        q0s = CATnWrapper(arrayfun(@(C) [1 zeros(1,size(e,1))]',1:size(s,1),'Un',0),2);
        q0_sym = sym(q0s);
        e_sym = sym(e);
        e1_sym = sym(e1);
        for j=1:n_sites
            q0_sym(2:(j+1),:) = q0_sym(2:(j+1),:) - repmat(e_sym(j,:),[j 1]).*q0_sym(1:j,:);
        end
        if isequal(sort(e_sym(imag(e_sym)>0)),sort(conj(e_sym(imag(e_sym)<0))))
            q0_sym = real(q0_sym);
        end
        q_sym = permute(CATnWrapper(arrayfun(@(target) flip_function(q0_sym(:,target),e1_sym(:,target)),1:n_genes,'Un',0),2),[2 1]);
        Px1_sym = fliplr(alphas_sym.*q_sym);%Probability distribution added in reverse order of coefficients
        Px1_normFactor_sym = repmat(squeeze(sum(Px1_sym,2)),[1 size(Px1_sym,2)]);
        Px1normalized_sym = Px1_sym./Px1_normFactor_sym;
        Px1 = double(Px1normalized_sym);
    else
        alphas_hpf = repmat(prod(hpf(full(p1)),2),[1 n_sites+1]);%is alpha needed if i just normalize it in the end
        q0s = CATnWrapper(arrayfun(@(C) [1 zeros(1,size(e,1))]',1:size(s,1),'Un',0),2);
        q0_hpf = hpf(q0s);
        e_hpf = hpf(e);
        e1_hpf = hpf(e1);
        for j=1:n_sites
            q0_hpf(2:(j+1),:) = q0_hpf(2:(j+1),:) - repmat(e_hpf(j,:),[j 1]).*q0_hpf(1:j,:);
        end
        q_hpf = permute(CATnWrapper(arrayfun(@(target) flip_function(q0_hpf(:,target),e1_hpf(:,target)),1:n_genes,'Un',0),2),[2 1]);
        Px1_hpf = fliplr(alphas_hpf.*q_hpf);%Probability distribution added in reverse order of coefficients
        Px1_normFactor_hpf = repmat(squeeze(sum(Px1_hpf,2)),[1 size(Px1_hpf,2)]);
        Px1normalized_hpf = Px1_hpf./Px1_normFactor_hpf;
        Px1 = double(Px1normalized_hpf);
    end
end
Qx1 = cumsum(Px1,2);%Cumulative distribution
end