function [Px1] = F_DiscretePossionBinomialMultiCell(p)
%% This function generates the poission binomial distribution for the odds
% that any number of targets binding sites are bound by probes given the
% probability that each binding site is bound by a probe. This distribution 
% is computed parallelized for every on/off-target using the possion binomials 
% probability generating function. The Poisson binomial is a distribution 
% of the sum of independent Bernoulli trials that can have distinct 
% success probability. This is for Probability Targets x Sites x Cell Lines.
n_genes = size(p,1);n_sites = size(p,2);n_cells = size(p,3);
p1 = p;
p1(p1==0) = 1;
alphas = permute(repmat(CATnWrapper(arrayfun(@(cell) prod(full(p1(:,:,cell)),2),1:n_cells,'Un',0),2),[1 1 n_sites+1]),[1 3 2]);%is alpha needed if i just normalize it in the end
s = -(1-p)./p;
e1 = permute(s,[2 3 1]);
e = permute(s,[2 3 1]);
% Strip out infinities
e(isnan(e)) = 0;
e(isinf(e)) = 0;
%diff between being empty no finite.
% Expand recursion formula
q0 = CATnWrapper(arrayfun(@(C) CATnWrapper(arrayfun(@(C) [1 zeros(1,size(e,1))]',1:size(s,3),'Un',0),2),1:size(s,1),'Un',0),3);
for j=1:n_sites
    q0(2:(j+1),:,:) = q0(2:(j+1),:,:) - e(j,:,:).*q0(1:j,:,:);
end
% if only one site is bound??
if isequal(sort(e(imag(e)>0)),sort(conj(e(imag(e)<0))))
    q0 = real(q0);
end
flip_function = @(Row,Status) flip([flip(Row(1:1+sum(~isinf(Status)))); Row(2+sum(~isinf(Status)):end)]); 
if (sum(isinf(q0)+isnan(q0),'all')+sum(alphas==0,'all')==0)
q = permute(CATnWrapper(arrayfun(@(target) CATnWrapper(arrayfun(@(cell) flip_function(q0(:,cell,target),e1(:,cell,target)),1:n_cells,'Un',0),2),1:n_genes,'Un',0),3),[3 1 2]);
%add step that checks if in all cells if there is a cell with zero probes bound for target to get zero bound distribution and not NaN.
Px1 = fliplr(alphas.*q);%Probability distribution added in reverse order of coefficients
Px1 = Px1./permute(repmat(squeeze(sum(Px1,2)),[1 1 size(Px1,2)]),[1 3 2]);
else
alphas_sym = permute(repmat(CATnWrapper(arrayfun(@(cell) prod(sym(full(p1(:,:,cell))),2),1:n_cells,'Un',0),2),[1 1 n_sites+1]),[1 3 2]);%is alpha needed if i just normalize it in the end    
q0s = CATnWrapper(arrayfun(@(C) CATnWrapper(arrayfun(@(C) [1 zeros(1,size(e,1))]',1:size(s,3),'Un',0),2),1:size(s,1),'Un',0),3);
q0_sym = sym(q0s);
e_sym = sym(e);
e1_sym = sym(e1);
for j=1:n_sites
    q0_sym(2:(j+1),:,:) = q0_sym(2:(j+1),:,:) - repmat(e_sym(j,:,:),[j 1 1]).*q0_sym(1:j,:,:);
end
if isequal(sort(e_sym(imag(e_sym)>0)),sort(conj(e_sym(imag(e_sym)<0))))
    q0_sym = real(q0_sym);
end
q_sym = permute(CATnWrapper(arrayfun(@(target) CATnWrapper(arrayfun(@(cell) flip_function(q0_sym(:,cell,target),e1_sym(:,cell,target)),1:n_cells,'Un',0),2),1:n_genes,'Un',0),3),[3 1 2]);
Px1_sym = fliplr(alphas_sym.*q_sym);%Probability distribution added in reverse order of coefficients
Px1_normFactor_sym =permute(repmat(squeeze(sum(Px1_sym,2)),[1 1 size(Px1_sym,2)]),[1 3 2]);
Px1normalized_sym = Px1_sym./Px1_normFactor_sym;
Px1 = double(Px1normalized_sym);
end
Qx1 = cumsum(Px1,2);%Cumulative distribution
end