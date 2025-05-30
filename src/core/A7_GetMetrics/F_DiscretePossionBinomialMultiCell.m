function [Px1] = F_DiscretePossionBinomialMultiCell(p)
n_genes = size(p,1);n_sites = size(p,2);n_cells = size(p,3);
s = -(1-p)./p;
e = permute(s,[2 3 1]);
e1 = permute(s,[2 3 1]);
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
q = permute(CATnWrapper(arrayfun(@(target) CATnWrapper(arrayfun(@(cell) flip_function(q0(:,cell,target),e1(:,cell,target)),1:n_cells,'Un',0),2),1:n_genes,'Un',0),3),[3 1 2]);
p1 = p;
p1(p1==0) = 1;
alphas = permute(repmat(CATnWrapper(arrayfun(@(cell) prod(full(p1(:,:,cell)),2),1:n_cells,'Un',0),2),[1 1 n_sites+1]),[1 3 2]);%is alpha needed if i just normalize it in the end
Px1 = fliplr(alphas.*q);%Probability distribution added in reverse order of coefficients
Px1 = Px1./permute(repmat(squeeze(sum(Px1,2)),[1 1 size(Px1,2)]),[1 3 2]);
Qx1 = cumsum(Px1,2);%Cumulative distribution
end