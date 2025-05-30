function [Px1] = F_DiscretePossionBinomialMultiCellMultiModel(p)
n_genes = size(p,1);n_sites = size(p,2);n_models = size(p,3);n_cells = size(p,4);
s = -(1-p)./p;
e = permute(s,[2 3 4 1]); %[Sites Models Cells Targets]
e1 = permute(s,[2 3 4 1]);%[Sites Models Cells Targets]
% Strip out infinities
e(isnan(e)) = 0;
e(isinf(e)) = 0;
%diff between being empty no finite.
% Expand recursion formula
q0 = CATnWrapper(arrayfun(@(C) CATnWrapper(arrayfun(@(C) CATnWrapper(arrayfun(@(C) [1 zeros(1,size(e,1))]',1:size(s,3),'Un',0),2),1:size(s,4),'Un',0),3),1:size(s,1),'Un',0),4);%[Sites Models Cell Target]
for j=1:n_sites
    q0(2:(j+1),:,:,:) = q0(2:(j+1),:,:,:) - e(j,:,:,:).*q0(1:j,:,:,:);
end
if isequal(sort(e(imag(e)>0)),sort(conj(e(imag(e)<0))))
    q0 = real(q0);
end
flip_function = @(Row,Status) flip([flip(Row(1:1+sum(~isinf(Status)))); Row(2+sum(~isinf(Status)):end)]); 
q = permute(CATnWrapper(arrayfun(@(target) CATnWrapper(arrayfun(@(cell) CATnWrapper(arrayfun(@(model) flip_function(q0(:,model,cell,target),e1(:,model,cell,target)),1:n_models,'Un',0),2),1:n_cells,'Un',0),3),1:n_genes,'Un',0),4),[4 1 2 3]);
%q [Target Sites Models Cell]
p1 = p;
p1(p1==0) = 1;
alphas = permute(repmat(CATnWrapper(arrayfun(@(cell) CATnWrapper(arrayfun(@(model) prod(full(squeeze(p1(:,:,model,cell))),2),1:n_models,'Un',0),2),1:n_cells,'Un',0),3),[1 1 1 n_sites+1]),[1 4 2 3]);%is alpha needed if i just normalize it in the end
Px1 = fliplr(alphas.*q);%Probability distribution added in reverse order of coefficients
Px1 = Px1./permute(repmat(squeeze(sum(Px1,2)),[1 1 1 size(Px1,2)]),[1 4 2 3]);
Qx1 = cumsum(Px1,2);%Cumulative distribution
end