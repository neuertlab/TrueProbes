function [Px1] = F_DiscretePossionBinomialStable_2D(p)
   C_Func = @(n) exp(2*pi*sqrt(-1)/(n+1)); 
   Sub_Func1 = @(P,n,L,M_Basis) prod(cell2mat(arrayfun(@(m)1+(C_Func(n)^L-1)*P(:,m),M_Basis,'Un',0)),2); 
   Sub_Func2 = @(P,n,k,M_Basis) sum(cell2mat(arrayfun(@(L) Sub_Func1(P,n,L,M_Basis)*C_Func(n)^(-L*k),0:n,'Un',0)),2)/(n+1);
   
   %break_up into regions num sites less than 20 then do convolution
   M_Basis = 1:28;
   Px1 = cell2mat(arrayfun(@(k) Sub_Func2(p,length(M_Basis),k,M_Basis),0:length(M_Basis)-1,'Un',0));
   
   
   Px1 = cell2mat(arrayfun(@(k) Sub_Func2(p,size(p,2),k,M_Basis),0:size(p,2),'Un',0));
   
   Px1 = real(Px1);
   Px1 = Px1./repmat(sum(Px1,2),[1 size(Px1,2)]);
   
   %problem when have ones that are zero
   %mean(sum(p)),  var =  sumP*(1-p)
%    sum(p,2);
%    sum(p*(1-p),2);
   p(1,p(1,:)>0)
   full(p(144,p(144,:)>0))
   %F_DiscretePossionBinomial(p(1,p(1,:)>0))
   ss = find(p(144,:)>0)
   F_DiscretePossionBinomial(p(144,ss(1:10)))
   %N = length(p) p>0
   %alpha prod
   % s = -(1-p)./p
   % alpha*poly(s);
   p = [0 1/2 1/2 1/2 0.99];
   x = F_DiscretePossionBinomial(p(2:end));
   0   0.125000000000000   0.375000000000000   0.375000000000000   0.125000000000000
   p1 = p;
   p1(p1==0) = 1;
   alpha = prod(p1)
   s = -(1-p)./p;
   alpha*poly(s);
   
end