function [Px1] = F_DiscretePossionBinomialMulti(p)
% This function generates the poission binomial distribution for the odds
% that any number of targets binding sites are bound by probes given the
% probability that each binding site is bound by a probe. This distribution 
% is computed parallelized for every on/off-target using the possion binomials 
% probability generating function. The Poisson binomial is a distribution 
% of the sum of independent Bernoulli trials that can have distinct 
% success probability.

% Note: this code removes cases where binding probability is zero as this 
% messes up the function by adding infinity terms into the transformation.

   s = -(1-p)./p;
   q = cell2mat(arrayfun(@(x) [zeros(1,sum(isinf(s(x,:)))) poly(s(x,:))]',1:size(p,1),'Un',0))';%poly(s)
   p1 = p;
   p1(p1==0) = 1;
   alphas = prod(full(p1),2);%is alpha needed if i just normalize it in the end
   Px1 = fliplr(alphas.*q);%Probability distribution added in reverse order of coefficients
   Px1 = Px1./repmat(sum(Px1,2),[1 size(Px1,2)]);
   Qx1 = cumsum(Px1,2);%Cumulative distribution
end