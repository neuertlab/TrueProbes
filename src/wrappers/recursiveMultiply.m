function E = recursiveMultiply(X)
E = X{1};
for k = 1:length(X)-1
   E = E.*X{k+1}; 
end
end
