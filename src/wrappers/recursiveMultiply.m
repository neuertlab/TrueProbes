function E = recursiveMultiply(X)
n = length(X);
E = X{1};
for k = 1:length(X(
   E = E.*X{k+1}; 
end

end
