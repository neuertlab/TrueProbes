function E = removeNaNvInfWrapper(A)
E = A;
E(isinf(A))=0;
E(isnan(A))=0;
end