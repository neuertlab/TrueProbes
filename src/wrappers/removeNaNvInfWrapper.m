<<<<<<< HEAD
<<<<<<< HEAD
function E = removeNaNvInfWrapper(A)
E = A;
E(isinf(A))=0;
E(isnan(A))=0;
=======
function E = removeNaNvInfWrapper(A)
E = A;
E(isinf(A))=0;
E(isnan(A))=0;
>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
=======
function E = removeNaNvInfWrapper(A)
E = A;
E(isinf(A))=0;
E(isnan(A))=0;
>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
end