<<<<<<< HEAD
<<<<<<< HEAD
function w_out = ind2sub_Wrapper3D(A,k)
[y,x,z] = ind2sub(size(A),find(A==k));
w_out.x = x;
w_out.y = y;
w_out.z = z;
=======
function w_out = ind2sub_Wrapper3D(A,k)
[y,x,z] = ind2sub(size(A),find(A==k));
w_out.x = x;
w_out.y = y;
w_out.z = z;
>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
=======
function w_out = ind2sub_Wrapper3D(A,k)
[y,x,z] = ind2sub(size(A),find(A==k));
w_out.x = x;
w_out.y = y;
w_out.z = z;
>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
end