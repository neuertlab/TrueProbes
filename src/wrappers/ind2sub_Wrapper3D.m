function w_out = ind2sub_Wrapper3D(A,k)
[y,x,z] = ind2sub(size(A),find(A==k));
w_out.x = x;
w_out.y = y;
w_out.z = z;
end