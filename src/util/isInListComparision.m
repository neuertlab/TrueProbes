<<<<<<< HEAD
<<<<<<< HEAD
function out = isInListComparision(x,C_GList_T,GList_T_ic)
try
   out = find(find(strcmp(x,C_GList_T))==GList_T_ic,1);
catch
   out = 0; 
end

=======
function out = isInListComparision(x,C_GList_T,GList_T_ic)
try
   out = find(find(strcmp(x,C_GList_T))==GList_T_ic,1);
catch
   out = 0; 
end

>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
=======
function out = isInListComparision(x,C_GList_T,GList_T_ic)
try
   out = find(find(strcmp(x,C_GList_T))==GList_T_ic,1);
catch
   out = 0; 
end

>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
end