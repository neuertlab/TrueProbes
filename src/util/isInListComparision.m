function out = isInListComparision(x,C_GList_T,GList_T_ic)
try
   out = find(find(strcmp(x,C_GList_T))==GList_T_ic,1);
catch
   out = 0; 
end

end