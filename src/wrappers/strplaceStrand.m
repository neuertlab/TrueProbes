function o = strplaceStrand(vec)
   if (prod(vec)>0)
       o = 'Plus/Plus';
   else
       o = 'Plus/Minus';
   end
end