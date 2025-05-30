function D = getD(seq,K_Context,K_values)
% This function parses meso-scale model triple nucleotide binding sequence 
% cases to get parameter values for Morse Potential Depth.

   Idx = find(strcmp(seq,K_Context));
   if (length(seq)==5&&isempty(Idx))
      Idx2 = find(strcmp(reverse(seq),K_Context));
      Idx = union(Idx,Idx2);
   end
   if (length(seq)==9&&isempty(Idx))
      seqr = strcat(seq(1),reverse(seq(6:8)),seq(5),reverse(seq(2:4)),seq(9));
      Idx2 = find(strcmp(seqr,K_Context));
      Idx = union(Idx,Idx2);
   end
   if (length(seq)==11&&isempty(Idx))
      seqr = strcat(seq(1),reverse(seq(7:10)),seq(6),reverse(seq(2:5)),seq(11));
      Idx2 = find(strcmp(seqr,K_Context));
      Idx = union(Idx,Idx2);
   end
   if (~isempty(Idx))
       D = mean(cell2mat(K_values(Idx)));
   else
       D = 0;
   end
end