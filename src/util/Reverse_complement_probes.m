function [rev_comp] = Reverse_complement_probes(probes)
%This takes a list of DNA sequences and makes the next column column containing the reverse complement 
%Import probes from excel as cell array, then change what probes is equal
% to below
rev_comp = {};
for i = 1:size(probes,1)
    seq = probes{i,1};
    forw = 1:size(seq,2);
    rev = size(seq,2) + 1 - forw;
    revcom = '';
    for j = forw
        if seq(rev(j)) == 'a' | seq(rev(j)) == 'A'
            revcom(j) = 't';
        elseif seq(rev(j)) == 't' | seq(rev(j)) == 'T'
            revcom(j) = 'a';
        elseif seq(rev(j)) == 'g' | seq(rev(j)) == 'G'
            revcom(j) = 'c';
        elseif seq(rev(j)) == 'c' | seq(rev(j)) == 'C'
            revcom(j) = 'g';
        end
    end
    rev_comp{i,1} = revcom;
=======
function [rev_comp] = Reverse_complement_probes(probes)
%This takes a list of DNA sequences and makes the next column column containing the reverse complement 
%Import probes from excel as cell array, then change what probes is equal
% to below
rev_comp = {};
for i = 1:size(probes,1)
    seq = probes{i,1};
    forw = 1:size(seq,2);
    rev = size(seq,2) + 1 - forw;
    revcom = '';
    for j = forw
        if seq(rev(j)) == 'a' | seq(rev(j)) == 'A'
            revcom(j) = 't';
        elseif seq(rev(j)) == 't' | seq(rev(j)) == 'T'
            revcom(j) = 'a';
        elseif seq(rev(j)) == 'g' | seq(rev(j)) == 'G'
            revcom(j) = 'c';
        elseif seq(rev(j)) == 'c' | seq(rev(j)) == 'C'
            revcom(j) = 'g';
        end
    end
    rev_comp{i,1} = revcom;
>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
=======
function [rev_comp] = Reverse_complement_probes(probes)
%This takes a list of DNA sequences and makes the next column column containing the reverse complement 
%Import probes from excel as cell array, then change what probes is equal
% to below
rev_comp = {};
for i = 1:size(probes,1)
    seq = probes{i,1};
    forw = 1:size(seq,2);
    rev = size(seq,2) + 1 - forw;
    revcom = '';
    for j = forw
        if seq(rev(j)) == 'a' | seq(rev(j)) == 'A'
            revcom(j) = 't';
        elseif seq(rev(j)) == 't' | seq(rev(j)) == 'T'
            revcom(j) = 'a';
        elseif seq(rev(j)) == 'g' | seq(rev(j)) == 'G'
            revcom(j) = 'c';
        elseif seq(rev(j)) == 'c' | seq(rev(j)) == 'C'
            revcom(j) = 'g';
        end
    end
    rev_comp{i,1} = revcom;
>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
end