


wA = 331.2;
wC = 307.2;
wG = 347.2;
wT = 322.2;
outGroup2 = [1 4 8 9 11];
for n = 1:length(outGroup2)
probes = ProbeDesignResults3.SoftwareResults(outGroup2(n)).probes;
Pset = ProbeDesignResults3.SoftwareResults(outGroup2(n)).Designed_Probes;
ProbeWeight = @(x) wA*count(seqrcomplement(probes{x,2}),'A')+...
                   wC*count(seqrcomplement(probes{x,2}),'C')+...
                   wG*count(seqrcomplement(probes{x,2}),'T')+...
                   wT*count(seqrcomplement(probes{x,2}),'G');
AverageWeight(n) = mean(arrayfun(@(x) ProbeWeight(x),Pset));
TotalWeight(n) = sum(arrayfun(@(x) ProbeWeight(x),Pset));
end