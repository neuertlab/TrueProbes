function [TargetSiteLocations,ProbeTargetSiteAlignment_NumMatchesAndMisMatches,ProbeTargetSiteAlignment_DanglingEnds,ProbeTargetSiteAlignment_GapLengthAndSegments,ProbeTargetSiteAlignment_Score,ProbeTargetSiteAlignment_BitScore,ProbeTargetSiteAlignment_Evalue] = A_AlignmentProperties_SiteMapFormatted(gene_table,settings,DoesProbeBindSite2,MolProbesAtEvents,Mol_ProbesAtEventsID)
gene_table = sortrows(gene_table,[7 6],'ascend');
gene_table = gene_table(gene_table.Match>=settings.MinHomologySearchTargetSize,:);
MinusStrandedHits = find(contains(gene_table.Strand,'Minus'));
gene_table_NamesZ = convertCharsToStrings(gene_table.Name);
contains_RNA = find(ismember(gene_table_NamesZ,settings.RNAdbParser));clear gene_table_NamesZ
RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);clear contains_RNA
gene_table = gene_table(setdiff(1:size(gene_table,1),RNA_MissedFilteredHits),:);clear RNA_MissedFilteredHits
gene_table.Ax = min(gene_table.SubjectIndices,[],2);
gene_table.Bx = max(gene_table.SubjectIndices,[],2);
gene_table = sortrows(gene_table,[7 13],'ascend');
gene_table_probe_lengths = cellfun(@length,gene_table.ProbeSequence);
gene_table_alignment_three_prime_probe_dangling_end = gene_table.QueryIndices(:,1)-1;
gene_table_alignment_five_prime_probe_dangling_end = gene_table_probe_lengths-gene_table.QueryIndices(:,2);
gene_table_alignment_n_total_gap_length = cellfun(@(x) sum(max(nt2int(x([1 3],:)),[],1)>4),gene_table.Alignment);
gene_table_alignment_n_gap_segments = cellfun(@(x) sum(diff([false, logical(max(nt2int(x([1 3],:)),[],1)>4), false])==1),gene_table.Alignment);
gene_table_alignment_n_subsitutions =  cellfun(@(x) sum(diff(double(nt2int(x([1 3],max(nt2int(x([1 3],:)),[],1)<5))))~=0),gene_table.Alignment);
gene_table_alignment_lengths = cellfun(@(x) size(x,2) ,gene_table.Alignment);
gene_table_alignment_n_canonical_matches = gene_table_alignment_lengths - gene_table_alignment_n_subsitutions - gene_table_alignment_n_total_gap_length;%also from identity
%alignment conversion frequency
%gene_table_alignment_n_watson_and_subsitution_conversions = CATnWrapper(cellfun(@(x) reshape(accumarray(double(nt2int(x([1 3],max(nt2int(x([1 3],:)),[],1)<5)))', 1, [4 4], @sum, 0), [], 1)', gene_table.Alignment(1:3), 'Un', 0),1)
%[int2nt(repelem((1:N)',N)),int2nt(repmat((1:N)',N,1)),int2str(reshape(accumarray(double(nt2int(x([1 3],max(nt2int(x([1 3],:)),[],1)<5)))', 1, [4 4], @sum, 0), [], 1))]
%rot90([int2nt(repelem((1:N)',N)),int2nt(repmat((1:N)',N,1)),int2str(reshape(accumarray(double(nt2int(x([13],max(nt2int(x([1 3],:)),[],1)<5)))', 1, [4 4], @sum, 0), [], 1))]) 
%flip([int2nt(repelem((1:N)',N)),int2nt(repmat((1:N)',N,1))]')'  (P->T) probe nucleotide to other nucleotide
 % targetMatch = arrayfun(@(x) strrep(gene_table.Alignment{x}(3,:),'-','N'),1:size(gene_table,1),'UniformOutput',false);%slow not so slow
 %    probeMatch = arrayfun(@(x) seqrcomplement(strrep(gene_table.Alignment{x}(1,:),'-','N')),1:size(gene_table,1),'UniformOutput',false);%slow
 % maybe also report bit score as well in probe checker general 10
Lambda = 1.374;
Kappa = 0.71;
BitScore_Function = @(L,K,S0) (L*S0-log(K))/log(2);
% settings.BlastParameters
%https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/algo/blast/core/blast_stat.c
%line 608 and lower
%more off-targets that pass when we do not account for identity length vs alignment length
%alignment length min of 15nt,
% 100*sum(gene_table_alignment_lengths>=settings.MinHomologySearchTargetSize)/length(gene_table_alignment_lengths);
% 100*sum(gene_table.Match>=settings.MinHomologySearchTargetSize)/length(gene_table.Match);
% 100*sum(gene_table_alignment_n_canonical_matches>=settings.MinHomologySearchTargetSize)/length(gene_table_alignment_n_canonical_matches);
TargetSiteLocations = ndSparse.build([length(DoesProbeBindSite2),size(DoesProbeBindSite2,2),size(DoesProbeBindSite2,3),2],0);%P T S L
ProbeTargetSiteAlignment_NumMatchesAndMisMatches = ndSparse.build([length(DoesProbeBindSite2),size(DoesProbeBindSite2,2),size(DoesProbeBindSite2,3),2],0);%P T S 
ProbeTargetSiteAlignment_DanglingEnds = ndSparse.build([length(DoesProbeBindSite2),size(DoesProbeBindSite2,2),size(DoesProbeBindSite2,3),2],0);%P T S L
ProbeTargetSiteAlignment_GapLengthAndSegments = ndSparse.build([length(DoesProbeBindSite2),size(DoesProbeBindSite2,2),size(DoesProbeBindSite2,3),2],0);%P T S L
ProbeTargetSiteAlignment_Score = ndSparse.build([length(DoesProbeBindSite2),size(DoesProbeBindSite2,2),size(DoesProbeBindSite2,3)],0);%P T S L
ProbeTargetSiteAlignment_BitScore = ndSparse.build([length(DoesProbeBindSite2),size(DoesProbeBindSite2,2),size(DoesProbeBindSite2,3)],0);%P T S 
ProbeTargetSiteAlignment_Evalue = ndSparse.build([length(DoesProbeBindSite2),size(DoesProbeBindSite2,2),size(DoesProbeBindSite2,3)],0);%P T S 
filtMolProbesAtEvents = cell(1,size(DoesProbeBindSite2,2)); 
        for T=1:size(DoesProbeBindSite2,2) %ALL SITES/PROBES IN IT
        %Filter check if Probe in DPS2 for target or in
            filtMolProbesAtEvents{T} = arrayfun(@(S) MolProbesAtEvents{T}{S}(DoesProbeBindSite2(MolProbesAtEvents{T}{S},T,S)==1),1:length(MolProbesAtEvents{T}),'Un',0);
            Sz = cell2mat(arrayfun(@(S) S*ones(1,length(Mol_ProbesAtEventsID{T}{S})),1:length(MolProbesAtEvents{T}),'Un',0));
            Pz = cell2mat(arrayfun(@(S) gene_table.ProbeNum(Mol_ProbesAtEventsID{T}{S})',1:length(MolProbesAtEvents{T}),'Un',0));
            Az = cell2mat(arrayfun(@(S) gene_table.Ax(Mol_ProbesAtEventsID{T}{S})',1:length(MolProbesAtEvents{T}),'Un',0));
            Bz = cell2mat(arrayfun(@(S) gene_table.Bx(Mol_ProbesAtEventsID{T}{S})',1:length(MolProbesAtEvents{T}),'Un',0));
            FivePrimeDanglingEndZ = cell2mat(arrayfun(@(S) gene_table_alignment_five_prime_probe_dangling_end(Mol_ProbesAtEventsID{T}{S})',1:length(MolProbesAtEvents{T}),'Un',0));
            ThreePrimeDanglingEndZ = cell2mat(arrayfun(@(S) gene_table_alignment_three_prime_probe_dangling_end(Mol_ProbesAtEventsID{T}{S})',1:length(MolProbesAtEvents{T}),'Un',0));
            GapLengthsZ = cell2mat(arrayfun(@(S) gene_table_alignment_n_total_gap_length(Mol_ProbesAtEventsID{T}{S})',1:length(MolProbesAtEvents{T}),'Un',0));
            GapSegmentsZ = cell2mat(arrayfun(@(S) gene_table_alignment_n_gap_segments(Mol_ProbesAtEventsID{T}{S})',1:length(MolProbesAtEvents{T}),'Un',0));
            NumMisMatchesZ = cell2mat(arrayfun(@(S) gene_table_alignment_n_subsitutions(Mol_ProbesAtEventsID{T}{S})',1:length(MolProbesAtEvents{T}),'Un',0));
            NumMatchesZ = cell2mat(arrayfun(@(S) gene_table_alignment_n_canonical_matches(Mol_ProbesAtEventsID{T}{S})',1:length(MolProbesAtEvents{T}),'Un',0));
            ScoreZ = cell2mat(arrayfun(@(S) gene_table.Score(Mol_ProbesAtEventsID{T}{S})',1:length(MolProbesAtEvents{T}),'Un',0));
            EvalueZ = cell2mat(arrayfun(@(S) gene_table.Expect(Mol_ProbesAtEventsID{T}{S})',1:length(MolProbesAtEvents{T}),'Un',0));
            BitScoreZ = BitScore_Function(Lambda,Kappa,ScoreZ);
            if (~isempty(Pz))
                vecA = sub2ind(size(TargetSiteLocations),Pz,T*ones(1,length(Sz)),Sz,ones(1,length(Sz)));
                vecB = sub2ind(size(TargetSiteLocations),Pz,T*ones(1,length(Sz)),Sz,2*ones(1,length(Sz)));
                TargetSiteLocations(vecA) = Az;
                TargetSiteLocations(vecB) = Bz;
                ProbeTargetSiteAlignment_DanglingEnds(vecA) = FivePrimeDanglingEndZ;
                ProbeTargetSiteAlignment_DanglingEnds(vecB) = ThreePrimeDanglingEndZ;
                ProbeTargetSiteAlignment_GapLengthAndSegments(vecA) = GapLengthsZ;
                ProbeTargetSiteAlignment_GapLengthAndSegments(vecB) = GapSegmentsZ;
                ProbeTargetSiteAlignment_NumMatchesAndMisMatches(vecA) = NumMatchesZ;
                ProbeTargetSiteAlignment_NumMatchesAndMisMatches(vecB) = NumMisMatchesZ;
                vecC = sub2ind(size(DoesProbeBindSite2),Pz,T*ones(1,length(Sz)),Sz);
                ProbeTargetSiteAlignment_Score(vecC) = ScoreZ;
                ProbeTargetSiteAlignment_BitScore(vecC) = BitScoreZ;
                ProbeTargetSiteAlignment_Evalue(vecC) = EvalueZ;
            end
        end  
end



