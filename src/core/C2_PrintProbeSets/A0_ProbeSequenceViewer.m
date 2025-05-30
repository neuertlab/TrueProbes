seqalignviewer


%seqviewer  ARF4 sequence
GeneSeqs = fastaread('OtherSoftwareProbeDesign/Oligostan/Targetseqs.fa');
dict = containers.Map({GeneSeqs.Header},arrayfun(@(x) strcat('ENSG',num2str(x)),1:length(GeneSeqs),'Un',0));
seqTarget = GeneSeqs(str2double(extractAfter(dict(ProbeDesignResults1.GeneTarget),'ENSG'))).Sequence;
ci = 1;
for n = 1:N_Software
    try
        probes = ProbeDesignResults3.SoftwareResults(outGroup(n)).probes;
        N_FullProbes = length(ProbeDesignResults3.SoftwareResults(outGroup(n)).Designed_ProbesOG);
        N_TrimProbes = length(ProbeDesignResults3.SoftwareResults(outGroup(n)).Designed_Probes);
        Pset0 = ProbeDesignResults3.SoftwareResults(outGroup(n)).Designed_ProbesOG;  
        Pset_SpecificitySorted = ProbeDesignResults3.SoftwareResults(outGroup(n)).Designed_Probes_ZigZagSelectionSorted;
        Pset_SpecificitySorted(ismember(Pset_SpecificitySorted,ProbeDesignResults3.SoftwareResults(outGroup1(n)).ProbesWithRibosomalHits))=[];  
        Pset_SpecificitySortedOG = ProbeDesignResults3.SoftwareResults(outGroup(n)).Designed_Probes_ZigZagSelectionSortedOG;
        Pset_SpecificitySortedOG(ismember(Pset_SpecificitySortedOG,ProbeDesignResults3.SoftwareResults(outGroup1(n)).ProbesWithRibosomalHits))=[]; 
        Pset_SpecificitySortedOG(ismember(Pset_SpecificitySortedOG,Pset_SpecificitySorted))=[];  
        Pset_SpecificitySorted = [Pset_SpecificitySorted Pset_SpecificitySortedOG];
        for p = 1:N_FullProbes
            if (p<=N_TrimProbes)
                num = num2str(p);  
            else
                num = strcat('(',num2str(p),')');
            end
            nam = strcat(orgID,Name,'-',IS_ID,'-',VS_Abreviation(outGroup(n)),'-',num);
            probe_name{ci} = nam{1};
            probe_seqs{ci} = seqrcomplement(probes{Pset_SpecificitySorted(p),2});%complement.
            rs = localalign(seqTarget,probes{Pset_SpecificitySorted(p),2});
            probe_locationFront{ci} = rs.Start(1);
            probe_locationEnd{ci} = rs.Stop(1);
            ci = ci+1;
        end
    catch
    end
end
T = table(probe_seqs',probe_name','VariableNames',["Sequence","Header"]);
T2 = table(probe_name',probe_seqs',probe_locationFront',probe_locationEnd','VariableNames',["Name","Sequence","StartLocation","EndLocation"]);
T3 = T2;
for nam = 1:5
vs = find(contains({T2.Name{:}},names{nam})); 
[~,idx] = sortrows(T3(vs,:),[3 4],'ascend');
T3(vs,:) = T3(vs(idx),:);
end

T3 = sorttable(T2,[3 4])


T1 = T;
nm = size(T1,1);

alignment_match = char();
alignment_binding = char();
alignment_match = repelem('-',6,length(ProbeDesignResults1.Gene)+length(seqTarget)+1);
alignment_binding = repelem('-',6,length(ProbeDesignResults1.Gene)+length(seqTarget)+1);
alignment_match(1,1:length(ProbeDesignResults1.Gene)+1) = strcat(ProbeDesignResults1.Gene,':');
alignment_binding(1,1:length(ProbeDesignResults1.Gene)+1) = strcat(ProbeDesignResults1.Gene,':');
alignment_match(1,length(ProbeDesignResults1.Gene)+2:length(ProbeDesignResults1.Gene)+length(seqTarget)+1) = seqTarget;
alignment_binding(1,length(ProbeDesignResults1.Gene)+2:length(ProbeDesignResults1.Gene)+length(seqTarget)+1) = seqTarget;
imBlock = zeros(5,length(seqTarget));
names = {'-TS','-SL','-OS','-PS','-MF'};
ck = 5;
for nam = 1:5
    vs = find(contains({T1.Header{:}},names{nam})); 
    alignment_match(nam+1,1:3) = names{nam};
    alignment_binding(nam+1,1:3) = names{nam};
    for v = vs
        rs = localalign(seqTarget,seqrcomplement(T1.Sequence{v}));
        alignment_match(nam+1,length(ProbeDesignResults1.Gene)+1+[rs.Start(1):rs.Stop(1)]) = seqrcomplement(T1.Sequence{v});
        alignment_binding(nam+1,length(ProbeDesignResults1.Gene)+1+[rs.Start(1):rs.Stop(1)]) = seqreverse(T1.Sequence{v}); 
        imBlock(nam,[rs.Start(1):rs.Stop(1)]) = ck;
        ck=ck+1;
        if (ck>10)
            ck = 1;
        end
        
    end
end
figure(2)
imagesc(imBlock);
ax = gca;
xlabel('Nucleotide base-pair position (bp)','FontWeight','bold','FontSize',20)
ylabel('Software','FontWeight','bold','FontSize',20)
title('ARF4 NM001660.4 RNA-FISH Probes Position on mRNA Transcript','FontWeight','bold','FontSize',20); 
ax.XTick = 0:50:1600;
ax.YTick = 1:5;
ax.YTickLabel = {'TS','SL','OS','PS','MF'};
ax.FontWeight = 'bold';
ax.FontSize = 20;
cz = parula;
cz(1,:) = [1 1 1];
colormap(cz);


%wrap if 
Nz = 50;
nwraps = ceil(length(seqTarget)/Nz);
Ri = mod(length(seqTarget),Nz);
spacel = 2*numel(num2str(length(seqTarget)))+1;
alignment_match_wrap = repelem('-',6*ceil((length(seqTarget))/Nz),spacel+length(ProbeDesignResults1.Gene)+1+Nz);
alignment_binding_wrap = repelem('-',6*ceil((length(seqTarget))/Nz),spacel+length(ProbeDesignResults1.Gene)+1+Nz);
for vi = 0:nwraps-1    
alignment_match_wrap(1+6*vi,spacel+[1:length(ProbeDesignResults1.Gene)+1]) = strcat(ProbeDesignResults1.Gene,':');
alignment_binding_wrap(1+6*vi,spacel+[1:length(ProbeDesignResults1.Gene)+1]) = strcat(ProbeDesignResults1.Gene,':');
alignment_match_wrap(6*vi+[2:6],1:spacel+length(ProbeDesignResults1.Gene)+1) = repelem(' ',5,spacel+length(ProbeDesignResults1.Gene)+1); 
alignment_binding_wrap(6*vi+[2:6],1:spacel+length(ProbeDesignResults1.Gene)+1) = repelem(' ',5,spacel+length(ProbeDesignResults1.Gene)+1); 
    for nam = 1:5 
        alignment_match_wrap(nam+1+6*vi,spacel+[1:3]) = names{nam};
        alignment_binding_wrap(nam+1+6*vi,spacel+[1:3]) = names{nam}; 
    end
    alignment_match_wrap(6*vi+[1:6],1:spacel) = repelem(' ',6,spacel); 
    alignment_binding_wrap(6*vi+[1:6],1:spacel) = repelem(' ',6,spacel); 
    if (vi<nwraps-1)
    alignment_match_wrap(6*vi+[1:6],1:length(strcat(num2str(1+Nz*vi),'-',num2str(Nz*(vi+1))))) = repelem(strcat(num2str(1+Nz*vi),'-',num2str(Nz*(vi+1))),6,1); 
    alignment_binding_wrap(6*vi+[1:6],1:length(strcat(num2str(1+Nz*vi),'-',num2str(Nz*(vi+1))))) = repelem(strcat(num2str(1+Nz*vi),'-',num2str(Nz*(vi+1))),6,1); 
    alignment_match_wrap(6*vi+[1:6],spacel+length(ProbeDesignResults1.Gene)+1+[1:Nz]) = alignment_match(1:6,length(ProbeDesignResults1.Gene)+1+[1+Nz*vi:Nz*(vi+1)]); 
    alignment_binding_wrap(6*vi+[1:6],spacel+length(ProbeDesignResults1.Gene)+1+[1:Nz]) = alignment_binding(1:6,length(ProbeDesignResults1.Gene)+1+[1+Nz*vi:Nz*(vi+1)]); 
    else
    alignment_match_wrap(6*vi+[1:6],1:length(strcat(num2str(1+Nz*vi),'-',num2str(length(seqTarget))))) = repelem(strcat(num2str(1+Nz*vi),'-',num2str(Nz*(vi+1))),6,1); 
    alignment_binding_wrap(6*vi+[1:6],1:length(strcat(num2str(1+Nz*vi),'-',num2str(length(seqTarget))))) = repelem(strcat(num2str(1+Nz*vi),'-',num2str(Nz*(vi+1))),6,1);
        alignment_match_wrap(6*vi+[1:6],spacel+length(ProbeDesignResults1.Gene)+1+[1:Ri]) = alignment_match(1:6,length(ProbeDesignResults1.Gene)+1+[1+Nz*vi:1+Nz*vi+Ri-1]); 
        alignment_binding_wrap(6*vi+[1:6],spacel+length(ProbeDesignResults1.Gene)+1+[1:Ri]) = alignment_binding(1:6,length(ProbeDesignResults1.Gene)+1+[1+Nz*vi:1+Nz*vi+Ri-1]); 
    end
end



seqalignviewer(alignment_match(:,length(ProbeDesignResults1.Gene)+2:end))
length(seqTarget)

filename_out = 'output2.txt';      
fid = fopen(filename_out, 'w' ); %// open file to writing
fprintf( fid, '%s\n', string(alignment_binding_wrap) ); %// print string to file
fclose( fid );



seqviewer(alignment_match(1,8:end))
distances = seqpdist(table2struct(T1),'method','jukes-cantor','indels','pair');
phylotree = seqlinkage(distances,'single',table2struct(T1));
seqsMultiAlign = multialign(table2struct(T1),phylotree,'ScoringMatrix','NUC44');%with reverse complement. order
seqalignviewer(seqsMultiAlign)
seqviewer({T.probe_sequence{find(contains([T.probe_name{:}],'-TS'))}})




