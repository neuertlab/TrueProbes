function T = A_generateProbeOrderTable(ProbeDesignResults3)
MaxProbesPerTable = 48;
outGroup = [1 4 8 9 11];
outGroup1 = [1 1 8 9 11];
VS_Abreviation = {'TS','SL','SL','SL','SL','SL','SL',...
                  'OS','PS','MF','MF'};
V_Cell = @(R) 4/3*pi*(R^3)/10^15;Rcell = 10;%um to L 
Noff_Pset_SpecificitySorted = cell(1,length(outGroup));
Noff_SortedF = cell(1,length(outGroup));
Noff_SortedCumulativeF = cell(1,length(outGroup));
Coff_Pset_SpecificitySorted = cell(1,length(outGroup2));
Coff_SortedF = cell(1,length(outGroup2));
Coff_SortedCumulativeF = cell(1,length(outGroup2));
RM3 = cell(1,length(LegMade));
RM3_P = cell(1,length(LegMade));
SpecSorted = cell(1,length(LegMade));
G_Func = @(F) gradient(log(F));
for v = 1:length(outGroup)
    Pset_SpecificitySorted = ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes_ZigZagSelectionSorted;
    Pset_SpecificitySorted(ismember(Pset_SpecificitySorted,ProbeDesignResults1.SoftwareResults(outGroup1(v)).ProbesWithRibosomalHits))=[];
    Noff_Pset_SpecificitySorted{v} = Pset_SpecificitySorted;
    Noff_SortedF{v} = arrayfun(@(n) length([ProbeDesignResults3.SoftwareResults(outGroup1(v)).Tvec_RNA{Pset_SpecificitySorted(n)}]),1:length(Pset_SpecificitySorted));
    Noff_SortedCumulativeF{v} = arrayfun(@(n) length([ProbeDesignResults3.SoftwareResults(outGroup1(v)).Tvec_RNA{Pset_SpecificitySorted(1:n)}]),1:length(Pset_SpecificitySorted));
end 
for v = 1:length(outGroup2)
    cTargetSites_Bound = ProbeDesignResults4.SoftwareResults(outGroup2(v)).ModelMetrics.cTargetSites_Bound;
    OFF_IDs = ProbeDesignResults1.SoftwareResults(outGroup2(v)).OFF_IDs;
    DoesProbeBindSite = ProbeDesignResults1.SoftwareResults(outGroup2(v)).DoesProbeBindSite;
    Js = @(x) find(sum(squeeze(sum(DoesProbeBindSite(x,:,:),1)),2)>0);
    Pset0 = ProbeDesignResults3.SoftwareResults(outGroup2(v)).Designed_Probes;  
    Sx = unique(cell2mat(arrayfun(@(x) find(sum(DoesProbeBindSite(Pset0,x,:),1)>0)',Js(Pset0),'Un',0)));
    TYs = squeeze(cTargetSites_Bound(:,OFF_IDs,Sx,m,t,d,c));%pts
    Py = sum(sum(TYs,2),3); 
    [~,Pset_sorted1] = sort(full(Py),'ascend');
    Pset_SpecificitySorted = Pset0(Pset_sorted1);
    Pset_SpecificitySorted(ismember(Pset_SpecificitySorted,ProbeDesignResults1.SoftwareResults(outGroup21(v)).ProbesWithRibosomalHits))=[]; 
    Pset_sorted = arrayfun(@(x) find(Pset0==x),Pset_SpecificitySorted);
    Coff_Pset_SpecificitySorted{v} = Pset_SpecificitySorted;
    ProbeDesignResults4.SoftwareResults(outGroup2(v)).CumulativeModelMetricsOG.Sorted_Coff(1:length(Pset_sorted),c,m,d,t) = ...
    arrayfun(@(n) full(sum(sum(sum(TYs(Pset_sorted(1:n),:,:),1,'omitnan'),2,'omitnan'),3,'omitnan')),1:length(Pset_sorted))*(V_Cell(Rcell)*6.022*10^23); 
    Coff_SortedF{v} = arrayfun(@(n) full(sum(sum(sum(TYs(Pset_sorted(n),:,:),1,'omitnan'),2,'omitnan'),3,'omitnan')),1:length(Pset_sorted))*(V_Cell(Rcell)*6.022*10^23);     
    Coff_SortedCumulativeF{v} = arrayfun(@(n) full(sum(sum(sum(TYs(Pset_sorted(1:n),:,:),1,'omitnan'),2,'omitnan'),3,'omitnan')),1:length(Pset_sorted))*(V_Cell(Rcell)*6.022*10^23);  
end
DecisionFunc = Noff_SortedCumulativeF;
for v = 1:length(outGroup)
    RM1 = max(find(islocalmin(G_Func(DecisionFunc{v}))))+2;
    RM2 = max(find(diff(diff(log(DecisionFunc{v})))>0))+2;
    RM3{outGroup(v)} = max(RM1,RM2):length(DecisionFunc{v});
    RM3_P{outGroup(v)} =  Noff_Pset_SpecificitySorted{v}(max(RM1,RM2):length(DecisionFunc{v}));
    SpecSorted{outGroup(v)} = Noff_Pset_SpecificitySorted{v};
end    
DecisionFunc = Coff_SortedCumulativeF;
for v = 1:length(outGroup2)
    RM1 = max(find(islocalmin(G_Func(DecisionFunc{v}))))+2;
    RM2 = max(find(diff(diff(log(DecisionFunc{v})))>0))+2;
    RM3{outGroup2(v)} = max(RM1,RM2):length(DecisionFunc{v});
    RM3_P{outGroup2(v)} = Coff_Pset_SpecificitySorted{v}(max(RM1,RM2):length(DecisionFunc{v})); 
    SpecSorted{outGroup2(v)} = Coff_Pset_SpecificitySorted{v};
end
for v = 1:length(outGroup)
    Pset0 = ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes;
    ProbeDesignResults2.SoftwareResults(outGroup(v)).Designed_ProbesOG = ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes;
    ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_ProbesOG = ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes;
    ProbeDesignResults4.SoftwareResults(outGroup(v)).Designed_ProbesOG = ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes;
    Pset0(ismember(Pset0,RM3_P{outGroup(v)}))=[];
    S1 = @(x) ProbeDesignResults3.SoftwareResults(outGroup1(v)).TPvec_RNA{x};
    S2 = @(x) ProbeDesignResults3.SoftwareResults(outGroup1(v)).TSvec_RNA{x};
    S3 = @(x) ProbeDesignResults3.SoftwareResults(outGroup1(v)).Nvec_RNAmulti(x);
% Expression T*ones(1,length(TPvec_RNA(x)) NonDNA_IDs
%     TP = cell2mat(ProbeDesignResults3.SoftwareResults(outGroup1(v)).TPvec_RNA);
%         try
%     TS = cell2mat(ProbeDesignResults3.SoftwareResults(outGroup1(v)).TSvec_RNA);
%         catch
%     TS = cell2mat(ProbeDesignResults3.SoftwareResults(outGroup1(v)).TSvec_RNA'); 
        Gi = @(Y,x) sum(arrayfun(@(y) double(sum(ismember(Y,ProbeDesignResults3.SoftwareResults(outGroup1(v)).TPvec_RNA{x}(S2(x)==y)))>0), unique(S2(x))));   
        G0 = zeros(1,length(ProbeDesignResults1.SoftwareResults(outGroup1(v)).OFF_RNAIDs));
        for i = 1:length(G0)
            G0(i) = Gi(Pset0,i);
        end 
        while (max(G0)>nTrimLim)
            if (max(G0)>nTrimLim)
                T1 = find(G0==max(G0));
                try
                TP1 = cell2mat(arrayfun(@(x) S1(x),T1,'Un',0));
                catch
                TP1 = cell2mat(arrayfun(@(x) S1(x)',T1,'Un',0));    
                end
                NS = arrayfun(@(x) sum(double(TP1==x)),Pset0);
                RM5 = find(NS==max(NS));
                RM = Pset0(RM5(find(S3(Pset0(RM5))==max(S3(Pset0(RM5))),1)));
                Pset0(Pset0==RM) = [];
                for i = 1:length(G0)
                    G0(i) = Gi(Pset0,i);
                end            
            end
        end
ProbeDesignResults2.SoftwareResults(outGroup(v)).Designed_Probes = Pset0;
ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes = Pset0;
ProbeDesignResults4.SoftwareResults(outGroup(v)).Designed_Probes = Pset0;
end
for v = 1:length(outGroup3)
    PsetOG = ProbeDesignResults2.SoftwareResults(outGroup3(v)).Designed_ProbesOG;
    PsetT = ProbeDesignResults2.SoftwareResults(outGroup3(v)).Designed_Probes;
    PsetT(ismember(PsetT,ProbeDesignResults1.SoftwareResults(outGroup31(v)).ProbesWithRibosomalHits))=[];
    RMT = setdiff(PsetOG,PsetT);
    ProbeDesignResults2.SoftwareResults(outGroup3(v)).Designed_Probes_ZigZagSelectionSortedOG = ProbeDesignResults3.SoftwareResults(outGroup3(v)).Designed_Probes_ZigZagSelectionSorted;
    ProbeDesignResults3.SoftwareResults(outGroup3(v)).Designed_Probes_ZigZagSelectionSortedOG = ProbeDesignResults3.SoftwareResults(outGroup3(v)).Designed_Probes_ZigZagSelectionSorted;
    ProbeDesignResults4.SoftwareResults(outGroup3(v)).Designed_Probes_ZigZagSelectionSortedOG = ProbeDesignResults3.SoftwareResults(outGroup3(v)).Designed_Probes_ZigZagSelectionSorted;
    Pset_TrimmedSpecificitySorted = ProbeDesignResults3.SoftwareResults(outGroup3(v)).Designed_Probes_ZigZagSelectionSorted;
    Pset_TrimmedSpecificitySorted(ismember(Pset_TrimmedSpecificitySorted,ProbeDesignResults1.SoftwareResults(outGroup31(v)).ProbesWithRibosomalHits))=[];
    Pset_TrimmedSpecificitySorted(ismember(Pset_TrimmedSpecificitySorted,RMT))=[];
    ProbeDesignResults2.SoftwareResults(outGroup3(v)).Designed_Probes_ZigZagSelectionSorted = Pset_TrimmedSpecificitySorted;
    ProbeDesignResults3.SoftwareResults(outGroup3(v)).Designed_Probes_ZigZagSelectionSorted = Pset_TrimmedSpecificitySorted;
    ProbeDesignResults4.SoftwareResults(outGroup3(v)).Designed_Probes_ZigZagSelectionSorted = Pset_TrimmedSpecificitySorted; 
end

N_Software = length(outGroup);
orgID = lower(ProbeDesignResults3.Organism(1));
Name = extractBetween(ProbeDesignResults3.Gene,'(',')');
RefName = ProbeDesignResults3.GeneTarget;
Isoform = strtrim(extractBetween(ProbeDesignResults1.Description,'transcript variant',','));
if (~isempty(Isoform))
    IS_ID = strcat('IS',Isoform);
else
    IS_ID = 'IS1';
end
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
            ci = ci+1;
        end
    catch
    end
end
T = table(probe_name',probe_seqs','VariableNames',["probe_name","probe_sequence"]);
outPutFile1 = strcat(orgID,Name,'-',IS_ID,'-',RefName);
filename = strcat(outPutFile1,'-ProbeTable.xlsx');
writetable(T,filename{:},'Sheet',1,'Range','A1'); 
%RNA plus/plus (probe) vs plus/minus  (seqrcomplement)
    %sequence complement plus/plus
end