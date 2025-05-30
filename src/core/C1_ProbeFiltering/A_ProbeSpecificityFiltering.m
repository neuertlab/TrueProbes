function [ProbeDesignResults1,ProbeDesignResults2,ProbeDesignResults3,ProbeDesignResults4] = A_ProbeSpecificityFiltering(vs,ProbeDesignResults1,ProbeDesignResults2,ProbeDesignResults3,ProbeDesignResults4)
V_Cell = @(R) 4/3*pi*(R^3)/10^15;%um to L 
LegMade = vs{1};
Rcell = vs{2};
outGroup = vs{3};
outGroup1 = vs{4};
outGroup2 = vs{5};
outGroup21 = vs{6};
outGroup3 = vs{7};
outGroup31 = vs{8};
CI = vs{9};
stateInput = vs{10};
nTrimLim = vs{11};

c = stateInput(1);
m = stateInput(2);
d = stateInput(3);
t = stateInput(4);
Has_Repeated_probes = zeros(1,length(outGroup));
%fix for repeated probes 
for v = 1:length(outGroup)
    if (~isempty(ProbeDesignResults3.SoftwareResults(outGroup(v)).settings)) 
        probes = ProbeDesignResults3.SoftwareResults(outGroup(v)).probes;  
        designed_probes = ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes; 
        designed_probes_sorted = ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes_ZigZagSelectionSorted;
        probePos = [probes{designed_probes,3}];
        %detected repeats need to be removed.
        if (length(unique(probePos))<length(designed_probes))
           [~,ia,~] = unique(probePos); 
           NonrepeatProbes = sort(designed_probes(ia)); 
           ProbeDesignResults2.SoftwareResults(outGroup(v)).Designed_Probes = NonrepeatProbes;
           ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes = NonrepeatProbes;
           ProbeDesignResults4.SoftwareResults(outGroup(v)).Designed_Probes = NonrepeatProbes;
           designed_Nonrepeat_probes_sorted = designed_probes_sorted(ismember(designed_probes_sorted,NonrepeatProbes));
           ProbeDesignResults2.SoftwareResults(outGroup(v)).Designed_Probes_ZigZagSelectionSorted = designed_Nonrepeat_probes_sorted;  
           ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes_ZigZagSelectionSorted = designed_Nonrepeat_probes_sorted;  
           ProbeDesignResults4.SoftwareResults(outGroup(v)).Designed_Probes_ZigZagSelectionSorted = designed_Nonrepeat_probes_sorted; 
           NonrepeatProbes = [];
           designed_Nonrepeat_probes_sorted = [];
           Has_Repeated_probes(v) = 1;
        end
    end
end


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
    if (~isempty(ProbeDesignResults3.SoftwareResults(outGroup(v)).settings))
        Pset_SpecificitySorted = ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes_ZigZagSelectionSorted;
        Pset_SpecificitySorted(ismember(Pset_SpecificitySorted,ProbeDesignResults1.SoftwareResults(outGroup1(v)).ProbesWithRibosomalHits))=[];
        DNA_IDs = ProbeDesignResults1.SoftwareResults(outGroup(v)).DNA_IDs;
        NonDNA_IDs = ProbeDesignResults1.SoftwareResults(outGroup(v)).NonDNA_IDs;
        ExpressionMatrix = ProbeDesignResults3.SoftwareResults(outGroup1(v)).ExpressionMatrix;
        NumNonUniformConditions = size(ExpressionMatrix,2);
        ExpressionMatrix(DNA_IDs,NumNonUniformConditions+1) = 2;
        ExpressionMatrix(NonDNA_IDs,NumNonUniformConditions+1) = 100;
        ExpressionMatrix = ExpressionMatrix';  
        Noff_Pset_SpecificitySorted{v} = Pset_SpecificitySorted;       
        Noff_SortedF{v} = arrayfun(@(n) sum(ExpressionMatrix(CI,[ProbeDesignResults3.SoftwareResults(outGroup1(v)).Tvec_RNA{Pset_SpecificitySorted(n)}])),1:length(Pset_SpecificitySorted));
        Noff_SortedCumulativeF{v} = arrayfun(@(n) sum(ExpressionMatrix(CI,[ProbeDesignResults3.SoftwareResults(outGroup1(v)).Tvec_RNA{Pset_SpecificitySorted(1:n)}])),1:length(Pset_SpecificitySorted));
    end
end 


for v = 1:length(outGroup2)
    if (~isempty(ProbeDesignResults3.SoftwareResults(outGroup2(v)).settings))
        OFF_IDs = ProbeDesignResults1.SoftwareResults(outGroup2(v)).OFF_IDs;
        DoesProbeBindSite = ProbeDesignResults1.SoftwareResults(outGroup2(v)).DoesProbeBindSite;
        Js = @(x) find(sum(squeeze(sum(DoesProbeBindSite(x,:,:),1)),2)>0);
        Pset0 = ProbeDesignResults3.SoftwareResults(outGroup2(v)).Designed_Probes;
        Sx = unique(cell2mat(arrayfun(@(x) find(sum(DoesProbeBindSite(Pset0,x,:),1)>0)',Js(Pset0),'Un',0)));
        if (CI==97)
            cTargetSites_Bound = ProbeDesignResults4.SoftwareResults(outGroup2(v)).ModelMetrics.cTargetSites_Bound;
        else 
            cTargetSites_Bound = ProbeDesignResults4.SoftwareResults(outGroup2(v)).c_TargetSites_Bound; 
        end
        TYs = squeeze(cTargetSites_Bound(1:length(Pset0),OFF_IDs,Sx,m,t,d,c));%pts 
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
end
DecisionFunc = Noff_SortedCumulativeF;
for v = 1:length(outGroup)
    if (~isempty(ProbeDesignResults3.SoftwareResults(outGroup(v)).settings))
        RM1 = max(find(islocalmin(G_Func(DecisionFunc{v}))))+2;
        RM2 = max(find(diff(diff(log(DecisionFunc{v})))>0))+2;
        RM3{outGroup(v)} = max(RM1,RM2):length(DecisionFunc{v});
        RM3_P{outGroup(v)} =  Noff_Pset_SpecificitySorted{v}(max(RM1,RM2):length(DecisionFunc{v}));
        SpecSorted{outGroup(v)} = Noff_Pset_SpecificitySorted{v};
    end
end    
DecisionFunc = Coff_SortedCumulativeF;
for v = 1:length(outGroup2)
    if (~isempty(ProbeDesignResults3.SoftwareResults(outGroup2(v)).settings))
        RM1 = max(find(islocalmin(G_Func(DecisionFunc{v}))))+2;
        RM2 = max(find(diff(diff(log(DecisionFunc{v})))>0))+2;
        RM3{outGroup2(v)} = max(RM1,RM2):length(DecisionFunc{v});
        RM3_P{outGroup2(v)} = Coff_Pset_SpecificitySorted{v}(max(RM1,RM2):length(DecisionFunc{v})); 
        SpecSorted{outGroup2(v)} = Coff_Pset_SpecificitySorted{v};
    end
end
for v = 1:length(outGroup)
    if (~isempty(ProbeDesignResults3.SoftwareResults(outGroup(v)).settings))
        DNA_IDs = ProbeDesignResults1.SoftwareResults(outGroup(v)).DNA_IDs;
        NonDNA_IDs = ProbeDesignResults1.SoftwareResults(outGroup(v)).NonDNA_IDs;
        ExpressionMatrix = ProbeDesignResults3.SoftwareResults(outGroup1(v)).ExpressionMatrix;
        NumNonUniformConditions = size(ExpressionMatrix,2);
        ExpressionMatrix(DNA_IDs,NumNonUniformConditions+1) = 2;
        ExpressionMatrix(NonDNA_IDs,NumNonUniformConditions+1) = 100;
        ExpressionMatrix = ExpressionMatrix';
        OFF_RNAIDs = ProbeDesignResults1.SoftwareResults(outGroup1(v)).OFF_RNAIDs;
        Pset0 = ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes;
        ProbeDesignResults2.SoftwareResults(outGroup(v)).Designed_ProbesOG = ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes;
        ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_ProbesOG = ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes;
        ProbeDesignResults4.SoftwareResults(outGroup(v)).Designed_ProbesOG = ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes;
        Pset0(ismember(Pset0,RM3_P{outGroup(v)}))=[];
        S1 = @(x) ProbeDesignResults3.SoftwareResults(outGroup1(v)).TPvec_RNA{x};
        S2 = @(x) ProbeDesignResults3.SoftwareResults(outGroup1(v)).TSvec_RNA{x};
        S3 = @(x) ProbeDesignResults3.SoftwareResults(outGroup1(v)).Nvec_RNAmulti(x);
        Gi = @(Y,x) double(ExpressionMatrix(CI,OFF_RNAIDs(x))>0)*sum(arrayfun(@(y) double(sum(ismember(Y,ProbeDesignResults3.SoftwareResults(outGroup1(v)).TPvec_RNA{x}(S2(x)==y)))>0), unique(S2(x))));   %probe binds target
        G0 = zeros(1,length(OFF_RNAIDs)); 
        for i = 1:length(G0)
            G0(i) = Gi(Pset0,i);
        end
        while (max(G0)>nTrimLim)
            if (max(G0)>nTrimLim)
                T1 = find(G0==max(G0));
                TP1 = cell2mat(arrayfun(@(x) S1(x),T1,'Un',0));
                EP = cell2mat(arrayfun(@(x) ExpressionMatrix(CI,OFF_RNAIDs(x))*ones(1,length(S1(x))),T1,'Un',0));
                NS = arrayfun(@(x) sum(EP.*double(TP1==x)),Pset0);
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
end
for v = 1:length(outGroup3)
    if (~isempty(ProbeDesignResults3.SoftwareResults(outGroup3(v)).settings))
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
end
end