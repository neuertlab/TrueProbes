function chosenProbes = A_ZigZagProbeSelection_V5(probes,gene_table,settings,addSelfProb,packOptimal,Kon,Nvec_RNAmulti,Off_Score,Specificity_Score,Tvec_RNA,Svec_RNA,TPvec_RNA,TSvec_RNA,TPvec_logKOFF_RNA,TPvec_logKOFFdivON_RNA,TPvec_logKONdivOFF_RNA,ExpressionMatrix,DoesProbeBindSite,Kb)
% This Function slects probes to design in a zig-zag fashion. 
% Sorting probes first on based on the number of off-targets,
% followed by the difference between on-target, off-target,
% and secondary structure binding affinity.
% Probe are grouped by those with and those without any off-targets
% For probes without off-targets probes are picked to maximally cover 
% the target with the strongest on-target binding affinity.
% Remaining probes are picked according to the zig-zag selection, 
% Accounted for with cross-dimerization potential with all probes already designed.

GeneName = settings.GeneName;
GeneTarget = settings.transcript_IDs;
gene_table = sortrows(gene_table,[7 6],'ascend');
gene_table = gene_table(gene_table.Match>=settings.MinHomologySearchTargetSize,:);
MinusStrandedHits = find(contains(gene_table.Strand,'Minus'));
if (strcmp(settings.referenceType,'RefSeq'))
RNA_IDs_1 = find(contains(gene_table.Name,'NM_'));
RNA_IDs_2 = find(contains(gene_table.Name,'NR_'));
RNA_IDs_3 = find(contains(gene_table.Name,'XM_'));
RNA_IDs_4 = find(contains(gene_table.Name,'XR_'));
contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);

else

end
RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);
gene_table = gene_table(setdiff(1:size(gene_table,1),RNA_MissedFilteredHits),:);
gene_table.Ax = min(gene_table.SubjectIndices,[],2);
gene_table.Bx = max(gene_table.SubjectIndices,[],2);
gene_table = sortrows(gene_table,[7 13],'ascend');
Names = unique(gene_table.Name);
Names = convertCharsToStrings(Names);
uniNames = extractBefore(Names,'.');


DNA_IDs_1 = find(contains(uniNames,'NC_'));%IDs
DNA_IDs_2 = find(contains(uniNames,'NT_'));%IDs
DNA_IDs_3 = find(contains(uniNames,'NW_'));%IDs
NonDNA_IDs_1 = find(~contains(uniNames,'NC_'));%IDs
NonDNA_IDs_2 = find(~contains(uniNames,'NT_'));%IDs
NonDNA_IDs_3 = find(~contains(uniNames,'NW_'));%


DNA_IDs =union(union(DNA_IDs_1,DNA_IDs_2),DNA_IDs_3).';
NonDNA_IDs = intersect(intersect(NonDNA_IDs_1,NonDNA_IDs_2),NonDNA_IDs_3).';
ON_RNAIDs = find(strcmp(uniNames,extractBefore(GeneTarget,'.')));
OFF_RNAIDs = setdiff(NonDNA_IDs,ON_RNAIDs);
ON_RNAIDs_Isos = find(contains(Names,GeneName));
Desired_Isoforms = find(contains(uniNames,extractBefore(GeneTarget,'.')));
UnDesired_Isoforms = setdiff(ON_RNAIDs_Isos,Desired_Isoforms);


OFF_RNAIDs_minusIsos = setdiff(OFF_RNAIDs,UnDesired_Isoforms);
ON_RNAIDs = ON_RNAIDs_Isos;
OFF_RNAIDs = OFF_RNAIDs_minusIsos;



% NumRibosomalHits = zeros(1,size(probes,1)); 
Iz = find(endsWith(gene_table.Name,', ribosomal RNA'));
ProbesWithRibosomalHits = unique(gene_table.ProbeNum(Iz));
% RiboHits = gene_table3.ProbeNum(Iz);
% NumRibosomalHits(ProbesWithRibosomalHits) = cell2mat(arrayfun(@(x)sum(ismember(RiboHits,x)),ProbesWithRibosomalHits,'Un',0));
%Allow Matrix
spacing = settings.ProbeSpacing;
ProbesWithoutRibosomalHits = setdiff(1:size(probes,1),ProbesWithRibosomalHits);
AllowableProbes = setdiff(1:size(probes,1),ProbesWithRibosomalHits);
AllowMatrix = zeros(length(AllowableProbes),length(AllowableProbes));
for u1 = 1:length(AllowableProbes)
    u = AllowableProbes(u1);
    for v1 = 1:length(AllowableProbes)
       v = AllowableProbes(v1);
       if (u~=v)
       AllowMatrix(u1,v1) = double(length(intersect([probes{u,3}-spacing:probes{u,3}+length(probes{u,2})-1+spacing],[probes{v,3}-spacing:probes{v,3}+length(probes{v,2})-1+spacing]))>spacing);   
       end
    end
end



%Threshold Specificity and Num Off-Targets
N_model = settings.N_model;

Lpmin = min(cell2mat(cellfun(@length,{probes{:,2}},'UniformOutput',false)));
Probes_WithNo_RNAOFF = AllowableProbes(Nvec_RNAmulti(AllowableProbes)==0);
Js = @(x) find(sum(squeeze(sum(DoesProbeBindSite(x,:,:),1)),2)>0);
Js_Sites = @(x) find(sum(sum(DoesProbeBindSite(x,Js(x),:),1),2)>0);
Js_RNA = @(x) NonDNA_IDs(ismember(NonDNA_IDs,Js(x)));
Js_OFFRNA = @(x) OFF_RNAIDs(ismember(OFF_RNAIDs,Js(x)));
Js_DNA = @(x) DNA_IDs(ismember(DNA_IDs,Js(x)));
Sx =@(x,Z) arrayfun(@(y) find(squeeze(DoesProbeBindSite(x,y,:))==1)',Z,'Un',0);

   %% Intial Conditions for Selecting Probes
   chosenProbes = [];
   probe_poses = zeros(size(probes,1),1);
   for pos1 = 1:size(probes,1)
       probe_poses(pos1) = probes{pos1,3};
   end
   spacing_req = settings.ProbeSpacing;
   spacing_matrix = zeros(1,spacing_req*2+Lpmin*2 + max(probe_poses(:,1)));  
   NumOffTargetOptions = unique(Nvec_RNAmulti(AllowableProbes));
   Probes_WithNOFF_targets = arrayfun(@(i) AllowableProbes(Nvec_RNAmulti(AllowableProbes)==NumOffTargetOptions(i)),1:length(NumOffTargetOptions),'Un',0);
   
%     for v = 1:length(outGroup)
%     Pset0 = ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes;
%     Pset0(ismember(Pset0,RM3_P{outGroup(v)}))=[];
%     Js = @(x) find(sum(squeeze(sum(ProbeDesignResults1.SoftwareResults(outGroup(v)).DoesProbeBindSite(x,:,:),1)),2)>0);
%     Sx = unique(cell2mat(arrayfun(@(x) find(sum(ProbeDesignResults1.SoftwareResults(outGroup(v)).DoesProbeBindSite(Pset0,x,:),1)>0)',Js(Pset0),'Un',0)));
%     Ti = find(squeeze(sum(sum(ProbeDesignResults1.SoftwareResults(outGroup(v)).DoesProbeBindSite(Pset0,:,:),1),3))>0);
%     T0 = setdiff(Ti,[ProbeDesignResults1.SoftwareResults(outGroup(v)).Desired_Isoforms; ProbeDesignResults1.SoftwareResults(outGroup(v)).UnDesired_Isoforms]);
%     %EXPRESSion matrix 
%     % Expression T*ones(1,length(TPvec_RNA(x)) NonDNA_IDs
%     TP = cell2mat(ProbeDesignResults3.SoftwareResults(outGroup1(v)).TPvec_RNA);
%         try
%     TS = cell2mat(ProbeDesignResults3.SoftwareResults(outGroup1(v)).TSvec_RNA);
%         catch
%     TS = cell2mat(ProbeDesignResults3.SoftwareResults(outGroup1(v)).TSvec_RNA');
%         end
%     S1 = @(x) ProbeDesignResults3.SoftwareResults(outGroup1(v)).TPvec_RNA{x};
%     S2 = @(x) ProbeDesignResults3.SoftwareResults(outGroup1(v)).TSvec_RNA{x};
%     S3 = @(x) ProbeDesignResults3.SoftwareResults(outGroup1(v)).Nvec_RNAmulti(x);
%     Gi = @(Y,x) sum(arrayfun(@(y) double(sum(ismember(Y,ProbeDesignResults3.SoftwareResults(outGroup1(v)).TPvec_RNA{x}(S2(x)==y)))>0), unique(S2(x))));   
%     G0 = zeros(1,length(ProbeDesignResults1.SoftwareResults(outGroup1(v)).OFF_RNAIDs));
%         for i = 1:length(G0)
%             G0(i) = Gi(Pset0,i);
%         end 
%         while (max(G0)>nTrimLim)
%             if (max(G0)>nTrimLim)
%                 T1 = find(G0==max(G0));
%                 try
%                 TP1 = cell2mat(arrayfun(@(x) S1(x),T1,'Un',0));
%                 catch
%                 TP1 = cell2mat(arrayfun(@(x) S1(x)',T1,'Un',0));    
%                 end
%                 NS = arrayfun(@(x) sum(double(TP1==x)),Pset0);
%                 RM5 = find(NS==max(NS));
%                 RM = Pset0(RM5(find(S3(Pset0(RM5))==max(S3(Pset0(RM5))),1)));
%                 Pset0(Pset0==RM) = [];
%                 for i = 1:length(G0)
%                     G0(i) = Gi(Pset0,i);
%                 end            
%             end
%         end
%     end
%    

try
probe_Ks = zeros(1,size(probes,1));
Tvec_sub_pick = arrayfun(@(x) Tvec_RNA{x},1:size(probes,1),'Un',0);
Svec_sub_pick = arrayfun(@(x) Svec_RNA{x},1:size(probes,1),'Un',0);
TvecID_sub_pick = cellfun(@(x) find(ismember(OFF_RNAIDs,x)),Tvec_sub_pick,'Un',0);
TTvec_sub_pick = cellfun(@(x) cell2mat(arrayfun(@(y) repmat(OFF_RNAIDs(y),[1 length(TPvec_RNA{y})]),x,'Un',0)),TvecID_sub_pick,'Un',0);
TPvec_sub_pick = cellfun(@(x) [TPvec_RNA{x}],TvecID_sub_pick,'Un',0);
TSvec_sub_pick = cellfun(@(x)cell2mat(arrayfun(@(y)TSvec_RNA{y}',x,'Un',0)),TvecID_sub_pick,'Un',0);
ProbesAtDifSites = arrayfun(@(x)cell2mat(arrayfun(@(y)TPvec_sub_pick{x}(TSvec_sub_pick{x}(TTvec_sub_pick{x}==Tvec_sub_pick{x}(y))~=Svec_sub_pick{x}(y)),1:length(Tvec_sub_pick{x}),'Un',0)),...
                            1:size(probes,1),'Un',0);                        
ProbesAtSameSites = arrayfun(@(x)cell2mat(arrayfun(@(y)TPvec_sub_pick{x}(TSvec_sub_pick{x}(TTvec_sub_pick{x}==Tvec_sub_pick{x}(y))==Svec_sub_pick{x}(y)),1:length(Tvec_sub_pick{x}),'Un',0)),...
                            1:size(probes,1),'Un',0);
ProbesAtDifSites_AllowedBySpacing = arrayfun(@(x) ...    
                            ProbesAtDifSites{x}(sum(AllowMatrix(ismember(AllowableProbes,[x chosenProbes]),ismember(AllowableProbes,ProbesAtDifSites{x})),1)==0),...
                            1:size(probes,1),'Un',0);    
ProbesAtSameSites_AllowedBySpacing = arrayfun(@(x) ...    
                            ProbesAtSameSites{x}(sum(AllowMatrix(ismember(AllowableProbes,[x chosenProbes]),ismember(AllowableProbes,ProbesAtSameSites{x})),1)==0),...
                            1:size(probes,1),'Un',0);                        
% TPvec_logKOFF_RNA_sub_pick = cellfun(@(x) [TPvec_logKOFF_RNA{x}],TvecID_sub_pick,'Un',0);
% TPvec_logKOFFdivON_RNA_sub_pick = cellfun(@(x) [TPvec_logKOFFdivON_RNA{x}],TvecID_sub_pick,'Un',0);
% TPvec_logKONdivOFF_RNA_sub_pick = cellfun(@(x) [TPvec_logKONdivOFF_RNA{x}],TvecID_sub_pick,'Un',0);
Num_Multi_Dif = cellfun(@length,ProbesAtDifSites);    
Num_Multi_Same = cellfun(@(x)length(unique(x)),ProbesAtSameSites);  
Num_MultiSite_Dif = cellfun(@(x)length(unique(x)),ProbesAtDifSites_AllowedBySpacing);    
Num_MultiSite_Same = cellfun(@length,ProbesAtSameSites_AllowedBySpacing);    
Num_Excluded_Primary = arrayfun(@(x)sum(AllowMatrix(ismember(AllowableProbes,x),ismember(AllowableProbes,x))),1:size(probes,1));
catch
end
    %probes off-targets find compares to other probes binding same target
    %in Kon, Koff or Koff/Kon.
  %
  
  %GC ON high, gc off low,  avoid self and cross
    
%check for a probes off-target what is the affinity of another probe
%binding their Koff or specificity min and max vs if their is no other.
%Probe Stats 
%   1. Ks, Kon 
%   2. For Each OFF-Target
%       2a. Koff 
%       2b. Koff/Kon
%       2c. Target ID
%       2d. Target Site
%       2e. Target Expression Min, 
%       2f. Target Expression Max,
%       2g. Num of other Probes that bind the OFF-Target
%       2g. Num of other Probes that bind the OFF-Target and do not overlap with this probe
%       2g. Koff Min For other probes bind on OFF-target
%       2h. Koff Max For other probes bind on OFF-target
%       2g. Koff/Kon Min For other probes bind on OFF-target
%       2h. Koff/Kon Max For other probes bind on OFF-target
%       Singular Probe Kon, Num Across Probes that bind any of current
%       probes OFF-targets
% Target_MatchID_table = cell(1,size(probes,1));
% Target_Single_table = cell(1,size(probes,1));
% TargetSingleSummary_Table = cell(1,size(probes,1));
% 
% probe_KS = zeros(1,size(probes,1));
% probe_KDi = zeros(1,size(probes,1));
% for p = 1:size(probes,1)
%     try
%     [KS,KCD] = A_JH_GenerateSecondaryStructureInfo_V2(probes,p,settings);
%     probe_KS(p) = full(KS(p));
%     probe_KDi(p) = 2*full(KCD(p,p));
%     catch
%     end 
%     if (Nvec_RNAmulti(p)>0)
%         T_uni = Tvec_RNA{p};%T_uni = unique(TTvec_sub_pick{p}(TPvec_sub_pick{p}~=p));
%         P_other = arrayfun(@(x) TPvec_sub_pick{p}(double(TPvec_sub_pick{p}~=p).*double(TTvec_sub_pick{p}==T_uni(x))==1),1:length(T_uni),'Un',0);
        %S_other = arrayfun(@(x) TSvec_sub_pick{p}(double(TPvec_sub_pick{p}~=p).*double(TTvec_sub_pick{p}==T_uni(x))==1),1:length(T_uni),'Un',0);
%         Num_P_other = cellfun(@length,P_other);
%         Num_P_otherAllowed = cell2mat(arrayfun(@(x) sum(sum(AllowMatrix(ismember(AllowableProbes,p),ismember(AllowableProbes,P_other{x}))==0)),1:length(T_uni),'Un',0));                                                            
%         P_inc = arrayfun(@(x) TPvec_sub_pick{p}(double(TTvec_sub_pick{p}==T_uni(x))==1),1:length(T_uni),'Un',0);
%         S_inc = arrayfun(@(x) TSvec_sub_pick{p}(double(TTvec_sub_pick{p}==T_uni(x))==1),1:length(T_uni),'Un',0);
%         T_KOFF = arrayfun(@(x) diag(full(squeeze(Kb(P_inc{x},T_uni(x),S_inc{x}))))',1:length(T_uni),'Un',0);
%         T_KOFFdivON = arrayfun(@(x) diag(full(squeeze(Kb(P_inc{x},T_uni(x),S_inc{x}))))'./Kon(P_inc{x}),1:length(T_uni),'Un',0);
%         T_KOFFmin = cell2mat(cellfun(@(x) min(x),T_KOFF,'Un',0));
%         T_KOFFmax = cell2mat(cellfun(@(x) max(x),T_KOFF,'Un',0));
%         T_KOFFdivONmin = cell2mat(cellfun(@(x) min(x(:)),T_KOFFdivON,'Un',0));
%         T_KOFFdivONmax = cell2mat(cellfun(@(x) max(x(:)),T_KOFFdivON,'Un',0));
%         rangeI = T_KOFFmax-T_KOFFmin;
%         T_K = diag(full(squeeze(Kb(p,Tvec_RNA{p},cell2mat(Sx(p,Js_OFFRNA(p)))))));
%         T_KOFF_DIF = diag(full(squeeze(Kb(p,Tvec_RNA{p},cell2mat(Sx(p,Js_OFFRNA(p)))))))-T_KOFFmin';
%         T_KOFF_LogRANGE = log10(rangeI+1)';
%         T_KOFF_RANK = cell2mat(arrayfun(@(x) find(ismember(unique(T_KOFF{x}),T_K(x)))-1,1:length(T_uni),'Un',0))';       
%         rangeI2 = T_KOFFdivONmax-T_KOFFdivONmin;   
%         T_KOFFdivON_DIF = diag(full(squeeze(Kb(p,Tvec_RNA{p},cell2mat(Sx(p,Js_OFFRNA(p)))))))/Kon(p)-T_KOFFdivONmin';
%         T_KOFFdivON_LogRANGE = log10(rangeI2+1)';   
%         Target_MatchID_table{p} = [p*ones(length(Tvec_RNA{p}),1) Tvec_RNA{p}' Svec_RNA{p}' ... %Target ID/Target Site
%             ...  %percentage of cell line with zero, log10 sum across cell line, mean median fano 
%         Num_P_other' Num_P_otherAllowed' ... %Num P overlap / overlap allowed  
%         log10(diag(full(squeeze(Kb(p,Tvec_RNA{p},cell2mat(Sx(p,Js_OFFRNA(p)))))))) ... %Koff
%         log10(T_KOFFmin)' log10(T_KOFFmax)' T_KOFF_LogRANGE ... %Koff_min Koff_max Koff_range
%         log10(T_KOFF_DIF+1) T_KOFF_RANK... %log10(Koff-Koff_min+1) Koff_Rank
%         log10(diag(full(squeeze(Kb(p,Tvec_RNA{p},cell2mat(Sx(p,Js_OFFRNA(p)))))))/Kon(p)) ... %Koff/Kon 
%         log10(T_KOFFdivONmin)' log10(T_KOFFdivONmax)' T_KOFFdivON_LogRANGE ... %Koff/Kon_min Koff/Kon_max log10(Koff/Kon-Koff/Kon_min+1)
% %         log10(T_KOFFdivON_DIF+1)];%log10(Koff/Kon - Koff/Kon_min+1) KoffdivON_rank            
%     end
% end

%probes fewer off-target binding affinity, use instead

%Avoid List

%dont be strict on ranks as some ranks have 1 or 2 probes only.
if (~isempty(Probes_WithNo_RNAOFF))
    Intersection_NoRNAOFF=AllowMatrix(ismember(AllowableProbes,Probes_WithNo_RNAOFF),ismember(AllowableProbes,Probes_WithNo_RNAOFF));
    IntersectionMap1 = triu(Intersection_NoRNAOFF)-diag(ones(1,length(Probes_WithNo_RNAOFF)));
    IntersectinMap_2 = [sum(abs(IntersectionMap1),1) 1];
    IsLands = find(IntersectinMap_2==1);%1 is start of new group. 
    MacroGroups = arrayfun(@(x) IsLands(x):IsLands(x+1)-1,1:length(IsLands)-1,'Un',0);
    Range_Front = 1:min(Probes_WithNo_RNAOFF)-1;
    Range_MacroGroups_logKon = cell(1,length(MacroGroups));
    Range_MacroGroups = cell(1,length(MacroGroups));
    Range_InterGroups = cell(1,length(MacroGroups));
    Range_BetweenGroups = cell(1,length(MacroGroups)-1);
    N_MacroGroups = length(MacroGroups);
    for i = 1:length(MacroGroups)
        Range_MacroGroups_logKon{i} = log10(Kon(Probes_WithNo_RNAOFF(MacroGroups{i})));
        Range_MacroGroups{i} = Probes_WithNo_RNAOFF(MacroGroups{i}); 
        Range_InterGroups{i} = setdiff(min(Probes_WithNo_RNAOFF(MacroGroups{i})):max(Probes_WithNo_RNAOFF(MacroGroups{i})),Probes_WithNo_RNAOFF(MacroGroups{i}));
    end
    for i = 1:length(MacroGroups)-1
        Range_BetweenGroups{i} = max(Probes_WithNo_RNAOFF(MacroGroups{i}))+1:min(Probes_WithNo_RNAOFF(MacroGroups{i+1}))-1;
    end
    Range_End = max(Probes_WithNo_RNAOFF)+1:size(probes,1);
    %% Pick Probes with Zero OFF-targets
    InterGroupSize = cellfun(@length,Range_InterGroups);
    GroupsWithoutSplit = find(InterGroupSize==0);%Groups with No split inside
    GroupsWithSplit = find(InterGroupSize>0);%Groups with No split inside
    GroupsWithSinglePick = find(cellfun(@(x) all(diff(IntersectinMap_2(x))),MacroGroups));%pack one no-off probes from group
    GroupsWithMultiPicks = find(cellfun(@(x) all(diff(IntersectinMap_2(x))),MacroGroups)==0);%pack multiple no-off probes from group
    Groups_s1a = GroupsWithSinglePick;
    for vs = 1:length(Groups_s1a)
        vi = Groups_s1a(vs);
        sub_pick_Probes = Range_MacroGroups{vi};
%         [K_S,K_CD,~,~] = A_JH_GenerateSecondaryStructureInfo(probes,[chosenProbes sub_pick_Probes],settings);
        [K_S,~,~,~,~,~,~,...
        ~,K_CD,~,~,~,~,~,~,~,...
        ~,~] = A_JH_GenerateSecondaryStructureInfo_V2(probes,[chosenProbes sub_pick_Probes],settings); 
        K_S = sum(squeeze(K_S(:,:,N_model)),2,'omitnan');
        K_CD = sum(squeeze(K_CD(:,:,:,N_model)),3,'omitnan');
        sub_pick_Probes_Kon = Range_MacroGroups_logKon{vi}';
        sub_pick_Probes_Ks = log10(K_S(sub_pick_Probes)+1)';
        sub_pick_Probes_Kcd = log10(K_CD(sub_pick_Probes,chosenProbes)+1)';%sub x Nset
        sub_pick_Probes_Kons = sub_pick_Probes_Kon-sub_pick_Probes_Ks-sum(sub_pick_Probes_Kcd,1);
        if (N_MacroGroups>1)
            switch vi
                case 1
                    RangeBefore = Range_Front;
                    RangeAfter = Range_BetweenGroups{vi};
                case length(MacroGroups)
                    RangeBefore = Range_BetweenGroups{vi-1};
                    RangeAfter = Range_End;
                otherwise
                    RangeBefore = Range_BetweenGroups{vi-1};
                    RangeAfter = Range_BetweenGroups{vi};
            end
            vk = 2*double(length(RangeBefore)<Lpmin)+3*double(length(RangeAfter)<Lpmin);
        else
            vk = 0;
        end
        if (addSelfProb)
             sub_pick_Probes_Kdec = sub_pick_Probes_Kons;
        else
             sub_pick_Probes_Kdec = sub_pick_Probes_Kon;
        end
        switch vk
            case 0
                sub_pick_Probes_Consider = sub_pick_Probes(sub_pick_Probes_Kdec==max(sub_pick_Probes_Kdec));
            case 2%Before Range is small for picking (might want to change for packing)
                sub_pick_Probes_Consider = sub_pick_Probes(sub_pick_Probes_Kdec==max(sub_pick_Probes_Kdec));
            case 3%After Range is small for picking (might want to change for packing)
                sub_pick_Probes_Consider = sub_pick_Probes(sub_pick_Probes_Kdec==max(sub_pick_Probes_Kdec));
            case 5%Before&After Range is small for picking (might want to change for packing)
                sub_pick_Probes_Consider = sub_pick_Probes(sub_pick_Probes_Kdec==max(sub_pick_Probes_Kdec));
        end       
        probe_num = sub_pick_Probes_Consider(1);Lp = length(probes{probe_num,2});
        if (spacing_matrix(probes{probe_num,3}+spacing_req+Lp-1) == 0 && ismember(probe_num,AllowableProbes)==1)%If the space for the probe is open
            chosenProbes = [chosenProbes sub_pick_Probes_Consider(1)];
            spacing_matrix(probes{probe_num,3}:probes{probe_num,3} + spacing_req*2+Lp+ Lpmin-1) = 1; 
        end
    end     
    
    
    
    
    
   %Stage 1b [Groups with no splits, multi pick]
    Groups_s1b = intersect(GroupsWithoutSplit,GroupsWithMultiPicks);
    for vs = 1:length(Groups_s1b)  %just occurs in long stretch 
        vi = Groups_s1b(vs);    
        Pi = Range_MacroGroups{vi};
        sub_pick_Probes = Range_MacroGroups{vi};
        [K_S,K_CD,~,~] = A_JH_GenerateSecondaryStructureInfo_V2(probes,[chosenProbes sub_pick_Probes],settings);
        K_S = sum(squeeze(K_S(:,:,N_model)),2,'omitnan');
        K_CD = sum(squeeze(K_CD(:,:,:,N_model)),3,'omitnan');
        sub_pick_Probes_Kon = Range_MacroGroups_logKon{vi};
        sub_pick_Probes_Ks = log10(K_S(sub_pick_Probes)+1)';
        sub_pick_Probes_Kcd = log10(K_CD(sub_pick_Probes,chosenProbes)+1)';%sub x Nset
        sub_pick_Probes_Kons = sub_pick_Probes_Kon-sub_pick_Probes_Ks-sum(sub_pick_Probes_Kcd,1);
        if (addSelfProb)
             sub_pick_Probes_Kdec = sub_pick_Probes_Kons;
        else
             sub_pick_Probes_Kdec = sub_pick_Probes_Kon;
        end
        if (N_MacroGroups>1)
            switch vi
                case 1
                    RangeBefore = Range_Front;
                    RangeAfter = Range_BetweenGroups{vi};
                case length(MacroGroups)
                    RangeBefore = Range_BetweenGroups{vi-1};
                    RangeAfter = Range_End;
                otherwise
                    RangeBefore = Range_BetweenGroups{vi-1};
                    RangeAfter = Range_BetweenGroups{vi};
            end
            vk = 2*double(length(RangeBefore)<Lpmin)+3*double(length(RangeAfter)<Lpmin);
        else
            vk = 0;
        end
        finishGroup = 0;
        while (finishGroup==0)
             if (packOptimal)
                 Num_Excluded = arrayfun(@(x)sum(AllowMatrix(ismember(AllowableProbes,Pi),ismember(AllowableProbes,Pi(x)))),1:length(Pi));
                 Gap_Between = diff(Pi)-1;
                 Priority_Pick = Pi(Num_Excluded==min(Num_Excluded));
             else%Pick best, Kon, and remove overlap and repeat
                 [~,temp] = sort(sub_pick_Probes_Kdec(ismember(sub_pick_Probes,Pi)),'descend');
                 Priority_Pick = Pi(temp(1));
             end
        probe_num = Priority_Pick;Lp = length(probes{probe_num,2});
             if (spacing_matrix(probes{probe_num,3}+spacing_req+Lp-1) == 0 && ismember(probe_num,AllowableProbes)==1)%If the space for the probe is open
                 chosenProbes = [chosenProbes Priority_Pick];
                 spacing_matrix(probes{probe_num,3}:probes{probe_num,3} + spacing_req*2+Lp+ Lpmin-1) = 1; 
                 Pi_Overlap= Pi(AllowMatrix(ismember(AllowableProbes,Pi),ismember(AllowableProbes,Priority_Pick))==1);
                 Pi = setdiff(Pi,[Priority_Pick Pi_Overlap]);     
             end
             if (isempty(Pi))
                 finishGroup = 1; 
             end
        end      
    end
    %Stage 1c [Groups with splits, multi pick]
    Groups_s1c = intersect(GroupsWithSplit,GroupsWithMultiPicks);
    for vs = 1:length(Groups_s1c) 
        vi = Groups_s1c(vs);
        Pi = Range_MacroGroups{vi};            
        sub_pick_Probes = Range_MacroGroups{vi};
        [K_S,K_CD,~,~] = A_JH_GenerateSecondaryStructureInfo_V2(probes,[chosenProbes sub_pick_Probes],settings);
        K_S = sum(squeeze(K_S(:,:,N_model)),2,'omitnan');
        K_CD = sum(squeeze(K_CD(:,:,:,N_model)),3,'omitnan');
        sub_pick_Probes_Kon = Range_MacroGroups_logKon{vi};
        sub_pick_Probes_Ks = log10(K_S(sub_pick_Probes)+1)';
        sub_pick_Probes_Kcd = log10(K_CD(sub_pick_Probes,chosenProbes)+1)';
        sub_pick_Probes_Kons = sub_pick_Probes_Kon-sub_pick_Probes_Ks-sum(sub_pick_Probes_Kcd,1);
        if (addSelfProb)
             sub_pick_Probes_Kdec = sub_pick_Probes_Kons;
        else
             sub_pick_Probes_Kdec = sub_pick_Probes_Kon;
        end
        if (N_MacroGroups>1)
        switch vi
             case 1
                 RangeBefore = Range_Front;
                 RangeAfter = Range_BetweenGroups{vi};
             case length(MacroGroups)
                 RangeBefore = Range_BetweenGroups{vi-1};
                 RangeAfter = Range_End;
            otherwise
                 RangeBefore = Range_BetweenGroups{vi-1};
                 RangeAfter = Range_BetweenGroups{vi};
        end
        RangeInter = Range_InterGroups{vi};
        vk = 2*double(length(RangeBefore)<Lpmin)+3*double(length(RangeAfter)<Lpmin);
        else
        vk = 0;    
        end
        finishGroup = 0;
        while (finishGroup==0)
             if (packOptimal)
                 Num_Excluded = arrayfun(@(x)sum(AllowMatrix(ismember(AllowableProbes,Pi),ismember(AllowableProbes,Pi(x)))),1:length(Pi));
                 Gap_Between = diff(Pi)-1;
                 Priority_Pick = Pi(Num_Excluded==min(Num_Excluded));
                 if (length(Priority_Pick)>1)
                     [~,temp] = sort(sub_pick_Probes_Kdec(ismember(sub_pick_Probes,Priority_Pick)),'descend');
                     Priority_Pick = Priority_Pick(temp(1));
                 end  
             else%Pick best, Kon, and remove overlap and repeat
                 [~,temp] = sort(sub_pick_Probes_Kdec(ismember(sub_pick_Probes,Pi)),'descend');
                 Priority_Pick = Pi(temp(1));
             end
             probe_num = Priority_Pick;Lp = length(probes{probe_num,2});
             if (spacing_matrix(probes{probe_num,3}+spacing_req+Lp-1) == 0 && ismember(probe_num,AllowableProbes)==1)%If the space for the probe is open
                 chosenProbes = [chosenProbes Priority_Pick];
                 spacing_matrix(probes{probe_num,3}:probes{probe_num,3} + spacing_req*2+Lp+ Lpmin-1) = 1; 
                 Pi_Overlap= Pi(AllowMatrix(ismember(AllowableProbes,Pi),ismember(AllowableProbes,Priority_Pick))==1);
                 Pi = setdiff(Pi,[Priority_Pick Pi_Overlap]);     
             end
             if (isempty(Pi))
                 finishGroup = 1; 
             end
        end       
    end
end   


%expression and affinity discounting
  %% Stage 2 [Pick Probes With OFF-Targets]
for vsk = 1:length(NumOffTargetOptions)
    %Tally Off-target T and Sites to exclude from repeating? option
    if (NumOffTargetOptions(vsk)>0)
        sub_pick_Probes = Probes_WithNOFF_targets{vsk};
        sub_pick_Probes = setdiff(sub_pick_Probes,ProbesWithRibosomalHits);
        vis = [];
        for probe_num = sub_pick_Probes
            Lp = length(probes{probe_num,2});
            if (spacing_matrix(probes{probe_num,3}+spacing_req+Lp-1) == 0 && ismember(probe_num,AllowableProbes)==1)%If the space for the probe is open
                vis = [vis probe_num];
            end
        end
           sub_pick_Probes_Pass = vis;
           %sub_pick_Probes_Pass = sub_pick_Probes(sum(AllowMatrix(ismember(AllowableProbes,chosenProbes),ismember(AllowableProbes,sub_pick_Probes)),1)==0);
        if (~isempty(sub_pick_Probes_Pass))
            sub_pick_Probes_Kon = log10(Kon(sub_pick_Probes_Pass));
            [K_S,~,~,~,~,~,~,~,K_CD,~,~,~,~,~,~,~,~,~] = A_JH_GenerateSecondaryStructureInfo_V2(probes,[chosenProbes sub_pick_Probes_Pass],settings);
            K_S = sum(squeeze(K_S(:,:,N_model)),2,'omitnan');
            K_CD = sum(squeeze(K_CD(:,:,:,N_model)),3,'omitnan');
            sub_pick_Probes_Ks = log10(K_S(sub_pick_Probes_Pass)+1)';
            sub_pick_Probes_Kcd = log10(K_CD(sub_pick_Probes_Pass,chosenProbes)+1)';
            sub_pick_Probes_Kons = sub_pick_Probes_Kon-sub_pick_Probes_Ks-sum(sub_pick_Probes_Kcd,1);
            sub_pick_Probes_SpecificityScore = Specificity_Score(sub_pick_Probes_Pass);
            sub_pick_Probes_OffScore = Off_Score(sub_pick_Probes_Pass);
            sub_pick_Probes_SpecificityScoreS = sub_pick_Probes_SpecificityScore +sub_pick_Probes_Ks;
            Pi = sub_pick_Probes_Pass;
            finishOptions = 0;
            while (finishOptions==0)        
                Tvec_sub_pick = arrayfun(@(x) Tvec_RNA{x},Pi,'Un',0);
                Svec_sub_pick = arrayfun(@(x) Svec_RNA{x},Pi,'Un',0);
                TvecID_sub_pick = cellfun(@(x) find(ismember(OFF_RNAIDs,x)),Tvec_sub_pick,'Un',0);
                TTvec_sub_pick = cellfun(@(x) cell2mat(arrayfun(@(y) repmat(OFF_RNAIDs(y),[1 length(TPvec_RNA{y})]),x,'Un',0)),TvecID_sub_pick,'Un',0);
                TPvec_sub_pick = cellfun(@(x) [TPvec_RNA{x}],TvecID_sub_pick,'Un',0);
                TSvec_sub_pick = cellfun(@(x)cell2mat(arrayfun(@(y)TSvec_RNA{y}',x,'Un',0)),TvecID_sub_pick,'Un',0);   
%                 TPvec_logKOFF_RNA_sub_pick = cellfun(@(x) [TPvec_logKOFF_RNA{x}],TvecID_sub_pick,'Un',0);
%                 TPvec_logKOFFdivON_RNA_sub_pick = cellfun(@(x) [TPvec_logKOFFdivON_RNA{x}],TvecID_sub_pick,'Un',0);
%                 TPvec_logKONdivOFF_RNA_sub_pick = cellfun(@(x) [TPvec_logKONdivOFF_RNA{x}],TvecID_sub_pick,'Un',0);
%                 %probes that bind same RNA at different site
                ProbesAtDifSites = arrayfun(@(x)cell2mat(arrayfun(@(y)TPvec_sub_pick{x}(TSvec_sub_pick{x}(TTvec_sub_pick{x}==Tvec_sub_pick{x}(y))~=Svec_sub_pick{x}(y)),1:length(Tvec_sub_pick{x}),'Un',0)),...
                    1:length(Pi),'Un',0);
                ProbesAtDifSites_AllowedBySpacing = arrayfun(@(x) ...    
                    ProbesAtDifSites{x}(sum(AllowMatrix(ismember(AllowableProbes,[Pi(x) chosenProbes]),ismember(AllowableProbes,ProbesAtDifSites{x})),1)==0),...
                    1:length(Pi),'Un',0);
                Num_MultiSite = cellfun(@length,ProbesAtDifSites_AllowedBySpacing);    
                Num_Excluded = arrayfun(@(x)sum(AllowMatrix(ismember(AllowableProbes,Pi),ismember(AllowableProbes,Pi(x)))),1:length(Pi));
                Gap_Between = diff(Pi)-1; 
                %restrict TPvec where TSvec is different for probe
                clear TSvec_sub_pick TPsub_sub_pick Tvec_sub_pick Svec_sub_pick TvecID_sub_pick                                         
                %Pick Probe 
                %a little more complex
                Priority_Pick = Pi(sub_pick_Probes_SpecificityScoreS(ismember(sub_pick_Probes_Pass,Pi))==min(sub_pick_Probes_SpecificityScoreS(ismember(sub_pick_Probes_Pass,Pi))));                
%               Toff_Target = Tvec_RNA{Priority_Pick};
%               Toff_Site = Svec_RNA{Priority_Pick};
%               Toff_KOFF = Tvec_logKOFF_RNA{Priority_Pick};
%               Toff_KOFFdivON = Tvec_logKOFFdivON_RNA{Priority_Pick};
%               multiplicityRestriction_Targets = [];
%               multiplicityRestriction_notSites = [];

%               ProbesAtDifSites  add penalty
                probe_num = Priority_Pick(1);Lp = length(probes{probe_num,2});
                if (spacing_matrix(probes{probe_num,3}+spacing_req+Lp-1) == 0 && ismember(probe_num,AllowableProbes)==1)%If the space for the probe is open
                    chosenProbes = [chosenProbes Priority_Pick(1)];
                    spacing_matrix(probes{probe_num,3}:probes{probe_num,3} + spacing_req*2+Lp+ Lpmin-1) = 1; 
                    Pi_Overlap= Pi(AllowMatrix(ismember(AllowableProbes,Pi),ismember(AllowableProbes,Priority_Pick(1)))==1);
                    Pi = setdiff(Pi,[Priority_Pick(1) Pi_Overlap]);     
                end
                if (isempty(Pi))
                    finishOptions = 1; 
                end
            end      
        end
    end
end      

% 
%Pset0 = ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes;  
%Pset = ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes;
%Pset(ismember(Pset,PProbesWithRibosomalHits))=[];
%Pset_SpecificitySorted = ProbeDesignResults3.SoftwareResults(outGroup(v)).Designed_Probes_ZigZagSelectionSorted;
%Pset_SpecificitySorted(ismember(Pset_SpecificitySorted,ProbesWithRibosomalHits))=[];
%   
%Noff_Pset_SpecificitySorted = Pset_SpecificitySorted;
%Noff_SortedF = arrayfun(@(n) length([Tvec_RNA{Pset_SpecificitySorted(n)}]),1:length(Pset_SpecificitySorted));
%Noff_SortedCumulativeF = arrayfun(@(n) length([Tvec_RNA{Pset_SpecificitySorted(1:n)}]),1:length(Pset_SpecificitySorted));
%G_Func = @(F) gradient(log(F));
%DecisionFunc = Noff_SortedCumulativeF;
%RM1 = max(find(islocalmin(G_Func(DecisionFunc))))+2;
%RM2 = max(find(diff(diff(log(DecisionFunc)))>0))+2;
%RM3{outGroup(v)} = max(RM1,RM2):length(DecisionFunc);
%RM3_P{outGroup(v)} =  Noff_Pset_SpecificitySorted(max(RM1,RM2):length(DecisionFunc));
%SpecSorted{outGroup(v)} = Noff_Pset_SpecificitySorted;  

end   