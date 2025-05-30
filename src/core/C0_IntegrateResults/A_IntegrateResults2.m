function A_IntegrateResults2(id)
%CD14 probe 41 does not appear in binding site map in CD14 v7Data
    inputs1 = {...  
    {'NM_000805.5'},{},{}, 'Human','(GAST)','(GAST)','17',{},1 ;...          %ENST00000329402.4 ,  465bp, 1-Iso 
    {'NM_005564.5'},{},{}, 'Human','(LCN2)','(LCN2)','9',{},1 ;...           %ENST00000277480.7 ,  820bp, 2-Isos
    {'NM_016307.4'},{},{}, 'Human','(PRRX2)','(PRRX2)','9',{},1 ;...         %ENST00000372469.6 , 1305bp, 2-Isos
    {'NM_000591.4'},{},{}, 'Human','(CD14)','(CD14)','5',{},1 ;...           %ENST00000302014.11, 1356bp, 4-Isos
    {'NM_000045.4'},{},{}, 'Human','(ARG1)','(ARG1)','6',{},1 ;...           %ENST00000368087.8 , 1447bp, 4-Isos
    {'NM_001305654.2'},{},{}, 'Human','(LIME1)','(LIME1)','20',{},1 ;...     %N/A               , 1559bp, 3-Isos closest match ENST00000632538.1	
    {'NM_001012414.3'},{},{}, 'Human','(TRIM61)','(TRIM61)','4',{},1 ;...    %N/A               , 1571bp, 2-Isos
    {'NM_018953.4'},{},{}, 'Human','(HOXC5)','(HOXC5)','12',{},1 ;...        %ENST00000312492.3 , 1611bp, 2-Isos
    {'NM_001165255.2'},{},{}, 'Human','(MNX1)','(MNX1)','7',{},1 ;...        %N/A               , 1619bp, 2-Isos
    {'NM_206963.2'},{},{}, 'Human','(RARRES1)','(RARRES1)','3',{},1 ;...     %ENST00000237696.10, 1713bp, 3-Isos
    {'NM_001286732.2'},{},{}, 'Human','(GRTP1)','(GRTP1)','13',{},1 ;...     %N/A               , 1831bp, 4-Isos
    {'NM_002220.3'},{},{}, 'Human','(ITPKA)','(ITPKA)','15',{},1 ;...        %ENST00000260386.7 , 1825bp, 2-Isos
    {'NM_001357731.1'},{},{}, 'Human','(EIF2S3B)','(EIF2S3B)','12',{},1 ;... %N/A               , 1866bp, 2-Isos
    {'NM_001375593.1'},{},{}, 'Human','(PRR7)','(PRR7)','5',{},1 ;...        %N/A               , 1897bp, 5-Isos
    {'NM_013357.2'},{},{}, 'Human','(PURG)','(PURG)','8',{},1 ;...           %N/A               , 1972bp, 4-Isos
    {'NM_001099773.2'},{},{}, 'Human','(CYP11A1)','(CYP11A1)','15',{},1 ;... %N/A               , 1986bp, 2-Isos
    {'NM_019065.3'},{},{}, 'Human','(NECAB2)','(NECAB2)','16',{},1 ;...      %ENST00000305202.9 , 1996bp, 4-Isos
    {''},{},{}, 'Mice','(XITE)','(XITE)','',{},1 ;... 
    {'NR_028380.1'},{},{}, 'Mice','(FTX)','(FTX)','X',{},1 ;... %4387 ENSMUST has 35 isoforms
    {'NR_028381.1'},{},{}, 'Mice','(FTX)','(FTX)','X',{},1 ;... %4272 ENSMUST00000130063.8	
    {''},{},{}, 'Mice','(RNF12)','(RNF12)','',{},1 ;... 
    {'NR_002844.2'},{},{}, 'Mice','(TSIX)','(TSIX)','X',{},1 ;... %4306
    {'NR_001570.2'},{},{}, 'Mice','(XIST)','(XIST)','X',{},1 ;... %12250
    {'NR_001463.3'},{},{}, 'Mice','(XIST)','(XIST)','X',{},1 ;... %17918   
    {'NR_025508.3'},{},{}, 'Mice','(JPX)','(JPX)','X',{},1 ;...%3802  ENSMUST00000181020.10	
    {''},{},{}, 'Mice','(FTX3)','(FTX3)','',{},1 ;... 
    {'NM_173740.3'},{},{}, 'Mice','(MAOA)','(MAOA)','X',{},1 ;... %4161 ENSMUST00000026013.6
    {'NR_038835.1'},{},{}, 'Human','(MNX1-AS1)','(MNX1-AS1)','7',{},1 ;... %992bp closest EMSEMBL 1041bp ENST00000480284.1	
    {'NR_147077.1'},{},{}, 'Human','(MNX1-AS2)','(MNX1-AS2)','7',{},1 ;... %362bp closest EMSEMBL 371bp ENST00000429228.1
    {'NR_120385.1'},{},{}, 'Human','(GRTP1-AS1)','(GRTP1-AS1)','13',{},1 ;... % 608bp closest EMSEMBL 465  ENST00000423246.1	
    {'NR_046541.1'},{},{}, 'Human','(GRTP1-AS1)','(GRTP1-AS1)','13',{},1 ;... % 1564bp closest EMSEMBL 1171 ENST00000669000.1	
    {'NR_038915.1'},{},{}, 'Human','(PRR7-AS1)','(PRR7-AS1)','5',{},1 ;... % 1600bp closest EMSEMBL 741bp ENST00000425316.3	
    {'NR_038916.1'},{},{}, 'Human','(PRR7-AS1)','(PRR7-AS1)','5',{},1 ;... % 1401bp closest EMSEMBL	741bp ENST00000425316.3	 
    {'NR_160935.1'},{},{}, 'Human','(FAM218A)','(TRIM61-AS1)','4',{},1 ;... % 2184bp closest EMSEMBL is 2175bp ENST00000648094.1
    {'NR_132102.1'},{},{}, 'Human','(PRRX2-AS1)','(PRRX2-AS1)','9',{},1 ;... % 502bp closest EMSEMBL	429bp ENST00000440413.1	
    {'NR_047498.1'},{},{}, 'Human','(LINC00853)','(PDZK1IP1-AS1)','1',{},1 ;... % 669bp closest EMSEMBL	661bp ENST00000429328.2
    {'NM_001318485.2'},{},{}, 'Human','(MRPL20)','(MRPL20)','1',{},1 ;... % 
    {'NM_002046.7'},{},{}, 'Human','(GAPDH)','(GAPDH)','12',{},1 ;... % 
    {'NM_001553.3'},{},{}, 'Human','(IGFBP7)','(IGFBP7)','4',{},1 ;... % 
    {'NM_015710.5'},{},{}, 'Human','(NOP53)','(NOP53)','19',{},1 ;... % 
    {'NM_014220.3'},{},{}, 'Human','(TM4SF1)','(TM4SF1)','3',{},1 ;... % 
    {'NM_001660.4'},{},{}, 'Human','(ARF4)','(ARF4)','3',{},1 ;... % 
    {'NM_004583.4'},{},{}, 'Human','(RAB5C)','(RAB5C)','17',{},1 ;... % 
    {'NM_004905.3'},{},{}, 'Human','(PRDX6)','(PRDX6)','1',{},1 ;... % 
    {'NM_003405.4'},{},{}, 'Human','(YWHAH)','(YWHAH)','22',{},1 ;... % 
    {'NM_201434.3'},{},{}, 'Human','(RAB5C)','(RAB5C)','17',{},1 ;... % 
    {'NM_001101.5'},{},{}, 'Human','(ACTB)','(ACTB)','7',{},1 ;... % 
    {'NM_001252039.2'},{},{}, 'Human','(RAB5C)','(RAB5C)','17',{},1 ;... % 
    {'NM_003380.5'},{},{}, 'Human','(VIM)','(VIM)','10',{},1 ;... % 
    {'NR_132370.1'},{},{}, 'Human','(ARF4-AS1)','(ARF4-AS1)','3',{},1 ;... % 
    {'NR_132382.1'},{},{}, 'Human','(NOP53-AS1)','(NOP53-AS1)','19',{},1 ;... % 
    {'NR_132369.1'},{},{}, 'Human','(ARF4-AS1)','(ARF4-AS1)','3',{},1 ;... % 
    {'NR_109810.1'},{},{}, 'Human','(TM4SF1-AS1)','(TM4SF1-AS1)','3',{},1 ;... % 
    {'NR_046650.1'},{},{}, 'Human','(TM4SF1-AS1)','(TM4SF1-AS1)','3',{},1 ;... % 
    {'NR_109809.1'},{},{}, 'Human','(TM4SF1-AS1)','(TM4SF1-AS1)','3',{},1 ;... % 
    {'NR_126507.1'},{},{}, 'Human','(YWHAH-AS1)','(YWHAH-AS1)','22',{},1 ;... % 
    {'NR_171018.1'},{},{}, 'Human','(YWHAH-AS1)','(YWHAH-AS1)','22',{},1 ;... % 
    {'NR_108060.1'},{},{}, 'Human','(VIM-AS1)','(VIM-AS1)','10',{},1 ;... % 
    {'NR_171020.1'},{},{}, 'Human','(YWHAH-AS1)','(YWHAH-AS1)','22',{},1 ;... % 
    {'NR_034081.1'},{},{}, 'Human','(IGFBP7-AS1)','(IGFBP7-AS1)','17',{},1 ;... % 
    {'NR_171019.1'},{},{}, 'Human','(YWHAH-AS1)','(YWHAH-AS1)','22',{},1 ;... % 
    {'XR_001752888.2'},{},{}, 'Human','(RAB5C-AS1)','(RAB5C-AS1)','17',{},1 ;... % 
    {'NR_015434.1'},{},{}, 'Human','(MRPL20-AS1)','(MRPL20-AS1)','1',{},1 ;... % 
    {'NR_125960.1'},{},{}, 'Human','(PRDX6-AS1)','(PRDX6-AS1)','1',{},1 ;... % 
    {'NR_108061.1'},{},{}, 'Human','(VIM-AS1)','(VIM-AS1)','10',{},1 ;... % 
    {'NM_001199893.2'},{},{}, 'Human','(ACTG2)','(ACTG2)','2',{},1 ;... %1372bp actin, gamma-enteric smooth muscle
    {'NM_005159.5'},{},{}, 'Human','(ACTC1)','(ACTC1)','15',{},1 ;... %1382bp actin, alpha cardiac muscle 1   
    {'NM_001406485.1'},{},{}, 'Human','(ACTC1)','(ACTC1)','15',{},1 ;... %1410bp actin, alpha cardiac muscle 1  
    {'NM_001406482.1'},{},{}, 'Human','(ACTC1)','(ACTC1)','15',{},1 ;... %1416bp actin, alpha cardiac muscle 1    
    {'NM_001406484.1'},{},{}, 'Human','(ACTC1)','(ACTC1)','15',{},1 ;... %1463bp actin, alpha cardiac muscle 1 
    {'NM_001100.4'},{},{}, 'Human','(ACTA1)','(ACTA1)','1',{},1 ;... % 1491bp actin alpha 1, skeletal striated muscle
    {'NM_001615.4'},{},{}, 'Human','(ACTG2)','(ACTG2)','2',{},1 ;... %1501bp actin, gamma-enteric smooth muscle
    {'NM_001024675.2'},{},{}, 'Human','(ACTL10)','(ACTL10)','20',{},1 ;... %1583bp actin-like protein 10
    {'NM_177989.4'},{},{}, 'Human','(ACTL6A)','(ACTL6A)','3',{},1 ;... %1710bp 	actin-like protein 6A   
    {'NM_001406483.1'},{},{}, 'Human','(ACTC1)','(ACTC1)','15',{},1 ;... %1730bp actin, alpha cardiac muscle 1   
    {'NR_037688.3'},{},{}, 'Human','(ACTG1)','(ACTG1)','17',{},1 ;... %1770bp    actin gamma 1 noncoding RNA
    {'NM_004301.5'},{},{}, 'Human','(ACTL6A)','(ACTL6A)','3',{},1 ;... %1854bp 	actin-like protein 6A    
    {'NM_178042.4'},{},{}, 'Human','(ACTL6A)','(ACTL6A)','3',{},1 ;... %1902bp 	actin-like protein 6A 
    {'NM_001614.5'},{},{}, 'Human','(ACTG1)','(ACTG1)','17',{},1 ;... %1919bp    actin gamma 1, cytoplasmic 2
    {'NM_001199954.3'},{},{}, 'Human','(ACTG1)','(ACTG1)','17',{},1 ;... %2038bp   actin gamma 1, cytoplasmic 2
    {'NM_024855.4'},{},{}, 'Human','(ACTR5)','(ACTR5)','20',{},1 ;... %2547bp  ,actin-related protein 5, Arp5, Ino80M
    {'NM_022899.5'},{},{}, 'Human','(ACTR8)','(ACTR8)','3',{},1 ;... %3579bp   	actin-related protein 8, ARP8, Ino80N      
    {'NM_001410774.1'},{},{}, 'Human','(ACTR8)','(ACTR8)','3',{},1 ;... %3579bp   	actin-related protein 8, ARP8, Ino80N   
    {'NM_001145442.1'},{},{}, 'Human','(POTEM)','(POTEM)','14',{},1 ;... %6666bp 	putative POTE ankyrin domain family member M
    };
PnNonZero = @(Pn) Pn(Pn~=0);
edgeNonZero = @(Edge,Pn) Edge(Pn~=0);
PnUp = @(Pn,k) Pn(Pn>=k);
edgeUp = @(Edge,Pn,k) Edge(Pn>=k);
PnUp2 = @(Pn,k) [zeros(1,2) Pn(Pn>=k) zeros(1,2)];
edgeUp2 = @(Edge,Pn,Edge0,Pn0,k) sort([min(Edge0(Pn0>=k)) min(Edge(Pn>=k))-1/10 Edge(Pn>=k) max(Edge(Pn>=k))+1/10 max(Edge0(Pn0>=k))]);
%Neuert Lab initial Probes,simulated annealing Probes, Old Designer Probes,
% Stellaris, Oligostan, MERFISH, PaintSHOP, Unspecific (design in reverse)


%slight problem being a probe dilution dependent effect on on/off-target binding.
%NameTile gs.Name(find(contains(gs.Name,'NM_001553.3'),1))

%% Settings Specification
saveRoot = '/nobackup/p_neuert_lab/Jason/2023-12-08 Probe Designer_JH/';
customBlastDatabase = 'N/A';
gene_num = id;
max_probes = 96;
minProbeSize = 20;
maxProbeSize = 20;
MininumProbeSpacing = 3;
cellPreset = 1;
updateLocation = 1;currLoc = 3;
HybridizationTemperature = 37;
SaltConcentration = 0.05;
RemoveMisMatches = 1;
SpecificityThreshold = 2;
DecisionAlgorithmFloorSize = 0.5;
RemoveRibosomalHits = 1;
RunOffline = 1;
DoAllGenesHaveSameExpression = 0;
UseGeneOverTranscLevelExpression = 0; 
nullRNAcopynumber = 100;
nullDNAcopynumber = 2;
%Save Settings
T_hybrid = HybridizationTemperature;
%%Load Data

%Neuert Lab initial Probes,simulated annealing Probes, Old Designer Probes,
% Stellaris, Oligostan, MERFISH, PaintSHOP, Unspecific (design in reverse)

%% Settings Specification
%Save Settings
%% Compile Results from  Supporting Output Files
settings.FolderRootName = strcat(saveRoot,inputs1{gene_num,5},'_',strjoin(inputs1{gene_num,1},'_'));
settings.rootName = strjoin(inputs1{gene_num,1},'_');


Nsoftware = 11;
designerName0 = '_NLPDS';
Nmodel = 4;RoundN = 1;Z = 1;iz = 1;addSelfProb = 1;packOptimal = 1;spacing = 0;
load([settings.FolderRootName '/' inputs1{gene_num,5} '_' num2str(HybridizationTemperature) '_DesignResults_CarryOver.mat'],'CarryOver','CarryOver2','CarryOver3')

for id = 1
ProbeDesignResults2.Gene = inputs1{gene_num,5};
ProbeDesignResults2.GeneTarget = settings.rootName;
ProbeDesignResults2.Organism = inputs1{gene_num,4};
ProbeDesignResults2.FolderName = settings.FolderRootName;
    for id2 = 1:Nsoftware
        switch id2
            case 1
                designerName = '_NLPDS';
            case 2
                designerName = '_Stellaris_L5';
            case 3
                designerName = '_Stellaris_L4'; 
            case 4
                designerName = '_Stellaris_L3';
            case 5
                designerName = '_Stellaris_L2';  
            case 6
                designerName = '_Stellaris_L1';
            case 7
                designerName = '_Stellaris_L0';  
            case 8
                designerName = '_OligoStanHT';
            case 9
                designerName = '_PaintSHOP';
            case 10
                designerName = '_ZZLMERFISH';
            case 11
                designerName = '_AIBSMERFISH';  
        end    
        ProbeDesignResults2.Software{id2} = designerName; 
        try
        if (~contains(designerName,'Stellaris'))    
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_probes' designerName '.mat'],'probes');
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_ExpressionInfo' designerName '.mat'],'ExpressionMatrix','settings');
        load([settings.FolderRootName '/' inputs1{gene_num,5}  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' designerName '.mat'],'Kb_mod')     
        load([settings.FolderRootName '/' inputs1{gene_num,5}  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' designerName '.mat'],'Kb_Complement')      
        load([settings.FolderRootName '/'  inputs1{gene_num,5} '_BindingMatrices' designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','dCp_mod','Tm_mod')
        load([settings.FolderRootName '/' inputs1{gene_num,5}  '_BindingMatrices2' designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','dCp_Complement')
        if (id2==1)
           load([settings.FolderRootName '/' inputs1{gene_num,5} '_' inputs1{gene_num,1}{:}  '_chosen.mat'],'chosenProbes');  
        else
           chosenProbes = 1:size(probes,1);
        end
        AllowableProbes = CarryOver.AllowableProbes{id2};
        ProbesWithRibosomalHits = CarryOver.ProbesWithRibosomalHits{id2};
        ON_RNAIDs = CarryOver.ON_RNAIDs{id2};
        OFF_RNAIDs = CarryOver.OFF_RNAIDs{id2};
        DNA_IDs = CarryOver2.DNA_IDs{id2};
        NonDNA_IDs = CarryOver2.NonDNA_IDs{id2};
        Pset_ON_IDs  = CarryOver3.Pset_ON_IDs{id2};
        Pset_OFF_IDs  = CarryOver3.Pset_OFF_IDs{id2};
        NumNonUniformConditions = size(ExpressionMatrix,2);%ADD null case in
        ExpressionMatrix(DNA_IDs,NumNonUniformConditions+1) = 2;
        ExpressionMatrix(NonDNA_IDs,NumNonUniformConditions+1) = 100; 
        Pset = setdiff(chosenProbes,ProbesWithRibosomalHits); 
        Kon = full(max(squeeze(Kb_mod(:,ON_RNAIDs,:,Nmodel)),[],2))';clear NumNonUniformConditions
        Koff = max(squeeze(Kb_mod(:,:,:,Nmodel)),[],3);
        Koff(:,ON_RNAIDs)=0;
        Koff = full(Koff);
        rangeKon = [unique(round(log10(Kon),RoundN)) max(round(log10(Kon),RoundN))+1];
        ProbeTarget_Kon_All = round(log10(Kon(AllowableProbes)),RoundN);
        ProbeTarget_Kon_S = round(log10(Kon(Pset)),RoundN);
        ProbeDesignResults2.SoftwareResults(id2).settings = settings;
        ProbeDesignResults2.SoftwareResults(id2).Kon = Kon;
        ProbeDesignResults2.SoftwareResults(id2).Koff = Koff;
        ProbeDesignResults2.SoftwareResults(id2).Kb = Kb_mod;clear Kb_mod
        ProbeDesignResults2.SoftwareResults(id2).Kb_Complement = Kb_Complement;clear Kb_Complement
        ProbeDesignResults2.SoftwareResults(id2).dHeq_mod = dHeq_mod; clear dHeq_mod
        ProbeDesignResults2.SoftwareResults(id2).dSeq_mod = dSeq_mod; clear dSeq_mod
        ProbeDesignResults2.SoftwareResults(id2).dHf_mod = dHf_mod; clear dHf_mod
        ProbeDesignResults2.SoftwareResults(id2).dSf_mod = dSf_mod; clear dSf_mod
        ProbeDesignResults2.SoftwareResults(id2).dHr_mod = dHr_mod; clear dHr_mod
        ProbeDesignResults2.SoftwareResults(id2).dSr_mod = dSr_mod; clear dSr_mod        
        ProbeDesignResults2.SoftwareResults(id2).dCp_mod = dCp_mod; clear dCp_mod
        ProbeDesignResults2.SoftwareResults(id2).Tm_mod = Tm_mod; clear Tm_mod
        ProbeDesignResults2.SoftwareResults(id2).dHeq_Complement = dHeq_Complement; clear dHeq_Complement
        ProbeDesignResults2.SoftwareResults(id2).dSeq_Complement = dSeq_Complement; clear dSeq_Complement
        ProbeDesignResults2.SoftwareResults(id2).dHf_Complement = dHf_Complement; clear dHf_Complement
        ProbeDesignResults2.SoftwareResults(id2).dSf_Complement = dSf_Complement; clear dSf_Complement
        ProbeDesignResults2.SoftwareResults(id2).dHr_Complement = dHr_Complement; clear dHr_Complement
        ProbeDesignResults2.SoftwareResults(id2).dSr_Complement = dSr_Complement; clear dSr_Complement        
        ProbeDesignResults2.SoftwareResults(id2).dCp_Complement = dCp_Complement; clear dCp_Complement
        Target_OnOff_Specificity =  zeros(size(probes,1),size(Koff,2));
        null_weighted_KON_KOFF_ratio = zeros(1,size(probes,1));
        exp_weighted_KON_KOFF_ratio_TissueSpecific = zeros(size(probes,1),size(ExpressionMatrix,2));
        exp_weighted_KON_KOFF_ratio_TissueSpecific_ONnorm = zeros(size(probes,1),size(ExpressionMatrix,2)); 
        for p=1:size(probes,1)
        Target_OnOff_Specificity(p,:) = Kon(p)./Koff(p,:);
        null_weighted_KON_KOFF_ratio(p) = Kon(p)/sum(Koff(p,OFF_RNAIDs));
        exp_weighted_KON_KOFF_ratio_TissueSpecific(p,:) = sum(Kon(p)*ExpressionMatrix(ON_RNAIDs,:))./((squeeze(Koff(p,OFF_RNAIDs))*ExpressionMatrix(OFF_RNAIDs,:))+1*iz);
        exp_weighted_KON_KOFF_ratio_TissueSpecific_ONnorm(p,:) = (sum(Kon(p)*ExpressionMatrix(ON_RNAIDs,:))./((squeeze(Koff(p,OFF_RNAIDs))*ExpressionMatrix(OFF_RNAIDs,:))+1*iz))./sum(ExpressionMatrix(ON_RNAIDs,:),1);
        end 
        clear Kon Koff
        Probe_Specificities.exp_weighted_KON_KOFF_ratio_TissueSpecific = exp_weighted_KON_KOFF_ratio_TissueSpecific;clear exp_weighted_KON_KOFF_ratio_TissueSpecific
        Probe_Specificities.exp_weighted_KON_KOFF_ratio_TissueSpecific_KONNorm = exp_weighted_KON_KOFF_ratio_TissueSpecific_ONnorm;clear exp_weighted_KON_KOFF_ratio_TissueSpecific_ONnorm
        Probe_Specificities.null_weighted_KON_KOFF_ratio = null_weighted_KON_KOFF_ratio;clear null_weighted_KON_KOFF_ratio
        Probe_Specificities.Target_OnOff_Specificity = Target_OnOff_Specificity;clear Target_OnOff_Specificity
        Probe_SpecificitySorted_CellLine = zeros(length(chosenProbes),size(ExpressionMatrix,2));
        Probe_SpecificitySorted_CellLine_KONNorm = zeros(length(chosenProbes),size(ExpressionMatrix,2));
        for v = 1:size(ExpressionMatrix,2)
        [~, expOrder] = sort(Probe_Specificities.exp_weighted_KON_KOFF_ratio_TissueSpecific(chosenProbes,v),'descend','MissingPlacement','last');
        Probe_SpecificitySorted_CellLine(:,v) = chosenProbes(expOrder);
        [~, expOrder] = sort(Probe_Specificities.exp_weighted_KON_KOFF_ratio_TissueSpecific_KONNorm(chosenProbes,v),'descend','MissingPlacement','last');
        Probe_SpecificitySorted_CellLine_KONNorm(:,v) = chosenProbes(expOrder);
        end
        [~, expOrder] = sort(Probe_Specificities.null_weighted_KON_KOFF_ratio(chosenProbes),'descend','MissingPlacement','last');
        Probe_SpecificitySorted_Null= chosenProbes(expOrder);clear expOrder
        Probe_Specificities.Probe_SpecificitySorted_CellLine = Probe_SpecificitySorted_CellLine;clear Probe_SpecificitySorted_CellLine
        Probe_Specificities.Probe_SpecificitySorted_CellLine_KONNorm = Probe_SpecificitySorted_CellLine_KONNorm;clear Probe_SpecificitySorted_CellLine_KONNorm
        Probe_Specificities.Probe_SpecificitySorted_Null = Probe_SpecificitySorted_Null;clear Probe_SpecificitySorted_Null
        ProbeDesignResults2.SoftwareResults(id2).Probe_Specificities = Probe_Specificities;clear Probe_Specificities
        clear ExpressionMatrix
        try
        [K_S,~,~,~,~,~,~,...
        ~,K_CD,~,~,~,~,~,~,~,...
        ~,~] = A_JH_GenerateSecondaryStructureInfo_V2(probes,Pset,settings); 
        Ks = full(sum(squeeze(K_S(Pset,:,Nmodel)),2,'omitnan'));clear K_S
        Kd = reshape(full(sum(squeeze(K_CD(Pset,Pset,:,Nmodel)),3,'omitnan')),[],1);clear K_CD
        catch
           Ks = [];
           Kd = [];
        end
        ProbeDesignResults2.SoftwareResults(id2).probes = probes; clear probes
        ProbeDesignResults2.SoftwareResults(id2).Designed_Probes = chosenProbes; clear chosenProbes
        ProbeDesignResults2.SoftwareResults(id2).Probe_Ks = Ks;clear Ks
        ProbeDesignResults2.SoftwareResults(id2).Probe_Kd = Kd;clear Kd
        ProbeTarget_Tm = CarryOver3.ProbeTarget_Tm{id2};
        ProbeTarget_GC = CarryOver3.ProbeTarget_GC{id2};
        rangeTm = [-1 0:5:100];
        rangeGC = [-1 0:5:100];
        ProbeTarget_GCOn_All = ProbeTarget_GC(Pset_ON_IDs(AllowableProbes));
        ProbeTarget_GCOff_All = ProbeTarget_GC(Pset_OFF_IDs(AllowableProbes));
        ProbeTarget_TmOn_All = round(ProbeTarget_Tm(Z,Pset_ON_IDs(AllowableProbes)));
        ProbeTarget_TmOff_All = round(ProbeTarget_Tm(Z,Pset_OFF_IDs(AllowableProbes)));
        ProbeTarget_GCOn_S = ProbeTarget_GC(Pset_ON_IDs(Pset));
        ProbeTarget_GCOff_S = ProbeTarget_GC(Pset_OFF_IDs(Pset));clear ProbeTargetGC
        ProbeTarget_TmOn_S = round(ProbeTarget_Tm(Z,Pset_ON_IDs(Pset)));
        ProbeTarget_TmOff_S = round(ProbeTarget_Tm(Z,Pset_OFF_IDs(Pset)));clear ProbeTarget_Tm
        [PnKon_All,edgeKon_All] = histcounts(ProbeTarget_Kon_All,'BinEdges',rangeKon,'Normalization','Probability');
        [PnGCOn_All,edgeGCOn_All] = histcounts(ProbeTarget_GCOn_All,'BinEdges',rangeGC,'Normalization','Probability');
        [PnGCOff_All,edgeGCOff_All] = histcounts(ProbeTarget_GCOff_All,'BinEdges',rangeGC,'Normalization','Probability');
        [PnTmOn_All,edgeTmOn_All] = histcounts(ProbeTarget_TmOn_All,'BinEdges',rangeTm,'Normalization','Probability');
        [PnTmOff_All,edgeTmOff_All] = histcounts(ProbeTarget_TmOff_All,'BinEdges',rangeTm,'Normalization','Probability');
        [PnKon_S,edgeKon_S] = histcounts(ProbeTarget_Kon_S,'BinEdges',rangeKon,'Normalization','Probability');hold on
        [PnGCOn_S,edgeGCOn_S] = histcounts(ProbeTarget_GCOn_S,'BinEdges',rangeGC,'Normalization','Probability');hold on
        [PnGCOff_S,edgeGCOff_S] = histcounts(ProbeTarget_GCOff_S,'BinEdges',rangeGC,'Normalization','Probability');hold on
        [PnTmOn_S,edgeTmOn_S] = histcounts(ProbeTarget_TmOn_S,'BinEdges',rangeTm,'Normalization','Probability');hold on
        [PnTmOff_S,edgeTmOff_S] = histcounts(ProbeTarget_TmOff_S,'BinEdges',rangeTm,'Normalization','Probability');hold on
        ProbeDesignResults2.SoftwareResults(id2).AllowableProbes = AllowableProbes;
        ProbeDesignResults2.SoftwareResults(id2).PnKon_S = PnKon_S;clear PnKon_S
        ProbeDesignResults2.SoftwareResults(id2).edgeKon_S = edgeKon_S;clear edgeKon_S
        ProbeDesignResults2.SoftwareResults(id2).PnGCOn_S = PnGCOn_S;clear PnGCOn_S
        ProbeDesignResults2.SoftwareResults(id2).edgeGCOn_S = edgeGCOn_S;clear edgeGCOn_S
        ProbeDesignResults2.SoftwareResults(id2).PnGCOff_S = PnGCOff_S;clear PnGCOff_S
        ProbeDesignResults2.SoftwareResults(id2).edgeGCOff_S = edgeGCOff_S;clear edgeGCOff_S
        ProbeDesignResults2.SoftwareResults(id2).PnTmOn_S = PnTmOn_S;clear PnTmOn_S
        ProbeDesignResults2.SoftwareResults(id2).edgeTmOn_S = edgeTmOn_S;clear edgeTmOn_S
        ProbeDesignResults2.SoftwareResults(id2).PnTmOff_S = PnTmOff_S;clear PnTmOff_S
        ProbeDesignResults2.SoftwareResults(id2).edgeTmOff_S = edgeTmOff_S;clear edgeTmOff_S;
        ProbeDesignResults2.SoftwareResults(id2).PnKon_All = PnKon_All;clear PnKon_All
        ProbeDesignResults2.SoftwareResults(id2).edgeKon_All = edgeKon_All;clear edgeKon_All
        ProbeDesignResults2.SoftwareResults(id2).PnGCOn_All = PnGCOn_All;clear PnGCOn_All
        ProbeDesignResults2.SoftwareResults(id2).edgeGCOn_All = edgeGCOn_All;clear edgeGCOn_All
        ProbeDesignResults2.SoftwareResults(id2).PnGCOff_All = PnGCOff_All;clear PnGCOff_All
        ProbeDesignResults2.SoftwareResults(id2).edgeGCOff_All = edgeGCOff_All;clear edgeGCOff_All
        ProbeDesignResults2.SoftwareResults(id2).PnTmOn_All = PnTmOn_All;clear PnTmOn_All;
        ProbeDesignResults2.SoftwareResults(id2).edgeTmOn_All = edgeTmOn_All;clear edgeTmOn_All
        ProbeDesignResults2.SoftwareResults(id2).PnTmOff_All = PnTmOff_All;clear PnTmOff_All
        ProbeDesignResults2.SoftwareResults(id2).edgeTmOff_All = edgeTmOff_All;clear edgeTmOff_All
        else
        %SL
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_probes' designerName0 '.mat'],'probes');
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_ExpressionInfo' designerName0 '.mat'],'ExpressionMatrix','settings');
        load([settings.FolderRootName '/' inputs1{gene_num,5}  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' designerName0 '.mat'],'Kb_mod')     
        load([settings.FolderRootName '/' inputs1{gene_num,5}  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' designerName0 '.mat'],'Kb_Complement') 
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_probes' designerName '.mat'],'chosenProbes');
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_BindingMatrices' designerName0 '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','dCp_mod','Tm_mod')
        load([settings.FolderRootName '/' inputs1{gene_num,5}  '_BindingMatrices2' designerName0 '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','dCp_Complement')
       
        AllowableProbes = CarryOver.AllowableProbes{id2};
        ProbesWithRibosomalHits = CarryOver.ProbesWithRibosomalHits{id2};
        ON_RNAIDs = CarryOver.ON_RNAIDs{id2};
        OFF_RNAIDs = CarryOver.OFF_RNAIDs{id2};
        DNA_IDs = CarryOver2.DNA_IDs{id2};
        NonDNA_IDs = CarryOver2.NonDNA_IDs{id2};
        Pset_ON_IDs  = CarryOver3.Pset_ON_IDs{id2};
        Pset_OFF_IDs  = CarryOver3.Pset_OFF_IDs{id2};
        NumNonUniformConditions = size(ExpressionMatrix,2);%ADD null case in
        ExpressionMatrix(DNA_IDs,NumNonUniformConditions+1) = 2;
        ExpressionMatrix(NonDNA_IDs,NumNonUniformConditions+1) = 100;
        Pset = setdiff(chosenProbes,ProbesWithRibosomalHits);
        Kon = full(max(squeeze(Kb_mod(:,ON_RNAIDs,:,4)),[],2))';clear NumNonUniformConditions
        Koff = max(squeeze(Kb_mod(:,:,:,4)),[],3);
        Koff(:,ON_RNAIDs)=0;
        Koff = full(Koff);
        rangeKon = [unique(round(log10(Kon),RoundN)) max(round(log10(Kon),RoundN))+1];
        ProbeTarget_Kon_All = round(log10(Kon(AllowableProbes)),RoundN);
        ProbeTarget_Kon_S = round(log10(Kon(Pset)),RoundN);
        
        ProbeDesignResults2.SoftwareResults(id2).settings = settings;
        ProbeDesignResults2.SoftwareResults(id2).Kon = Kon;
        ProbeDesignResults2.SoftwareResults(id2).Koff = Koff;
        ProbeDesignResults2.SoftwareResults(id2).Kb = Kb_mod;clear Kb_mod
        ProbeDesignResults2.SoftwareResults(id2).Kb_Complement = Kb_Complement;clear Kb_Complement
        ProbeDesignResults2.SoftwareResults(id2).dHeq_mod = dHeq_mod; clear dHeq_mod
        ProbeDesignResults2.SoftwareResults(id2).dSeq_mod = dSeq_mod; clear dSeq_mod
        ProbeDesignResults2.SoftwareResults(id2).dHf_mod = dHf_mod; clear dHf_mod
        ProbeDesignResults2.SoftwareResults(id2).dSf_mod = dSf_mod; clear dSf_mod
        ProbeDesignResults2.SoftwareResults(id2).dHr_mod = dHr_mod; clear dHr_mod
        ProbeDesignResults2.SoftwareResults(id2).dSr_mod = dSr_mod; clear dSr_mod        
        ProbeDesignResults2.SoftwareResults(id2).dCp_mod = dCp_mod; clear dCp_mod
        ProbeDesignResults2.SoftwareResults(id2).Tm_mod = Tm_mod; clear Tm_mod
        ProbeDesignResults2.SoftwareResults(id2).dHeq_Complement = dHeq_Complement; clear dHeq_Complement
        ProbeDesignResults2.SoftwareResults(id2).dSeq_Complement = dSeq_Complement; clear dSeq_Complement
        ProbeDesignResults2.SoftwareResults(id2).dHf_Complement = dHf_Complement; clear dHf_Complement
        ProbeDesignResults2.SoftwareResults(id2).dSf_Complement = dSf_Complement; clear dSf_Complement
        ProbeDesignResults2.SoftwareResults(id2).dHr_Complement = dHr_Complement; clear dHr_Complement
        ProbeDesignResults2.SoftwareResults(id2).dSr_Complement = dSr_Complement; clear dSr_Complement        
        ProbeDesignResults2.SoftwareResults(id2).dCp_Complement = dCp_Complement; clear dCp_Complement
        Target_OnOff_Specificity =  zeros(size(probes,1),size(Koff,2));
        null_weighted_KON_KOFF_ratio = zeros(1,size(probes,1));
        exp_weighted_KON_KOFF_ratio_TissueSpecific = zeros(size(probes,1),size(ExpressionMatrix,2));
        exp_weighted_KON_KOFF_ratio_TissueSpecific_ONnorm = zeros(size(probes,1),size(ExpressionMatrix,2)); 
        for p=1:size(probes,1)
        Target_OnOff_Specificity(p,:) = Kon(p)./Koff(p,:);
        null_weighted_KON_KOFF_ratio(p) = Kon(p)/sum(Koff(p,OFF_RNAIDs));
        exp_weighted_KON_KOFF_ratio_TissueSpecific(p,:) = sum(Kon(p)*ExpressionMatrix(ON_RNAIDs,:))./((squeeze(Koff(p,OFF_RNAIDs))*ExpressionMatrix(OFF_RNAIDs,:))+1*iz);
        exp_weighted_KON_KOFF_ratio_TissueSpecific_ONnorm(p,:) = (sum(Kon(p)*ExpressionMatrix(ON_RNAIDs,:))./((squeeze(Koff(p,OFF_RNAIDs))*ExpressionMatrix(OFF_RNAIDs,:))+1*iz))./sum(ExpressionMatrix(ON_RNAIDs,:),1);
        end 
        clear Kon Koff
        Probe_Specificities.exp_weighted_KON_KOFF_ratio_TissueSpecific = exp_weighted_KON_KOFF_ratio_TissueSpecific;clear exp_weighted_KON_KOFF_ratio_TissueSpecific
        Probe_Specificities.exp_weighted_KON_KOFF_ratio_TissueSpecific_KONNorm = exp_weighted_KON_KOFF_ratio_TissueSpecific_ONnorm;clear exp_weighted_KON_KOFF_ratio_TissueSpecific_ONnorm
        Probe_Specificities.null_weighted_KON_KOFF_ratio = null_weighted_KON_KOFF_ratio;clear null_weighted_KON_KOFF_ratio
        Probe_Specificities.Target_OnOff_Specificity = Target_OnOff_Specificity;clear Target_OnOff_Specificity
        Probe_SpecificitySorted_CellLine = zeros(length(chosenProbes),size(ExpressionMatrix,2));
        Probe_SpecificitySorted_CellLine_KONNorm = zeros(length(chosenProbes),size(ExpressionMatrix,2));
        for v = 1:size(ExpressionMatrix,2)
        [~, expOrder] = sort(Probe_Specificities.exp_weighted_KON_KOFF_ratio_TissueSpecific(chosenProbes,v),'descend','MissingPlacement','last');
        Probe_SpecificitySorted_CellLine(:,v) = chosenProbes(expOrder);
        [~, expOrder] = sort(Probe_Specificities.exp_weighted_KON_KOFF_ratio_TissueSpecific_KONNorm(chosenProbes,v),'descend','MissingPlacement','last');
        Probe_SpecificitySorted_CellLine_KONNorm(:,v) = chosenProbes(expOrder);
        end
        [~, expOrder] = sort(Probe_Specificities.null_weighted_KON_KOFF_ratio(chosenProbes),'descend','MissingPlacement','last');
        Probe_SpecificitySorted_Null= chosenProbes(expOrder);clear expOrder
        Probe_Specificities.Probe_SpecificitySorted_CellLine = Probe_SpecificitySorted_CellLine;clear Probe_SpecificitySorted_CellLine
        Probe_Specificities.Probe_SpecificitySorted_CellLine_KONNorm = Probe_SpecificitySorted_CellLine_KONNorm;clear Probe_SpecificitySorted_CellLine_KONNorm
        Probe_Specificities.Probe_SpecificitySorted_Null = Probe_SpecificitySorted_Null;clear Probe_SpecificitySorted_Null
        ProbeDesignResults2.SoftwareResults(id2).Probe_Specificities = Probe_Specificities;clear Probe_Specificities
        clear ExpressionMatrix
        try
        [K_S,~,~,~,~,~,~,...
        ~,K_CD,~,~,~,~,~,~,~,...
        ~,~] = A_JH_GenerateSecondaryStructureInfo_V2(probes,Pset,settings); 
        Ks = full(sum(squeeze(K_S(Pset,:,Nmodel)),2,'omitnan'));clear K_S
        Kd = reshape(full(sum(squeeze(K_CD(Pset,Pset,:,Nmodel)),3,'omitnan')),[],1);clear K_CD
        catch
           Ks = [];
           Kd = [];
        end
        ProbeDesignResults2.SoftwareResults(id2).probes = probes; clear probes
        ProbeDesignResults2.SoftwareResults(id2).Designed_Probes = chosenProbes; clear chosenProbes
        ProbeDesignResults2.SoftwareResults(id2).Probe_Ks = Ks;clear Ks
        ProbeDesignResults2.SoftwareResults(id2).Probe_Kd = Kd;clear Kd
        ProbeTarget_Tm = CarryOver3.ProbeTarget_Tm{id2};
        ProbeTarget_GC = CarryOver3.ProbeTarget_GC{id2};
        rangeTm = [-1 0:5:100];
        rangeGC = [-1 0:5:100];
        ProbeTarget_GCOn_All = ProbeTarget_GC(Pset_ON_IDs(AllowableProbes));
        ProbeTarget_GCOff_All = ProbeTarget_GC(Pset_OFF_IDs(AllowableProbes));
        ProbeTarget_TmOn_All = round(ProbeTarget_Tm(Z,Pset_ON_IDs(AllowableProbes)));
        ProbeTarget_TmOff_All = round(ProbeTarget_Tm(Z,Pset_OFF_IDs(AllowableProbes)));
        ProbeTarget_GCOn_S = ProbeTarget_GC(Pset_ON_IDs(Pset));
        ProbeTarget_GCOff_S = ProbeTarget_GC(Pset_OFF_IDs(Pset));clear ProbeTargetGC
        ProbeTarget_TmOn_S = round(ProbeTarget_Tm(Z,Pset_ON_IDs(Pset)));
        ProbeTarget_TmOff_S = round(ProbeTarget_Tm(Z,Pset_OFF_IDs(Pset)));clear ProbeTarget_Tm
        [PnKon_All,edgeKon_All] = histcounts(ProbeTarget_Kon_All,'BinEdges',rangeKon,'Normalization','Probability');
        [PnGCOn_All,edgeGCOn_All] = histcounts(ProbeTarget_GCOn_All,'BinEdges',rangeGC,'Normalization','Probability');
        [PnGCOff_All,edgeGCOff_All] = histcounts(ProbeTarget_GCOff_All,'BinEdges',rangeGC,'Normalization','Probability');
        [PnTmOn_All,edgeTmOn_All] = histcounts(ProbeTarget_TmOn_All,'BinEdges',rangeTm,'Normalization','Probability');
        [PnTmOff_All,edgeTmOff_All] = histcounts(ProbeTarget_TmOff_All,'BinEdges',rangeTm,'Normalization','Probability');
        [PnKon_S,edgeKon_S] = histcounts(ProbeTarget_Kon_S,'BinEdges',rangeKon,'Normalization','Probability');hold on
        [PnGCOn_S,edgeGCOn_S] = histcounts(ProbeTarget_GCOn_S,'BinEdges',rangeGC,'Normalization','Probability');hold on
        [PnGCOff_S,edgeGCOff_S] = histcounts(ProbeTarget_GCOff_S,'BinEdges',rangeGC,'Normalization','Probability');hold on
        [PnTmOn_S,edgeTmOn_S] = histcounts(ProbeTarget_TmOn_S,'BinEdges',rangeTm,'Normalization','Probability');hold on
        [PnTmOff_S,edgeTmOff_S] = histcounts(ProbeTarget_TmOff_S,'BinEdges',rangeTm,'Normalization','Probability');hold on
        ProbeDesignResults2.SoftwareResults(id2).AllowableProbes = AllowableProbes;
        ProbeDesignResults2.SoftwareResults(id2).PnKon_S = PnKon_S;clear PnKon_S
        ProbeDesignResults2.SoftwareResults(id2).edgeKon_S = edgeKon_S;clear edgeKon_S
        ProbeDesignResults2.SoftwareResults(id2).PnGCOn_S = PnGCOn_S;clear PnGCOn_S
        ProbeDesignResults2.SoftwareResults(id2).edgeGCOn_S = edgeGCOn_S;clear edgeGCOn_S
        ProbeDesignResults2.SoftwareResults(id2).PnGCOff_S = PnGCOff_S;clear PnGCOff_S
        ProbeDesignResults2.SoftwareResults(id2).edgeGCOff_S = edgeGCOff_S;clear edgeGCOff_S
        ProbeDesignResults2.SoftwareResults(id2).PnTmOn_S = PnTmOn_S;clear PnTmOn_S
        ProbeDesignResults2.SoftwareResults(id2).edgeTmOn_S = edgeTmOn_S;clear edgeTmOn_S
        ProbeDesignResults2.SoftwareResults(id2).PnTmOff_S = PnTmOff_S;clear PnTmOff_S
        ProbeDesignResults2.SoftwareResults(id2).edgeTmOff_S = edgeTmOff_S;clear edgeTmOff_S;
        ProbeDesignResults2.SoftwareResults(id2).PnKon_All = PnKon_All;clear PnKon_All
        ProbeDesignResults2.SoftwareResults(id2).edgeKon_All = edgeKon_All;clear edgeKon_All
        ProbeDesignResults2.SoftwareResults(id2).PnGCOn_All = PnGCOn_All;clear PnGCOn_All
        ProbeDesignResults2.SoftwareResults(id2).edgeGCOn_All = edgeGCOn_All;clear edgeGCOn_All
        ProbeDesignResults2.SoftwareResults(id2).PnGCOff_All = PnGCOff_All;clear PnGCOff_All
        ProbeDesignResults2.SoftwareResults(id2).edgeGCOff_All = edgeGCOff_All;clear edgeGCOff_All
        ProbeDesignResults2.SoftwareResults(id2).PnTmOn_All = PnTmOn_All;clear PnTmOn_All;
        ProbeDesignResults2.SoftwareResults(id2).edgeTmOn_All = edgeTmOn_All;clear edgeTmOn_All
        ProbeDesignResults2.SoftwareResults(id2).PnTmOff_All = PnTmOff_All;clear PnTmOff_All
        ProbeDesignResults2.SoftwareResults(id2).edgeTmOff_All = edgeTmOff_All;clear edgeTmOff_All
        end
        catch
        end
    end
    save([settings.FolderRootName '/' inputs1{gene_num,5} '_' num2str(HybridizationTemperature) '_DesignResults2.mat'],'ProbeDesignResults2','-v7.3')
end


end


            