function A_IntegrateResults1(id)
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
ProbeDesignResults1.Gene = inputs1{gene_num,5};
ProbeDesignResults1.GeneTarget = settings.rootName;
ProbeDesignResults1.Organism = inputs1{gene_num,4};
ProbeDesignResults1.FolderName = settings.FolderRootName;
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
        ProbeDesignResults1.Software{id2} = designerName; 
        try
        if (~contains(designerName,'Stellaris'))
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_probes' designerName '.mat'],'probes');
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_hits_table' designerName '.mat'],'gene_table')  
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_ExpressionInfo' designerName '.mat'],'settings');
        load([settings.FolderRootName '/' settings.GeneName '_binding_hits_map' designerName '.mat'],'DoesProbeBindSite2'); 
        GeneName = settings.GeneName;
        GeneTarget = settings.transcript_IDs;
        gene_table = sortrows(gene_table,[7 6],'ascend');
        gene_table2 = gene_table(gene_table.Match>=15,:);
        MinusStrandedHits = find(contains(gene_table2.Strand,'Minus'));clear gene_table
        RNA_IDs_1 = find(contains(gene_table2.Name,'NM_'));
        RNA_IDs_2 = find(contains(gene_table2.Name,'NR_'));
        RNA_IDs_3 = find(contains(gene_table2.Name,'XM_'));
        RNA_IDs_4 = find(contains(gene_table2.Name,'XR_'));
        contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
        RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);clear RNA_IDs_*
        gene_table2 = gene_table2(setdiff(1:size(gene_table2,1),RNA_MissedFilteredHits),:);
        gene_table2.Ax = min(gene_table2.SubjectIndices,[],2);clear MinusStrandedHits
        gene_table2.Bx = max(gene_table2.SubjectIndices,[],2);clear RNA_MissedFilteredHits
        gene_table3 = sortrows(gene_table2,[7 13],'ascend');
        Iz = find(endsWith(gene_table3.Name,', ribosomal RNA'));clear gene_table2
        ProbesWithRibosomalHits = unique(gene_table3.ProbeNum(Iz));clear Iz
        AllowableProbes = setdiff(1:size(probes,1),ProbesWithRibosomalHits);
        ProbeDesignResults1.Description = gene_table3.Name(find(contains(gene_table3.Name,settings.rootName),1));
        ProbeDesignResults1.SoftwareResults(id2).DoesProbeBindSite = DoesProbeBindSite2;clear DoesProbeBindSite2
        Names = unique(gene_table3.Name);
        Names = convertCharsToStrings(Names);
        Isoforms = find(contains(Names,GeneName));
        uniNames = extractBefore(Names,'.');
        try
            Desired_Isoforms = find(contains(uniNames,extractBefore(GeneTarget{:},'.')));
        catch
            Desired_Isoforms = find(contains(uniNames,extractBefore(GeneTarget,'.')));  
        end 
        UnDesired_Isoforms = setdiff(Isoforms,Desired_Isoforms);
        ON_RNAIDs = find(strcmp(uniNames,extractBefore(GeneTarget,'.')));
        NonDNA_IDs_1 = find(~contains(uniNames,'NC_'));%IDs
        NonDNA_IDs_2 = find(~contains(uniNames,'NT_'));%IDs
        NonDNA_IDs_3 = find(~contains(uniNames,'NW_'));%
        NonDNA_IDs = intersect(intersect(NonDNA_IDs_1,NonDNA_IDs_2),NonDNA_IDs_3).';clear NonDNA_IDs_*
        Names_NonDNA = Names(NonDNA_IDs);
        PN = gene_table3.ProbeNum(NonDNA_IDs);
        AS = arrayfun(@(x)seqrcomplement(gene_table3.Alignment{NonDNA_IDs(x)}(3,:)),1:length(NonDNA_IDs),'Un',0);
        AS = cellfun(@(x) erase(x,'-'),AS,'Un',0);clear gene_table3 
        ProbeTarget_GC = cell2mat(arrayfun(@(x) oligoprop(AS{x}).GC,1:length(NonDNA_IDs),'Un',0));
        ProbeTarget_Tm = cell2mat(arrayfun(@(x) oligoprop(AS{x},'Salt',settings.SaltConcentration,'Temp',settings.HybridizationTemperature).Tm',1:length(NonDNA_IDs),'Un',0));
        clear AS
        Hit_ON_RNAIDs_Isos = find(contains(Names_NonDNA,GeneName));clear Names_NonDNA
        Hit_OFF_RNAIDs_minusIsos = setdiff(1:length(NonDNA_IDs),Hit_ON_RNAIDs_Isos);
        Hit_ON_IDs = Hit_ON_RNAIDs_Isos;
        Hit_OFF_IDs = Hit_OFF_RNAIDs_minusIsos;
        ProbeDesignResultsC.SoftwareResults(id2).Pset_ON_IDs = @(Q) setdiff(find(ismember(PN,Q)),Hit_ON_IDs);
        ProbeDesignResultsC.SoftwareResults(id2).Pset_OFF_IDs = @(Q) setdiff(find(ismember(PN,Q)),Hit_OFF_IDs);
        OFF_RNAIDs = setdiff(NonDNA_IDs,ON_RNAIDs);
        OFF_IDs = OFF_RNAIDs;
        OFF_RNAIDs_minusIsos = setdiff(OFF_RNAIDs,UnDesired_Isoforms);
        ON_IDs = Desired_Isoforms;
        OFF_IDs = OFF_RNAIDs_minusIsos;
        DNA_IDs_1 = find(contains(uniNames,'NC_'));
        DNA_IDs_2 = find(contains(uniNames,'NT_'));%IDs
        DNA_IDs_3 = find(contains(uniNames,'NW_'));%IDs
        DNA_IDs =union(union(DNA_IDs_1,DNA_IDs_2),DNA_IDs_3).';clear DNA_IDs_*
        ProbeDesignResults1.SoftwareResults(id2).uniNames = uniNames;clear uniNames
        ProbeDesignResults1.SoftwareResults(id2).Names = Names;clear Names
        ProbeDesignResults1.SoftwareResults(id2).ON_IDs = ON_IDs;
        ProbeDesignResults1.SoftwareResults(id2).ON_RNAIDs = ON_RNAIDs;
        ProbeDesignResults1.SoftwareResults(id2).OFF_IDs = OFF_IDs;
        ProbeDesignResults1.SoftwareResults(id2).OFF_RNAIDs = OFF_RNAIDs;
        ProbeDesignResults1.SoftwareResults(id2).OFF_RNAIDs_minusIsos = OFF_RNAIDs_minusIsos;
        ProbeDesignResults1.SoftwareResults(id2).NonDNA_IDs = NonDNA_IDs;
        ProbeDesignResults1.SoftwareResults(id2).DNA_IDs = DNA_IDs; 
        ProbeDesignResults1.SoftwareResults(id2).Desired_Isoforms = Desired_Isoforms; 
        ProbeDesignResults1.SoftwareResults(id2).UnDesired_Isoforms = UnDesired_Isoforms; 
        ProbeDesignResults1.SoftwareResults(id2).ProbesWithRibosomalHits = ProbesWithRibosomalHits; 
        ProbeDesignResults1.SoftwareResults(id2).AllowableProbes = AllowableProbes;
        ProbeDesignResultsC.SoftwareResults(id2).ProbeTarget_GC = ProbeTarget_GC;clear ProbeTarget_GC
        ProbeDesignResultsC.SoftwareResults(id2).ProbeTarget_Tm = ProbeTarget_Tm;clear ProbeTarget_Tm
        clear UnDesired_Isoforms Desired_Isoforms DNA_IDs NonDNA_IDs RNA_IDs OFF_RNAIDs_minusIsos OFF_RNAIDs OFF_IDs ON_IDs ON_RNAIDs
        else
        %SL  
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_probes' designerName0 '.mat'],'probes');
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_hits_table' designerName0 '.mat'],'gene_table')  
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_ExpressionInfo' designerName0 '.mat'],'settings');
        load([settings.FolderRootName '/' settings.GeneName '_binding_hits_map' designerName0 '.mat'],'DoesProbeBindSite2'); 
        GeneName = settings.GeneName;
        GeneTarget = settings.transcript_IDs;
        gene_table = sortrows(gene_table,[7 6],'ascend');
        gene_table2 = gene_table(gene_table.Match>=15,:);
        MinusStrandedHits = find(contains(gene_table2.Strand,'Minus'));clear gene_table
        RNA_IDs_1 = find(contains(gene_table2.Name,'NM_'));
        RNA_IDs_2 = find(contains(gene_table2.Name,'NR_'));
        RNA_IDs_3 = find(contains(gene_table2.Name,'XM_'));
        RNA_IDs_4 = find(contains(gene_table2.Name,'XR_'));
        contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
        RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);clear RNA_IDs_*
        gene_table2 = gene_table2(setdiff(1:size(gene_table2,1),RNA_MissedFilteredHits),:);
        gene_table2.Ax = min(gene_table2.SubjectIndices,[],2);clear MinusStrandedHits
        gene_table2.Bx = max(gene_table2.SubjectIndices,[],2);clear RNA_MissedFilteredHits
        gene_table3 = sortrows(gene_table2,[7 13],'ascend');
        Iz = find(endsWith(gene_table3.Name,', ribosomal RNA'));clear gene_table2
        ProbesWithRibosomalHits = unique(gene_table3.ProbeNum(Iz));clear Iz
        AllowableProbes = setdiff(1:size(probes,1),ProbesWithRibosomalHits);
        ProbeDesignResults1.Description = gene_table3.Name(find(contains(gene_table3.Name,settings.rootName),1));
        ProbeDesignResults1.SoftwareResults(id2).DoesProbeBindSite = DoesProbeBindSite2;clear DoesProbeBindSite2
        Names = unique(gene_table3.Name);
        Names = convertCharsToStrings(Names);
        Isoforms = find(contains(Names,GeneName));
        uniNames = extractBefore(Names,'.');
        try
            Desired_Isoforms = find(contains(uniNames,extractBefore(GeneTarget{:},'.')));
        catch
            Desired_Isoforms = find(contains(uniNames,extractBefore(GeneTarget,'.')));  
        end 
        UnDesired_Isoforms = setdiff(Isoforms,Desired_Isoforms);
        ON_RNAIDs = find(strcmp(uniNames,extractBefore(GeneTarget,'.')));
        NonDNA_IDs_1 = find(~contains(uniNames,'NC_'));%IDs
        NonDNA_IDs_2 = find(~contains(uniNames,'NT_'));%IDs
        NonDNA_IDs_3 = find(~contains(uniNames,'NW_'));%
        NonDNA_IDs = intersect(intersect(NonDNA_IDs_1,NonDNA_IDs_2),NonDNA_IDs_3).';clear NonDNA_IDs_*
        Names_NonDNA = Names(NonDNA_IDs);
        PN = gene_table3.ProbeNum(NonDNA_IDs);
        AS = arrayfun(@(x)seqrcomplement(gene_table3.Alignment{NonDNA_IDs(x)}(3,:)),1:length(NonDNA_IDs),'Un',0);
        AS = cellfun(@(x) erase(x,'-'),AS,'Un',0);clear gene_table3 
        ProbeTarget_GC = cell2mat(arrayfun(@(x) oligoprop(AS{x}).GC,1:length(NonDNA_IDs),'Un',0));
        ProbeTarget_Tm = cell2mat(arrayfun(@(x) oligoprop(AS{x},'Salt',settings.SaltConcentration,'Temp',settings.HybridizationTemperature).Tm',1:length(NonDNA_IDs),'Un',0));
        clear AS
        Hit_ON_RNAIDs_Isos = find(contains(Names_NonDNA,GeneName));clear Names_NonDNA
        Hit_OFF_RNAIDs_minusIsos = setdiff(1:length(NonDNA_IDs),Hit_ON_RNAIDs_Isos);
        Hit_ON_IDs = Hit_ON_RNAIDs_Isos;
        Hit_OFF_IDs = Hit_OFF_RNAIDs_minusIsos;
        ProbeDesignResultsC.SoftwareResults(id2).Pset_ON_IDs = @(Q) setdiff(find(ismember(PN,Q)),Hit_ON_IDs);
        ProbeDesignResultsC.SoftwareResults(id2).Pset_OFF_IDs = @(Q) setdiff(find(ismember(PN,Q)),Hit_OFF_IDs);
        OFF_RNAIDs = setdiff(NonDNA_IDs,ON_RNAIDs);
        OFF_IDs = OFF_RNAIDs;
        OFF_RNAIDs_minusIsos = setdiff(OFF_RNAIDs,UnDesired_Isoforms);
        ON_IDs = Desired_Isoforms;
        OFF_IDs = OFF_RNAIDs_minusIsos;
        DNA_IDs_1 = find(contains(uniNames,'NC_'));
        DNA_IDs_2 = find(contains(uniNames,'NT_'));%IDs
        DNA_IDs_3 = find(contains(uniNames,'NW_'));%IDs
        DNA_IDs =union(union(DNA_IDs_1,DNA_IDs_2),DNA_IDs_3).';clear DNA_IDs_*
        ProbeDesignResults1.SoftwareResults(id2).uniNames = uniNames;clear uniNames
        ProbeDesignResults1.SoftwareResults(id2).Names = Names;clear Names
        ProbeDesignResults1.SoftwareResults(id2).ON_IDs = ON_IDs;
        ProbeDesignResults1.SoftwareResults(id2).ON_RNAIDs = ON_RNAIDs;
        ProbeDesignResults1.SoftwareResults(id2).OFF_IDs = OFF_IDs;
        ProbeDesignResults1.SoftwareResults(id2).OFF_RNAIDs = OFF_RNAIDs;
        ProbeDesignResults1.SoftwareResults(id2).OFF_RNAIDs_minusIsos = OFF_RNAIDs_minusIsos;
        ProbeDesignResults1.SoftwareResults(id2).NonDNA_IDs = NonDNA_IDs;
        ProbeDesignResults1.SoftwareResults(id2).DNA_IDs = DNA_IDs; 
        ProbeDesignResults1.SoftwareResults(id2).Desired_Isoforms = Desired_Isoforms; 
        ProbeDesignResults1.SoftwareResults(id2).UnDesired_Isoforms = UnDesired_Isoforms; 
        ProbeDesignResults1.SoftwareResults(id2).ProbesWithRibosomalHits = ProbesWithRibosomalHits; 
        ProbeDesignResults1.SoftwareResults(id2).AllowableProbes = AllowableProbes;
        ProbeDesignResultsC.SoftwareResults(id2).ProbeTarget_GC = ProbeTarget_GC;clear ProbeTarget_GC
        ProbeDesignResultsC.SoftwareResults(id2).ProbeTarget_Tm = ProbeTarget_Tm;clear ProbeTarget_Tm
        clear UnDesired_Isoforms Desired_Isoforms DNA_IDs NonDNA_IDs RNA_IDs OFF_RNAIDs_minusIsos OFF_RNAIDs OFF_IDs ON_IDs ON_RNAIDs        
        end
        catch    
        end
    end
    save([settings.FolderRootName '/' inputs1{gene_num,5} '_' num2str(HybridizationTemperature) '_DesignResults1.mat'],'ProbeDesignResults1','-v7.3')
CarryOver.ON_RNAIDs = {ProbeDesignResults1.SoftwareResults(:).ON_RNAIDs};
CarryOver.OFF_RNAIDs = {ProbeDesignResults1.SoftwareResults(:).OFF_RNAIDs};
CarryOver.AllowableProbes = {ProbeDesignResults1.SoftwareResults(:).AllowableProbes};
CarryOver.ProbesWithRibosomalHits = {ProbeDesignResults1.SoftwareResults(:).ProbesWithRibosomalHits};
CarryOver2.ON_IDs = {ProbeDesignResults1.SoftwareResults(:).ON_IDs};
CarryOver2.OFF_IDs = {ProbeDesignResults1.SoftwareResults(:).OFF_IDs};
CarryOver2.DNA_IDs = {ProbeDesignResults1.SoftwareResults(:).DNA_IDs};
CarryOver2.NonDNA_IDs = {ProbeDesignResults1.SoftwareResults(:).NonDNA_IDs};
CarryOver2.UnDesired_Isoforms = {ProbeDesignResults1.SoftwareResults(:).UnDesired_Isoforms};
CarryOver3.ProbeTarget_Tm = {ProbeDesignResultsC.SoftwareResults(:).ProbeTarget_Tm};
CarryOver3.ProbeTarget_GC = {ProbeDesignResultsC.SoftwareResults(:).ProbeTarget_GC};
CarryOver3.Pset_ON_IDs = {ProbeDesignResultsC.SoftwareResults(:).Pset_ON_IDs};
CarryOver3.Pset_OFF_IDs = {ProbeDesignResultsC.SoftwareResults(:).Pset_OFF_IDs};
save([settings.FolderRootName '/' inputs1{gene_num,5} '_' num2str(HybridizationTemperature) '_DesignResults_CarryOver.mat'],'CarryOver','CarryOver2','CarryOver3','-v7.3')




end


