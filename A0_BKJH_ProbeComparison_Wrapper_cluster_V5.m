function A0_BKJH_ProbeComparison_Wrapper_cluster_V5(id,cluster,id2)
% Add Yeast Expression, update to use GTF
% Update that way more species can be used and alignment of database files
% and GTF with reference genome, and add custom gene expression file
%Update so that more species can be included in list for their blast database
%Update for using Python BLAST Local for MATLAB 2021+
%stellaris should probability be updated using reversesequence match when

%Run for predefined set of probe sequences
%Specify Genes, [Expr], [Tfixed, Topt, Tgradient] & Application
%Modes [Single Isoform, Multiple Isoforms, Multile Gene Isoforms]
%Specify Data
%RNA/Protein Binding Data [mask model]
%CHIP/TT-seq data
%Genome/Transcriptome
%SNP Data
%Expr Data
%GTF/Annotation Files [exons,introns, name, position, etc.]
%Generates Primary Probes
%BKJH_Probe_Generator([Lmin Lmax],inputs1{gene_num,1},inputs1{gene_num,2},inputs1{gene_num,3},inputs1{gene_num,8},inputs1{gene_num,9},settings.isOffline,settings.SEQdbRoot);
%BKJH_MultiplexProbe_Generator([Lmin Lmax],inputs1{gene_num,1},inputs1{gene_num,2},inputs1{gene_num,3},inputs1{gene_num,8},inputs1{gene_num,9},settings.isOffline,settings.SEQdbRoot);
%Generate Secondary Probes
%BLAST Primary Probes
%BLAST Secondary Probes
%Get Expression Data
%Get Thermodynamic EQ Info [T,Distinct ON/OFF-Targets, Models]
%Get Thermodynamic kf/kfr Info [T,Distinct ON/OFF-Targets, Models]
%Get Binding Site Map and Location Mapping
%Get Landscape Properties
%Get Nascent RNA Molecules
%Get RNA Secondary Structure
%Get Optimization Function /Simulator
%Probe Designer
%Get Probe Secondary Structure
%Get Metrics
%EQ Model and Kinetic Model
%Generate Predictions
%Predict Spot Detection
%Print Final Probes
%Output Predictions/Metrics
%Output Figures
%Aggregate Results Across Probe Sets
%Cross-Software Figures
%Aggregate Results Across Genes/Isoforms
%Cross-Gene Cross-Software Results
%outline of changes
%Update FolderName to not reference temperature
%size of every variable in byyrd
%update t3emp variation using two point encoding. just need dH, dS, dCp to
%define dG for all T. for different models.
%code that stores location info in 3d map [SiteLoc]
%code that re-assigns probe and target ID in 3d-4d maps
% code that combines maps
% code that designs probes [mix designer with V7_Out codes]
% solves eq/kinetics/metrics with multi-transcript output
% all forked with/without temp changes or gradients.
%update to have higher dimensions for versions of code/ cell array
%might turn 4D multiTargetLocations into two 3d matricies and change
%combineMaps code to for those
%% Genes To Design Probes For
inputs1 = {...
    {'NM_001660.4'},{},{}, 'Human','(ARF4)','(ARF4)','3',{},1 ;... %
    {'NM_000805.5'},{},{}, 'Human','(GAST)','(GAST)','17',{},1 ;...          %ENST00000329402.4 ,  465bp, 1-Iso
    {'NR_002844.2'},{},{}, 'Mice','(TSIX)','(TSIX)','X',{},1 ;... %4306
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

%% Settings Specification
%Neuert Lab initial Probes,simulated annealing Probes, Old Designer Probes,
% Stellaris, Oligostan, MERFISH, PaintSHOP,
switch id2
    case 1
        designerName = '_oldNeuertLab';
        minProbeSize = 20;
        maxProbeSize = 20;
    case 2
        designerName = '_Stellaris';
        minProbeSize = 20;%only plot stellaris groups with distinct sets of probes designed (group together ones with same)
        maxProbeSize = 20;%only L1-L5, then first with 25+ from L5 to L1. do they pick diff probes or just subset of L1
    case 3
        designerName = '_OligoStanHT';
        minProbeSize = 23;
        maxProbeSize = 38;
    case 4
        designerName = '_PaintSHOP';
        minProbeSize = 30;
        maxProbeSize = 37;
        PaintSHOPfile = 'hg38_refseq_newBalance.tsv';%iso-resolved
        %PaintSHOPfile = 'hg38_iso_refseq_newBalance.tsv';%iso-flattened
    case 5
        designerName = '_ZZLMERFISH';%diff levels
        minProbeSize = 30;
        maxProbeSize = 30;
    case 6
        designerName = '_AIBSMERFISH';%diff levels
        minProbeSize = 30;
        maxProbeSize = 30;
    case 7
        designerName = '_LengthOptimized';%diff levels
        minProbeSize = 20;
        maxProbeSize = 38;
end
referenceTypes = {'EMBL','RefSeq'};
saveRoot = strcat('output',filesep);
BLASTpath = 'src/thirdparty/blast-2.8.1+/bin/blastn';
customBlastDatabase_DNA = 'N/A';
customBlastDatabase_RNA = 'N/A';
SingleOverMultiplex = 1;
AllIsoforms = 0;
gene_num = id;
referenceType = 'RefSeq';
max_probes = 96;
MininumProbeSpacing = 3;
Nmodel = 4;
cellPreset = 1;
HybridizationTemperature = 37;
SaltConcentration = 0.05;
RemoveMisMatches = 1;
SpecificityThreshold = 2;
DecisionAlgorithmFloorSize = 0.5;
RemoveRibosomalHits = 1;
RunOffline = 1;
withNascentTranscripts = 0;
MinHomologySearchTargetSize = 15;
UseGeneOverTranscLevelExpression = 0;
nullRNAcopynumber = 100;
nullDNAcopynumber = 2;
BLASTrna = 1;
BLASTdna = 0;
batchSize = 20;
batchSize = 10;
targetBatchSize = 200;
ParsingPreference = 1;
addSelfProb = 1;
packOptimal = 1;
targetTypes = [1 0];
removeUndesiredIsos=1;

if (strcmp(inputs1{gene_num,4},'Yeast'))
    DoAllGenesHaveSameExpression = 1;%same expr.
else
    DoAllGenesHaveSameExpression = 0;%diff expr
end


%Save Settings
settings.FolderRootName = strcat(saveRoot,inputs1{gene_num,5},'_',strjoin(inputs1{gene_num,1},'_'));
settings.rootName = strjoin(inputs1{gene_num,1},'_');
settings.designerName = designerName;
addSelfProb = 1;packOptimal = 1;
targetTypes = [1 0];removeUndesiredIsos=1;
CellTypeID = settings.CellType_ExprID;
expValType = settings.expressionValType;% 1-4 (expCounts,expValues,mean(CellTypeExpValues),one cell types CellTypeExpValues)

%Design Target Info
settings.Organism = inputs1{gene_num,4};
settings.GeneName = inputs1{gene_num,5};
settings.ChrNum = inputs1{gene_num,7};
settings.GeneChr = strcat('chr',inputs1{gene_num,7});
settings.transcript_IDs = inputs1{gene_num,1};
settings.ProbeSpacing = MininumProbeSpacing;
settings.RemoveMisMatches = RemoveMisMatches;
settings.SaltConcentration = SaltConcentration;
settings.HybridizationTemperature = HybridizationTemperature;
settings.withNascent = withNascentTranscripts;
settings.BLASTdna = BLASTdna;
settings.BLASTrna = BLASTrna;
settings.BLASTbatchSize = batchSize;
settings.BLASTsimultaneousParsingOverSequentialParsing = ParsingPreference;
settings.BLASTpath = BLASTpath;
if (length(inputs1{gene_num,1})==1)
        %single
        refInfo = inputs1{gene_num,1}{1}(1:2);
        if (ismember(refInfo,{'NR','XR','NM','XM'}))
        settings.referenceType = 'RefSeq';
        elseif (strcmp(refInfo,'EN'))
        settings.referenceType = 'ENSEMBL';
        else
        end
else
        % multiple

end
%Type of Databases     
if (strcmp(settings.referenceType,'RefSeq'))
settings.hLocRoot = 'data/DatabaseData/Blast_Databases/Human/NCBI_RefSeq/';
settings.hLoc = 'data/DatabaseData/Blast_Databases/Human/NCBI_RefSeq/Human_NCBI_genomic';
settings.hLoc2 = 'data/DatabaseData/Blast_Databases/Human/NCBI_RefSeq/Human_NCBI_transcript';
settings.mLocRoot =  'data/DatabaseData/Blast_Databases/Mouse/NCBI_RefSeq/';
settings.mLoc = 'data/DatabaseData/Blast_Databases/Mouse/NCBI_RefSeq/Mouse_NCBI_genomic';
settings.mLoc2 = 'data/DatabaseData/Blast_Databases/Mouse/NCBI_RefSeq/Mouse_NCBI_transcript';
settings.yLocRoot =  'data/DatabaseData/Blast_Databases/Yeast/NCBI_RefSeq/';
settings.yLoc =  'data/DatabaseData/Blast_Databases/Yeast/NCBI_RefSeq/Yeast_NCBI_genomic';
settings.yLoc2 =  'data/DatabaseData/Blast_Databases/Yeast/NCBI_RefSeq/Yeast_NCBI_transcript';
elseif (strcmp(settings.referenceType,'ENSEMBL'))
settings.hLocRoot = 'data/DatabaseData/Blast_Databases/Human/EMBL_EBI/';
settings.hLoc = 'data/DatabaseData/Blast_Databases/Human/EMBL_EBI/Human_ENSEMBL_genomic';%chrom number, MT, X, Y and scaffold KI and GL
settings.hLoc2 = 'data/DatabaseData/Blast_Databases/Human/EMBL_EBI/Human_ENSEMBL_transcript';%cdna and ncrna
settings.mLocRoot =  'data/DatabaseData/Blast_Databases/EMBL_EBI/';
settings.mLoc = 'data/DatabaseData/Blast_Databases/Mouse/EMBL_EBI/Mouse_ENSEMBL_genomic';
settings.mLoc2 = 'data/DatabaseData/Blast_Databases/Mouse/EMBL_EBI/Mouse_ENSEMBL_transcript';%ENMUST, cdna and ncRNA
settings.yLocRoot =  'data/DatabaseData/Blast_Databases/Yeast/EBML_EBI/';%SGD
settings.yLoc =  'data/DatabaseData/Blast_Databases/Yeast/EMBL_EBI/Yeast_ENSEMBL_genomic';%chromsome roman numerals and Mito
settings.yLoc2 =  'data/DatabaseData/Blast_Databases/Yeast/EMBL_EBI/Yeast_ENSEMBL_transcript';%SGD cdna and ncRNA
elseif (strcmp(settings.referenceType,'Gencode'))
settings.hLocRoot = 'data/DatabaseData/Blast_Databases/Human/GENCODE/';%
settings.hLoc = 'data/DatabaseData/Blast_Databases/Human/GENCODE/Human_GENCODE_genomic';%chrM chrX chrY
settings.hLoc2 = 'data/DatabaseData/Blast_Databases/Human/GENCODE/Human_GENCODE_transcript';%ENST and ENSG
settings.mLocRoot =  'data/DatabaseData/Blast_Databases/Mouse/GENCODE/';%ENSMUST... | ENSMUSG ....|OTTOMUST OTTOMUSG  |Name protein coding]
settings.mLoc = 'data/DatabaseData/Blast_Databases/Mouse/GENCODE/Mouse_GENCODE_genomic';%chr1,... chrX, chrY, chrM
settings.mLoc2 = 'data/DatabaseData/Blast_Databases/Mouse/GENCODE/Mouse_GENCODE_transcript';
settings.yLocRoot =  'data/DatabaseData/Blast_Databases/Yeast/YGD/';    
settings.yLoc =  'data/DatabaseData/Blast_Databases/Yeast/YGD/Yeast_SGD_genomic';
settings.yLoc2 =  'data/DatabaseData/Blast_Databases/Yeast/YGD/Yeast_SGD_transcript';
end
settings.TargetBatchSize = targetBatchSize;
settings.SingleOrMulti = SingleOverMultiplex;
settings.AllIsoforms = AllIsoforms;
settings.MinProbeSize = minProbeSize;
settings.MaxProbeSize = maxProbeSize;
settings.MinHomologySearchTargetSize = MinHomologySearchTargetSize;

%Gene Expression Parameters
settings.DoAllGenesHaveSameExpression = DoAllGenesHaveSameExpression;
settings.HumanSpecific.HumanExpGeneOrTransc = UseGeneOverTranscLevelExpression; % 1 (Gene/EMBL GENEID) , 0 (Transcript/EMBL Transcript ID)
settings.UseRegularDNAExpression = 1;%0 use DNA expression from gene expression track in expression data, 1 set expression to 2 for DNA.
settings.UniformRNAExpression = nullRNAcopynumber;%if assuming no differences in gene's expression sets level.
settings.DNAPloidy = nullDNAcopynumber;
settings.otherBlastDatabase = customBlastDatabase_DNA;
settings.otherBlastDatabase2 = customBlastDatabase_RNA;

%Cluster Parameters
settings.clusterStatus = cluster;
settings.isOffline = RunOffline;

%Selecting Probes Specifications
settings.maxProbes = max_probes;
settings.SpecificityThreshold = SpecificityThreshold;
settings.ChoosingFloorStepSize = DecisionAlgorithmFloorSize;
settings.RemoveProbesWithRibosomalHits = RemoveRibosomalHits;


switch cellPreset
    case 1%Average Tissue
        settings.expressionValType = 1;% 1-4 (expCounts,expValues,mean(CellTypeExpValues),one cell types CellTypeExpValues)
        settings.HumanSpecific.Ontology = 'Normal';%Cancer or Normal
        settings.HumanSpecific.TissueOrTissueAndCellType = 0; % 0 (TissueOnly) 1 (Tissue/Cell Type)
        settings.HumanSpecific.SCOutputTrack = 1;
        settings.HumanSpecific.SCTracks = 1;
        settings.CellType_ExprID = 1;
    case 2%Jurkat/T lymphocyte (NK cells, CD4 T cells, B cells, white blood cell)
        settings.expressionValType = 1;
        settings.HumanSpecific.Ontology = 'Normal';%Cancer or Normal
        settings.HumanSpecific.TissueOrTissueAndCellType = 1; % 0 (TissueOnly) 1 (Tissue/Cell Type)
        settings.HumanSpecific.SCOutputTrack = 13;%13 Kidney, 4 PBMC
        settings.HumanSpecific.SCTracks = 13;
        settings.CellType_ExprID = 3;%3 CD4 T cell, 2 PBMC CD4 T cell
    case 3%THP1 Kidney
        settings.expressionValType = 1;
        settings.HumanSpecific.Ontology = 'Normal';%Cancer or Normal
        settings.HumanSpecific.TissueOrTissueAndCellType = 1; % 0 (TissueOnly) 1 (Tissue/Cell Type)
        settings.HumanSpecific.SCOutputTrack = 13;%13 Kidney, 4 PBMC
        settings.HumanSpecific.SCTracks = 13;
        settings.CellType_ExprID = 11;%11 kidney mononuclear phagocyte, 5 PBMC monocyte
    case 4%PBMC
        settings.expressionValType = 4;
        settings.HumanSpecific.Ontology = 'Normal';%Cancer or Normal
        settings.HumanSpecific.TissueOrTissueAndCellType = 1; % 0 (TissueOnly) 1 (Tissue/Cell Type)
        settings.HumanSpecific.SCOutputTrack = 4;
        settings.HumanSpecific.SCTracks = 4;
        %settings.CellType_ExprID = 1;
    case 5%Normal Colon
        settings.expressionValType = 4;
        settings.HumanSpecific.Ontology = 'Normal';%Cancer or Normal
        settings.HumanSpecific.TissueOrTissueAndCellType = 1; % 0 (TissueOnly) 1 (Tissue/Cell Type)
        settings.HumanSpecific.SCOutputTrack = 1;
        settings.HumanSpecific.SCTracks = 1;
        %settings.CellType_ExprID = 1;
    case 6%Colon Cancer Cells
        settings.expressionValType = 4;
        settings.HumanSpecific.Ontology = 'Cancer';%Cancer or Normal
        settings.CellType_ExprID = 8;
        %settings.HumanSpecific.TissueOrTissueAndCellType = 0; % 0 (TissueOnly) 1 (Tissue/Cell Type)
        %settings.HumanSpecific.SCOutputTrack = 8;
        %settings.HumanSpecific.SCTracks = 8;
    case 7%Pancreatic Islet Alpha Cell
        settings.expressionValType = 4;
        settings.HumanSpecific.Ontology = 'Normal';%Cancer or Normal
        settings.HumanSpecific.TissueOrTissueAndCellType = 1; % 0 (TissueOnly) 1 (Tissue/Cell Type)
        settings.HumanSpecific.SCOutputTrack = 22;
        settings.HumanSpecific.SCTracks = 22;
        settings.CellType_ExprID = 3;
    case 8%Pancreatic Islet Beta Cell
        settings.expressionValType = 4;
        settings.HumanSpecific.Ontology = 'Normal';%Cancer or Normal
        settings.HumanSpecific.TissueOrTissueAndCellType = 1; % 0 (TissueOnly) 1 (Tissue/Cell Type)
        settings.HumanSpecific.SCOutputTrack = 22;
        settings.HumanSpecific.SCTracks = 22;
        settings.CellType_ExprID = 4;
    case 9%Pancreatic Islet Delta Cell
        settings.expressionValType = 4;
        settings.HumanSpecific.Ontology = 'Normal';%Cancer or Normal
        settings.HumanSpecific.TissueOrTissueAndCellType = 1; % 0 (TissueOnly) 1 (Tissue/Cell Type)
        settings.HumanSpecific.SCOutputTrack = 22;
        settings.HumanSpecific.SCTracks = 22;
        settings.CellType_ExprID = 5;
end
settings.MouseExpressionFile = 'DatabaseData/tabulamuris_barChart.bed';
settings.ExpressionVariableNames = {'chrom','chromStart','chromEnd','name',...
    'score','strand','name2','expCount','expScores','dataOffset','dataLen'};
settings.ExpressionVariableNames2 = {'chrom','chromStart','chromEnd','name',...
    'score','strand','geneId','geneType','expCount','expScores'};
settings.HumanHPAcellLine_TranscExpressionFile = 'data/DatabaseData/rna_celline.tsv';
settings.HumanTCGA_TranscExpressionFile = 'data/DatabaseData/tcgaTranscExpr.bed';
settings.HumanTCGA_GeneExpressionFile = 'data/DatabaseData/tcgaGeneExpr.bed';
settings.HumanGTEX_TranscExpressionFile = 'data/DatabaseData/gtexTranscExpr.bed';
settings.HumanGTEX_GeneExpressionFile = 'data/DatabaseData/gtexGeneV8.txt';
settings.Human_wgEncodeGencodeRefSeqFile = 'data/DatabaseData/wgEncodeGencodeRefSeqV44.txt';
settings.Human_wgEncodeGencodeAttributesFile = 'data/DatabaseData/wgEncodeGencodeAttrsV44.txt';
settings.Human_wgEncodeGencodeCompFile = 'data/DatabaseData/wgEncodeGencodeCompV44.txt';
settings.wgEncodeGencodeAttributesVariableNames = {'geneId','geneName','geneType','geneStatus',...
    'transcriptId','transcriptName','transcriptType','transcriptStatus',...
    'havanaGeneId','havanaTranscriptId','ccdsId','level','transcriptClass','proteinId'};
settings.wgEncodeGencodeCompVariableNames = {'nx','transcriptId','chrom','strand','txStart',...
    'txEnd','cdsStart','cdsEnd','exonCount',...
    'exonStarts','exonEnds','score','name2','cdsStartStat','cdsEndStat','exonFrames'};
settings.Human_GencodeRefSeqMetadataFile = 'data/DatabaseData/gencode.v44.metadata.RefSeq';
settings.Mouse_wgEncodeGencodeRefSeqFile = 'data/DatabaseData/wgEncodeGencodeRefSeqVM25.txt';
settings.Mouse_GencodeRefSeqMetadataFile = 'data/DatabaseData/gencode.vM30.metadata.RefSeq';
settings.Mouse_wgEncodeGencodeAttributesFile = 'data/DatabaseData/wgEncodeGencodeAttrsVM25.txt';
settings.Mouse_wgEncodeGencodeCompFile = 'data/DatabaseData/wgEncodeGencodeCompVM25.txt';
settings.Custom_TranscExpressionFile = 'N/A';
settings.CustomTranscExpressionVariableNames ='N/A';
settings.Custom_GeneExpressionFile = 'N/A';
settings.CustomGeneExpressionVariableNames = 'N/A';
settings.CustomExpGeneOrTransc = 1;
settings.Custom_wgEncodeGencodeRefSeqFile = 'N/A';
settings.Custom_GencodeRefSeqMetadataFile = 'N/A';
settings.Custom_wgEncodeGencodeAttributesFile = 'N/A';
settings.Custom_wgEncodeGencodeCompFile = 'N/A';
settings.YeastExpressionFile = 'N/A';
settings.Yeast_wgEncodeGencodeRefSeqFile = 'N/A';
settings.Yeast_GencodeRefSeqMetadataFile = 'N/A';
settings.Yeast_wgEncodeGencodeAttributesFile = 'N/A';
settings.Yeast_wgEncodeGencodeCompFile = 'N/A';
settings.scTracksBedFiles = {'tabulasapiens_tissue_cell_type','colonWang_cell_type','ileumWang_cell_type','rectumWang_cell_type'};
settings.HumanGenomeAssemblyReportFile = 'data/Metadata/Human/GCF_000001405.39_GRCh38.p13_assembly_report.txt';
settings.MouseGenomeAssemblyReportFile = 'data/Metadata/Mouse/GCA_000001635.9_GRCm39_assembly_report.txt';
settings.YeastGenomeAssemblyReportFile = 'data/Metadata/Yeast/S288C/GCF_000146045.2_R64_assembly_report.txt';
settings.CustomGenomeAssemblyReportFile = 'N/A';
settings.Custom_X_ChromNumber = 1;
settings.Custom_Y_ChromNumber = 1;
settings.Custom_MT_ChromNumber = 1;
settings.scTracks = {'colonWangCellType'};  %
'tabulasapiens_tissue_cell_type';
%Save Settings
settings.FolderRootName = strcat(saveRoot,inputs1{gene_num,5},'_',strjoin(inputs1{gene_num,1},'_'));
settings.rootName = strjoin(inputs1{gene_num,1},'_');

%Probe Properties Specifications
settings.AddOneOrZeroInNormalization = 1;%In GetProbeProperties

%RNA Secondary Structure Specification
settings.SolveStructure = 0;
settings.SecondaryStructureFileRoot = '/data/DatabaseData/dbnFiles/';

%Selecting Probes Specifications
settings.maxProbes = max_probes;
settings.SpecificityThreshold = SpecificityThreshold;
settings.ChoosingFloorStepSize = DecisionAlgorithmFloorSize;
settings.RemoveProbesWithRibosomalHits = RemoveRibosomalHits;
settings.designerName = designerName;
if (strcmp(settings.Organism,'Human'))
    settings.SEQdbRoot = settings.hLocRoot;
elseif (strcmp(settings.Organism,'Mouse'))
    settings.SEQdbRoot = settings.mLocRoot;
elseif (strcmp(settings.Organism,'Yeast'))
    settings.SEQdbRoot = settings.yLocRoot;
else
    settings.SEQdbRoot = settings.otherBlastDatabase;
end
T_hybrid = HybridizationTemperature;
Lmin = minProbeSize; Lmax = maxProbeSize;

%% Load Annotation File
if (strcmp(settings.Organism,'Human'))
    optsUCSC = detectImportOptions(settings.Human_wgEncodeGencodeRefSeqFile);
    optsUCSC.VariableNames = {'Var1','Var2','Var3'};
    Gencode_db = readtable(settings.Human_GencodeRefSeqMetadataFile,optsUCSC);
    optsUCSC2 = detectImportOptions(settings.Human_wgEncodeGencodeAttributesFile);
    optsUCSC2.VariableNames = settings.wgEncodeGencodeAttributesVariableNames;
    EMBLAttrAlign_db = readtable(settings.Human_wgEncodeGencodeAttributesFile,optsUCSC2);
    optsUCSC3 = detectImportOptions(settings.Human_wgEncodeGencodeCompFile);
    optsUCSC3.VariableNames = settings.wgEncodeGencodeCompVariableNames;
    EMBLComp_db = readtable(settings.Human_wgEncodeGencodeCompFile,optsUCSC3);
elseif (strcmp(settings.Organism,'Mouse'))
    optsUCSC = detectImportOptions(settings.Mouse_wgEncodeGencodeRefSeqFile);
    optsUCSC.VariableNames = {'Var1','Var2','Var3'};
    Gencode_db = readtable(settings.Mouse_GencodeRefSeqMetadataFile,optsUCSC);
    optsUCSC2 = detectImportOptions(settings.Mouse_wgEncodeGencodeAttributesFile);
    optsUCSC2.VariableNames =settings.wgEncodeGencodeAttributesVariableNames;
    EMBLAttrAlign_db = readtable(settings.Mouse_wgEncodeGencodeCompFile,optsUCSC2);
    optsUCSC3 = detectImportOptions(settings.Mouse_wgEncodeGencodeCompFile);
    optsUCSC3.VariableNames = settings.wgEncodeGencodeCompVariableNames;
    EMBLComp_db = readtable(settings.Mouse_wgEncodeGencodeCompFile,optsUCSC3);
elseif (strcmp(settings.Organism,'Yeast'))
    optsUCSC = detectImportOptions(settings.Yeast_wgEncodeGencodeRefSeqFile);
    optsUCSC.VariableNames = {'Var1','Var2','Var3'};
    Gencode_db = readtable(settings.Yeast_GencodeRefSeqMetadataFile,optsUCSC);
    optsUCSC2 = detectImportOptions(settings.Yeast_wgEncodeGencodeAttributesFile);
    optsUCSC2.VariableNames = settings.wgEncodeGencodeAttributesVariableNames;
    EMBLAttrAlign_db = readtable(settings.Yeast_wgEncodeGencodeAttributesFile,optsUCSC2);
    optsUCSC3 = detectImportOptions(settings.Yeast_wgEncodeGencodeCompFile);
    optsUCSC3.VariableNames = settings.wgEncodeGencodeCompVariableNames;
    EMBLComp_db = readtable(settings.Yeast_wgEncodeGencodeCompFile,optsUCSC3);
else  %custom organism
    optsUCSC = detectImportOptions(settings.Custom_wgEncodeGencodeRefSeqFile);
    optsUCSC.VariableNames = {'Var1','Var2','Var3'};
    Gencode_db = readtable(settings.Custom_GencodeRefSeqMetadataFile,optsUCSC);
    optsUCSC2 = detectImportOptions(settings.Custom_wgEncodeGencodeAttributesFile);
    optsUCSC2.VariableNames = settings.wgEncodeGencodeAttributesVariableNames;
    EMBLAttrAlign_db = readtable(settings.Custom_wgEncodeGencodeAttributesFile,optsUCSC2);
    optsUCSC3 = detectImportOptions(settings.Custom_wgEncodeGencodeCompFile);
    optsUCSC3.VariableNames = settings.wgEncodeGencodeCompVariableNames;
    EMBLComp_db = readtable(settings.Custom_wgEncodeGencodeCompFile,optsUCSC3);
end

%% Generate Folder
if (settings.SingleOrMulti==1&&settings.AllIsoforms == 0)%One Gene/One Isoform
    if (not(isfolder([inputs1{gene_num,5} '_' settings.rootName])))
        mkdir([inputs1{gene_num,5} '_' settings.rootName])
    end
    settings.FolderName = [inputs1{gene_num,5} '_' settings.rootName];
elseif (settings.SingleOrMulti==1&&settings.AllIsoforms == 1)%One Gene/All Isoform
    ENST_ID = Gencode_db.Var1{strcmp(inputs1{gene_num,1}{1},Gencode_db.Var2)};
    ENST_GeneName = EMBLAttrAlign_db.geneName{strcmp(ENST_ID,EMBLAttrAlign_db.transcriptId)};
    ENST_GeneIDs = EMBLAttrAlign_db.transcriptId(strcmp(ENST_GeneName,EMBLAttrAlign_db.geneName));
    EMBL_ChrNums = arrayfun(@(x) extractAfter(EMBLComp_db.chrom{strcmp(ENST_GeneIDs{x},EMBLComp_db.transcriptId)},'chr'),1:length(ENST_GeneIDs),'Un',0);
    ENST_RefSeqIDs = cell(1,size(ENST_GeneIDs,1));
    for k = 1:size(ENST_GeneIDs,1)
        if (sum(strcmp(ENST_GeneIDs{k},Gencode_db.Var1))>0)
            ENST_RefSeqIDs{k} = Gencode_db.Var2(find(strcmp(ENST_GeneIDs{k},Gencode_db.Var1)));
        end
    end
    ENST_RefSeqIDz = unique(vertcat(ENST_RefSeqIDs{:}));
    Exists_Match = find(~cellfun(@isempty,ENST_RefSeqIDs));
    match_Found = arrayfun(@(x) find(arrayfun(@(y) sum(contains(ENST_RefSeqIDs{y},ENST_RefSeqIDz{x})),find(~cellfun(@isempty,ENST_RefSeqIDs))),1),1:length(ENST_RefSeqIDz));
    Isoform_chrom = {EMBL_ChrNums{Exists_Match(match_Found)}};
    FoldName = cell(1,size(ENST_RefSeqIDz,1));
    for k = 1:size(ENST_RefSeqIDz,1)
        if (not(isfolder([inputs1{gene_num,5} '_' ENST_RefSeqIDz{k}])))
            mkdir([inputs1{gene_num,5} '_' ENST_RefSeqIDz{k}])
        end
        FoldName{k} = [inputs1{gene_num,5} '_' ENST_RefSeqIDz{k}];
    end
    settings.rootName = strjoin(ENST_RefSeqIDz,'_');
    if (not(isfolder([inputs1{gene_num,5} '_' settings.rootName])))
        mkdir([inputs1{gene_num,5} '_' settings.rootName])
    end
    All_chrom = Isoform_chrom;
    settings.FolderName = [inputs1{gene_num,5} '_' settings.rootName];
elseif (settings.SingleOrMulti==0&&settings.AllIsoforms == 0)%Multi Gene/One Isoform
    %multiplex folder name
    Multi_GeneName = cell(1,size(inputs1{gene_num,1},2));
    Multi_GeneENST_ID = cell(1,size(inputs1{gene_num,1},2));
    FoldName = cell(1,size(inputs1{gene_num,1},2));
    for k = 1:size(inputs1{gene_num,1},2)
        ENST_ID = Gencode_db.Var1{strcmp(inputs1{gene_num,1}{k},Gencode_db.Var2)};
        Multi_GeneName{k} = EMBLAttrAlign_db.geneName{strcmp(ENST_ID,EMBLAttrAlign_db.transcriptId)};
        Multi_GeneENST_ID{k} = ENST_ID;
    end
    Transcript_chrom = arrayfun(@(x) extractAfter(EMBLComp_db.chrom{strcmp(Multi_GeneENST_ID{x},EMBLComp_db.transcriptId)},'chr'),1:length(Multi_GeneENST_ID),'Un',0);
    for k = 1:size(Multi_GeneName,2)
        if (not(isfolder(['(' Multi_GeneName{k} ')_' inputs1{gene_num,1}{k}])))
            mkdir(['(' Multi_GeneName{k} ')_' inputs1{gene_num,1}{k}])
        end
        FoldName{k} = ['(' Multi_GeneName{k} ')_' inputs1{gene_num,1}{k}];
    end
    if (not(isfolder(strjoin(FoldName,'_'))))
        mkdir(strjoin(FoldName,'_'))
    end
    All_chrom = Transcript_chrom;
    settings.FolderName = strjoin(FoldName,'_');
elseif (settings.SingleOrMulti==0&&settings.AllIsoforms == 1)%Multi Gene/All Isoform
    Multi_GeneName = cell(1,size(inputs1{gene_num,1},2));
    Multi_ENST_RefSeqIDz = cell(1,size(inputs1{gene_num,1},2));
    TranscriptIsoform_chrom = cell(1,size(inputs1{gene_num,1},2));
    for k = 1:size(inputs1{gene_num,1},2)
        ENST_ID = Gencode_db.Var1{strcmp(inputs1{gene_num,1}{k},Gencode_db.Var2)};
        Multi_GeneName{k} = EMBLAttrAlign_db.geneName{strcmp(ENST_ID,EMBLAttrAlign_db.transcriptId)};
        ENST_GeneName = EMBLAttrAlign_db.geneName{strcmp(ENST_ID,EMBLAttrAlign_db.transcriptId)};
        ENST_GeneIDs = EMBLAttrAlign_db.transcriptId(strcmp(ENST_GeneName,EMBLAttrAlign_db.geneName));
        EMBL_ChrNums = arrayfun(@(x) extractAfter(EMBLComp_db.chrom{strcmp(ENST_GeneIDs{x},EMBLComp_db.transcriptId)},'chr'),1:length(ENST_GeneIDs),'Un',0);
        ENST_RefSeqIDs = cell(1,size(ENST_GeneIDs,1));
        for k2 = 1:size(ENST_GeneIDs,1)
            if (sum(strcmp(ENST_GeneIDs{k2},Gencode_db.Var1))>0)
                ENST_RefSeqIDs{k2} = Gencode_db.Var2(find(strcmp(ENST_GeneIDs{k2},Gencode_db.Var1)));
            end
        end
        Multi_ENST_RefSeqIDz{k} = unique(vertcat(ENST_RefSeqIDs{:}));
        Exists_Match = find(~cellfun(@isempty,ENST_RefSeqIDs));
        match_Found = arrayfun(@(x) find(arrayfun(@(y) sum(contains(ENST_RefSeqIDs{y},Multi_ENST_RefSeqIDz{k})),find(~cellfun(@isempty,ENST_RefSeqIDs))),1),1:length(Multi_ENST_RefSeqIDz{k}));
        TranscriptIsoform_chrom{k} = {EMBL_ChrNums{Exists_Match(match_Found)}};
    end
    All_ENST_RefSeqIDz = unique(vertcat(Multi_ENST_RefSeqIDz{:}));
    kc = zeros(size(Multi_ENST_RefSeqIDz,2),max(cellfun(@length,Multi_ENST_RefSeqIDz)));
    for k = 1:size(Multi_ENST_RefSeqIDz,2)
        for k2 = 1:size(Multi_ENST_RefSeqIDz{k},1)
            kc(k,k2) = find(strcmp(All_ENST_RefSeqIDz,Multi_ENST_RefSeqIDz{k}{k2}));
        end
    end
    FoldName = cell(1,size(Multi_ENST_RefSeqIDz,2));
    tempName = cell(1,size(inputs1{gene_num,1},2));
    for k = 1:size(inputs1{gene_num,1},2)
        for k2 = 1:size(Multi_ENST_RefSeqIDz{k},1)
            if (not(isfolder(['(' Multi_GeneName{k} ')_' Multi_ENST_RefSeqIDz{k}{k2}])))
                mkdir(['(' Multi_GeneName{k} ')_' Multi_ENST_RefSeqIDz{k}{k2}])
            end
            FoldName{kc(k,k2)} = ['(' Multi_GeneName{k} ')_' Multi_ENST_RefSeqIDz{k}{k2}];
        end
        tempName{k} = ['(' Multi_GeneName{k} ')_' strjoin(Multi_ENST_RefSeqIDz{k},'_')];
    end
    if (not(isfolder(strjoin(tempName,'_'))))
        mkdir(strjoin(tempName,'_'))
    end
    All_chrom = TranscriptIsoform_chrom;
    settings.FolderName = strjoin(tempName,'_');
end
settings.FolderRootName = strcat(saveRoot,settings.FolderName);
settings.rootName = strjoin(inputs1{gene_num,1},'_');


designerName0 = '_NLPDS';
GeneName = inputs1{gene_num,5};
%% Generate Probes
%Update for multiple targets probe has associated list of on-target genes
if (settings.SingleOrMulti==1&&settings.AllIsoforms == 0)%One Gene/One Isoform
    try
        load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_probes' designerName '.mat'],'probes')
    catch
        [init_probes,~,~,~,~] = ...
            BKJH_Probe_Generator([Lmin Lmax],inputs1{gene_num,1},inputs1{gene_num,2},inputs1{gene_num,3},inputs1{gene_num,8},inputs1{gene_num,9},settings.isOffline,settings.SEQdbRoot);
        probes = init_probes;
        %% Get Probes from other software
        Lpmin = min(cell2mat(cellfun(@length,{probes{:,2}},'UniformOutput',false)));
        Lpmax = max(cell2mat(cellfun(@length,{probes{:,2}},'UniformOutput',false)));
        switch id2
            case 1
                ONP = readtable(['Probe Designer ' GeneName(2:end-1) ' probes'],'Range','B2:B100','ReadVariableNames',false);
                OldNeuertProbes = ONP.Var1.';
                OldNeuertProbes = convertCharsToStrings(OldNeuertProbes);
                for m = 1:length(OldNeuertProbes)
                    try
                        OldProbes(m) = find(cell2mat(cellfun(@(x) strcmpi(x,seqrcomplement(OldNeuertProbes{m})),{probes{:,2}},'UniformOutput',false)));
                    catch
                        OldProbes(m) = NaN;
                    end
                end
                chosenProbes = OldProbes(~isnan(OldProbes));
            case 2
                AllStellarisProbes = [];
                StellarisProbesLi = cell(1,6);
                for v = 0:5
                    Li = strcat('L',num2str(v));
                    SP = readtable(strcat('OtherSoftwareProbeDesign/Stellaris/', GeneName, '_', inputs1{gene_num,1}{:}, '_Stellaris.xlsx'),'Range','D2:D100','Sheet',Li,'ReadVariableNames',false);
                    try
                        StellarisProbes = SP.Var1.';
                        StellarisProbesLi{v+1} = StellarisProbes(~isnan(StellarisProbes));
                        AllStellarisProbes = [AllStellarisProbes StellarisProbes(~isnan(StellarisProbes))];
                    catch
                        StellarisProbesLi{v+1}  = [];
                    end
                end
                load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_probes' designerName0 '.mat'],'probes')
                temp = probes;
                MaxTile = size(probes,1);
                for v = 0:5
                    designerNameLi = strcat('_Stellaris_L',num2str(v));
                    chosenProbes = StellarisProbesLi{v+1};
                    chosenProbes(chosenProbes>MaxTile) = [];
                    StellarisProbesLi{v+1} = chosenProbes;
                    probes = temp(chosenProbes,:);
                    save([saveRoot filesep settings.FolderName filesep settings.FolderName '_probes' designerNameLi '.mat'],'chosenProbes','probes','-v7.3')
                end
                probes = temp;
                chosenProbes = unique(AllStellarisProbes);
                chosenProbes(chosenProbes>MaxTile) = [];
            case 3
                transcript_ID = inputs1{gene_num,1};
                oligoStan_probeseqs = cell(1,size(transcript_ID,1));
                GeneSeqs = fastaread('OtherSoftwareProbeDesign/Oligostan/Targetseqs.fa');
                dict = containers.Map({GeneSeqs.Header},arrayfun(@(x) strcat('ENSG',num2str(x)),1:length(GeneSeqs),'Un',0));
                clear GeneSeqs
                if (gene_num<66)
                    OligoStan_oligos = readtable(['OtherSoftwareProbeDesign/Oligostan/Distance1_Design/Merged_design.txt'],'FileType','delimitedtext');
                else
                    OligoStan_oligos = readtable(['OtherSoftwareProbeDesign/Oligostan/Design_ACTs/Merged_design.txt'],'FileType','delimitedtext');
                end
                for z = 1:size(transcript_ID,1)
                    oligoStan_probeseqs{z} = OligoStan_oligos.Seq(strcmp(OligoStan_oligos.ENSG,dict(transcript_ID{z})));
                end
                oligoStan_probeseqs = cat(2,oligoStan_probeseqs{:});
                clear dict OligoStan_oligos
                for m = 1:length(oligoStan_probeseqs)
                    %remove trailing empty space
                    try
                        OligoStanProbes(m) = find(cell2mat(cellfun(@(x) strcmpi(x,seqrcomplement(oligoStan_probeseqs{m})),{probes{:,2}},'UniformOutput',false)));
                    catch
                        OligoStanProbes(m) = NaN;
                    end
                end
                chosenProbes = OligoStanProbes(~isnan(OligoStanProbes));
            case 4
                transcript_ID = inputs1{gene_num,1};
                paintOptions = detectImportOptions(['OtherSoftwareProbeDesign/PaintSHOP/' PaintSHOPfile],'FileType','delimitedtext');
                paintOptions.VariableNames = {'chrom','start','stop','sequence',...
                    'Tm','on_score','off_score','repeat','prob','max_kmer','strand','refseq','transcript_id','gene_id'};
                paint_oligos = tdfread(['OtherSoftwareProbeDesign/PaintSHOP/' PaintSHOPfile],'\t',paintOptions);
                paintOligos2 = cell2struct(struct2cell(paint_oligos),paintOptions.VariableNames);
                clear paintOptions paint_oligos
                paintOligos3 = struct2table(paintOligos2);clear paintOligos2
                paintOligos4 = table2struct(paintOligos3);clear paintOligos3
                for z = 1:size(transcript_ID,1)
                    B = extractBefore(transcript_ID{z},'.');
                    probe_rows{z} = find(cell2mat(arrayfun(@(x) strcmp(strtrim(paintOligos4(x).refseq),B),1:size(paintOligos4,1),'UniformOutput',false)));
                end
                probe_rows = cat(2,probe_rows{:});
                probe_rows2 = find(cell2mat(arrayfun(@(x) strcmpi(strtrim(paintOligos4(x).gene_id),GeneName(2:end-1)),1:size(paintOligos4,1),'UniformOutput',false)));
                probe_rows = intersect(probe_rows,probe_rows2);
                PaintSHOP_oligos = strtrim({paintOligos4(probe_rows).sequence});%why not just PaintSHOP_oligos
                clear paintOligos4
                for m = 1:length(PaintSHOP_oligos)
                    %remove trailing empty space
                    try
                        PaintShopProbes(m) = find(cell2mat(cellfun(@(x) strcmpi(x,seqrcomplement(PaintSHOP_oligos{m})),{probes{:,2}},'UniformOutput',false)));
                    catch
                        PaintShopProbes(m) = NaN;
                    end
                end
                chosenProbes = PaintShopProbes(~isnan(PaintShopProbes));
                clear PaintSHOP_oligos
            case 5
                MERFISH_oligos = fastaread('OtherSoftwareProbeDesign/MERFISH/L1R1_oligos.fasta');%trDesigner inputs might need to be customized
                probe_rows = find(cell2mat(arrayfun(@(x) contains(MERFISH_oligos(x).Header,inputs1{id,1}),1:size(MERFISH_oligos,1),'UniformOutput',false)));
                MERFISH_barcoded_oligos = arrayfun(@(x) strsplit(MERFISH_oligos(x).Sequence," "),probe_rows,'UniformOutput',false); % array of array 8. 4th entry
                MERFISH_seqs = cell(1,length(MERFISH_barcoded_oligos));
                for m = 1:length(MERFISH_barcoded_oligos)
                    lengths_sq = cellfun(@length,MERFISH_barcoded_oligos{m});
                    lens1 = find(lengths_sq==1);
                    MERFISH_seqs{m} = MERFISH_barcoded_oligos{m}{lens1(2)-1};
                end
                for m = 1:length(MERFISH_barcoded_oligos)
                    try
                        MERFISHProbes(m) = find(cell2mat(cellfun(@(x) strcmpi(x,seqrcomplement(MERFISH_seqs{m})),{probes{:,2}},'UniformOutput',false)));
                    catch
                        MERFISHProbes(m) = NaN;
                    end
                end
                chosenProbes = MERFISHProbes(~isnan(MERFISHProbes));
                probePos = [probes{chosenProbes,3}];
                %detected repeats need to be removed.
                if (length(unique(probePos))<length(chosen_probes))
                    [~,ia,~] = unique(probePos);
                    chosenProbes = sort(chosenProbes(ia));
                end
            case 6
                if (id<66)
                    MERFISH_oligos = fastaread('OtherSoftwareProbeDesign/MERFISH/L1R2_oligos.fasta');%trDesigner inputs might need to be customized
                else
                    MERFISH_oligos = fastaread('OtherSoftwareProbeDesign/MERFISH/L1R2_actins_oligos.fasta');%trDesigner inputs might need to be customized
                end
                probe_rows = find(cell2mat(arrayfun(@(x) contains(MERFISH_oligos(x).Header,inputs1{id,1}),1:size(MERFISH_oligos,1),'UniformOutput',false)));
                MERFISH_barcoded_oligos = arrayfun(@(x) strsplit(MERFISH_oligos(x).Sequence," "),probe_rows,'UniformOutput',false); % array of array 8. 4th entry
                MERFISH_seqs = cell(1,length(MERFISH_barcoded_oligos));
                for m = 1:length(MERFISH_barcoded_oligos)
                    lengths_sq = cellfun(@length,MERFISH_barcoded_oligos{m});
                    lens1 = find(lengths_sq==1);
                    MERFISH_seqs{m} = MERFISH_barcoded_oligos{m}{lens1(2)-1};
                end
                for m = 1:length(MERFISH_barcoded_oligos)
                    try
                        MERFISHProbes(m) = find(cell2mat(cellfun(@(x) strcmpi(x,seqrcomplement(MERFISH_seqs{m})),{probes{:,2}},'UniformOutput',false)));
                    catch
                        MERFISHProbes(m) = NaN;
                    end
                end
                chosenProbes = MERFISHProbes(~isnan(MERFISHProbes));
                probePos = [probes{chosenProbes,3}];
                %detected repeats need to be removed.
                if (length(unique(probePos))<length(chosenProbes))
                    [~,ia,~] = unique(probePos);
                    chosenProbes = sort(chosenProbes(ia));
                end
            case 7
                try
                    load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_probes' designerName '.mat'],'probes')
                catch
                    [init_probes,~,~,~,~] = ...
                        BKJH_Probe_Generator([Lmin Lmax],inputs1{gene_num,1},inputs1{gene_num,2},inputs1{gene_num,3},inputs1{gene_num,8},inputs1{gene_num,9},settings.isOffline,settings.SEQdbRoot);
                    probes = init_probes;
                    save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_probes' designerName '.mat'],'probes','-v7.3')
                end
            otherwise
                ProbeSets = randomProbeSets(probes,[Lpmin Lpmax],settings.ProbeSpacing,1,max_probes);
                chosenProbes = ProbeSets{1};
        end
        if (id2<7)
            probes = probes(chosenProbes,:);
            save([saveRoot filesep settings.FolderName filesep settings.FolderName '_probes' designerName '.mat'],'probes','-v7.3')
        end
    end
    if (id2==2)
        load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_probes' designerName0 '.mat'],'probes')
        MaxTile = size(probes,1);
        AllStellarisProbes = [];
        StellarisProbesLi = cell(1,6);
        for v = 0:5
            Li = strcat('L',num2str(v));
            SP = readtable(strcat('OtherSoftwareProbeDesign/Stellaris/', GeneName, '_', inputs1{gene_num,1}{:}, '_Stellaris.xlsx'),'Range','D2:D100','Sheet',Li,'ReadVariableNames',false);
            try
                StellarisProbes = SP.Var1.';
                StellarisProbesLi{v+1} = StellarisProbes(~isnan(StellarisProbes));
                chosenProbes = StellarisProbesLi{v+1};
                chosenProbes(chosenProbes>MaxTile) = [];
                StellarisProbesLi{v+1} = chosenProbes;
                AllStellarisProbes = [AllStellarisProbes StellarisProbes(~isnan(StellarisProbes))];
            catch
                StellarisProbesLi{v+1}  = [];
            end
        end
        temp = probes;
        for v = 0:5
            designerNameLi = strcat('_Stellaris_L',num2str(v));
            chosenProbes = StellarisProbesLi{v+1};
            probes = temp(chosenProbes,:);
            save([saveRoot filesep settings.FolderName filesep settings.FolderName '_probes' designerNameLi '.mat'],'chosenProbes','probes','-v7.3')
        end
        probes = temp;
        chosenProbes = unique(AllStellarisProbes);
        chosenProbes(chosenProbes>MaxTile) = [];
        probes = probes(chosenProbes,:);
        save([saveRoot filesep settings.FolderName filesep settings.FolderName '_probes' designerName '.mat'],'probes','-v7.3')
        probes = temp;
    end
elseif (settings.SingleOrMulti==1&&settings.AllIsoforms == 1)%One Gene/All Isoform
    Isoform_probes = cell(1,size(ENST_RefSeqIDz,1));
    for k = 1:size(ENST_RefSeqIDz,1)
        try
            load([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_probes' designerName '.mat'],'probes');
        catch
            [init_probes,~,~,~,~] = ...
                BKJH_Probe_Generator([Lmin Lmax],ENST_RefSeqIDz(k),inputs1{gene_num,2},inputs1{gene_num,3},inputs1{gene_num,8},inputs1{gene_num,9},settings.isOffline,settings.SEQdbRoot);
            probes = init_probes;
        end
        Isoform_probes{k} = probes;
    end
    try
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_probes' designerName '.mat'],'probes');
    catch
        [init_probes,~,~,~,~] = ... %multi-gene/isoform
            BKJH_MultiplexProbe_Generator([Lmin Lmax],ENST_RefSeqIDz',inputs1{gene_num,2},inputs1{gene_num,3},inputs1{gene_num,8},inputs1{gene_num,9},settings.isOffline,settings.SEQdbRoot);
        probes = init_probes;
    end
    AllIsoforms_probes = probes;
elseif (settings.SingleOrMulti==0&&settings.AllIsoforms == 0)%Multi Gene/One Isoform
    Transcript_probes = cell(1,size(inputs1{gene_num,1},2));
    for k = 1:size(inputs1{gene_num,1},2)
        try
            load([saveRoot filesep FoldName{k} filesep FoldName{k} '_probes' designerName '.mat'],'probes');
        catch
            [init_probes,~,~,~,~] = ...
                BKJH_Probe_Generator([Lmin Lmax],inputs1{gene_num,1}(k),inputs1{gene_num,2},inputs1{gene_num,3},inputs1{gene_num,8},inputs1{gene_num,9},settings.isOffline,settings.SEQdbRoot);
            probes = init_probes;
        end
        Transcript_probes{k} = probes;
    end
    try
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_probes' designerName '.mat'],'probes');
    catch
        [init_probes,~,~,~,~] = ... %multi-gene/isoform
            BKJH_MultiplexProbe_Generator([Lmin Lmax],inputs1{gene_num,1},inputs1{gene_num,2},inputs1{gene_num,3},inputs1{gene_num,8},inputs1{gene_num,9},settings.isOffline,settings.SEQdbRoot);
        probes = init_probes;
    end
    AllTranscripts_probes = probes;
elseif (settings.SingleOrMulti==0&&settings.AllIsoforms == 1)%Multi Gene/All Isoform
    TranscriptIsoform_probes = cell(1,size(All_ENST_RefSeqIDz,1));
    for k = 1:size(inputs1{gene_num,1},2)
        for k2 = 1:length(Multi_ENST_RefSeqIDz{k})
            try
                load([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_probes' designerName '.mat'],'probes');
            catch
                [init_probes,~,~,~,~] = ...
                    BKJH_Probe_Generator([Lmin Lmax],Multi_ENST_RefSeqIDz{k}(k2),inputs1{gene_num,2},inputs1{gene_num,3},inputs1{gene_num,8},inputs1{gene_num,9},settings.isOffline,settings.SEQdbRoot);
                probes = init_probes;
            end
            TranscriptIsoform_probes{kc(k,k2)} = probes;
        end
    end
    try
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_probes' designerName '.mat'],'probes');
    catch
        [init_probes,~,~,~,~] = ... %multi-gene/isoform
            BKJH_MultiplexProbe_Generator([Lmin Lmax],All_ENST_RefSeqIDz',inputs1{gene_num,2},inputs1{gene_num,3},inputs1{gene_num,8},inputs1{gene_num,9},settings.isOffline,settings.SEQdbRoot);
        probes = init_probes;
    end
    AllTranscriptIsoforms_probes = probes;
end


%% BLAST Probes
needToCombine = 1;
if (settings.SingleOrMulti==1&&settings.AllIsoforms == 0)%One Gene/One Isoform
    if (id2~=2)
        try
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_hits_table' designerName '.mat'],'gene_table')
        catch                                                                  %Organism           %(Gene Name)        %(Gene Name)        %ChrNum
            [~,gene_table] = Probe_checker_general_JH10(probes,inputs1{gene_num,4},inputs1{gene_num,5},inputs1{gene_num,6},inputs1{gene_num,7},settings);
            save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_hits_table' designerName '.mat'],'-mat','gene_table','-v7.3');
        end
    else
        try
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_hits_table' designerName0 '.mat'],'gene_table')
        catch                                                                  %Organism           %(Gene Name)        %(Gene Name)        %ChrNum
            [~,gene_table] = Probe_checker_general_JH10(probes,inputs1{gene_num,4},inputs1{gene_num,5},inputs1{gene_num,6},inputs1{gene_num,7},settings);
            save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_hits_table' designerName0 '.mat'],'-mat','gene_table','-v7.3');
        end
    end
elseif (settings.SingleOrMulti==1&&settings.AllIsoforms == 1)%One Gene/All Isoform
    try
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_hits_table' designerName '.mat'],'gene_table')
        needToCombine = 0;
    catch
        Isoform_gene_table = cell(1,size(ENST_RefSeqIDz,1));
        for k = 1:size(ENST_RefSeqIDz,1)
            try
                load([saveRoot '/(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '/(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_hits_table' designerName '.mat'],'gene_table')
            catch
                [~,gene_table] = Probe_checker_general_JH10(Isoform_probes{k},inputs1{gene_num,4},inputs1{gene_num,5},inputs1{gene_num,6},Isoform_chrom{k},settings);
                save([saveRoot '/(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '/(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_hits_table' designerName '.mat'],'-mat','gene_table','-v7.3');
            end
            Isoform_gene_table{k} = gene_table;
        end
    end
elseif (settings.SingleOrMulti==0&&settings.AllIsoforms == 0)%Multi Gene/One Isoform
    try
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_hits_table' designerName '.mat'],'gene_table')
        needToCombine = 0;
    catch
        Transcript_gene_table =  cell(1,size(inputs1{gene_num,1},2));
        for k = 1:size(inputs1{gene_num,1},2)
            try
                load([saveRoot filesep FoldName{k} filesep FoldName{k}  '_hits_table' designerName '.mat'],'gene_table')
            catch
                [~,gene_table] = Probe_checker_general_JH10(Transcript_probes{k},inputs1{gene_num,4},['(' Multi_GeneName{k} ')'],['(' Multi_GeneName{k} ')'],Transcript_chrom{k},settings);
                save([saveRoot filesep FoldName{k} filesep FoldName{k}  '_hits_table' designerName '.mat'],'-mat','gene_table','-v7.3');
            end
            Transcript_gene_table{k} = gene_table;
        end
    end
elseif (settings.SingleOrMulti==0&&settings.AllIsoforms == 1)%Multi Gene/All Isoform
    try
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_hits_table' designerName '.mat'],'gene_table')
        needToCombine = 0;
    catch
        TranscriptIsoform_gene_table =  cell(1,size(All_ENST_RefSeqIDz,1));
        for k = 1:size(inputs1{gene_num,1},2)
            for k2 = 1:length(Multi_ENST_RefSeqIDz{k})
                try
                    load([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_hits_table' designerName '.mat'],'gene_table')
                catch                                                                                                                                                                    %affected
                    [~,gene_table] = Probe_checker_general_JH10(TranscriptIsoform_probes{k},inputs1{gene_num,4},['(' Multi_GeneName{k} ')'],['(' Multi_GeneName{k} ')'],TranscriptIsoform_chrom{k}{k2},settings);
                    save([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_hits_table' designerName '.mat'],'-mat','gene_table','-v7.3');
                end
                TranscriptIsoform_gene_table{kc(k,k2)} = gene_table;
            end
        end
    end
end
%% Add saving combined hits and hits_table
if (needToCombine)
    if (settings.SingleOrMulti==1&&settings.AllIsoforms == 1)%One Gene/All Isoform
        UpdateProbeNumAssociation = cell(1,size(ENST_RefSeqIDz,1));
        UpdateTargetNumAssociation = cell(1,size(ENST_RefSeqIDz,1));
        Isoform_DPS_Names = cell(1,size(ENST_RefSeqIDz,1));
        Isoform_gene_table3 = cell(1,size(ENST_RefSeqIDz,1));
        modIsoform_gene_table = Isoform_gene_table;
        NE_Numbers = [AllIsoforms_probes{:,1}];
        NE_Positions = vertcat(AllIsoforms_probes{:,3});
        for k = 1:size(ENST_RefSeqIDz,1)
            OGP_Numbers = [Isoform_probes{k}{:,1}];
            OG_Position = [Isoform_probes{k}{:,3}];
            UpdateProbeNumAssociation{k} = [OGP_Numbers;...
                arrayfun(@(x) NE_Numbers(NE_Positions(:,k)==OG_Position(x)),1:length(OG_Position))]';
        end
        for k = 1:size(ENST_RefSeqIDz,1)
            old_Nums = Isoform_gene_table{k}.ProbeNum;
            dict = containers.Map(UpdateProbeNumAssociation{k}(:,1),UpdateProbeNumAssociation{k}(:,2));
            new_Nums = arrayfun(@(x) dict(old_Nums(x)),1:length(old_Nums))';
            modIsoform_gene_table{k}.ProbeNum = new_Nums;
        end
        AllIsoforms_gene_table = vertcat(modIsoform_gene_table{:});
        gene_table = AllIsoforms_gene_table;
        %Target numbers Association
        for k =  1:size(ENST_RefSeqIDz,1)
            temp_gene_table = modIsoform_gene_table{k};
            temp_gene_table = sortrows(temp_gene_table,[7 6],'ascend');
            temp_gene_table = temp_gene_table(temp_gene_table.Match>=settings.MinHomologySearchTargetSize,:);
            MinusStrandedHits = find(contains(temp_gene_table.Strand,'Minus'));
            if (strcmp(settings.referenceType,'RefSeq'))
                RNA_IDs_1 = find(contains(temp_gene_table.Name,'NM_'));
                RNA_IDs_2 = find(contains(temp_gene_table.Name,'NR_'));
                RNA_IDs_3 = find(contains(temp_gene_table.Name,'XM_'));
                RNA_IDs_4 = find(contains(temp_gene_table.Name,'XR_'));
                contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
            elseif (strcmp(settings.referenceType,'ENSEMBL'))
                RNA_IDs_1 = find(contains(temp_gene_table.Name,'EN'));
                RNA_IDs_2 = find(contains(temp_gene_table.Name,'EN'));
                RNA_IDs_3 = find(contains(temp_gene_table.Name,'EN'));
                RNA_IDs_4 = find(contains(temp_gene_table.Name,'EN'));
                contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
            else
                ss = 1;
                % RNA_IDs_1 = find(contains(gene_table.Name,'EN'));
                % RNA_IDs_2 = find(contains(gene_table.Name,'EN'));
                % RNA_IDs_3 = find(contains(gene_table.Name,'EN'));
                % RNA_IDs_4 = find(contains(gene_table.Name,'EN'));
                % contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
            end
            RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);
            temp_gene_table = temp_gene_table(setdiff(1:size(temp_gene_table,1),RNA_MissedFilteredHits),:);
            temp_gene_table.Ax = min(temp_gene_table.SubjectIndices,[],2);
            temp_gene_table.Bx = max(temp_gene_table.SubjectIndices,[],2);
            temp_gene_table = sortrows(temp_gene_table,[7 13],'ascend');
            Isoform_gene_table3{k} = temp_gene_table;
            Names = unique(temp_gene_table.Name);
            Names = convertCharsToStrings(Names);
            Isoform_DPS_Names{k} = Names;
        end
        temp_gene_table = AllIsoforms_gene_table;
        temp_gene_table = sortrows(temp_gene_table,[7 6],'ascend');
        temp_gene_table = temp_gene_table(temp_gene_table.Match>=settings.MinHomologySearchTargetSize,:);
        MinusStrandedHits = find(contains(temp_gene_table.Strand,'Minus'));
        if (strcmp(settings.referenceType,'RefSeq'))
            RNA_IDs_1 = find(contains(temp_gene_table.Name,'NM_'));
            RNA_IDs_2 = find(contains(temp_gene_table.Name,'NR_'));
            RNA_IDs_3 = find(contains(temp_gene_table.Name,'XM_'));
            RNA_IDs_4 = find(contains(temp_gene_table.Name,'XR_'));
            contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
        elseif (strcmp(settings.referenceType,'ENSEMBL'))
            RNA_IDs_1 = find(contains(temp_gene_table.Name,'EN'));
            RNA_IDs_2 = find(contains(temp_gene_table.Name,'EN'));
            RNA_IDs_3 = find(contains(temp_gene_table.Name,'EN'));
            RNA_IDs_4 = find(contains(temp_gene_table.Name,'EN'));
            contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
        else
            ss = 1;
            % RNA_IDs_1 = find(contains(gene_table.Name,'EN'));
            % RNA_IDs_2 = find(contains(gene_table.Name,'EN'));
            % RNA_IDs_3 = find(contains(gene_table.Name,'EN'));
            % RNA_IDs_4 = find(contains(gene_table.Name,'EN'));
            % contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
        end
        RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);
        temp_gene_table = temp_gene_table(setdiff(1:size(temp_gene_table,1),RNA_MissedFilteredHits),:);
        temp_gene_table.Ax = min(temp_gene_table.SubjectIndices,[],2);
        temp_gene_table.Bx = max(temp_gene_table.SubjectIndices,[],2);
        temp_gene_table = sortrows(temp_gene_table,[7 13],'ascend');
        Names = unique(temp_gene_table.Name);
        Names = convertCharsToStrings(Names);
        AllIsoforms_DPS_Names = Names;
        NE_TargetNum = 1:length(AllIsoforms_DPS_Names);
        for k = 1:size(ENST_RefSeqIDz,1)
            OGT_Numbers = 1:length(Isoform_DPS_Names{k});
            UpdateTargetNumAssociation{k} = [OGT_Numbers;...
                cellfun(@(x) NE_TargetNum(strcmp(x,AllIsoforms_DPS_Names)),Isoform_DPS_Names{k})']';
        end
    elseif (settings.SingleOrMulti==0&&settings.AllIsoforms == 0)%Multi Gene/One Isoform
        UpdateProbeNumAssociation = cell(1,size(inputs1{gene_num,1},2));
        UpdateTargetNumAssociation = cell(1,size(inputs1{gene_num,1},2));
        Transcript_DPS_Names = cell(1,size(inputs1{gene_num,1},2));
        Transcript_gene_table3 = cell(1,size(inputs1{gene_num,1},2));
        modTranscript_gene_table = Transcript_gene_table;
        NE_Numbers = [AllTranscripts_probes{:,1}];
        NE_Positions = vertcat(AllTranscripts_probes{:,3});
        for k = 1:size(inputs1{gene_num,1},2)
            OGP_Numbers = [Transcript_probes{k}{:,1}];
            OG_Position = [Transcript_probes{k}{:,3}];
            UpdateProbeNumAssociation{k} = [OGP_Numbers;...
                arrayfun(@(x) NE_Numbers(NE_Positions(:,k)==OG_Position(x)),1:length(OG_Position))]';
        end
        for k = 1:size(inputs1{gene_num,1},2)
            old_Nums = Transcript_gene_table{k}.ProbeNum;
            dict = containers.Map(UpdateProbeNumAssociation{k}(:,1),UpdateProbeNumAssociation{k}(:,2));
            new_Nums = arrayfun(@(x) dict(old_Nums(x)),1:length(old_Nums))';
            modTranscript_gene_table{k}.ProbeNum = new_Nums;
        end
        AllTranscript_gene_table = vertcat(modTranscript_gene_table{:});
        gene_table = AllTranscript_gene_table;
        %Target numbers Association
        for k =  1:size(inputs1{gene_num,1},2)
            temp_gene_table = modTranscript_gene_table{k};
            temp_gene_table = sortrows(temp_gene_table,[7 6],'ascend');
            temp_gene_table = temp_gene_table(temp_gene_table.Match>=settings.MinHomologySearchTargetSize,:);
            MinusStrandedHits = find(contains(temp_gene_table.Strand,'Minus'));
            if (strcmp(settings.referenceType,'RefSeq'))
                RNA_IDs_1 = find(contains(temp_gene_table.Name,'NM_'));
                RNA_IDs_2 = find(contains(temp_gene_table.Name,'NR_'));
                RNA_IDs_3 = find(contains(temp_gene_table.Name,'XM_'));
                RNA_IDs_4 = find(contains(temp_gene_table.Name,'XR_'));
                contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
            elseif (strcmp(settings.referenceType,'ENSEMBL'))
                RNA_IDs_1 = find(contains(temp_gene_table.Name,'EN'));
                RNA_IDs_2 = find(contains(temp_gene_table.Name,'EN'));
                RNA_IDs_3 = find(contains(temp_gene_table.Name,'EN'));
                RNA_IDs_4 = find(contains(temp_gene_table.Name,'EN'));
                contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
            else
                ss = 1;
                % RNA_IDs_1 = find(contains(gene_table.Name,'EN'));
                % RNA_IDs_2 = find(contains(gene_table.Name,'EN'));
                % RNA_IDs_3 = find(contains(gene_table.Name,'EN'));
                % RNA_IDs_4 = find(contains(gene_table.Name,'EN'));
                % contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
            end
            RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);
            temp_gene_table = temp_gene_table(setdiff(1:size(temp_gene_table,1),RNA_MissedFilteredHits),:);
            temp_gene_table.Ax = min(temp_gene_table.SubjectIndices,[],2);
            temp_gene_table.Bx = max(temp_gene_table.SubjectIndices,[],2);
            temp_gene_table = sortrows(temp_gene_table,[7 13],'ascend');
            Transcript_gene_table3{k} = temp_gene_table;
            Names = unique(temp_gene_table.Name);
            Names = convertCharsToStrings(Names);
            Transcript_DPS_Names{k} = Names;
        end
        temp_gene_table = AllTranscript_gene_table;
        temp_gene_table = sortrows(temp_gene_table,[7 6],'ascend');
        temp_gene_table = temp_gene_table(temp_gene_table.Match>=settings.MinHomologySearchTargetSize,:);
        MinusStrandedHits = find(contains(temp_gene_table.Strand,'Minus'));
        if (strcmp(settings.referenceType,'RefSeq'))
            RNA_IDs_1 = find(contains(temp_gene_table.Name,'NM_'));
            RNA_IDs_2 = find(contains(temp_gene_table.Name,'NR_'));
            RNA_IDs_3 = find(contains(temp_gene_table.Name,'XM_'));
            RNA_IDs_4 = find(contains(temp_gene_table.Name,'XR_'));
            contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
        elseif (strcmp(settings.referenceType,'ENSEMBL'))
            RNA_IDs_1 = find(contains(temp_gene_table.Name,'EN'));
            RNA_IDs_2 = find(contains(temp_gene_table.Name,'EN'));
            RNA_IDs_3 = find(contains(temp_gene_table.Name,'EN'));
            RNA_IDs_4 = find(contains(temp_gene_table.Name,'EN'));
            contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
        else
            ss = 1;
            % RNA_IDs_1 = find(contains(gene_table.Name,'EN'));
            % RNA_IDs_2 = find(contains(gene_table.Name,'EN'));
            % RNA_IDs_3 = find(contains(gene_table.Name,'EN'));
            % RNA_IDs_4 = find(contains(gene_table.Name,'EN'));
            % contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
        end
        RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);
        temp_gene_table = temp_gene_table(setdiff(1:size(temp_gene_table,1),RNA_MissedFilteredHits),:);
        temp_gene_table.Ax = min(temp_gene_table.SubjectIndices,[],2);
        temp_gene_table.Bx = max(temp_gene_table.SubjectIndices,[],2);
        temp_gene_table = sortrows(temp_gene_table,[7 13],'ascend');
        Names = unique(temp_gene_table.Name);
        Names = convertCharsToStrings(Names);
        AllTranscript_DPS_Names = Names;
        NE_TargetNum = 1:length(AllTranscript_DPS_Names);
        for k = 1:size(inputs1{gene_num,1},2)
            OGT_Numbers = 1:length(Transcript_DPS_Names{k});
            UpdateTargetNumAssociation{k} = [OGT_Numbers;...
                cellfun(@(x) NE_TargetNum(strcmp(x,AllTranscript_DPS_Names)),Transcript_DPS_Names{k})']';
        end
    elseif (settings.SingleOrMulti==0&&settings.AllIsoforms == 1)%Multi Gene/All Isoform
        UpdateProbeNumAssociation = cell(1,size(All_ENST_RefSeqIDz,1));
        UpdateTargetNumAssociation = cell(1,size(All_ENST_RefSeqIDz,1));
        TranscriptIsoform_DPS_Names = cell(1,size(All_ENST_RefSeqIDz,1));
        TranscriptIsoform_gene_table3 = cell(1,size(All_ENST_RefSeqIDz,1));
        modTranscriptIsoform_gene_table = TranscriptIsoform_gene_table;
        NE_Numbers = [AllTranscriptIsoforms_probes{:,1}];
        NE_Positions = vertcat(AllTranscriptIsoforms_probes{:,3});
        for k = 1:size(All_ENST_RefSeqIDz,1)
            OGP_Numbers = [TranscriptIsoform_probes{k}{:,1}];
            OG_Position = [TranscriptIsoform_probes{k}{:,3}];
            UpdateProbeNumAssociation{k} = [OGP_Numbers;...
                arrayfun(@(x) NE_Numbers(NE_Positions(:,k)==OG_Position(x)),1:length(OG_Position))]';
        end
        for k = 1:size(All_ENST_RefSeqIDz,1)
            old_Nums = TranscriptIsoform_gene_table{k}.ProbeNum;
            dict = containers.Map(UpdateProbeNumAssociation{k}(:,1),UpdateProbeNumAssociation{k}(:,2));
            new_Nums = arrayfun(@(x) dict(old_Nums(x)),1:length(old_Nums))';
            modTranscriptIsoform_gene_table{k}.ProbeNum = new_Nums;
        end
        AllTranscriptIsoform_gene_table = vertcat(modTranscriptIsoform_gene_table{:});
        gene_table = AllTranscriptIsoform_gene_table;
        %Target numbers Association
        for k =  1:size(inputs1{gene_num,1},2)
            temp_gene_table = modTranscriptIsoform_gene_table{k};
            temp_gene_table = sortrows(temp_gene_table,[7 6],'ascend');
            temp_gene_table = temp_gene_table(temp_gene_table.Match>=settings.MinHomologySearchTargetSize,:);
            MinusStrandedHits = find(contains(temp_gene_table.Strand,'Minus'));
            if (strcmp(settings.referenceType,'RefSeq'))
                RNA_IDs_1 = find(contains(temp_gene_table.Name,'NM_'));
                RNA_IDs_2 = find(contains(temp_gene_table.Name,'NR_'));
                RNA_IDs_3 = find(contains(temp_gene_table.Name,'XM_'));
                RNA_IDs_4 = find(contains(temp_gene_table.Name,'XR_'));
                contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
            elseif (strcmp(settings.referenceType,'ENSEMBL'))
                RNA_IDs_1 = find(contains(temp_gene_table.Name,'EN'));
                RNA_IDs_2 = find(contains(temp_gene_table.Name,'EN'));
                RNA_IDs_3 = find(contains(temp_gene_table.Name,'EN'));
                RNA_IDs_4 = find(contains(temp_gene_table.Name,'EN'));
                contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
            else
                ss = 1;
                % RNA_IDs_1 = find(contains(gene_table.Name,'EN'));
                % RNA_IDs_2 = find(contains(gene_table.Name,'EN'));
                % RNA_IDs_3 = find(contains(gene_table.Name,'EN'));
                % RNA_IDs_4 = find(contains(gene_table.Name,'EN'));
                % contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
            end
            RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);
            temp_gene_table = temp_gene_table(setdiff(1:size(temp_gene_table,1),RNA_MissedFilteredHits),:);
            temp_gene_table.Ax = min(temp_gene_table.SubjectIndices,[],2);
            temp_gene_table.Bx = max(temp_gene_table.SubjectIndices,[],2);
            temp_gene_table = sortrows(temp_gene_table,[7 13],'ascend');
            TranscriptIsoform_gene_table3{k} = temp_gene_table;
            Names = unique(temp_gene_table.Name);
            Names = convertCharsToStrings(Names);
            TranscriptIsoform_DPS_Names{k} = Names;
        end
        temp_gene_table = AllTranscriptIsoform_gene_table;
        temp_gene_table = sortrows(temp_gene_table,[7 6],'ascend');
        temp_gene_table = temp_gene_table(temp_gene_table.Match>=settings.MinHomologySearchTargetSize,:);
        MinusStrandedHits = find(contains(temp_gene_table.Strand,'Minus'));
            if (strcmp(settings.referenceType,'RefSeq'))
                RNA_IDs_1 = find(contains(temp_gene_table.Name,'NM_'));
                RNA_IDs_2 = find(contains(temp_gene_table.Name,'NR_'));
                RNA_IDs_3 = find(contains(temp_gene_table.Name,'XM_'));
                RNA_IDs_4 = find(contains(temp_gene_table.Name,'XR_'));
                contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
            elseif (strcmp(settings.referenceType,'ENSEMBL'))
                RNA_IDs_1 = find(contains(temp_gene_table.Name,'EN'));
                RNA_IDs_2 = find(contains(temp_gene_table.Name,'EN'));
                RNA_IDs_3 = find(contains(temp_gene_table.Name,'EN'));
                RNA_IDs_4 = find(contains(temp_gene_table.Name,'EN'));
                contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
            else
                ss = 1;
                % RNA_IDs_1 = find(contains(gene_table.Name,'EN'));
                % RNA_IDs_2 = find(contains(gene_table.Name,'EN'));
                % RNA_IDs_3 = find(contains(gene_table.Name,'EN'));
                % RNA_IDs_4 = find(contains(gene_table.Name,'EN'));
                % contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
            end
        RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);
        temp_gene_table = temp_gene_table(setdiff(1:size(temp_gene_table,1),RNA_MissedFilteredHits),:);
        temp_gene_table.Ax = min(temp_gene_table.SubjectIndices,[],2);
        temp_gene_table.Bx = max(temp_gene_table.SubjectIndices,[],2);
        temp_gene_table = sortrows(temp_gene_table,[7 13],'ascend');
        Names = unique(temp_gene_table.Name);
        Names = convertCharsToStrings(Names);
        AllTranscriptIsoform_DPS_Names = Names;
        NE_TargetNum = 1:length(AllTranscriptIsoform_DPS_Names);
        for k = 1:size(All_ENST_RefSeqIDz,1)
            OGT_Numbers = 1:length(TranscriptIsoform_DPS_Names{k});
            UpdateTargetNumAssociation{k} = [OGT_Numbers;...
                cellfun(@(x) NE_TargetNum(strcmp(x,AllTranscriptIsoform_DPS_Names)),TranscriptIsoform_DPS_Names{k})']';
        end
    end
end
%% Get Gene Expression Information

if (settings.SingleOrMulti==1&&settings.AllIsoforms == 0)%One Gene/One Isoform
    if (id2~=2)
        try
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_ExpressionInfo' designerName '.mat'],'ExpressionMatrix');
        catch
            [ExpressionMatrix,get_expression_time] = A_JH_GetExpressionInfo_V2(gene_table,settings);
            save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_ExpressionInfo' designerName '.mat'],'ExpressionMatrix','get_expression_time','settings','-v7.3');
        end
    else
        try
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_ExpressionInfo' designerName0 '.mat'],'ExpressionMatrix');
        catch
            [ExpressionMatrix,get_expression_time] = A_JH_GetExpressionInfo_V2(gene_table,settings);
            save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_ExpressionInfo' designerName0 '.mat'],'ExpressionMatrix','get_expression_time','settings','-v7.3');
        end

    end
else %Other with combined gene_table
    if (id2~=2)
        try
            load([saveRoot filesep settings.FolderName filesep settings.FolderName '_ExpressionInfo' designerName '.mat'],'ExpressionMatrix');
        catch
            [ExpressionMatrix,get_expression_time] = A_JH_GetExpressionInfo_V2(gene_table,settings);
            save([saveRoot filesep settings.FolderName filesep settings.FolderName '_ExpressionInfo' designerName '.mat'],'ExpressionMatrix','get_expression_time','settings','-v7.3');
        end
    else
        try
            load([saveRoot filesep settings.FolderName filesep settings.FolderName '_ExpressionInfo' designerName0 '.mat'],'ExpressionMatrix');
        catch
            [ExpressionMatrix,get_expression_time] = A_JH_GetExpressionInfo_V2(gene_table,settings);
            save([saveRoot filesep settings.FolderName filesep settings.FolderName '_ExpressionInfo' designerName0 '.mat'],'ExpressionMatrix','get_expression_time','settings','-v7.3');
        end
    end
end


settings.saveRoot = saveRoot;
%% Get Thermodynamic Information (On-Target,Off-Target)
if (settings.SingleOrMulti==1&&settings.AllIsoforms == 0)%One Gene/One Isoform
    if (id2~=2)
        try
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_Tm' num2str(settings.HybridizationTemperature) '_OnOffThermoInfo' designerName '.mat'],'Kon','Koff','Kb_Match');
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dHInfo' designerName '.mat'],'-mat','dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match');
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dSInfo' designerName '.mat'],'-mat','dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match');
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_TmInfo' designerName '.mat'],'-mat','Tm_on','Tm_Match');
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dCpInfo' designerName '.mat'],'-mat','dCpon_eq','dCpeq_Match');
        catch
            [Kb_Match,Kon,Koff,dHeq_Match,dSeq_Match,dHf_Match,dSf_Match,dHr_Match,dSr_Match,dCpeq_Match,dHon_eq,dSon_eq,dHon_f,dSon_f,dHon_r,dSon_r,dCpon_eq,Tm_on,Tm_Match] = ...
                A_JH_GenerateThermoInfo_V6(probes,gene_table,inputs1{gene_num,5},settings);%add Kon Koff
            save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_Tm' num2str(settings.HybridizationTemperature) '_OnOffThermoInfo' designerName '.mat'],'-mat','Kon','Koff','Kb_Match','-v7.3');
            save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dHInfo' designerName '.mat'],'-mat','dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match','-v7.3');
            save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dSInfo' designerName '.mat'],'-mat','dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match','-v7.3');
            save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_TmInfo' designerName '.mat'],'-mat','Tm_on','Tm_Match','-v7.3');
            save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dCpInfo' designerName '.mat'],'-mat','dCpon_eq','dCpeq_Match','-v7.3');
        end
    else
        designerName0 = '_NLPDS';
        try
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_Tm' num2str(settings.HybridizationTemperature) '_OnOffThermoInfo' designerName0 '.mat'],'Kon','Koff','Kb_Match');
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dHInfo' designerName0 '.mat'],'-mat','dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match');
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dSInfo' designerName0 '.mat'],'-mat','dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match');
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_TmInfo' designerName0 '.mat'],'-mat','Tm_on','Tm_Match');
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dCpInfo' designerName0 '.mat'],'-mat','dCpon_eq','dCpeq_Match');
        catch
            [Kb_Match,Kon,Koff,dHeq_Match,dSeq_Match,dHf_Match,dSf_Match,dHr_Match,dSr_Match,dCpeq_Match,dHon_eq,dSon_eq,dHon_f,dSon_f,dHon_r,dSon_r,dCpon_eq,Tm_on,Tm_Match] = ...
                A_JH_GenerateThermoInfo_V5(probes,gene_table,inputs1{gene_num,5},settings);%add Kon Koff
            save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_Tm' num2str(settings.HybridizationTemperature) '_OnOffThermoInfo' designerName0 '.mat'],'-mat','Kon','Koff','Kb_Match','-v7.3');
            save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dHInfo' designerName0 '.mat'],'-mat','dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match','-v7.3');
            save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dSInfo' designerName0 '.mat'],'-mat','dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match','-v7.3');
            save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_TmInfo' designerName0 '.mat'],'-mat','Tm_on','Tm_Match','-v7.3');
            save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dCpInfo' designerName0 '.mat'],'-mat','dCpon_eq','dCpeq_Match','-v7.3');
        end
    end

elseif(settings.SingleOrMulti==1&&settings.AllIsoforms == 1)%One Gene/All Isoform
    try
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_OnOffThermoInfo' designerName '.mat'],'Kon','Koff');
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_BindingMatrices' designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod','dCp_mod')
        if (settings.BLASTdna)
            load([saveRoot filesep settings.FolderName filesep settings.FolderName '_BindingMatrices2' designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','dCp_Complement')
        end
        needToCombine = 0;
    catch
        Multi_Kb_Match = cell(1,size(ENST_RefSeqIDz,1));
        Multi_Kon = cell(1,size(ENST_RefSeqIDz,1));
        Multi_Koff = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dHeq_Match = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dSeq_Match = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dHf_Match = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dSf_Match = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dHr_Match = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dSr_Match = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dCpeq_Match = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dHon_eq = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dSon_eq = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dHon_f = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dSon_f = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dHon_r = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dSon_r = cell(1,size(ENST_RefSeqIDz,1));
        Multi_Tm_on = cell(1,size(ENST_RefSeqIDz,1));
        Multi_Tm_Match = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dCpon_eq = cell(1,size(ENST_RefSeqIDz,1));
        for k = 1:size(ENST_RefSeqIDz,1)
            settings.GeneName = extractBefore(FoldName{k},'_');
            settings.transcript_IDs = extractAfter(FoldName{k},'_');
            settings.ChrNum = All_chrom{k};
            settings.GeneChr = strcat('chr',All_chrom{k});
            try
                load([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_OnOffThermoInfo' designerName '.mat'],'Kon','Koff','Kb_Match');
                load([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_dHInfo' designerName '.mat'],'-mat','dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match');
                load([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_dSInfo' designerName '.mat'],'-mat','dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match');
                load([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_TmInfo' designerName '.mat'],'-mat','Tm_on','Tm_Match');
                load([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_dCInfo' designerName '.mat'],'-mat','dCpon_eq','dCpeq_Match');
            catch
                [Kb_Match,Kon,Koff,dHeq_Match,dSeq_Match,dHf_Match,dSf_Match,dHr_Match,dSr_Match,dCpeq_Match,dHon_eq,dSon_eq,dHon_f,dSon_f,dHon_r,dSon_r,dCpon_eq,Tm_on,Tm_Match] = ...
                    A_JH_GenerateThermoInfo_V5(Isoform_probes{k},Isoform_gene_table{k},['(' ENST_GeneName ')_' ENST_RefSeqIDz{k}],settings);
                save([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_OnOffThermoInfo' designerName '.mat'],'-mat','Kon','Koff','Kb_Match','-v7.3');
                save([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_dHInfo' designerName '.mat'],'-mat','dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match','-v7.3');
                save([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_dSInfo' designerName '.mat'],'-mat','dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match','-v7.3');
                save([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_TmInfo' designerName '.mat'],'-mat','Tm_on','Tm_Match','-v7.3');
                save([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_dCpInfo' designerName '.mat'],'-mat','dCpon_eq','dCpeq_Match','-v7.3');
            end
            Multi_Kb_Match{k} = Kb_Match;
            Multi_Kon{k} = Kon;
            Multi_Koff{k} = Koff;
            Multi_dHeq_Match{k} = dHeq_Match;
            Multi_dSeq_Match{k} = dSeq_Match;
            Multi_dHf_Match{k} = dHf_Match;
            Multi_dSf_Match{k} = dSf_Match;
            Multi_dHr_Match{k} = dHr_Match;
            Multi_dSr_Match{k} = dSr_Match;
            Multi_dHon_eq{k} = dHon_eq;
            Multi_dSon_eq{k} = dSon_eq;
            Multi_dHon_f{k} = dHon_f;
            Multi_dSon_f{k} = dSon_f;
            Multi_dHon_r{k} = dHon_r;
            Multi_dSon_r{k} = dSon_r;
            Multi_Tm_on{k} = Tm_on;
            Multi_Tm_Match{k} = Tm_Match;
            Multi_dCpon_eq{k} = dCpon_eq;
            Multi_dCpeq_Match{k} = dCpeq_Match;
        end
    end
elseif (settings.SingleOrMulti==0&&settings.AllIsoforms == 0)%Multi Gene/One Isoform
    try
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_OnOffThermoInfo' designerName '.mat'],'Kon','Koff');
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_BindingMatrices' designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod','dCp_mod')
        if (settings.BLASTdna)
            load([saveRoot filesep settings.FolderName filesep settings.FolderName '_BindingMatrices2' designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','dCp_Complement')
        end
        needToCombine = 0;
    catch
        Multi_Kb_Match = cell(1,size(inputs1{gene_num,1},2));
        Multi_Kon = cell(1,size(inputs1{gene_num,1},2));
        Multi_Koff = cell(1,size(inputs1{gene_num,1},2));
        Multi_dHeq_Match = cell(1,size(inputs1{gene_num,1},2));
        Multi_dSeq_Match = cell(1,size(inputs1{gene_num,1},2));
        Multi_dHf_Match = cell(1,size(inputs1{gene_num,1},2));
        Multi_dSf_Match = cell(1,size(inputs1{gene_num,1},2));
        Multi_dHr_Match = cell(1,size(inputs1{gene_num,1},2));
        Multi_dSr_Match = cell(1,size(inputs1{gene_num,1},2));
        Multi_dCpeq_Match = cell(1,size(inputs1{gene_num,1},2));
        Multi_dHon_eq = cell(1,size(inputs1{gene_num,1},2));
        Multi_dSon_eq = cell(1,size(inputs1{gene_num,1},2));
        Multi_dHon_f = cell(1,size(inputs1{gene_num,1},2));
        Multi_dSon_f = cell(1,size(inputs1{gene_num,1},2));
        Multi_dHon_r = cell(1,size(inputs1{gene_num,1},2));
        Multi_dSon_r = cell(1,size(inputs1{gene_num,1},2));
        Multi_dCpon_eq = cell(1,size(inputs1{gene_num,1},2));
        Multi_Tm_on = cell(1,size(inputs1{gene_num,1},2));
        Multi_Tm_Match = cell(1,size(inputs1{gene_num,1},2));
        for k = 1:size(inputs1{gene_num,1},2)
            settings.GeneName = extractBefore(FoldName{k},'_');
            settings.transcript_IDs = extractAfter(FoldName{k},'_');
            settings.ChrNum = All_chrom{k};
            settings.GeneChr = strcat('chr',All_chrom{k});
            try
                load([saveRoot filesep FoldName{k} filesep FoldName{k} '_Tm' num2str(settings.HybridizationTemperature) '_OnOffThermoInfo' designerName '.mat'],'Kon','Koff','Kb_Match');
                load([saveRoot filesep FoldName{k} filesep FoldName{k} '_dHInfo' designerName '.mat'],'-mat','dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match');
                load([saveRoot filesep FoldName{k} filesep FoldName{k} '_dSInfo' designerName '.mat'],'-mat','dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match');
                load([saveRoot filesep FoldName{k} filesep FoldName{k} '_TmInfo' designerName '.mat'],'-mat','Tm_on','Tm_Match');
                load([saveRoot filesep FoldName{k} filesep FoldName{k} '_dCpInfo' designerName '.mat'],'-mat','dCpon_eq','dCpeq_Match');
            catch
                [Kb_Match,Kon,Koff,dHeq_Match,dSeq_Match,dHf_Match,dSf_Match,dHr_Match,dSr_Match,dCpeq_Match,dHon_eq,dSon_eq,dHon_f,dSon_f,dHon_r,dSon_r,dCpon_eq,Tm_on,Tm_Match] = ...
                    A_JH_GenerateThermoInfo_V5(Transcript_probes{k},Transcript_gene_table{k},FoldName{k},settings);
                save([saveRoot filesep FoldName{k} filesep FoldName{k} '_Tm' num2str(settings.HybridizationTemperature) '_OnOffThermoInfo' designerName '.mat'],'-mat','Kon','Koff','Kb_Match','-v7.3');
                save([saveRoot filesep FoldName{k} filesep FoldName{k} '_dHInfo' designerName '.mat'],'-mat','dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match','-v7.3');
                save([saveRoot filesep FoldName{k} filesep FoldName{k} '_dSInfo' designerName '.mat'],'-mat','dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match','-v7.3');
                save([saveRoot filesep FoldName{k} filesep FoldName{k} '_TmInfo' designerName '.mat'],'-mat','Tm_on','Tm_Match','-v7.3');
                save([saveRoot filesep FoldName{k} filesep FoldName{k} '_dCpInfo' designerName '.mat'],'-mat','dCpon_eq','dCpeq_Match','-v7.3');
            end
            Multi_Kb_Match{k} = Kb_Match;
            Multi_Kon{k} = Kon;
            Multi_Koff{k} = Koff;
            Multi_dHeq_Match{k} = dHeq_Match;
            Multi_dSeq_Match{k} = dSeq_Match;
            Multi_dHf_Match{k} = dHf_Match;
            Multi_dSf_Match{k} = dSf_Match;
            Multi_dHr_Match{k} = dHr_Match;
            Multi_dSr_Match{k} = dSr_Match;
            Multi_dHon_eq{k} = dHon_eq;
            Multi_dSon_eq{k} = dSon_eq;
            Multi_dHon_f{k} = dHon_f;
            Multi_dSon_f{k} = dSon_f;
            Multi_dHon_r{k} = dHon_r;
            Multi_dSon_r{k} = dSon_r;
            Multi_Tm_on{k} = Tm_on;
            Multi_Tm_Match{k} = Tm_Match;
            Multi_dCpon_eq{k} = dCpon_eq;
            Multi_dCpeq_Match{k} = dCpeq_Match;
        end
    end
elseif (settings.SingleOrMulti==0&&settings.AllIsoforms == 1)%Multi Gene/All Isoform
    try
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_OnOffThermoInfo' designerName '.mat'],'Kon','Koff');
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_BindingMatrices' designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod','dCp_mod')
        if (settings.BLASTdna)
            load([saveRoot filesep settings.FolderName filesep settings.FolderName '_BindingMatrices2' designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','dCp_Complement')
        end
        needToCombine = 0;
    catch
        Multi_Kb_Match = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_Kon = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_Koff = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dHeq_Match = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dSeq_Match = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dHf_Match = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dSf_Match = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dHr_Match = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dSr_Match = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dCpeq_Match = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dHon_eq = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dSon_eq = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dHon_f = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dSon_f = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dHon_r = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dSon_r = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dCpon_eq = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_Tm_on = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_Tm_Match = cell(1,size(All_ENST_RefSeqIDz,1));
        for k = 1:size(inputs1{gene_num,1},2)
            for k2 = 1:length(Multi_ENST_RefSeqIDz{k})
                settings.GeneName = extractBefore(FoldName{kc(k,k2)},'_');
                settings.transcript_IDs = extractAfter(FoldName{kc(k,k2)},'_');
                settings.ChrNum = All_chrom{k}{k2};
                settings.GeneChr = strcat('chr',All_chrom{k}{k2});
                try
                    load([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_OnOffThermoInfo' designerName '.mat'],'Kon','Koff','Kb_Match');
                    load([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_dHInfo' designerName '.mat'],'-mat','dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match');
                    load([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_dSInfo' designerName '.mat'],'-mat','dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match');
                    load([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_TmInfo' designerName '.mat'],'-mat','Tm_on','Tm_Match');
                    load([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_dCpInfo' designerName '.mat'],'-mat','dCpon_eq','dCpeq_Match');
                catch                                                                                                                                                                    %affected
                    [Kb_Match,Kon,Koff,dHeq_Match,dSeq_Match,dHf_Match,dSf_Match,dHr_Match,dSr_Match,dCpeq_Match,dHon_eq,dSon_eq,dHon_f,dSon_f,dHon_r,dSon_r,dCpon_eq,Tm_on,Tm_Match] = ...
                        A_JH_GenerateThermoInfo_V5(TranscriptIsoform_probes{kc(k,k2)},TranscriptIsoform_gene_table{kc(k,k2)},FoldName{kc(k,k2)},settings);
                    save([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_Tm' num2str(settings.HybridizationTemperature) '_OnOffThermoInfo' designerName '.mat'],'-mat','Kon','Koff','Kb_Match','-v7.3');
                    save([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_dHInfo' designerName '.mat'],'-mat','dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match','-v7.3');
                    save([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_dSInfo' designerName '.mat'],'-mat','dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match','-v7.3');
                    save([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_TmInfo' designerName '.mat'],'-mat','Tm_on','Tm_Match','-v7.3');
                    save([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_CpInfo' designerName '.mat'],'-mat','dCpon_eq','dCpeq_Match','-v7.3');
                end
                Multi_Kb_Match{kc(k,k2)} = Kb_Match;
                Multi_Kon{kc(k,k2)} = Kon;
                Multi_Koff{kc(k,k2)} = Koff;
                Multi_dHeq_Match{kc(k,k2)} = dHeq_Match;
                Multi_dSeq_Match{kc(k,k2)} = dSeq_Match;
                Multi_dHf_Match{kc(k,k2)} = dHf_Match;
                Multi_dSf_Match{kc(k,k2)} = dSf_Match;
                Multi_dHr_Match{kc(k,k2)} = dHr_Match;
                Multi_dSr_Match{kc(k,k2)} = dSr_Match;
                Multi_dHon_eq{kc(k,k2)} = dHon_eq;
                Multi_dSon_eq{kc(k,k2)} = dSon_eq;
                Multi_dHon_f{kc(k,k2)} = dHon_f;
                Multi_dSon_f{kc(k,k2)} = dSon_f;
                Multi_dHon_r{kc(k,k2)} = dHon_r;
                Multi_dSon_r{kc(k,k2)} = dSon_r;
                Multi_Tm_on{kc(k,k2)} = Tm_on;
                Multi_Tm_Match{kc(k,k2)} = Tm_Match;
                Multi_dCpon_eq{kc(k,k2)} = dCpon_eq;
                Multi_dCpeq_Match{kc(k,k2)} = dCpeq_Match;
            end
        end
    end
end
%% Get Binding Site Mapping and Energy [Edit inside of binding site mapping code for save directory info
%add nascent expression by updating capability of getNascentBindingInfo_JH2
if (settings.withNascent)
    if (strcmp(settings.Organism,'Human'))
        Expression_Null = settings.DoAllGenesHaveSameExpression; %do genes have same expression
        CellTypeID = settings.CellType_ExprID;
        HumanOntology = settings.HumanSpecific.Ontology;%Normal/Cancer
        HumanTissueOnly = settings.HumanSpecific.TissueOrTissueAndCellType; % 0 (TissueOnssue/Cell Type)
        HumanExpGeneOrTransc = settings.HumanSpecific.HumanExpGeneOrTransc; % 1 (Gene/EMBL GENEID) , 0 (Transcript/EMBL Transcript ID)
        Human_SelTrack = settings.HumanSpecific.SCOutputTrack;
        expValType = settings.expressionValType;% 1-4 (expCounts,expValues,mean(CellTypeExpValues),one cell types CellTypeExpValues)
        NascentRecordID = strcat(HumanOntology,num2str(Expression_Null),num2str(CellTypeID),num2str(HumanTissueOnly),...
            num2str(HumanExpGeneOrTransc),num2str(Human_SelTrack),num2str(expValType));
    else
        NascentRecordID = '1';
    end
else
    nascentInfo = [];
end
%% Get Binding Site Mapping and Energy [Edit inside of binding site mapping code for save directory info
%add nascent expression by updating capability of getNascentBindingInfo_JH2
if (settings.withNascent)
    if (strcmp(settings.Organism,'Human'))
        Expression_Null = settings.DoAllGenesHaveSameExpression; %do genes have same expression
        CellTypeID = settings.CellType_ExprID;
        HumanOntology = settings.HumanSpecific.Ontology;%Normal/Cancer
        HumanTissueOnly = settings.HumanSpecific.TissueOrTissueAndCellType; % 0 (TissueOnssue/Cell Type)
        HumanExpGeneOrTransc = settings.HumanSpecific.HumanExpGeneOrTransc; % 1 (Gene/EMBL GENEID) , 0 (Transcript/EMBL Transcript ID)
        Human_SelTrack = settings.HumanSpecific.SCOutputTrack;
        expValType = settings.expressionValType;% 1-4 (expCounts,expValues,mean(CellTypeExpValues),one cell types CellTypeExpValues)
        NascentRecordID = strcat(HumanOntology,num2str(Expression_Null),num2str(CellTypeID),num2str(HumanTissueOnly),...
            num2str(HumanExpGeneOrTransc),num2str(Human_SelTrack),num2str(expValType));
    else
        NascentRecordID = '1';
    end
else
    nascentInfo = [];
end
if (settings.SingleOrMulti==1&&settings.AllIsoforms == 0)%One Gene/One Isoform
    if (id2~=2)
        try
            load([settings.FolderRootName filesep settings.GeneName '_binding_hits_map' designerName '.mat'],'DoesProbeBindSite2','Num_of_Molecule_Sites')
            load([settings.FolderRootName filesep settings.GeneName  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' designerName '.mat'],'Kb_Complement')
            load([settings.FolderRootName filesep settings.GeneName  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' designerName '.mat'],'Kb_mod')
            load([settings.FolderRootName filesep settings.GeneName '_BindingMatrices' designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod','dCp_mod')
            if (settings.BLASTdna)
                load([settings.FolderRootName filesep settings.GeneName '_BindingMatrices2' designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','dCp_Complement')
            end
            if (settings.withNascent)
                load([settings.FolderRootName filesep settings.GeneName '_Tm' num2str(T_hybrid) '_NascentInfo' NascentRecordID designerName '.mat'],'nascentInfo')
            else
                nascentInfo = [];
            end
        catch
            [Kb_mod,Kb_Complement,DoesProbeBindSite2,Num_of_Molecule_Sites,MolProbesAtEvents,nascentInfo,...
                dHeq_mod,dSeq_mod,dHf_mod,dSf_mod,dHr_mod,dSr_mod,Tm_mod,dCp_mod,...
                dHeq_Complement,dSeq_Complement,dHf_Complement,dSf_Complement,dHr_Complement,dSr_Complement,dCp_Complement] = ...
                A_JH_GetSiteMapping_NonTile_V6(probes,settings,gene_table,Kb_Match,dHeq_Match,dSeq_Match,dHf_Match,dSf_Match,dHr_Match,dSr_Match,dCpeq_Match,Tm_Match);
        end
    else
        try
            load([settings.FolderRootName filesep settings.GeneName '_binding_hits_map' designerName0 '.mat'],'DoesProbeBindSite2','Num_of_Molecule_Sites')
            load([settings.FolderRootName filesep settings.GeneName  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' designerName0 '.mat'],'Kb_Complement')
            load([settings.FolderRootName filesep settings.GeneName  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' designerName0 '.mat'],'Kb_mod')
            load([settings.FolderRootName filesep settings.GeneName '_BindingMatrices' designerName0 '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod','dCp_mod')
            if (settings.BLASTdna)
                load([settings.FolderRootName filesep settings.GeneName '_BindingMatrices2' designerName0 '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','dCp_Complement')
            end
            if (settings.withNascent)
                load([settings.FolderRootName filesep settings.GeneName '_Tm' num2str(T_hybrid) '_NascentInfo' NascentRecordID designerName0 '.mat'],'nascentInfo')
            else
                nascentInfo = [];
            end
        catch
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_probes' designerName0 '.mat'],'probes')
            [Kb_mod,Kb_Complement,DoesProbeBindSite2,Num_of_Molecule_Sites,MolProbesAtEvents,nascentInfo,...
                dHeq_mod,dSeq_mod,dHf_mod,dSf_mod,dHr_mod,dSr_mod,Tm_mod,dCp_mod,...
                dHeq_Complement,dSeq_Complement,dHf_Complement,dSf_Complement,dHr_Complement,dSr_Complement,dCp_Complement] = ...
                A_JH_GetSiteMapping_NonTile_V6(probes,settings,gene_table,Kb_Match,dHeq_Match,dSeq_Match,dHf_Match,dSf_Match,dHr_Match,dSr_Match,dCpeq_Match,Tm_Match);
        end
    end


elseif(settings.SingleOrMulti==1&&settings.AllIsoforms == 1)%One Gene/All Isoform
    try
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_binding_hits_map' designerName '.mat'],'DoesProbeBindSite2','Num_of_Molecule_Sites')
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' designerName '.mat'],'Kb_mod')
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' designerName '.mat'],'Kb_Complement')
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_BindingMatrices' designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod','dCp_mod')
        if (settings.BLASTdna)
            load([saveRoot filesep settings.FolderName filesep settings.FolderName '_BindingMatrices2' designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','dCp_Complement')
        end
        if (settings.withNascent)
            load([saveRoot filesep settings.FolderName filesep settings.FolderName '_Tm' num2str(T_hybrid) '_NascentInfo' NascentRecordID designerName '.mat'],'nascentInfo')
        else
            nascentInfo = [];
        end
        needToCombine = 0;
    catch
        Multi_DPS2 = cell(1,size(ENST_RefSeqIDz,1));
        Multi_KbMod = cell(1,size(ENST_RefSeqIDz,1));
        Multi_NMS = cell(1,size(ENST_RefSeqIDz,1));
        Multi_MolProbesAtEvents = cell(1,size(ENST_RefSeqIDz,1));
        Multi_Mol_ProbesAtEventsID = cell(1,size(ENST_RefSeqIDz,1));
        Multi_nascentInfo = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dHeqMod = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dSeqMod = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dHfMod = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dSfMod = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dHrMod = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dSrMod = cell(1,size(ENST_RefSeqIDz,1));
        Multi_dCpMod = cell(1,size(ENST_RefSeqIDz,1));
        Multi_TmMod = cell(1,size(ENST_RefSeqIDz,1));
        if (settings.BLASTdna)
            Multi_KbComp = cell(1,size(ENST_RefSeqIDz,1));
            Multi_dHeq_Complement = cell(1,size(ENST_RefSeqIDz,1));
            Multi_dSeq_Complement = cell(1,size(ENST_RefSeqIDz,1));
            Multi_dHf_Complement = cell(1,size(ENST_RefSeqIDz,1));
            Multi_dSf_Complement = cell(1,size(ENST_RefSeqIDz,1));
            Multi_dHr_Complement = cell(1,size(ENST_RefSeqIDz,1));
            Multi_dSr_Complement = cell(1,size(ENST_RefSeqIDz,1));
            Multi_dCp_Complement = cell(1,size(ENST_RefSeqIDz,1));
        end
        for k = 1:size(ENST_RefSeqIDz,1)
            try
                load([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_binding_hits_map' designerName '.mat'],'DoesProbeBindSite2','Num_of_Molecule_Sites','MolProbesAtEvents','Mol_ProbesAtEventsID')
                load([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' designerName '.mat'],'Kb_mod')
                load([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_BindingMatrices' designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod','dCp_mod')
                load([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' designerName '.mat'],'Kb_Complement')
                if (settings.BLASTdna)
                    load([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k}  '_BindingMatrices2' designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','dCp_Complement')
                end
                if (settings.withNascent)
                    load([saveRoot filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} filesep '(' ENST_GeneName ')_' ENST_RefSeqIDz{k} '_Tm' num2str(T_hybrid) '_NascentInfo' NascentRecordID designerName '.mat'],'nascentInfo')
                else
                    nascentInfo = [];
                end
            catch
                settings.GeneName = extractBefore(FoldName{k},'_');
                settings.transcript_IDs = extractAfter(FoldName{k},'_');
                settings.ChrNum = All_chrom{k};
                settings.GeneChr = strcat('chr',All_chrom{k});
                [Kb_mod,Kb_Complement,DoesProbeBindSite2,Num_of_Molecule_Sites,MolProbesAtEvents,Mol_ProbesAtEventsID,nascentInfo,...
                    dHeq_mod,dSeq_mod,dHf_mod,dSf_mod,dHr_mod,dSr_mod,Tm_mod,dCp_mod,...
                    dHeq_Complement,dSeq_Complement,dHf_Complement,dSf_Complement,dHr_Complement,dSr_Complement,dCp_Complement] = ...
                    A_JH_GetSiteMapping_NonTile_V6(Isoform_probes{k},settings,Isoform_gene_table{k},Multi_Kb_Match{k},Multi_dHeq_Match{k},Multi_dSeq_Match{k},Multi_dHf_Match{k},Multi_dSf_Match{k},Multi_dHr_Match{k},Multi_dSr_Match{k},Multi_Tm_Match{k});
            end
            Multi_DPS2{k} = DoesProbeBindSite2;
            Multi_KbMod{k} = Kb_mod;
            Multi_KbComp{k} = Kb_Complement;
            Multi_NMS{k} = Num_of_Molecule_Sites;
            Multi_MolProbesAtEvents{k} = MolProbesAtEvents;
            Multi_Mol_ProbesAtEventsID{k} = Mol_ProbesAtEventsID;
            Multi_nascentInfo{k} = nascentInfo;
            Multi_dHeqMod{k} = dHeq_mod;
            Multi_dSeqMod{k} = dSeq_mod;
            Multi_dHfMod{k} = dHf_mod;
            Multi_dSfMod{k} = dSf_mod;
            Multi_dHrMod{k} = dHr_mod;
            Multi_dSrMod{k} = dSr_mod;
            Multi_dCpMod{k} = dCp_mod;
            Multi_TmMod{k} = Tm_mod;
            if (settings.BLASTdna)
                Multi_dHeq_Complement{k} = dHeq_Complement;
                Multi_dSeq_Complement{k} = dSeq_Complement;
                Multi_dHf_Complement{k} = dHf_Complement;
                Multi_dSf_Complement{k} = dSf_Complement;
                Multi_dHr_Complement{k} = dHr_Complement;
                Multi_dSr_Complement{k} = dSr_Complement;
                Multi_dCp_Complement{k} = dCp_Complement;
            end
        end
    end
elseif (settings.SingleOrMulti==0&&settings.AllIsoforms == 0)%Multi Gene/One Isoform
    try
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_binding_hits_map' designerName '.mat'],'DoesProbeBindSite2','Num_of_Molecule_Sites')
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' designerName '.mat'],'Kb_mod')
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_BindingMatrices' designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod','dCp_mod')
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' designerName '.mat'],'Kb_Complement')
        if (settings.BLASTdna)
            load([saveRoot filesep settings.FolderName filesep settings.FolderName '_BindingMatrices2' designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','dCp_Complement')
        end
        if (settings.withNascent)
            load([saveRoot filesep settings.FolderName filesep settings.FolderName '_Tm' num2str(T_hybrid) '_NascentInfo' NascentRecordID designerName '.mat'],'nascentInfo')
        else
            nascentInfo = [];
        end
        needToCombine = 0;
    catch
        Multi_DPS2 = cell(1,size(inputs1{gene_num,1},2));
        Multi_KbMod = cell(1,size(inputs1{gene_num,1},2));
        Multi_NMS = cell(1,size(inputs1{gene_num,1},2));
        Multi_MolProbesAtEvents = cell(1,size(inputs1{gene_num,1},2));
        Multi_Mol_ProbesAtEventsID = cell(1,size(inputs1{gene_num,1},2));
        Multi_nascentInfo = cell(1,size(inputs1{gene_num,1},2));
        Multi_dHeqMod = cell(1,size(inputs1{gene_num,1},2));
        Multi_dSeqMod = cell(1,size(inputs1{gene_num,1},2));
        Multi_dHfMod = cell(1,size(inputs1{gene_num,1},2));
        Multi_dSfMod = cell(1,size(inputs1{gene_num,1},2));
        Multi_dHrMod = cell(1,size(inputs1{gene_num,1},2));
        Multi_dSrMod = cell(1,size(inputs1{gene_num,1},2));
        Multi_dCpMod = cell(1,size(inputs1{gene_num,1},2));
        Multi_TmMod = cell(1,size(inputs1{gene_num,1},2));
        if (settings.BLASTdna)
            Multi_KbComp = cell(1,size(inputs1{gene_num,1},2));
            Multi_dHeq_Complement = cell(1,size(inputs1{gene_num,1},2));
            Multi_dSeq_Complement = cell(1,size(inputs1{gene_num,1},2));
            Multi_dHf_Complement = cell(1,size(inputs1{gene_num,1},2));
            Multi_dSf_Complement = cell(1,size(inputs1{gene_num,1},2));
            Multi_dHr_Complement = cell(1,size(inputs1{gene_num,1},2));
            Multi_dSr_Complement = cell(1,size(inputs1{gene_num,1},2));
            Multi_dCp_Complement = cell(1,size(inputs1{gene_num,1},2));
        end
        for k = 1:size(inputs1{gene_num,1},2)
            try
                load([saveRoot filesep FoldName{k} filesep extractBefore(FoldName{k},'_') '_binding_hits_map' designerName '.mat'],'DoesProbeBindSite2','Num_of_Molecule_Sites','MolProbesAtEvents','Mol_ProbesAtEventsID')
                load([saveRoot filesep FoldName{k} filesep extractBefore(FoldName{k},'_')  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' designerName '.mat'],'Kb_mod')
                load([saveRoot filesep FoldName{k} filesep extractBefore(FoldName{k},'_') '_BindingMatrices' designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod','dCp_mod')
                load([saveRoot filesep FoldName{k} filesep extractBefore(FoldName{k},'_')  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' designerName '.mat'],'Kb_Complement')
                if (settings.BLASTdna)
                    load([saveRoot filesep FoldName{k} filesep extractBefore(FoldName{k},'_') '_BindingMatrices2' designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','dCp_Complement')
                end
                if (settings.withNascent)
                    load([saveRoot filesep FoldName{k} filesep extractBefore(FoldName{k},'_') '_Tm' num2str(T_hybrid) '_NascentInfo' NascentRecordID designerName '.mat'],'nascentInfo')
                else
                    nascentInfo = [];
                end
            catch
                settings.GeneName = extractBefore(FoldName{k},'_');
                settings.transcript_IDs = extractAfter(FoldName{k},'_');
                settings.ChrNum = All_chrom{k};
                settings.GeneChr = strcat('chr',All_chrom{k});
                [Kb_mod,Kb_Complement,DoesProbeBindSite2,Num_of_Molecule_Sites,MolProbesAtEvents,Mol_ProbesAtEventsID,nascentInfo,...
                    dHeq_mod,dSeq_mod,dHf_mod,dSf_mod,dHr_mod,dSr_mod,Tm_mod,dCp_mod,...
                    dHeq_Complement,dSeq_Complement,dHf_Complement,dSf_Complement,dHr_Complement,dSr_Complement,dCp_Complement] = ...
                    A_JH_GetSiteMapping_NonTile_V6(Transcript_probes{k},settings,Transcript_gene_table{k},Multi_Kb_Match{k},Multi_dHeq_Match{k},Multi_dSeq_Match{k},Multi_dHf_Match{k},Multi_dSf_Match{k},Multi_dHr_Match{k},Multi_dSr_Match{k},Multi_Tm_Match{k});
            end
            Multi_DPS2{k} = DoesProbeBindSite2;
            Multi_KbMod{k} = Kb_mod;
            Multi_NMS{k} = Num_of_Molecule_Sites;
            Multi_MolProbesAtEvents{k} = MolProbesAtEvents;
            Multi_Mol_ProbesAtEventsID{k} = Mol_ProbesAtEventsID;
            Multi_nascentInfo{k} = nascentInfo;
            Multi_dHeqMod{k} = dHeq_mod;
            Multi_dSeqMod{k} = dSeq_mod;
            Multi_dHfMod{k} = dHf_mod;
            Multi_dSfMod{k} = dSf_mod;
            Multi_dHrMod{k} = dHr_mod;
            Multi_dSrMod{k} = dSr_mod;
            Multi_dCpMod{k} = dCp_mod;
            Multi_TmMod{k} = Tm_mod;
            if (settings.BLASTdna)
                Multi_KbComp{k} = Kb_Complement;
                Multi_dHeq_Complement{k} = dHeq_Complement;
                Multi_dSeq_Complement{k} = dSeq_Complement;
                Multi_dHf_Complement{k} = dHf_Complement;
                Multi_dSf_Complement{k} = dSf_Complement;
                Multi_dHr_Complement{k} = dHr_Complement;
                Multi_dSr_Complement{k} = dSr_Complement;
                Multi_dCp_Complement{k} = dCp_Complement;
            end
        end
    end
elseif (settings.SingleOrMulti==0&&settings.AllIsoforms == 1)%Multi Gene/All Isoform
    try
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_binding_hits_map' designerName '.mat'],'DoesProbeBindSite2','Num_of_Molecule_Sites')
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' designerName '.mat'],'Kb_mod')
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_BindingMatrices' designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod','dCp_mod')
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' designerName '.mat'],'Kb_Complement')
        if (settings.BLASTdna)
            load([saveRoot filesep settings.FolderName filesep settings.FolderName  '_BindingMatrices2' designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','dCp_Complement')
        end
        if (settings.withNascent)
            load([saveRoot filesep settings.FolderName filesep settings.FolderName '_Tm' num2str(T_hybrid) '_NascentInfo' NascentRecordID designerName '.mat'],'nascentInfo')
        else
            nascentInfo = [];
        end
        needToCombine = 0;
    catch
        Multi_DPS2 = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_KbMod = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_NMS = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_MolProbesAtEvents = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_Mol_ProbesAtEventsID{k} = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_nascentInfo = cell(1,size(inputs1{gene_num,1},2));
        Multi_dHeqMod = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dSeqMod = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dHfMod = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dSfMod = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dHrMod = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dSrMod = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_dCpMod = cell(1,size(All_ENST_RefSeqIDz,1));
        Multi_TmMod = cell(1,size(All_ENST_RefSeqIDz,1));
        if (settings.BLASTdna)
            Multi_KbComp = cell(1,size(All_ENST_RefSeqIDz,1));
            Multi_dHeq_Complement = cell(1,size(All_ENST_RefSeqIDz,1));
            Multi_dSeq_Complement = cell(1,size(All_ENST_RefSeqIDz,1));
            Multi_dHf_Complement = cell(1,size(All_ENST_RefSeqIDz,1));
            Multi_dSf_Complement = cell(1,size(All_ENST_RefSeqIDz,1));
            Multi_dHr_Complement = cell(1,size(All_ENST_RefSeqIDz,1));
            Multi_dSr_Complement = cell(1,size(All_ENST_RefSeqIDz,1));
            Multi_dCp_Complement = cell(1,size(All_ENST_RefSeqIDz,1));
        end
        for k = 1:size(inputs1{gene_num,1},2)
            for k2 = 1:length(Multi_ENST_RefSeqIDz{k})
                try
                    load([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_binding_hits_map' designerName '.mat'],'DoesProbeBindSite2','Num_of_Molecule_Sites','MolProbesAtEvents')
                    load([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' designerName '.mat'],'Kb_mod')
                    load([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_BindingMatrices' designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod','dCp_mod')
                    load([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' designerName '.mat'],'Kb_Complement')
                    if (settings.BLASTdna)
                        load([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_BindingMatrices2' designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','dCp_Complement')
                    end
                    if (settings.withNascent)
                        load([saveRoot filesep FoldName{kc(k,k2)} filesep FoldName{kc(k,k2)} '_Tm' num2str(T_hybrid) '_NascentInfo' NascentRecordID designerName '.mat'],'nascentInfo')
                    else
                        nascentInfo = [];
                    end
                catch
                    settings.GeneName = extractBefore(FoldName{kc(k,k2)},'_');
                    settings.transcript_IDs = extractAfter(FoldName{kc(k,k2)},'_');
                    settings.ChrNum = All_chrom{k}{k2};
                    settings.GeneChr = strcat('chr',All_chrom{k}{k2});
                    [Kb_mod,Kb_Complement,DoesProbeBindSite2,Num_of_Molecule_Sites,MolProbesAtEvents,Mol_ProbesAtEventsID,nascentInfo,...
                        dHeq_mod,dSeq_mod,dHf_mod,dSf_mod,dHr_mod,dSr_mod,Tm_mod,dCp_mod,...
                        dHeq_Complement,dSeq_Complement,dHf_Complement,dSf_Complement,dHr_Complement,dSr_Complement,dCp_Complement] = ...
                        A_JH_GetSiteMapping_NonTile_V6(TranscriptIsoform_probes{kc(k,k2)},settings,TranscriptIsoform_gene_table{kc(k,k2)},Multi_Kb_Match{kc(k,k2)},Multi_dHeq_Match{kc(k,k2)},Multi_dSeq_Match{kc(k,k2)},Multi_dHf_Match{kc(k,k2)},Multi_dSf_Match{kc(k,k2)},Multi_dHr_Match{kc(k,k2)},Multi_dSr_Match{kc(k,k2)},Multi_Tm_Match{kc(k,k2)});
                end
                Multi_DPS2{kc(k,k2)} = DoesProbeBindSite2;
                Multi_KbMod{kc(k,k2)} = Kb_mod;
                Multi_NMS{kc(k,k2)} = Num_of_Molecule_Sites;
                Multi_MolProbesAtEvents{kc(k,k2)} = MolProbesAtEvents;
                Multi_Mol_ProbesAtEventsID{kc(k,k2)} = Mol_ProbesAtEventsID;
                Multi_nascentInfo{kc(k,k2)} = nascentInfo;
                Multi_dHeqMod{kc(k,k2)} = dHeq_mod;
                Multi_dSeqMod{kc(k,k2)} = dSeq_mod;
                Multi_dHfMod{kc(k,k2)} = dHf_mod;
                Multi_dSfMod{kc(k,k2)} = dSf_mod;
                Multi_dHrMod{kc(k,k2)} = dHr_mod;
                Multi_dSrMod{kc(k,k2)} = dSr_mod;
                Multi_dCpMod{kc(k,k2)} = dCp_mod;
                Multi_TmMod{kc(k,k2)} = Tm_mod;
                if (settings.BLASTdna)
                    Multi_KbComp{kc(k,k2)} = Kb_Complement;
                    Multi_dHeq_Complement{kc(k,k2)} = dHeq_Complement;
                    Multi_dSeq_Complement{kc(k,k2)} = dSeq_Complement;
                    Multi_dHf_Complement{kc(k,k2)} = dHf_Complement;
                    Multi_dSf_Complement{kc(k,k2)} = dSf_Complement;
                    Multi_dHr_Complement{kc(k,k2)} = dHr_Complement;
                    Multi_dSr_Complement{kc(k,k2)} = dSr_Complement;
                    Multi_dCp_Complement{kc(k,k2)} = dCp_Complement;
                end
            end
        end
    end
end
%code that stores location info in 3d map [SiteLoc]
%code that re-assigns probe and target ID in 3d-4d maps
% code that combines maps
% code that designs probes [mix designer with V7_Out codes]
% solves eq/kinetics/metrics with multi-transcript output
% all forked with/without temp changes or gradients.

%might turn 4D multiTargetLocations into two 3d matricies and change
%combineMaps code to for those
if (settings.SingleOrMulti==1&&settings.AllIsoforms == 0)%One Gene/One Isoform
elseif (needToCombine)
    if(settings.SingleOrMulti==1&&settings.AllIsoforms == 1)%One Gene/All Isoform
        Multi_gene_table3 = Isoform_gene_table3;
        All_DPS_Names = AllIsoforms_DPS_Names;
        All_probes = AllIsoforms_probes;
    elseif (settings.SingleOrMulti==0&&settings.AllIsoforms == 0)%Multi Gene/One Isoform
        Multi_gene_table3 = Transcript_gene_table3;
        All_DPS_Names = AllTranscript_DPS_Names;
        All_probes = AllTranscripts_probes;
    elseif (settings.SingleOrMulti==0&&settings.AllIsoforms == 1)%Multi Gene/All Isoform
        Multi_gene_table3 = TranscriptIsoform_gene_table3;
        All_DPS_Names = AllTranscriptIsoform_DPS_Names;
        All_probes = AllTranscriptIsoforms_probes;
    end
    Multi_TargetSiteLocations = cell(1,length(Multi_MolProbesAtEvents));
    for k = 1:length(Multi_MolProbesAtEvents)
        Multi_TargetSiteLocations{k} = ndSparse.build([length(All_probes),size(Multi_DPS2{k},2),size(Multi_DPS2{k},3),2],0);%P T S L
        filtMolProbesAtEvents = cell(1,size(Multi_DPS2{k},2));
        for T=1:size(Multi_DPS2{k},2) %ALL SITES/PROBES IN IT
            %Filter check if Probe in DPS2 for target or in
            filtMolProbesAtEvents{T} = arrayfun(@(S) Multi_MolProbesAtEvents{k}{T}{S}(Multi_DPS2{k}(Multi_MolProbesAtEvents{k}{T}{S},T,S)==1),1:length(Multi_MolProbesAtEvents{k}{T}),'Un',0);
            Sz = cell2mat(arrayfun(@(S) S*ones(1,length(Multi_Mol_ProbesAtEventsID{k}{T}{S})),1:length(Multi_MolProbesAtEvents{k}{T}),'Un',0));
            Pz = cell2mat(arrayfun(@(S) Multi_gene_table3{k}.ProbeNum(Multi_Mol_ProbesAtEventsID{k}{T}{S})',1:length(Multi_MolProbesAtEvents{k}{T}),'Un',0));
            Az = cell2mat(arrayfun(@(S) Multi_gene_table3{k}.Ax(Multi_Mol_ProbesAtEventsID{k}{T}{S})',1:length(Multi_MolProbesAtEvents{k}{T}),'Un',0));
            Bz = cell2mat(arrayfun(@(S) Multi_gene_table3{k}.Bx(Multi_Mol_ProbesAtEventsID{k}{T}{S})',1:length(Multi_MolProbesAtEvents{k}{T}),'Un',0));
            if (~isempty(Pz))
                vecA = sub2ind(size(Multi_TargetSiteLocations{k}),Pz,T*ones(1,length(Sz)),Sz,ones(1,length(Sz)));
                vecB = sub2ind(size(Multi_TargetSiteLocations{k}),Pz,T*ones(1,length(Sz)),Sz,2*ones(1,length(Sz)));
                Multi_TargetSiteLocations{k}(vecA) = Az;
                Multi_TargetSiteLocations{k}(vecB) = Bz;
            end
        end
    end
    updated_Multi_DPS2 = cell(1,length(Multi_MolProbesAtEvents)); %get index of all non-zero elements.
    updated_Multi_KbMod = cell(1,length(Multi_MolProbesAtEvents));%index -> [P,T] -> [P',T'] ->index'
    updated_Multi_TargetSiteLocations = cell(1,length(Multi_MolProbesAtEvents));
    updated_Multi_Kon = cell(1,length(Multi_MolProbesAtEvents));
    updated_Multi_Koff = cell(1,length(Multi_MolProbesAtEvents));
    updated_Multi_dHeqMod = cell(1,length(Multi_MolProbesAtEvents));
    updated_Multi_dSeqMod = cell(1,length(Multi_MolProbesAtEvents));
    updated_Multi_dHfMod = cell(1,length(Multi_MolProbesAtEvents));
    updated_Multi_dSfMod = cell(1,length(Multi_MolProbesAtEvents));
    updated_Multi_dHrMod = cell(1,length(Multi_MolProbesAtEvents));
    updated_Multi_dSrMod = cell(1,length(Multi_MolProbesAtEvents));
    updated_Multi_dCpMod = cell(1,length(Multi_MolProbesAtEvents));
    updated_Multi_TmMod = cell(1,length(Multi_MolProbesAtEvents));
    if (settings.BLASTdna)
        updated_Multi_KbComp = cell(1,length(Multi_MolProbesAtEvents));        %subsasign(obj,S,rhs)
        updated_Multi_dHeq_Complement = cell(1,length(Multi_MolProbesAtEvents));
        updated_Multi_dSeq_Complement = cell(1,length(Multi_MolProbesAtEvents));
        updated_Multi_dHf_Complement = cell(1,length(Multi_MolProbesAtEvents));
        updated_Multi_dSf_Complement = cell(1,length(Multi_MolProbesAtEvents));
        updated_Multi_dHr_Complement = cell(1,length(Multi_MolProbesAtEvents));
        updated_Multi_dSr_Complement = cell(1,length(Multi_MolProbesAtEvents));
        updated_Multi_dCp_Complement = cell(1,length(Multi_MolProbesAtEvents));
    end
    for k = 1:size(inputs1{gene_num,1},2)
        dictP = containers.Map(UpdateProbeNumAssociation{k}(:,1),UpdateProbeNumAssociation{k}(:,2));
        dictT = containers.Map(UpdateTargetNumAssociation{k}(:,1),UpdateTargetNumAssociation{k}(:,2));
        updated_Multi_Kon{k} = ndSparse.build([size(All_probes,1),N_methods],0);
        updated_Multi_Koff{k} = ndSparse.build([size(All_probes,1),length(All_DPS_Names),N_methods],0);
        updated_Multi_DPS2{k} = ndSparse.build([size(All_probes,1),length(All_DPS_Names),size(Multi_DPS2{k},3)],0);
        updated_Multi_TargetSiteLocations{k} = ndSparse.build([size(All_probes,1),length(All_DPS_Names),max(cellfun(@(x) size(x,3),Multi_TargetSiteLocations)),2],0);
        if (settings.BLASTdna)
            updated_Multi_KbComp{k} = ndSparse.build([length(All_DPS_Names), max(cellfun(@(x) size(x,2),Multi_KbComp)),N_methods],0);
            updated_Multi_dHeq_Complement{k} = ndSparse.build([length(All_DPS_Names),size(Multi_dHeq_Complement{k},2),N_methods],0);
            updated_Multi_dSeq_Complement{k} = ndSparse.build([length(All_DPS_Names),size(Multi_dSeq_Complement{k},2),N_methods],0);
            updated_Multi_dCp_Complement{k} = ndSparse.build([length(All_DPS_Names),size(Multi_dCp_Complement{k},2),N_methods],0);
            updated_Multi_dHf_Complement{k} = ndSparse.build([length(All_DPS_Names),size(Multi_dHf_Complement{k},2),N_methods2],0);
            updated_Multi_dSf_Complement{k} = ndSparse.build([length(All_DPS_Names),size(Multi_dSf_Complement{k},2),N_methods2],0);
            updated_Multi_dHr_Complement{k} = ndSparse.build([length(All_DPS_Names),size(Multi_dHr_Complement{k},2),N_methods2],0);
            updated_Multi_dSr_Complement{k} = ndSparse.build([length(All_DPS_Names),size(Multi_dSr_Complement{k},2),N_methods2],0);
        end
        updated_Multi_KbMod{k} = ndSparse.build([size(All_probes,1),length(All_DPS_Names),size(Multi_KbMod{k},3),N_methods],0);
        updated_Multi_TmMod{k} = ndSparse.build([size(All_probes,1),length(All_DPS_Names),size(Multi_TmMod{k},3),N_methods+1],0);
        updated_Multi_dHeqMod{k} = ndSparse.build([size(All_probes,1),length(All_DPS_Names),size(Multi_dHeqMod{k},3),N_methods],0);
        updated_Multi_dSeqMod{k} = ndSparse.build([size(All_probes,1),length(All_DPS_Names),size(Multi_dSeqMod{k},3),N_methods],0);
        updated_Multi_dCpMod{k} = ndSparse.build([size(All_probes,1),length(All_DPS_Names),size(Multi_dCpMod{k},3),N_methods],0);
        updated_Multi_dHfMod{k} = ndSparse.build([size(All_probes,1),length(All_DPS_Names),size(Multi_dHfMod{k},3),N_methods2],0);
        updated_Multi_dSfMod{k} = ndSparse.build([size(All_probes,1),length(All_DPS_Names),size(Multi_dSfMod{k},3),N_methods2],0);
        updated_Multi_dHrMod{k} = ndSparse.build([size(All_probes,1),length(All_DPS_Names),size(Multi_dHrMod{k},3),N_methods2],0);
        updated_Multi_dSrMod{k} = ndSparse.build([size(All_probes,1),length(All_DPS_Names),size(Multi_dSrMod{k},3),N_methods2],0);

        updated_Multi_Kon{k}(UpdateProbeNumAssociation{k}(:,2),:) = Multi_Kon{k}(UpdateProbeNumAssociation{k}(:,1),:);
        linVec0 = find(Multi_Koff{k});
        [pVec0,tVec0,mVec0] = ind2sub(size(Multi_Koff{k}),linVec0);
        pVecPrime0 = arrayfun(@(x) dictP(pVec0(x)),1:length(pVec0));
        tVecPrime0 = arrayfun(@(x) dictT(tVec0(x)),1:length(tVec0));
        nonVec0 = sub2ind(size(updated_Multi_Koff{k}),pVecPrime0,tVecPrime0,mVec0');
        updated_Multi_Koff{k}(nonVec0) = Multi_Koff{k}(linVec0);

        linVec = find(Multi_DPS2{k});
        [pVec,tVec,sVec] = ind2sub(size(Multi_DPS2{k}),linVec);
        pVecPrime = arrayfun(@(x) dictP(pVec(x)),1:length(pVec));
        tVecPrime = arrayfun(@(x) dictT(tVec(x)),1:length(tVec));
        nonVec = sub2ind(size(updated_Multi_DPS2{k}),pVecPrime,tVecPrime,sVec');
        updated_Multi_DPS2{k}(nonVec) = Multi_DPS2{k}(linVec);

        linVec1 = find(Multi_TargetSiteLocations{k});
        [pVec1,tVec1,sVec1i,sVec1o] = ind2sub(size(Multi_TargetSiteLocations{k}),linVec1);
        tVecPrime1 = arrayfun(@(x) dictT(tVec1(x)),1:length(tVec1));
        nonVec1 = sub2ind(size(updated_Multi_TargetSiteLocations{k}),pVec1,tVecPrime1',sVec1i,sVec1o)';
        updated_Multi_TargetSiteLocations{k}(nonVec1) = Multi_TargetSiteLocations{k}(linVec1);%slow?

        linVec2 = find(Multi_KbMod{k});
        [pVec2,tVec2,sVec2,mVec2] = ind2sub(size(Multi_KbMod{k}),linVec2);
        pVecPrime2 = arrayfun(@(x) dictP(pVec2(x)),1:length(pVec2));
        tVecPrime2 = arrayfun(@(x) dictT(tVec2(x)),1:length(tVec2));
        nonVec2 = sub2ind(size(updated_Multi_KbMod{k}),pVecPrime2,tVecPrime2,sVec2,mVec2);
        updated_Multi_KbMod{k}(nonVec2) = Multi_KbMod{k}(linVec2);
        linVec3 = find(Multi_TmMod{k});
        [pVec3,tVec3,sVec3,mVec3] = ind2sub(size(Multi_TmMod{k}),linVec3);
        pVecPrime3 = arrayfun(@(x) dictP(pVec3(x)),1:length(pVec3));
        tVecPrime3 = arrayfun(@(x) dictT(tVec3(x)),1:length(tVec3));
        nonVec3 = sub2ind(size(updated_Multi_TmMod{k}),pVecPrime3,tVecPrime3,sVec3,mVec3);
        updated_Multi_TmMod{k}(nonVec3) = Multi_TmMod{k}(linVec3);

        if (settings.BLASTdna)
            linVec4 = find(Multi_KbComp{k});
            [tVec4,sVec4,mVec4] = ind2sub(size(Multi_KbComp{k}),linVec4);
            tVecPrime4 = arrayfun(@(x) dictT(tVec4(x)),1:length(tVec4));
            nonVec4 = sub2ind(size(updated_Multi_KbComp{k}),tVecPrime4,sVec4',mVec4');%slow?
            updated_Multi_KbComp{k}(nonVec4) = Multi_KbComp{k}(linVec4);
        end
        linVec5A = find(Multi_dHeqMod{k});
        linVec5B = find(Multi_dSeqMod{k});
        linVec5C = find(Multi_dCpMod{k});
        linVec5 = unique([linVec5A linVec5B linVec5C]);
        [pVec5,tVec5,sVec5,mVec5] = ind2sub(size(Multi_dHeqMod{k}),linVec5);
        pVecPrime5 = arrayfun(@(x) dictP(pVec5(x)),1:length(pVec5));
        tVecPrime5 = arrayfun(@(x) dictT(tVec5(x)),1:length(tVec5));
        nonVec5 = sub2ind(size(updated_Multi_dHeqMod{k}),pVecPrime5,tVecPrime5,sVec5,mVec5);
        updated_Multi_dHeqMod{k}(nonVec5) = Multi_dHeqMod{k}(linVec5);
        updated_Multi_dSeqMod{k}(nonVec5) = Multi_dSeqMod{k}(linVec5);
        updated_Multi_dCpMod{k}(nonVec5) = Multi_dCpMod{k}(linVec5);
        linVec6A = find(Multi_dHfMod{k});
        linVec6B = find(Multi_dSfMod{k});
        linVec6C = find(Multi_dHrMod{k});
        linVec6D = find(Multi_dSrMod{k});
        linVec6 = unique([linVec6A linVec6B linVec6C linVec6D]);
        [pVec6,tVec6,sVec6,mVec6] = ind2sub(size(Multi_dHfMod{k}),linVec6);
        pVecPrime6 = arrayfun(@(x) dictP(pVec6(x)),1:length(pVec6));
        tVecPrime6 = arrayfun(@(x) dictT(tVec6(x)),1:length(tVec6));
        nonVec6 = sub2ind(size(updated_Multi_dHfMod{k}),pVecPrime6,tVecPrime6,sVec6,mVec6);
        updated_Multi_dHfMod{k}(nonVec6) = Multi_dHfMod{k}(linVec6);
        updated_Multi_dSfMod{k}(nonVec6) = Multi_dSfMod{k}(linVec6);
        updated_Multi_dHrMod{k}(nonVec6) = Multi_dHrMod{k}(linVec6);
        updated_Multi_dSrMod{k}(nonVec6) = Multi_dSrMod{k}(linVec6);
        if (settings.BLASTdna)
            linVec7A = find(Multi_dHeq_Complement{k});
            linVec7B = find(Multi_dSeqMod{k});
            linVec7C = find(Multi_dCpMod{k});
            linVec7 = unique([linVec7A linVec7B linVec7C]);
            [tVec7,sVec7,mVec7] = ind2sub(size(Multi_dHeq_Complement{k}),linVec7);
            tVecPrime7 = arrayfun(@(x) dictT(tVec7(x)),1:length(tVec7));
            nonVec7 = sub2ind(size(updated_Multi_dHeq_Complement{k}),tVecPrime7,sVec7,mVec7);
            updated_Multi_dHeq_Complement{k}(nonVec7) = Multi_dHeq_Complement{k}(linVec7);
            updated_Multi_dSeq_Complement{k}(nonVec7) = Multi_dSeq_Complement{k}(linVec7);
            updated_Multi_dCp_Complement{k}(nonVec7) = Multi_dCp_Complement{k}(linVec7);
            linVec8A = find(Multi_dHf_Complement{k});
            linVec8B = find(Multi_dSf_Complement{k});
            linVec8C = find(Multi_dHr_Complement{k});
            linVec8D = find(Multi_dSr_Complement{k});
            linVec8 = unique([linVec8A linVec8B linVec8C linVec8D]);
            [tVec8,sVec8,mVec8] = ind2sub(size(Multi_dHf_Complement{k}),linVec8);
            tVecPrime8 = arrayfun(@(x) dictT(tVec8(x)),1:length(tVec8));
            nonVec8 = sub2ind(size(updated_Multi_dHfMod{k}),tVecPrime8,sVec8,mVec8);
            updated_Multi_dHf_Complement{k}(nonVec8) = Multi_dHf_Complement{k}(linVec8);
            updated_Multi_dSf_Complement{k}(nonVec8) = Multi_dSf_Complement{k}(linVec8);
            updated_Multi_dHr_Complement{k}(nonVec8) = Multi_dHr_Complement{k}(linVec8);
            updated_Multi_dSr_Complement{k}(nonVec8) = Multi_dSr_Complement{k}(linVec8);
        end
    end
end
if (settings.SingleOrMulti==1&&settings.AllIsoforms == 0)%One Gene/One Isoform
elseif (needToCombine)
    MapMultiVars = cell(1,length(Multi_gene_table3));
    for k = 1:length(Multi_gene_table3)
        MapMultiVars{k}{1} = updated_Multi_KbMod{k};
        MapMultiVars{k}{2} = updated_Multi_TmMod{k};
        MapMultiVars{k}{3} = updated_Multi_dHeqMod{k};
        MapMultiVars{k}{4} = updated_Multi_dSeqMod{k};
        MapMultiVars{k}{5} = updated_Multi_dCpMod{k};
        MapMultiVars{k}{6} = updated_Multi_dHfMod{k};
        MapMultiVars{k}{7} = updated_Multi_dSfMod{k};
        MapMultiVars{k}{8} = updated_Multi_dHrMod{k};
        MapMultiVars{k}{9} = updated_Multi_dSrMod{k};
        if (settings.BLASTdna)
            MapMultiVars{k}{10} = updated_Multi_KbComp{k};
            MapMultiVars{k}{11} = updated_Multi_dHeq_Complement{k};
            MapMultiVars{k}{12} = updated_Multi_dSeq_Complement{k};
            MapMultiVars{k}{13} = updated_Multi_dCp_Complement{k};
            MapMultiVars{k}{14} = updated_Multi_dHf_Complement{k};
            MapMultiVars{k}{15} = updated_Multi_dSf_Complement{k};
            MapMultiVars{k}{16} = updated_Multi_dHr_Complement{k};
            MapMultiVars{k}{17} = updated_Multi_dSr_Complement{k};
        end
    end
    [All_DPS, All_TargetSiteLocations,All_Vars,All_Num_of_Molecule_Sites] = combineMultipleBindingSiteMaps_V2(updated_Multi_DPS2, updated_Multi_TargetSiteLocations,MapMultiVars);
    All_KbMod = All_Vars{1};
    All_TmMod =  All_Vars{2};
    All_dHeqMod = All_Vars{3};
    All_dSeqMod = All_Vars{4};
    All_dCpMod = All_Vars{5};
    All_dHfMod = All_Vars{6};%ee
    All_dSfMod = All_Vars{7};
    All_dHrMod = All_Vars{8};
    All_dSrMod = All_Vars{9};
    if (settings.BLASTdna)
        All_KbComp = All_Vars{10};
        All_dHeq_Complement = All_Vars{11};
        All_dSeq_Complement = All_Vars{12};
        All_dCp_Complement = All_Vars{13};
        All_dHf_Complement = All_Vars{14};
        All_dSf_Complement = All_Vars{15};
        All_dHr_Complement = All_Vars{16};
        All_dSr_Complement = All_Vars{17};
    end
    Kon = max(horzcat(updated_Multi_Kon{:}),[],2);
    Koff = ndSparse.build([size(All_probes,1),length(All_DPS_Names),length(Multi_gene_table3)],0);
    for ki = 1:length(Multi_gene_table3)
        Koff(:,:,ki) = updated_Multi_Koff{ki};
    end
    Koff = max(Koff,[],3);
    DoesProbeBindSite2 = All_DPS;
    Num_of_Molecule_Sites = All_Num_of_Molecule_Sites;
    Kb_mod = All_KbMod;
    Kb_Complement = All_KbComp;
    TargetSiteLocations = All_TargetSiteLocations;
    if (settings.BLASTdna)
        dHeq_Complement = All_dHeq_Complement;
        dSeq_Complement = All_dSeq_Complement;
        dCp_Complement = All_dCp_Complement;
        dHf_Complement = All_dHf_Complement;
        dSf_Complement = All_dSf_Complement;
        dHr_Complement = All_dHr_Complement;
        dSr_Complement = All_dSr_Complement;
    end
    dHeq_mod = All_dHeqMod;
    dSeq_mod = All_dSeqMod;
    dCp_mod = All_dCpMod;
    dHf_mod = All_dHfMod;
    dSf_mod = All_dSfMod;
    dHr_mod = All_dHrMod;
    dSr_mod = All_dSrMod;
    Tm_mod = All_TmMod;
    save([saveRoot filesep settings.FolderName filesep settings.FolderName '_hits_table' designerName '.mat'],'gene_table','-v7.3')
    save([saveRoot filesep settings.FolderName filesep settings.FolderName '_binding_hits_map' designerName '.mat'],'DoesProbeBindSite2','Num_of_Molecule_Sites','-v7.3')
    save([saveRoot filesep settings.FolderName filesep settings.FolderName '_TargetSiteLocations' designerName '.mat'],'TargetSiteLocations','-v7.3')
    save([saveRoot filesep settings.FolderName filesep settings.FolderName '_BindingMatrices' designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod','dCp_mod','-v7.3')
    save([saveRoot filesep settings.FolderName filesep settings.FolderName '_BindingMatrices2' designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','dCp_Complement','-v7.3')
    save([saveRoot filesep settings.FolderName filesep settings.FolderName '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' designerName '.mat'],'Kb_Complement','-v7.3')

    save([saveRoot filesep settings.FolderName filesep settings.FolderName '_OnOffThermoInfo' designerName '.mat'],'Kon','Koff',-'v7.3');
    save([saveRoot filesep settings.FolderName filesep settings.FolderName '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' designerName '.mat'],'Kb_mod','-v7.3')
end
%% Get RNA Secondary Structure
%Update to load from database if precalculated
%Check if function fails and if can remove catch on getting Binding energy
%for RNA secondary structure
if (settings.SingleOrMulti==1&&settings.AllIsoforms == 0)%One Gene/One Isoform
    try
        load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_Tm' num2str(T_hybrid) '_RNASecondaryStructure' designerName '.mat'],'Ks_TjLi','IsSiteInLoop')
    catch
        try
            [IsSiteInLoop,Ks_TjLi,dHs_TjLi,dSs_TjLi,dGs_TjLi_eq,dHs_TjLi_f,dSs_TjLi_f,dGs_TjLi_f,dHs_TjLi_r,dSs_TjLi_r,dGs_TjLi_r,dCp_TjLi_eq] = ...
                A_JH_GetRNASecondaryStructures_V2(settings,gene_table,DoesProbeBindSite2,probes)
            save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_Tm' num2str(T_hybrid) '_RNASecondaryStructure' designerName '.mat'],'Ks_TjLi','IsSiteInLoop','-v7.3')
        catch
            Ks_TjLi = [];
            IsSiteInLoop = [];
        end
    end
else
    try
        load([saveRoot filesep settings.FolderName filesep settings.FolderName '_Tm' num2str(T_hybrid) '_RNASecondaryStructure' designerName '.mat'],'Ks_TjLi','IsSiteInLoop')
    catch
        try
            [IsSiteInLoop,Ks_TjLi,dHs_TjLi,dSs_TjLi,dGs_TjLi_eq,dHs_TjLi_f,dSs_TjLi_f,dGs_TjLi_f,dHs_TjLi_r,dSs_TjLi_r,dGs_TjLi_r,dCp_TjLi_eq] = ...
                A_JH_GetRNASecondaryStructures_V2(settings,gene_table,DoesProbeBindSite2,probes)
            save([saveRoot filesep settings.FolderName filesep settings.FolderName '_Tm' num2str(T_hybrid) '_RNASecondaryStructure' designerName '.mat'],'Ks_TjLi','IsSiteInLoop','-v7.3')
        catch
            Ks_TjLi = [];
            IsSiteInLoop = [];
        end
    end
end



ExprLevels_Mean = mean(ExpressionMatrix,2);
ExprLevels_Null = ones(size(ExpressionMatrix,1),1);
EKernel = ExprLevels_Null;
if (settings.SingleOrMulti==1&&settings.AllIsoforms == 0)%One Gene/One Isoform
    FoldName = [];
end
Kon = squeeze(Kon(:,Nmodel));
Koff = squeeze(Koff(:,:,Nmodel));
Kb_mod = squeeze(Kb_mod(:,:,:,Nmodel));
DoesProbeBindSite = DoesProbeBindSite2;
C_var{1} = dCp_mod;
C_var{2} = dHeq_mod;
C_var{3} = dSeq_mod;
C_var{4} = dHf_mod;
C_var{5} = dSf_mod;
C_var{6} = dHr_mod;
C_var{7} = dSr_mod;
if (settings.BLASTdna)
    C_var{8} = dHeq_Complement;
    C_var{9} = dCp_Complement;
    C_var{10} = dSeq_Complement;
    C_var{11} = dHf_Complement;
    C_var{12} = dSf_Complement;
    C_var{13} = dHr_Complement;
    C_var{14} = dSr_Complement;
    Kb_Complement = squeeze(Kb_Complement(:,:,Nmodel));
else
    dHeq_Complement = [];
    dCp_Complement = [];
    dSeq_Complement = [];
    dHf_Complement = [];
    dSf_Complement = [];
    dHr_Complement = [];
    dSr_Complement = [];
    Kb_Complement = [];
end

if (settings.SingleOrMulti==1&&settings.AllIsoforms == 0)%One Gene/One Isoform
    if (id2~=2)
        try
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_Tm' num2str(T_hybrid) '_BasicDesignerStats' designerName '.mat'],'Tvec_RNA','Svec_RNA','TPvec_RNA','TSvec_RNA','TPvec_logKOFF_RNA','TPvec_logKOFFdivON_RNA','TPvec_logKONdivOFF_RNA',...
                'Nvec_RNAmulti','Tvec_DNA','Svec_DNA','TPvec_DNA','TSvec_DNA','TPvec_logKOFF_DNA','TPvec_logKOFFdivON_DNA','TPvec_logKONdivOFF_DNA','TPvec_logKOFFdivCOMP_DNA','TPvec_logKCOMPdivOFF_DNA',...
                'Off_Score','Specificity_Score','NumRNAOffTargetOptions','Probes_WithNRNAOFF','NumDNAOffTargetOptions','Probes_WithNDNAOFF')
        catch
            try
                [Nvec_RNAmulti,RNAOFF_Score,RNASpecificity_Score,NumRNAOffTargetOptions,Probes_WithNRNAOFF,DNAOFF_Score,DNASpecificity_Score,NumDNAOffTargetOptions,Probes_WithNDNAOFF,Cout] = ...
                    A0_BasicDesignerStats(targetTypes,removeUndesiredIsos,gene_table,settings,FoldName,DoesProbeBindSite,Kon,Kb_mod,Kb_Complement,EKernel);
                Tvec_RNA = Cout{1}{1};Svec_RNA = Cout{1}{2};TPvec_RNA = Cout{1}{3};TSvec_RNA = Cout{1}{4};
                TPvec_logKOFF_RNA = Cout{1}{5};TPvec_logKOFFdivON_RNA = Cout{1}{6};TPvec_logKONdivOFF_RNA = Cout{1}{7};
                Tvec_DNA = Cout{2}{1};Svec_DNA = Cout{2}{2};TPvec_DNA = Cout{2}{3};TSvec_DNA = Cout{2}{4};
                TPvec_logKOFF_DNA = Cout{2}{5};TPvec_logKOFFdivON_DNA = Cout{2}{6};TPvec_logKONdivOFF_DNA = Cout{2}{7};
                TPvec_logKOFFdivCOMP_DNA = Cout{2}{8};TPvec_logKCOMPdivOFF_DNA = Cout{2}{9};
                Off_Score = RNAOFF_Score+DNAOFF_Score;
                Specificity_Score = RNASpecificity_Score+DNASpecificity_Score;
                save([saveRoot filesep settings.FolderName filesep settings.FolderName '_Tm' num2str(T_hybrid) '_BasicDesignerStats' designerName '.mat'],...
                    'Tvec_RNA','Svec_RNA','TPvec_RNA','TSvec_RNA','TPvec_logKOFF_RNA','TPvec_logKOFFdivON_RNA','TPvec_logKONdivOFF_RNA',...
                    'Tvec_DNA','Svec_DNA','TPvec_DNA','TSvec_DNA','TPvec_logKOFF_DNA','TPvec_logKOFFdivON_DNA','TPvec_logKONdivOFF_DNA','TPvec_logKOFFdivCOMP_DNA','TPvec_logKCOMPdivOFF_DNA',...
                    'Nvec_RNAmulti','Off_Score','Specificity_Score','NumRNAOffTargetOptions','Probes_WithNRNAOFF','NumDNAOffTargetOptions','Probes_WithNDNAOFF','-v7.3')
            catch
            end
        end
    else
        try
            load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_Tm' num2str(T_hybrid) '_BasicDesignerStats' designerName0 '.mat'],'Tvec_RNA','Svec_RNA','TPvec_RNA','TSvec_RNA','TPvec_logKOFF_RNA','TPvec_logKOFFdivON_RNA','TPvec_logKONdivOFF_RNA',...
                'Nvec_RNAmulti','Tvec_DNA','Svec_DNA','TPvec_DNA','TSvec_DNA','TPvec_logKOFF_DNA','TPvec_logKOFFdivON_DNA','TPvec_logKONdivOFF_DNA','TPvec_logKOFFdivCOMP_DNA','TPvec_logKCOMPdivOFF_DNA',...
                'Off_Score','Specificity_Score','NumRNAOffTargetOptions','Probes_WithNRNAOFF','NumDNAOffTargetOptions','Probes_WithNDNAOFF')
        catch
            try
                [Nvec_RNAmulti,RNAOFF_Score,RNASpecificity_Score,NumRNAOffTargetOptions,Probes_WithNRNAOFF,DNAOFF_Score,DNASpecificity_Score,NumDNAOffTargetOptions,Probes_WithNDNAOFF,Cout] = ...
                    A0_BasicDesignerStats(targetTypes,removeUndesiredIsos,gene_table,settings,FoldName,DoesProbeBindSite,Kon,Kb_mod,Kb_Complement,EKernel);
                Tvec_RNA = Cout{1}{1};Svec_RNA = Cout{1}{2};TPvec_RNA = Cout{1}{3};TSvec_RNA = Cout{1}{4};
                TPvec_logKOFF_RNA = Cout{1}{5};TPvec_logKOFFdivON_RNA = Cout{1}{6};TPvec_logKONdivOFF_RNA = Cout{1}{7};
                Tvec_DNA = Cout{2}{1};Svec_DNA = Cout{2}{2};TPvec_DNA = Cout{2}{3};TSvec_DNA = Cout{2}{4};
                TPvec_logKOFF_DNA = Cout{2}{5};TPvec_logKOFFdivON_DNA = Cout{2}{6};TPvec_logKONdivOFF_DNA = Cout{2}{7};
                TPvec_logKOFFdivCOMP_DNA = Cout{2}{8};TPvec_logKCOMPdivOFF_DNA = Cout{2}{9};
                Off_Score = RNAOFF_Score+DNAOFF_Score;
                Specificity_Score = RNASpecificity_Score+DNASpecificity_Score;
                save([saveRoot filesep settings.FolderName filesep settings.FolderName '_Tm' num2str(T_hybrid) '_BasicDesignerStats' designerName0 '.mat'],...
                    'Tvec_RNA','Svec_RNA','TPvec_RNA','TSvec_RNA','TPvec_logKOFF_RNA','TPvec_logKOFFdivON_RNA','TPvec_logKONdivOFF_RNA',...
                    'Tvec_DNA','Svec_DNA','TPvec_DNA','TSvec_DNA','TPvec_logKOFF_DNA','TPvec_logKOFFdivON_DNA','TPvec_logKONdivOFF_DNA','TPvec_logKOFFdivCOMP_DNA','TPvec_logKCOMPdivOFF_DNA',...
                    'Nvec_RNAmulti','Off_Score','Specificity_Score','NumRNAOffTargetOptions','Probes_WithNRNAOFF','NumDNAOffTargetOptions','Probes_WithNDNAOFF','-v7.3')
            catch
            end
        end
    end
else
end

if (id2~=2)
    %% Get Metric Information (Probes, Final Probe Set)
    try
        load([settings.FolderRootName filesep inputs1{gene_num,5} '_Tm' num2str(T_hybrid) '_ModelMetrics' designerName '.mat'],'ModelMetrics')
    catch%add changes to make fully work (add isoformflattened)
        chosenProbes = 1:size(probes,1);
        if (~isempty(chosenProbes))
            ModelMetrics = ...
                RNAsolver_JH(1:size(probes,1),settings,probes,gene_table,ExpressionMatrix,DoesProbeBindSite2,dHeq_mod,dSeq_mod,dCp_mod,dHeq_Complement,dSeq_Complement,dCp_Complement)
            save([settings.FolderRootName filesep inputs1{gene_num,5} '_Tm' num2str(T_hybrid) '_ModelMetrics' designerName '.mat'],'ModelMetrics','chosenProbes','settings','-v7.3')
        end
    end
else
    for v = 5:-1:0
        load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_probes' designerName0 '.mat'],'probes')
        designerName = strcat('_Stellaris_L',num2str(v));
        chosenProbes = StellarisProbesLi{v+1};
        %% Get Metric Information (Probes, Final Probe Set)
        if (~isempty(chosenProbes))
            try
                load([settings.FolderRootName filesep inputs1{gene_num,5} '_Tm' num2str(T_hybrid) '_ModelMetrics' designerName '.mat'],'ModelMetrics')
            catch%add changes to make fully work (add isoformflattened)                                                                             Num_of_Molecule_Sites
                ModelMetrics = ...
                    RNAsolver_JH2(chosenProbes,settings,probes,gene_table,ExpressionMatrix,DoesProbeBindSite2,dHeq_mod,dSeq_mod,dCp_mod,dHeq_Complement,dSeq_Complement,dCp_Complement)
                save([settings.FolderRootName filesep inputs1{gene_num,5} '_Tm' num2str(T_hybrid) '_ModelMetrics' designerName '.mat'],'ModelMetrics','chosenProbes','settings','-v7.3')
            end
        end
    end
end
%% Make/Save Comparison Figures & Figure Objects
Fin = 1
end
