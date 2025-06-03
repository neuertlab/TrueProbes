function A0_BKJH_ProbeDesign_Wrapper_cluster_V5(id,cluster)
% Main FIle for designing TrueProbes RNA-FISH probes
% Input Argument 1: id
% id is an integer and is the row of the input design table at the top of the script to run and design probes against. 
% 
% Input Argument 2: cluster
% cluster is an integer and determines the software parallelization pool between local or remote servers when running the script via Slurm. 
% The only difference between running on a cluster is the number of cores in the Slurm file.
% cluster = (0,if run locally and 1 if run on cluster)
% 
% Inputs1 Table
% Table describing all the targets to design probes against. Each row is a different desired target design. 
% Input table columns:
% 	1. Included target Accession IDs. Designs probes shared across all target accession number(s)
% 	2. Text Sequence Files to Include (files). Default is Empty
% 	3. Text Sequence Files to Exclude (files). Default is Empty
% 	4. Organism you are trying to design probes in.
% 	5. Gene Name 1. First name of gene in parenthesis
%     6. Gene Name 2. Second potential gene name to use instead of the first  in parenthesis
% 	7. Chromosome. The chromosome number or letter symbol
% 	8. Excluded target accession IDs. Removes probes in any exclusion accession number(s). Default is Empty
%    9. Strand. Which strand to design probes against for RNA default is  plus .

%% Genes To Design Probes For
inputs1 = {...
    {'NM_001660.4'},{},{}, 'Human','(ARF4)','(ARF4)','3',{},1 ;... %
    {'NM_181604.2'},{},{}, 'Human','(KRTAP6-2)','(KRTAP6-2)','21',{},1 ;...           %ENST00000277480.7 ,  820bp, 2-Isos
    {'NM_181602.2'},{},{}, 'Human','(KRTAP6-1)','(KRTAP6-1)','21',{},1 ;...           %ENST00000277480.7 ,  820bp, 2-Isos
    {'NR_002844.2'},{},{}, 'Mouse','(TSIX)','(TSIX)','X',{},1 ;... %4306
    {'NM_001179952.1'},{},{}, 'Yeast','(HSP12)','(HSP12)','VI',{},1 ;...           %ENST00000368087.8 , 1447bp, 4-Isos
    {'ENST00000818900.1'},{},{}, 'Human','(MNX1-AS1)','(MNX1-AS1)','7',{},1 ;...     %N/A               , 1559bp, 3-Isos closest match ENST00000632538.1
    {'ENSMUST00000181020.10'},{},{}, 'Mouse','(JPX)','(JPX)','X',{},1 ;...    %N/A               , 1571bp, 2-Isos
    {''},{},{}, 'Yeast','()','()','12',{},1 ;...        %ENST00000312492.3 , 1611bp, 2-Isos
    {''},{},{}, 'Human','()','()','7',{},1 ;...        %N/A               , 1619bp, 2-Isos
    {''},{},{}, 'Mouse','()','()','3',{},1 ;...     %ENST00000237696.10, 1713bp, 3-Isos
    {''},{},{}, 'Yeast','()','()','13',{},1 ;...     %N/A               , 1831bp, 4-Isos
    {'ENST00000260386.7'},{},{}, 'Human','(ITPKA)','(ITPKA)','15',{},1 ;...        %ENST00000260386.7 , 1825bp, 2-Isos
    {'NM_001357731.1'},{},{}, 'Human','(EIF2S3B)','(EIF2S3B)','12',{},1 ;... %N/A               , 1866bp, 2-Isos
    {'NM_001375593.1'},{},{}, 'Human','(PRR7)','(PRR7)','5',{},1 ;...        %N/A               , 1897bp, 5-Isos
    {'NM_013357.2'},{},{}, 'Human','(PURG)','(PURG)','8',{},1 ;...           %N/A               , 1972bp, 4-Isos
    {'NM_001099773.2'},{},{}, 'Human','(CYP11A1)','(CYP11A1)','15',{},1 ;... %N/A               , 1986bp, 2-Isos
    {'NM_019065.3'},{},{}, 'Human','(NECAB2)','(NECAB2)','16',{},1 ;...      %ENST00000305202.9 , 1996bp, 4-Isos
    {'ENST00000305202.9'},{},{}, 'Human','(NECAB2)','(NECAB2)','16',{},1 ;...
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
baseFolder = pwd;
addpath(genpath(baseFolder));
%% Main Probe Design Settings (Usually migt change)
max_probes = 96; % Max number of probes to design
minProbeSize = 20;% Min Probe Size
maxProbeSize = 20;% Max Probe Size
MininumProbeSpacing = 3; %Minimum Spacing Between Adjacent Probes
BLASTrna = 1; %Decides if you will BLAST the reference transcriptome (RNA targets)
BLASTdna = 0; %Decides if you will BLAST the reference genome (DNA targets)
customBlastDatabase_DNA = 'N/A'; % Location of custom DNA database
customBlastDatabase_RNA = 'N/A'; % Location of custom RNA database
HybridizationTemperature = 37; % Temperature for Evaluating Probe Design and Simulations
saveRoot = strcat('output',filesep);
designerName = '_TrueProbes';

%% Secondary Parameters (You Usually will not change)
Nmodel = 4; %Which Hybridization model to use for probe design and evaluation
RunOffline = 1;%Is Design Run Offline using databases or looks online to get genbank record
ParsingPreference = 1; % Sequentially Parse or Parallel parsing of blast hits together into hit table (1 sequentially, 0 parallel)
RemoveRibosomalHits = 1;% Filter out probes with targets hits to ribosomal proteins
SaltConcentration = 0.05; %Concentration of Salt in thermodynamic calculations mol/L
MinHomologySearchTargetSize = 15; % minimum off-target match size
probeBatchSize = 20;%batch size for parallelizing probe evaluations in probe design
targetBatchSize = 200;%batch size for parallelizing target evaluations in probe design
BLASTpath_Windows = strcat('src',filesep,'thirdparty',filesep,'ncbi-blast-2.8.1+-x64-win64',filesep,'blast-2.8.1+',filesep,'bin',filesep,'blastn'); %version of blast database used and path to blast
BLASTpath_Mac = strcat('src',filesep,'thirdparty',filesep,'ncbi-blast-2.8.1+-x64-macosx',filesep,'blast-2.8.1+',filesep,'bin',filesep,'blastn'); %version of blast database used and path to blast
BLASTpath_Linux = strcat('src',filesep,'thirdparty',filesep,'ncbi-blast-2.8.1+-x64-linux',filesep,'blast-2.8.1+',filesep,'bin',filesep,'blastn'); %version of blast database used and path to blast
%% BLAST Parameters (You Usually will not change)
outfmt = 5;
reward = 1;
penalty = -3;
word_size = 7;
gapopen = 5;
gapextend = 2;
evalue = 1000;
num_alignments = 1000;
dust = 'no';






SingleOverMultiplex = 1;
AllIsoforms = 0;

gene_num = id;
addSelfProb = 1;
packOptimal = 1;
targetTypes = [1 0];
removeUndesiredIsos=1;
UseGeneOverTranscLevelExpression = 0;
nullRNAcopynumber = 100;
nullDNAcopynumber = 2;
RemoveMisMatches = 1;
SpecificityThreshold = 2;
DecisionAlgorithmFloorSize = 0.5;
withNascentTranscripts = 0;
if (minProbeSize>maxProbeSize)
    msg = 'Error. Minimum probe size must be less than or equal to max probe size';
    error(msg)
end
if (MininumProbeSpacing<0)
    msg = 'Error. Minimum probe spacing must be greater than or equal to zero';
    error(msg)
end
if (max_probes<0)
    msg = 'Error. Maximum number of probes must be greater than zero';
    error(msg)
end


%Update so that more species can be included in list for their blast database
% Add Yeast Expression, update to use GTF
% Update that way more species can be used and alignment of database files
% and GTF with reference genome, and add custom gene expression file


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
%code that stores location info in 3d                ID in 3d-4d maps
% code that combines maps
% code that designs probes [mix designer with V7_Out codes]
% solves eq/kinetics/metrics with multi-transcript output
% all forked with/wit                                    hout temp changes or gradients.
%update to have higher dimensions for versions of code/ cell array
%might turn 4D multiTargetLocations into two 3d matricies and change
%combineMaps code to for those
probes = [];%  SNP allele specificity  (similar to yeast)
startup
%Banerjee RNA/DNA Improved nearest-neighbor parameters for the stability of RNA/DNA hybrids under a physiological condition


%find ncbi gene name alias's  (since its on ncbi page)


%37-65
settings.BUILD_STRING = '2024.11.12.00';
settings.VERSION_STRING = 'v1.1.1';


%read(seqBioIFobj,inputs1{x}{1}).Sequence

%% Settings Specification




if (strcmp(inputs1{gene_num,4},'Yeast'))
    DoAllGenesHaveSameExpression = 1;%same expr.
else
    DoAllGenesHaveSameExpression = 0;%diff expr
end


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






%Save Settings
settings.FolderRootName = strcat(saveRoot,inputs1{gene_num,5},'_',strjoin(inputs1{gene_num,1},'_'));
settings.rootName = strjoin(inputs1{gene_num,1},'_');
settings.designerName = designerName;
settings.saveRoot = saveRoot;
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
settings.BLASTbatchSize = probeBatchSize;
settings.BLASTsimultaneousParsingOverSequentialParsing = ParsingPreference;
settings.BLASTpath_Windows = BLASTpath_Windows;
settings.BLASTpath_Mac = BLASTpath_Mac;
settings.BLASTpath_Linux = BLASTpath_Linux;
settings.BlastParameters.outfmt = outfmt;
settings.BlastParameters.reward = reward;
settings.BlastParameters.penalty = penalty;
settings.BlastParameters.wordsize = word_size;
settings.BlastParameters.gapopen = gapopen;
settings.BlastParameters.gapextend = gapextend;
settings.BlastParameters.evalue = evalue;
settings.BlastParameters.num_alignments = num_alignments;
settings.BlastParameters.dust = dust;



settings.TargetBatchSize = targetBatchSize;
settings.SingleOrMulti = SingleOverMultiplex;
settings.AllIsoforms = AllIsoforms;
settings.MinProbeSize = minProbeSize;
settings.MaxProbeSize = maxProbeSize;
settings.MinHomologySearchTargetSize = MinHomologySearchTargetSize;
settings.N_model = Nmodel;

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

%ENSEMBL
%dna chromosome
%RNA ENST... chromosome, gene:ENSG...

%% Update Location of Databases & Needed Files
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


%% GTF and GFF Databases
if (strcmp(settings.referenceType,'RefSeq'))
    settings.hLocGTF = 'data/DatabaseData/GTF_Databases/Human/NCBI_RefSeq/GCF_000001405.40_GRCh38.p14_genomic.gtf';
    settings.hLocGFF = 'data/DatabaseData/GFF3_Databases/Human/NCBI_RefSeq/GCF_000001405.40_GRCh38.p14_genomic.gff';
    settings.mLocGTF = 'data/DatabaseData/GTF_Databases/Mouse/NCBI_RefSeq/GCF_000001635.27_GRCm39_genomic.gtf';
    settings.mLocGFF = 'data/DatabaseData/GFF3_Databases/Mouse/NCBI_RefSeq/GCF_000001635.27_GRCm39_genomic.gff';
    settings.yLocGTF =  'data/DatabaseData/GTF_Databases/Yeast/NCBI_RefSeq/GCF_000146045.2_R64_genomic.gtf';
    settings.yLocGFF =  'data/DatabaseData/GFF3_Databases/Yeast/NCBI_RefSeq/GCF_000146045.2_R64_genomic.gff';
elseif (strcmp(settings.referenceType,'ENSEMBL'))
    settings.hLocGTF = 'data/DatabaseData/GTF_Databases/Human/EMBL_EBI/Homo_sapiens.GRCh38.114.gtf';%chrom number, MT, X, Y and scaffold KI and GL
    settings.hLocGFF = 'data/DatabaseData/GFF3_Databases/Human/EMBL_EBI/Homo_sapiens.GRCh38.114.gff3';%cdna and ncrna
    settings.mLocGTF = 'data/DatabaseData/GTF_Databases/Mouse/EMBL_EBI/Mus_musculus.GRCm39.114.gtf';
    settings.mLocGFF = 'data/DatabaseData/GFF3_Databases/Mouse/EMBL_EBI/Mus_musculus.GRCm39.114.gff3';%ENMUST, cdna and ncRNA
    settings.yLocGTF =  'data/DatabaseData/GTF_Databases/Yeast/EMBL_EBI/Saccharomyces_cerevisiae.R64-1-1.1.114.gtf';%chromsome roman numerals and Mito
    settings.yLocGFF =  'data/DatabaseData/GFF3_Databases/Yeast/EMBL_EBI/Saccharomyces_cerevisiae.R64-1-1.1.114.gff3';%SGD cdna and ncRNA
elseif (strcmp(settings.referenceType,'Gencode'))
    settings.hLocGTF = 'data/DatabaseData/GTF_Databases/Human/GENCODE/gencode.v48.basic.annotation.gtf';%chrM chrX chrY
    settings.hLocGFF = 'data/DatabaseData/GFF3_Databases/Human/GENCODE/gencode.v48.basic.annotation.gff3';%ENST and ENSG
    settings.mLocGTF = 'data/DatabaseData/GTF_Databases/Mouse/GENCODE/gencode.vM37.basic.annotation.gtf';%chr1,... chrX, chrY, chrM
    settings.mLocGFF = 'data/DatabaseData/GFF3_Databases/Mouse/GENCODE/gencode.vM37.basic.annotation.gff3';
    settings.yLocGTF =  'data/DatabaseData/GTF_Databases/Yeast/YGD/saccharomyces_cerevisiae.20240529.gtf';
    settings.yLocGFF =  'data/DatabaseData/GFF3_Databases/Yeast/saccharomyces_cerevisiae.20240529.gff';
end



cellPreset = 1;
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

settings.MouseExpressionFile = 'data/DatabaseData/tabulamuris_barChart.bed';
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
settings.scTracks = {'colonWangCellType','tabulasapiens_tissue_cell_type'};


%installedToolbox = matlab.addons.toolbox.installToolbox('BLASTPlus.Support.Package.for.Bioinformatics.Toolbox.mltbx')

%Save Settings
settings.FolderRootName = strcat(saveRoot,inputs1{gene_num,5},'_',strjoin(inputs1{gene_num,1},'_'));
settings.rootName = strjoin(inputs1{gene_num,1},'_');
settings.designerName = designerName;

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
if (size(inputs1{id,1},2)>1)
    settings.SingleOrMulti = 1;
end



%% Load Annotation File
fprintf('\n')
fprintf("Loading genome and transcriptome annotation files")
tic
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
    % optsUCSC2 = detectImportOptions(settings.Yeast_wgEncodeGencodeAttributesFile);
    % optsUCSC2.VariableNames = settings.wgEncodeGencodeAttributesVariableNames;
    % EMBLAttrAlign_db = readtable(settings.Yeast_wgEncodeGencodeAttributesFile,optsUCSC2);
    % optsUCSC3 = detectImportOptions(settings.Yeast_wgEncodeGencodeCompFile);
    % optsUCSC3.VariableNames = settings.wgEncodeGencodeCompVariableNames;
    % EMBLComp_db = readtable(settings.Yeast_wgEncodeGencodeCompFile,optsUCSC3);
else  %custom organism
    optsUCSC2 = detectImportOptions(settings.Custom_wgEncodeGencodeAttributesFile);
    optsUCSC2.VariableNames = settings.wgEncodeGencodeAttributesVariableNames;
    EMBLAttrAlign_db = readtable(settings.Custom_wgEncodeGencodeAttributesFile,optsUCSC2);
    optsUCSC3 = detectImportOptions(settings.Custom_wgEncodeGencodeCompFile);
    optsUCSC3.VariableNames = settings.wgEncodeGencodeCompVariableNames;
    EMBLComp_db = readtable(settings.Custom_wgEncodeGencodeCompFile,optsUCSC3);
end
tEnd = toc;
fprintf('\n')
fprinf("Time elapsed to load annotation files %g seconds",round(tEnd,3,"significant"))

%% Generate Folder
if (not(isfolder([saveRoot])))
    mkdir([saveRoot])
end
fprintf('\n')
fprintf("Making probe target design output folder")
tic
if (settings.SingleOrMulti==1&&settings.AllIsoforms == 0)%One Gene/One Isoform
    if (not(isfolder([saveRoot inputs1{gene_num,5} '_' settings.rootName])))
        mkdir([saveRoot inputs1{gene_num,5} '_' settings.rootName])
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
        if (not(isfolder([saveRoot filesep inputs1{gene_num,5} '_' ENST_RefSeqIDz{k}])))
            mkdir([saveRoot filesep inputs1{gene_num,5} '_' ENST_RefSeqIDz{k}])
        end
        FoldName{k} = [inputs1{gene_num,5} '_' ENST_RefSeqIDz{k}];
    end
    settings.rootName = strjoin(ENST_RefSeqIDz,'_');
    if (not(isfolder([saveRoot filesep inputs1{gene_num,5} '_' settings.rootName])))
        mkdir([saveRoot filesep inputs1{gene_num,5} '_' settings.rootName])
    end
    All_chrom = Isoform_chrom;
    settings.FolderName = [inputs1{gene_num,5} '_' settings.rootName];
end
settings.FolderRootName = strcat('output',filesep,settings.FolderName);
settings.rootName = strjoin(inputs1{gene_num,1},'_');
settings.TargetLists = strjoin(inputs1{gene_num,1},', ');


%% Generate Probes
fprintf('\n')
fprintf(strcat("Designing probes for"," ",inputs1{gene_num,5},"Transcript IDs:"," ",settings.TargetLists))
fprintf('\n')
try
    load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_probes' designerName '.mat'],'probes')
    fprintf("Loading probe tile sequences")
catch
    fprintf("Generating probe tile sequences")
    try
        tic
        [init_probes,~,~,~,~] = ...
            BKJH_Probe_Generator([Lmin Lmax],inputs1{gene_num,1},inputs1{gene_num,2},inputs1{gene_num,3},inputs1{gene_num,8},inputs1{gene_num,9},settings.isOffline,settings.SEQdbRoot);
        probes = init_probes;
        save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_probes' designerName '.mat'],'probes','-v7.3')
        tEnd = toc;
        fprintf('\n')
        fprinf("Time elapsed to tile probes %g seconds",round(tEnd,3,"significant"))
    catch e
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end
%% BLAST Probes
fprintf('\n')
try
    load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_hits_table' designerName '.mat'],'gene_table')
    fprintf("Loading probe BLAST results table")
catch
    fprintf("BLASTING Probes in batches")
    try
        tic                                                                                                           %Organism           %(Gene Name)        %(Gene Name)        %ChrNum
        [~,gene_table] = Probe_checker_general_JH10(probes,inputs1{gene_num,4},inputs1{gene_num,5},inputs1{gene_num,6},inputs1{gene_num,7},settings);
        save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_hits_table' designerName '.mat'],'-mat','gene_table','-v7.3');
        tEnd = toc;fprintf('\n')
        fprinf("Time elapsed to generate probe BLAST results table %g seconds",round(tEnd,3,"significant"))
    catch e
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end
%% Get Gene Expression Information
fprintf('\n')
try
    load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_ExpressionInfo' designerName '.mat'],'ExpressionMatrix');
    fprintf("Loading BLAST hits gene expression information")
catch
    fprintf("Getting BLAST hits gene expression information")
    try
        tic
        [ExpressionMatrix,get_expression_time] = A_JH_GetExpressionInfo_V2(gene_table,settings);
        save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_ExpressionInfo' designerName '.mat'],'ExpressionMatrix','get_expression_time','settings','-v7.3');
        tEnd = toc;fprintf('\n')
        fprinf("Time elapsed to generate probe BLAST hits gene expression information %g seconds",round(tEnd,3,"significant"))
    catch e
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end
settings.saveRoot = saveRoot;

%% Get Thermodynamic Information (On-Target,Off-Target)
fprintf('\n')
try
    load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_Tm' num2str(settings.HybridizationTemperature) '_OnOffThermoInfo' designerName '.mat'],'Kon','Koff','Kb_Match');
    load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dHInfo' designerName '.mat'],'-mat','dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match');
    load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dSInfo' designerName '.mat'],'-mat','dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match');
    load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_TmInfo' designerName '.mat'],'-mat','Tm_on','Tm_Match');
    load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dCpInfo' designerName '.mat'],'-mat','dCpon_eq','dCpeq_Match');
    fprintf("Loading probe target thermodynamic information")
catch
    fprintf("Computing probe target thermodynamic information")
    tic
    try
        [Kb_Match,Kon,Koff,dHeq_Match,dSeq_Match,dHf_Match,dSf_Match,dHr_Match,dSr_Match,dCpeq_Match,dHon_eq,dSon_eq,dHon_f,dSon_f,dHon_r,dSon_r,dCpon_eq,Tm_on,Tm_Match] = ...
            A_JH_GenerateThermoInfo_V5(probes,gene_table,inputs1{gene_num,5},settings);%add Kon Koff
        save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_Tm' num2str(settings.HybridizationTemperature) '_OnOffThermoInfo' designerName '.mat'],'-mat','Kon','Koff','Kb_Match','-v7.3');
        save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dHInfo' designerName '.mat'],'-mat','dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match','-v7.3');
        save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dSInfo' designerName '.mat'],'-mat','dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match','-v7.3');
        save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_TmInfo' designerName '.mat'],'-mat','Tm_on','Tm_Match','-v7.3');
        save([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_dCpInfo' designerName '.mat'],'-mat','dCpon_eq','dCpeq_Match','-v7.3');
        toc
        tEnd = toc;fprintf('\n')
        fprinf("Time elapsed to compute probe target thermodynamic information %g seconds",round(tEnd,3,"significant"))
    catch e
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end
%% Get Binding Site Mapping and Energy
fprintf('\n')
try
    load([settings.FolderRootName filesep settings.GeneName '_binding_hits_map' designerName '.mat'],'DoesProbeBindSite2','Num_of_Molecule_Sites')
    load([settings.FolderRootName filesep settings.GeneName  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' designerName '.mat'],'Kb_Complement')
    load([settings.FolderRootName filesep settings.GeneName  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' designerName '.mat'],'Kb_mod')
    load([settings.FolderRootName filesep settings.GeneName '_BindingMatrices' designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod','dCp_mod')
    if (settings.BLASTdna)
        load([settings.FolderRootName filesep settings.GeneName '_BindingMatrices2' designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','dCp_Complement')
    end
    fprintf("Loading probe target binding site maps")
catch
    fprintf("Generating probe target binding site maps")
    try
        tic
        [Kb_mod,Kb_Complement,DoesProbeBindSite2,~,~,~,...
            dHeq_mod,dSeq_mod,dHf_mod,dSf_mod,dHr_mod,dSr_mod,~,dCp_mod,...
            dHeq_Complement,dSeq_Complement,dHf_Complement,dSf_Complement,dHr_Complement,dSr_Complement,dCp_Complement] = ...
            A_JH_GetSiteMapping_V6(probes,settings,gene_table,Kb_Match,dHeq_Match,dSeq_Match,dHf_Match,dSf_Match,dHr_Match,dSr_Match,dCpeq_Match,Tm_Match);
        tEnd = toc;fprintf('\n')
        fprinf("Time elapsed to create probe target binding site maps %g seconds",round(tEnd,3,"significant"))
    catch e
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end

%% Basic Stats for Designing Probes Function
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

fprintf('\n')
try
    load([settings.FolderRootName filesep inputs1{gene_num,5} '_' settings.rootName '_Tm' num2str(T_hybrid) '_BasicDesignerStats' designerName '.mat'],'Tvec_RNA','Svec_RNA','TPvec_RNA','TSvec_RNA','TPvec_logKOFF_RNA','TPvec_logKOFFdivON_RNA','TPvec_logKONdivOFF_RNA',...
        'Nvec_RNAmulti','Tvec_DNA','Svec_DNA','TPvec_DNA','TSvec_DNA','TPvec_logKOFF_DNA','TPvec_logKOFFdivON_DNA','TPvec_logKONdivOFF_DNA','TPvec_logKOFFdivCOMP_DNA','TPvec_logKCOMPdivOFF_DNA',...
        'Off_Score','Specificity_Score','NumRNAOffTargetOptions','Probes_WithNRNAOFF','NumDNAOffTargetOptions','Probes_WithNDNAOFF')
    fprintf("Loading probe target statistics information")
catch
    try
        fprintf("Computing probe target statistics information")
        tic
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
        tEnd = toc;fprintf('\n')
        fprinf("Time elapsed to compute probe target statistics %g seconds",round(tEnd,3,"significant"))
    catch e
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end

%% Selection of TrueProbes Probes
fprintf('\n')
try
    load([saveRoot filesep settings.FolderName filesep settings.FolderName '_chosen.mat'],'chosenProbes')
    fprintf("Loading TrueProbes designed probes")
catch
    try
        fprintf("Designing TrueProbes probe set")
        tic
        chosenProbes = A_ZigZagProbeSelection_V5(probes,gene_table,settings,addSelfProb,packOptimal,Kon,Nvec_RNAmulti,Off_Score,Specificity_Score,Tvec_RNA,Svec_RNA,TPvec_RNA,TSvec_RNA,TPvec_logKOFF_RNA,TPvec_logKOFFdivON_RNA,TPvec_logKONdivOFF_RNA,ExpressionMatrix,DoesProbeBindSite,Kb_mod);
        save([saveRoot filesep settings.FolderName filesep settings.FolderName '_chosen.mat'],'chosenProbes','-v7.3')
        tEnd = toc;fprintf('\n')
        fprinf("Time elapsed to select TrueProbes probes %g seconds",round(tEnd,3,"significant"))
    catch e
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end

%% Print Excel Spreedsheet of Probes
fprintf('\n')
fprintf("Printing TrueProbes probes to Excel spreadsheet")
Lpmin = min(cell2mat(cellfun(@length,{probes{:,2}},'UniformOutput',false)));
Lpmax = max(cell2mat(cellfun(@length,{probes{:,2}},'UniformOutput',false)));
TL_hits = zeros(length(min(15,Lpmin):Lpmax),size(probes,1));
for L = min(15,Lpmin):Lpmax
    TL_hits(L-14,:) = cell2mat(arrayfun(@(z) length(find((gene_table.ProbeNum==z).*(gene_table.Match == L)>0)),1:size(probes,1),'Un',0));
end
for v = 1:length(chosenProbes)
    probe_num = chosenProbes(v);
    Lp = length(probes{chosenProbes(v),2});
    final_probe_info{v,1} = settings.transcript_IDs;
    final_probe_info{v,2} = probes{chosenProbes(v),3};
    final_probe_info{v,3} = probes{chosenProbes(v),3}+Lp-1;
    final_probe_info{v,4} = probes{chosenProbes(v),1};
    [rev_comp] = Reverse_complement_probes({probes{chosenProbes(v),2}});  %find reverse complement
    final_probe_info{v,5} = rev_comp{1};
    final_probe_info{v,6} = length(probes{chosenProbes(v),2});
    final_probe_info{v,7} = oligoprop(rev_comp{1},'Salt',SaltConcentration,'Temp',T_hybrid).GC;
    final_probe_info{v,8} = oligoprop(rev_comp{1},'Salt',SaltConcentration,'Temp',T_hybrid).Tm(5);
    final_probe_info{v,9} = sum(TL_hits(6:end,chosenProbes(v)));
    final_probe_info{v,10} = TL_hits(5,probe_num);
    final_probe_info{v,11} = TL_hits(4,probe_num);
    final_probe_info{v,12} = TL_hits(3,probe_num);
    final_probe_info{v,13} = TL_hits(2,probe_num);
    final_probe_info{v,14} = TL_hits(1,probe_num);
end
T = cell2table(final_probe_info,...
    'VariableNames',{'Probe Target' 'Base_Position-start' 'Base_Position-stop' 'Probe ID' 'Probe Sequence' 'Probe Length' 'GC_fraction' 'Probe Tm' '20+bp matches' '19bp matches' '18bp matches' '17bp matches' '16bp matches' '15bp matches'});
filename = [saveRoot filesep settings.FolderName filesep settings.FolderName '_probes_final_' num2str(max_probes) 'max.xlsx'];
if exist(filename,'file')        %delete if already exists
    delete(filename)
end
writetable(T,filename,'Sheet',1,'Range','A1')
size(final_probe_info,1)

%% Get Metric Information (Probes, Final Probe Set)
%Using Concentrations Solve For Equilibrium and Get Distributions and probe set metrics for detection
fprintf('\n')
try
    load([settings.FolderRootName filesep inputs1{gene_num,5} '_Tm' num2str(T_hybrid) '_ModelMetrics' designerName '.mat'],'ModelMetrics')
    fprintf("Loading TrueProbes probe set kinetic-model simulation results and RNA-FISH probe design metrics")
catch
    fprintf("Computing TrueProbes probe set kinetic-model simulation results and RNA-FISH probe design metrics")
    try
        tic
        ModelMetrics = ...
            RNAsolver_JH(chosenProbes,settings,probes,gene_table,ExpressionMatrix,DoesProbeBindSite2,dHeq_mod,dSeq_mod,dCp_mod,dHeq_Complement,dSeq_Complement,dCp_Complement)
        save([settings.FolderRootName filesep inputs1{gene_num,5} '_Tm' num2str(T_hybrid) '_ModelMetrics' designerName '.mat'],'ModelMetrics','chosenProbes','settings','-v7.3')
        tEnd = toc;fprintf('\n')
        fprinf("Time elapsed to select TrueProbes probes %g seconds",round(tEnd,3,"significant"))
    catch e
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end
%% Get Output File with summary of analysis on designed probes

%% Filter out probes below specific specficity

%% Make/Save Output Figures & Figure Objects
%Output plots and excel or csv table with statistics

end
