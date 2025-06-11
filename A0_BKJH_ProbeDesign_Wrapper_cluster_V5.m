function A0_BKJH_ProbeDesign_Wrapper_cluster_V5(id,cluster)
input_file = 'TrueProbes_DesignTargets.csv';
input_parameters = 'TrueProbes_ParameterSettings.xml';
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
% 	1. Organism you are trying to design probes in.
% 	2. Included target Accession IDs. Designs probes shared across all target accession number(s)
% 	3. Excluded target accession IDs. Removes probes in any exclusion accession number(s). Default is Empty
% 	4. Text Sequence Files to Include (files). Default is Empty
% 	5. Text Sequence Files to Exclude (files). Default is Empty

%% Clean up if necessary
try ProgressBar.deleteAllTimers(); end;

%% Genes To Design Probes For
core_folders = genpath(strcat(pwd,filesep,'src',filesep,'core'));
modified_matlab_function_folders = genpath(strcat(pwd,filesep,'src',filesep,'modified_matlab_functions'));
util_folders = genpath(strcat(pwd,filesep,'src',filesep,'util'));
wrappers_folders = genpath(strcat(pwd,filesep,'src',filesep,'wrappers'));
mtp_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'modified_thirdparty'));
nd_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'ndSparse_G4_2025_01_19'));
multiprod_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'Multiprod_2009'));
xsum_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'XSum_2014_06_16'));
xml_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'xml2struct-master'));
trueSpot_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'TrueSpot-main'));
prog_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'MatlabProgressBar'));
addpath(core_folders);
addpath(modified_matlab_function_folders);
addpath(util_folders);
addpath(wrappers_folders);
addpath(nd_folders);
addpath(mtp_folders);
addpath(multiprod_folders);
addpath(xsum_folders);
addpath(xml_folders);
addpath(trueSpot_folders);
addpath(prog_folders);


gene_num = id;
saveRoot = strcat('output',filesep); 
input_file_opts = detectImportOptions(input_file);
inputs1 = readmatrix(input_file,input_file_opts);

Organism = inputs1{gene_num,1};
IncludeAccessionNumbers = split(inputs1{gene_num,2},',');
ExcludeAccessionNumbers = split(inputs1{gene_num,3},',');
InclusionSequenceFiles = split(inputs1{gene_num,4},',');
ExcludeSequenceFiles = split(inputs1{gene_num,5},',');
if (isempty(IncludeAccessionNumbers{:}))
      IncludeAccessionNumbers = {};
end
if (isempty(ExcludeAccessionNumbers{:}))
      ExcludeAccessionNumbers = {};
end
if (isempty(ExcludeSequenceFiles{:}))
      ExcludeSequenceFiles = {};
end
if (isempty(InclusionSequenceFiles{:}))
      InclusionSequenceFiles= {};
end
RunOffline = 1;%Is Design Run Offline using databases or looks online to get genbank record

inputsParameterSettings = readstruct(input_parameters);
settings.otherLocRoot =  inputsParameterSettings.customOrganism_Database_Location.customRoot_FASTA;
settings.otherBlastDatabase = inputsParameterSettings.customOrganism_Database_Location.customBlastDatabase_DNA; % Location of custom DNA database
settings.otherBlastDatabase2 = inputsParameterSettings.customOrganism_Database_Location.customBlastDatabase_RNA; % Location of custom RNA database
settings.otherGFF =  inputsParameterSettings.customOrganism_Database_Location.custom_GFF;
settings.otherGTF =  inputsParameterSettings.customOrganism_Database_Location.custom_GTF;

%% Main Probe Design Settings (Usually might change)
minProbeSize = inputsParameterSettings.MainProbe_Settings.minProbeSize;% Min Probe Size
maxProbeSize = inputsParameterSettings.MainProbe_Settings.maxProbeSize;% Max Probe Size
minProbeSpacing = inputsParameterSettings.MainProbe_Settings.minProbeSpacing; %Minimum Spacing Between Adjacent Probes
maxNumberOfProbes = inputsParameterSettings.MainProbe_Settings.maxNumberOfProbes; % Max number of probes to design
MinHomologySearchTargetSize =  inputsParameterSettings.MainProbe_Settings.MinHomologySearchTargetSize; % minimum off-target match size
BLASTrna = inputsParameterSettings.MainProbe_Settings.BLAST_RNA; %Decides if you will BLAST the reference transcriptome (RNA targets)
BLASTdna = inputsParameterSettings.MainProbe_Settings.BLAST_DNA; %Decides if you will BLAST the reference genome (DNA targets)
ExpressionReferenceForProbeDesign = inputsParameterSettings.MainProbe_Settings.ExpressionReferenceForProbeDesign; % which expression reference to use for designing probes, default 0 meaning equal expression.
targetStrand = inputsParameterSettings.MainProbe_Settings.targetStrand; %Target strand to design probes against. 1 for RNA

%% Design Filtering Settings (You Usually will not change)
RemoveProbesBindingRibosomalHits = inputsParameterSettings.DesignFiltering_Settings.RemoveProbesBindingRibosomalHits;% Filter out probes with targets hits to ribosomal proteins
packOptimal_ProbesWithNoOffTargets = inputsParameterSettings.DesignFiltering_Settings.packOptimal_ProbesWithNoOffTargets;%when designing probes without off-targets do optimal packing to get the most or normal selection not considering packing efficiency
IncludeSelfHybridizationInProbeSelection = inputsParameterSettings.DesignFiltering_Settings.IncludeSelfHybridizationInProbeSelection;%when designing probes consider probe self-hybridization when ranking probes on binding affinity

%% Thermodynamic Settings
Gibbs_Model = inputsParameterSettings.Thermodynamic_Settings.Gibbs_Model; %Which Hybridization model to use for probe design and evaluation
SaltConcentration = inputsParameterSettings.Thermodynamic_Settings.SaltConcentration; %Concentration of Salt in thermodynamic calculations mol/L
HybridizationTemperatureCelsius = inputsParameterSettings.Thermodynamic_Settings.HybridizationTemperature_Celsius; % Temperature for Evaluating Probe Design and Simulations
Tref = inputsParameterSettings.Thermodynamic_Settings.HeatCapacityReferenceTemperature_Celsius;
RemoveMisMatches = inputsParameterSettings.Thermodynamic_Settings.RemoveMisMatches;%remove mismatches from nearest neighbor quantification of probe binding

%% Paralellization Parameters (You Usually will not change unless the gene is very long)
Parallelization_probeBatchSize = inputsParameterSettings.Parallelization_Settings.Parallelization_probeBatchSize;%batch size for parallelizing probe evaluations in probe design
Parallelization_targetBatchSize = inputsParameterSettings.Parallelization_Settings.Parallelization_targetBatchSize;%batch size for parallelizing target evaluations in probe design
ParsingPreference = inputsParameterSettings.Parallelization_Settings.ParsingPreference; % Sequentially Parse or Parallel parsing of blast hits together into hit table (1 sequentially, 0 parallel)

%% Gene Expression Parameters (You Usually will not change)
DoAllGenesHaveSameExpression =  inputsParameterSettings.GeneExpression_Settings.DoAllGenesHaveSameExpression;%do all genes have the same expression level when designing probes and evaluating probe set metrics
UseGeneOverTranscLevelExpression = inputsParameterSettings.GeneExpression_Settings.UseGeneOverTranscLevelExpression;%Use gene level over transcript level expression information when designing probes  % 1 (Gene/EMBL GENEID) , 0 (Transcript/EMBL Transcript ID)
UseRegularDNAExpression = inputsParameterSettings.GeneExpression_Settings.UseRegularDNAExpression;%0 use DNA expression from gene expression track in expression data, 1 set expression to 2 for DNA.
nullRNAcopynumber = inputsParameterSettings.GeneExpression_Settings.nullRNAcopynumber;%average cell RNA copy number when evaluating probes without expression
nullDNAcopynumber = inputsParameterSettings.GeneExpression_Settings.nullDNAcopynumber;%average cell DNA copy number when evaluating probes without expression
logRatioTrim = inputsParameterSettings.GeneExpression_Settings.TMM_LogRatioTrim;
sumTrim = inputsParameterSettings.GeneExpression_Settings.TMM_SumTrim;
Acutoff = inputsParameterSettings.GeneExpression_Settings.TMM_Acutoff;
doWeighting = inputsParameterSettings.GeneExpression_Settings.TMM_doWeighting;

%% BLAST Parameters (You Usually will not change)
evalue =  inputsParameterSettings.BLASTN_Settings.evalue; %BLAST expectation value cutoff
dust =  convertStringsToChars(inputsParameterSettings.BLASTN_Settings.dust); % filter BLAST query sequences with DUST
gapextend = inputsParameterSettings.BLASTN_Settings.gapextend; %BLAST cost to extend a gap
gapopen = inputsParameterSettings.BLASTN_Settings.gapopen; %BLAST cost to open a gap
num_alignments = inputsParameterSettings.BLASTN_Settings.numAlignments; %Number of database sequneces to show num_alignments for
penalty = inputsParameterSettings.BLASTN_Settings.penalty; % BLAST penalty for nucleotide mismatch
reward = inputsParameterSettings.BLASTN_Settings.reward; % BLAST reward for a nucleotide match
word_size = inputsParameterSettings.BLASTN_Settings.wordsize; %BLAST word size for wordfinder algorithm 

%% Model Prediction Settings (You Usually will not change)
removeUndesiredIsoformsFromPredictionOffTargets=inputsParameterSettings.ModelSimulation_Settings.removeUndesiredIsoformsFromPrediction;%Removes other isoforms from the computation of probe statistics used in designing probes
AutoBackground_Mean = inputsParameterSettings.ModelSimulation_Settings.AutoFluoresenceBackground_MEAN;
AutoBackground_STD = inputsParameterSettings.ModelSimulation_Settings.AutoFluoresenceBackground_STD;
NumReferenceProbes = inputsParameterSettings.ModelSimulation_Settings.NumberOfProbesInReferenceSpots;
SpotIntensity_Mean = inputsParameterSettings.ModelSimulation_Settings.ReferenceSpotIntensity_MEAN;
SpotIntensity_STD = inputsParameterSettings.ModelSimulation_Settings.ReferenceSpotIntensity_STD;
Rspot = inputsParameterSettings.ModelSimulation_Settings.SpotRadius_Pixels;
Mean_Diameter = inputsParameterSettings.ModelSimulation_Settings.CellDiameter_Pixels;
Nstacks = inputsParameterSettings.ModelSimulation_Settings.NumberOfReferenceZStacks;
Rcell = inputsParameterSettings.ModelSimulation_Settings.CellRadius_Micron;%um to L
errThreshold = inputsParameterSettings.ModelSimulation_Settings.SolutionErrorTolerance;
MaxIter = inputsParameterSettings.ModelSimulation_Settings.MaxRecursiveEquilibriiumIterations;
InitialGuessConc = inputsParameterSettings.ModelSimulation_Settings.InitialFreeSolutionGuessConcentration_MicroMolar;
ProbeConcentration = inputsParameterSettings.ModelSimulation_Settings.ProbeConcentration_MicroMolar;%5uM
Dilution_Vector =  double(split(inputsParameterSettings.ModelSimulation_Settings.Dilution_Vector,','));
Temperature_Celsius_Model_Vector =  double(split(inputsParameterSettings.ModelSimulation_Settings.Dilution_Vector,','));
Gibbs_Model_Vector =  double(split(inputsParameterSettings.ModelSimulation_Settings.Gibbs_Model_Vector ,','));
Signal_StepSize = inputsParameterSettings.ModelSimulation_Settings.SignalStepSize;
Signal_MaxValue = inputsParameterSettings.ModelSimulation_Settings.SignalMaxValue;

if (isempty(Organism))
    msg = 'Error. There must be a organism specified with the design in order to quantify target hits in the reference genome or transcriptome.';
    error(msg)
end
if (isempty(IncludeAccessionNumbers)*isempty(InclusionSequenceFiles)==1)
    msg = 'Error. Must include at least one target Accession Number or target inclusion sequence file.';
    error(msg)
end
if (minProbeSize>maxProbeSize)
    msg = 'Error. Minimum probe size must be less than or equal to max probe size.';
    error(msg)
end
if (minProbeSpacing<0)
    msg = 'Error. Minimum probe spacing must be greater than or equal to zero.';
    error(msg)
end
if (maxNumberOfProbes<0)
    msg = 'Error. Maximum number of probes must be greater than zero.';
    error(msg)
end
if (sum(mod([gapextend gapopen num_alignments penalty reward word_size],1))>0)
    msg = 'Error. BLAST settings must be integers.';
    error(msg)
end
if (penalty>0)
    msg = 'Error. BLAST nucleotide mismatch penalty must be less than or equal to zero.';
    error(msg)
end
if (reward<0)
    msg = 'Error. BLAST nucleotide match reward must be greater than or equal to zero.';
    error(msg)
end
if (word_size<2)
    msg = 'Error. BLAST nucleotide match reward must be greater than or equal to two.';
    error(msg)
end
if (isMATLABReleaseOlderThan("R2024a"))
    msg = 'Error. \n MATLAB must Be version 2024a or higher for blastplus to work.';
    error(msg)
end
addons = matlab.addons.installedAddons;
if (~ismember("Bioinformatics Toolbox",addons.Name))
    msg = 'Error. You need to have the Bioinformatics Toolbox installed for the software to work properly (https://www.mathworks.com/products/bioinfo.html)';
    error(msg)
end
if (~ismember("Curve Fitting Toolbox",addons.Name))
    msg = 'Error. The Curve Fitting Toolbox must be installed for the software to work properly: (https://www.mathworks.com/products/curvefitting.html)';
       error(msg)
end
if (~ismember("Parallel Computing Toolbox",addons.Name))
    msg = 'Error. The Parallel Computing Toolbox must be installed for the software to work properly: (https://www.mathworks.com/products/parallel-computing.html)';
     error(msg)
end
if (~ismember("Signal Processing Toolbox",addons.Name))
    msg = 'Error. The Signal Processing Toolbox must be installed for the software to work properly: (https://www.mathworks.com/products/signal.html)';
     error(msg)
end
if (~ismember("Statistics and Machine Learning Toolbox",addons.Name))
    msg = 'Error. The Statistics and Machine Learning Toolbox must be installed for the software to work properly: (https://www.mathworks.com/products/statistics.html)';
       error(msg)   
end
if (~ismember("Symbolic Math Toolbox",addons.Name))
    msg = 'Error. The Symbolic Math Toolbox installed must be installed for the software to work properly: (https://www.mathworks.com/products/symbolic.html)';
   error(msg)
end
if (~ismember("BLAST+ Support Package for Bioinformatics Toolbox",addons.Name))
    msg = 'Error. The BLAST+ Support Package for Bioinformatics Toolbox must be installed for the software to work properly: (https://www.mathworks.com/matlabcentral/fileexchange/156414-blast-support-package-for-bioinformatics-toolbox)';
   error(msg) 
end
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
if (~strcmp(pwd,extractBefore(mfilePath,strcat(filesep,'A0'))))
    msg = 'Error. The script must be run in the TrueProbes main folder for the code to work properly';
   error(msg) 
end




designerName = '_TrueProbes';
probes = [];
startup;
settings.BUILD_STRING = '2025.06.04.00';
settings.VERSION_STRING = 'v1.1.1';
refInfo = IncludeAccessionNumbers{1}(1:2);
if (ismember(refInfo,{'NR','XR','NM','XM'}))
    settings.referenceType = 'RefSeq';
else
    settings.referenceType = 'ENSEMBL';
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
%% Update Location of Databases & Needed Files (You usually will not change)
if (ispc)
    settings.blastpath = strcat('src',filesep,'blast',filesep,'ncbi-blast-2.16.0+-x64-win64',filesep,'ncbi-blast-2.16.0+',filesep,'bin');
end
if (isunix)
    settings.blastpath  = strcat('src',filesep,'blast',filesep,'ncbi-blast-2.16.0+-x64-linux',filesep,'ncbi-blast-2.16.0+',filesep,'bin');
end
if (ismac)
    settings.blastpath = strcat('src',filesep,'blast',filesep,'ncbi-blast-2.16.0+-universal-macosx',filesep,'ncbi-blast-2.16.0+',filesep,'bin');
end
OrganismList = {'Human';'Mouse';'Rat';'Zebrafish';'Fruitfly';'Worm';'Yeast'};
DatabaseVars = {'Organism','Root_FASTA','BLASTDB_DNA','BLASTDB_RNA','GTF','GFF'};
ENSEMBL_ParserRNA = {'ENST','ENSMUST','ENSRNOT','ENSDART','','',''};
%DNA and NonDNA not works for non-ensembl format, need to use gtf to find strings in assembly as list to match
if (strcmp(settings.referenceType,'RefSeq'))
    LocDatabases = {'data/DatabaseData/Blast_Databases/Human/NCBI_RefSeq/','data/DatabaseData/Blast_Databases/Human/NCBI_RefSeq/Human_NCBI_genomic','data/DatabaseData/Blast_Databases/Human/NCBI_RefSeq/Human_NCBI_transcript',...
                                'data/DatabaseData/GTF_Databases/Human/NCBI_RefSeq/GCF_000001405.40_GRCh38.p14_genomic.gtf','data/DatabaseData/GFF3_Databases/Human/NCBI_RefSeq/GCF_000001405.40_GRCh38.p14_genomic.gff';
                                'data/DatabaseData/Blast_Databases/Mouse/NCBI_RefSeq/','data/DatabaseData/Blast_Databases/Mouse/NCBI_RefSeq/Mouse_NCBI_genomic','data/DatabaseData/Blast_Databases/Mouse/NCBI_RefSeq/Mouse_NCBI_transcript',...
                                'data/DatabaseData/GTF_Databases/Mouse/NCBI_RefSeq/GCF_000001635.27_GRCm39_genomic.gtf','data/DatabaseData/GFF3_Databases/Mouse/NCBI_RefSeq/GCF_000001635.27_GRCm39_genomic.gff';...
                                'data/DatabaseData/Blast_Databases/Rat/NCBI_RefSeq/','data/DatabaseData/Blast_Databases/Rat/NCBI_RefSeq/Rat_NCBI_genomic','data/DatabaseData/Blast_Databases/Rat/NCBI_RefSeq/Rat_NCBI_transcript',...
                                'data/DatabaseData/GTF_Databases/Rat/NCBI_RefSeq/_genomic.gtf','data/DatabaseData/GFF3_Databases/Rat/NCBI_RefSeq/_genomic.gff';...
                                'data/DatabaseData/Blast_Databases/Zebrafish/NCBI_RefSeq/','data/DatabaseData/Blast_Databases/Zebrafish/NCBI_RefSeq/Zebrafish_NCBI_genomic','data/DatabaseData/Blast_Databases/Zebrafish/NCBI_RefSeq/Zebrafish_NCBI_transcript',...
                                'data/DatabaseData/GTF_Databases/Zebrafish/NCBI_RefSeq/GCF_049306965.1_GRCz12tu_genomic.gtf','data/DatabaseData/GFF3_Databases/Zebrafish/NCBI_RefSeq/GCF_049306965.1_GRCz12tu_genomic.gff';
                                'data/DatabaseData/Blast_Databases/Fruitfly/NCBI_RefSeq/','data/DatabaseData/Blast_Databases/Fruitfly/NCBI_RefSeq/Fruitfly_NCBI_genomic','data/DatabaseData/Blast_Databases/Fruitfly/NCBI_RefSeq/Fruitfly_NCBI_transcript',...
                                'data/DatabaseData/GTF_Databases/Fruitfly/NCBI_RefSeq/_genomic.gtf','data/DatabaseData/GFF3_Databases/Fruitfly/NCBI_RefSeq/_genomic.gff';...
                                'data/DatabaseData/Blast_Databases/Worm/NCBI_RefSeq/','data/DatabaseData/Blast_Databases/Worm/NCBI_RefSeq/Worm_NCBI_genomic','data/DatabaseData/Blast_Databases/Worm/NCBI_RefSeq/Worm_NCBI_transcript',...
                                'data/DatabaseData/GTF_Databases/Worm/NCBI_RefSeq/_genomic.gtf','data/DatabaseData/GFF3_Databases/Worm/NCBI_RefSeq/_genomic.gff';...
                                'data/DatabaseData/Blast_Databases/Yeast/NCBI_RefSeq/','data/DatabaseData/Blast_Databases/Yeast/NCBI_RefSeq/Yeast_NCBI_genomic','data/DatabaseData/Blast_Databases/Yeast/NCBI_RefSeq/Yeast_NCBI_transcript',...
                                'data/DatabaseData/GTF_Databases/Yeast/NCBI_RefSeq/GCF_000146045.2_R64_genomic.gtf','data/DatabaseData/GFF3_Databases/Yeast/NCBI_RefSeq/GCF_000146045.2_R64_genomic.gff';...
                                };                
elseif (strcmp(settings.referenceType,'ENSEMBL'))
          LocDatabases = {'data/DatabaseData/Blast_Databases/Human/EMBL_EBI/','data/DatabaseData/Blast_Databases/Human/EMBL_EBI/Human_ENSEMBL_genomic','data/DatabaseData/Blast_Databases/Human/EMBL_EBI/Human_ENSEMBL_transcript',...
                                'data/DatabaseData/GTF_Databases/Human/EMBL_EBI/Homo_sapiens.GRCh38.114.gtf','data/DatabaseData/GFF3_Databases/Human/EMBL_EBI/Homo_sapiens.GRCh38.114.gff3';
                                'data/DatabaseData/Blast_Databases/Mouse/EMBL_EBI/','data/DatabaseData/Blast_Databases/Mouse/EMBL_EBI/Mouse_ENSEMBL_genomic','data/DatabaseData/Blast_Databases/Mouse/EMBL_EBI/Mouse_ENSEMBL_transcript',...
                                'data/DatabaseData/GTF_Databases/Mouse/EMBL_EBI/Mus_musculus.GRCm39.114.gtf','data/DatabaseData/GFF3_Databases/Mouse/EMBL_EBI/Mus_musculus.GRCm39.114.gff3';...
                                'data/DatabaseData/Blast_Databases/Rat/EMBL_EBI/','data/DatabaseData/Blast_Databases/Rat/EMBL_EBI/Rat_ENSEMBL_genomic','data/DatabaseData/Blast_Databases/Rat/EMBL_EBI/Rat_ENSEMBL_transcript',...
                                'data/DatabaseData/GTF_Databases/Rat/EMBL_EBI/.114.gtf','data/DatabaseData/GFF3_Databases/Rat/EMBL_EBI/.114.gff3';...
                                'data/DatabaseData/Blast_Databases/Zebrafish/EMBL_EBI/','data/DatabaseData/Blast_Databases/Zebrafish/EMBL_EBI/Zebrafish_ENSEMBL_genomic','data/DatabaseData/Blast_Databases/Zebrafish/EMBL_EBI/Zebrafish_ENSEMBL_transcript',...
                                'data/DatabaseData/GTF_Databases/Zebrafish/EMBL_EBI/Danio_rerio.GRCz11.114.gtf','data/DatabaseData/GFF3_Databases/Zebrafish/EMBL_EBI/Danio_rerio.GRCz11.114.gff3'; 
                                'data/DatabaseData/Blast_Databases/Fruitfly/EMBL_EBI/','data/DatabaseData/Blast_Databases/Fruitfly/EMBL_EBI/Fruitfly_ENSEMBL_genomic','data/DatabaseData/Blast_Databases/Fruitfly/EMBL_EBI/Fruitfly_ENSEMBL_transcript',...
                                'data/DatabaseData/GTF_Databases/Fruitfly/EMBL_EBI/.114.gtf','data/DatabaseData/GFF3_Databases/Fruitfly/EMBL_EBI/.114.gff3';...
                                'data/DatabaseData/Blast_Databases/Worm/EMBL_EBI/','data/DatabaseData/Blast_Databases/Worm/EMBL_EBI/Worm_ENSEMBL_genomic','data/DatabaseData/Blast_Databases/Worm/EMBL_EBI/Worm_ENSEMBL_transcript',...
                                'data/DatabaseData/GTF_Databases/Worm/EMBL_EBI/.114.gtf','data/DatabaseData/GFF3_Databases/Worm/EMBL_EBI/.114.gff3';...
                                'data/DatabaseData/Blast_Databases/Yeast/EBML_EBI/','data/DatabaseData/Blast_Databases/Yeast/EMBL_EBI/Yeast_ENSEMBL_genomic','data/DatabaseData/Blast_Databases/Yeast/EMBL_EBI/Yeast_ENSEMBL_transcript',...
                                'data/DatabaseData/GTF_Databases/Yeast/EMBL_EBI/Saccharomyces_cerevisiae.R64-1-1.1.114.gtf','data/DatabaseData/GFF3_Databases/Yeast/EMBL_EBI/Saccharomyces_cerevisiae.R64-1-1.1.114.gff3';...
                                };
end
DatabaseLocations = cell2table([OrganismList LocDatabases]);
DatabaseLocations.Properties.VariableNames = DatabaseVars;
for k = 2:length(DatabaseVars)
settings.(strcat("Loc",DatabaseLocations.Properties.VariableNames{k})) = dictionary(convertCharsToStrings(DatabaseLocations.Organism),convertCharsToStrings(DatabaseLocations.(DatabaseLocations.Properties.VariableNames{k})));
end
settings.EMBL_RNAparser = dictionary(convertCharsToStrings(DatabaseLocations.Organism),convertCharsToStrings(ENSEMBL_ParserRNA)');
%Set BLAST and Reference Genome/Transcriptome Location
if (isKey(settings.LocRoot_FASTA,Organism))
    settings.SEQdbRoot = char(settings.LocRoot_FASTA(Organism));
else
    settings.SEQdbRoot = char(settings.otherLocRoot);
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
settings.CustomGenomeAssemblyReportFile = 'N/A';
settings.Custom_X_ChromNumber = 1;
settings.Custom_Y_ChromNumber = 1;
settings.Custom_MT_ChromNumber = 1;
settings.scTracks = {'colonWangCellType','tabulasapiens_tissue_cell_type'};


%% Load Annotation File
fprintf('\n')
fprintf("Loading genome and transcriptome annotation files")
fprintf('\n')
tic
if (isKey(settings.LocRoot_FASTA,Organism))
    if isfile([settings.LocGFF(Organism)])
        GFFobj = GFFAnnotation(settings.LocGFF(Organism));
    else
         msg = strcat('Error. Need to provide GFF Location for'," ",Organism);
        error(msg)
    end
else
    if isfile([settings.otherGFF])
        GFFobj = GFFAnnotation(settings.otherGFF);
    else
        msg = 'Error. Need to provide GFF Location for the custom organism.';
        error(msg)
    end
end
    try
        if (isKey(settings.LocRoot_FASTA,Organism))
             if isfile([settings.LocGTF(Organism)])
             GTFobj = GTFAnnotation(settings.LocGTF(Organism));
             else
                 msg = strcat('Error. Need to provide GTF Location for'," ",Organism);
                 error(msg)
             end
        else
            if isfile([settings.otherGTF])
                GTFobj = GFFAnnotation(settings.otherGTF);
            else
                msg = 'Error. Need to provide GTF Location for the custom organism.';
                error(msg)
            end
        end
    catch
        if (isKey(settings.LocRoot_FASTA,Organism))
            fixedEmptyTranscriptIdEntriesFile = strcat(extractBefore(settings.LocGTF(Organism),'.gtf'),'_fixed.gtf');
        else
            fixedEmptyTranscriptIdEntriesFile = strcat(extractBefore(settings.otherGTF,'.gtf'),'_fixed.gtf');
        end
        if (isfile(fixedEmptyTranscriptIdEntriesFile))
            GTFobj = GTFAnnotation(fixedEmptyTranscriptIdEntriesFile);
        else
            if (isKey(settings.LocRoot_FASTA,Organism))
            originalFile = settings.LocGTF(Organism);
            else
             originalFile = settings.otherGTF;  
            end
            % Read in GTF file (skip last empty line)
            gtfLines = readlines(originalFile, EmptyLineRule="skip");
            % Find entries that do not contain "transcript_id" attribute (exclude header)
            invalidEntryIdx = ~contains(gtfLines, ("transcript_id" | textBoundary + "#"));
            % Add empty "transcript_id" attribute to the invalid entries
            gtfLines(invalidEntryIdx) = gtfLines(invalidEntryIdx) + ' transcript_id "";';
            % Write to a new GTF file that has the invalid entries fixed
            writelines(gtfLines, fixedEmptyTranscriptIdEntriesFile);
            GTFobj = GTFAnnotation(fixedEmptyTranscriptIdEntriesFile);
        end
    end
tEnd = toc;
fprintf('\n')
fprintf("Time elapsed to load annotation files %g seconds",round(tEnd,3,"significant"))
ids_gtf = GTFobj.Transcript(strcmp(GTFobj.Feature,'transcript'));
btypes = extractBetween(string(regexp(GTFobj.Attributes(strcmp(GTFobj.Feature,'transcript')),'(?<=transcript_biotype )"\w*"(?=;)','match')),'"','"');
Ribosomal_IDs = unique(ids_gtf(contains(btypes,'rRNA')));
if (strcmp(settings.referenceType,'RefSeq'))
    geneInfo = getGenes(GTFobj,"Transcript",IncludeAccessionNumbers);
    geneNames = geneInfo.GeneID{:};
    ids_gff_loc = find(strcmp(GFFobj.Feature,'region').*contains(GFFobj.Attributes,'NC'));
    chr_list = regexp(GFFobj.Attributes(ids_gff_loc),'(?<=Name=)\w*(?=;)','match');
    chrom_list = string(convertCharsToStrings(chr_list));
    ref_list = regexp(GFFobj.Attributes(ids_gff_loc),'(?<=ID=)\w*(?=.)','match');
    reference_list = string(convertCharsToStrings(ref_list));
    ChromosomeToReference = dictionary(reference_list,chrom_list);
    ReferenceToChromosome = dictionary(reference_list,chrom_list);
    geneChrNum = ReferenceToChromosome(extractBefore(string(geneInfo.Reference),'.'));
    geneReference_ID = extractBefore(string(geneInfo.Reference),'.');
elseif (strcmp(settings.referenceType,'ENSEMBL'))
    geneInfo = getGenes(GTFobj,"Transcript",extractBefore(IncludeAccessionNumbers,'.'));
    geneNames = geneInfo.GeneName{:};
    ids_gff_loc = find(strcmp(GFFobj.Feature,'chromosome'));
    chr_list = regexp(GFFobj.Attributes(ids_gff_loc),'(?<=chromosome:).*(?=;)','match');
    chrom_list = string(convertCharsToStrings(chr_list));
    ref_list = regexp(GFFobj.Attributes(ids_gff_loc),'(?<=Alias=).*(?=,)','match');
    reference_list = string(convertCharsToStrings(ref_list));
    refparts = split(reference_list,',');
    paired_ID_chrom = sort(refparts,2);
    ReferenceToChromosome = dictionary(paired_ID_chrom(:,1),chrom_list);
    ChromosomeToReference = dictionary(chrom_list,paired_ID_chrom(:,1));
    geneChrNum = string(geneInfo.Reference);
    geneReference_ID = extractBefore(ChromosomeToReference(string(geneInfo.Reference)),'.');
end
if (strcmp(settings.referenceType,'RefSeq'))
geneInfo_Table = getGenes(GTFobj,"Transcript",IncludeAccessionNumbers);
transcriptInfo_Table = getTranscripts(GTFobj,"Transcript",IncludeAccessionNumbers);
all_isoform_transcriptInfo_Table= getTranscripts(GTFobj,"Gene",getGenes(GTFobj,"Transcript",IncludeAccessionNumbers).GeneID);
Transcript_All_Isoform_IDs = cellfun(@(x) all_isoform_transcriptInfo_Table.Transcript(strcmp(all_isoform_transcriptInfo_Table.GeneID,x)), transcriptInfo_Table.GeneID,'Un',0); 
Transcript_UnDesired_Isoform_IDs = arrayfun(@(x) setdiff(Transcript_All_Isoform_IDs{x},IncludeAccessionNumbers{x}),1:length(IncludeAccessionNumbers),'Un',0);
Transcript_Joint_UnDesired_Isoforms_IDs = setdiff(all_isoform_transcriptInfo_Table.Transcript,transcriptInfo_Table.Transcript);
elseif (strcmp(settings.referenceType,'ENSEMBL'))
geneInfo_Table = getGenes(GTFobj,"Transcript",extractBefore(IncludeAccessionNumbers,'.'));
transcriptInfo_Table = getTranscripts(GTFobj,"Transcript",extractBefore(IncludeAccessionNumbers,'.'));
all_isoform_transcriptInfo_Table= getTranscripts(GTFobj,"Gene",getGenes(GTFobj,"Transcript",extractBefore(IncludeAccessionNumbers,'.')).GeneID);
Transcript_All_Isoform_IDs = cellfun(@(x) all_isoform_transcriptInfo_Table.Transcript(strcmp(all_isoform_transcriptInfo_Table.GeneID,x)), transcriptInfo_Table.GeneID,'Un',0); 
Transcript_UnDesired_Isoform_IDs = arrayfun(@(x) setdiff(Transcript_All_Isoform_IDs{x},extractBefore(IncludeAccessionNumbers{x},'.')),1:length(IncludeAccessionNumbers),'Un',0);
Transcript_Joint_UnDesired_Isoforms_IDs = setdiff(all_isoform_transcriptInfo_Table.Transcript,transcriptInfo_Table.Transcript);    
end
% for individual_transcripts = 1:length(IncludeAccessionNumbers)
%     exonInfo_Table{individual_transcripts} = getExons(GTFobj,"Transcript",IncludeAccessionNumbers{individual_transcripts});
%     segmentInfo_Table{individual_transcripts} = getSegments(GTFobj,"Transcript",IncludeAccessionNumbers{individual_transcripts});
%     all_isoform_exonInfo_Table{individual_transcripts} = getExons(GTFobj,"Gene",getGenes(GTFobj,"Transcript",IncludeAccessionNumbers{individual_transcripts}).GeneID{1});
%     all_isoform_segmentInfo_Table{individual_transcripts}  = getSegments(GTFobj,"Gene",getGenes(GTFobj,"Transcript",IncludeAccessionNumbers{individual_transcripts}).GeneID{1});
% end

%% Save Settings
settings.FolderRootName = strcat(saveRoot,'(',geneNames,')','_',strjoin(IncludeAccessionNumbers,'_'));
settings.rootName = strjoin(IncludeAccessionNumbers,'_');
settings.designerName = designerName;
settings.saveRoot = saveRoot;

%% Design Target Info
settings.Organism = Organism;
settings.GeneName = geneNames;
settings.ChrNum = geneChrNum;
settings.GeneChr = strcat('chr',geneChrNum);
settings.chromosome_IDs = {char(geneInfo.Reference)};
settings.transcript_IDs = IncludeAccessionNumbers;
settings.ribosomal_IDs = Ribosomal_IDs;
settings.transcript_IDs_desired = Transcript_All_Isoform_IDs;
settings.transcript_IDs_undesired = Transcript_UnDesired_Isoform_IDs;
settings.transcript_IDs_joint_undesired = Transcript_Joint_UnDesired_Isoforms_IDs ;
settings.ProbeSpacing = minProbeSpacing;
settings.RemoveMisMatches = RemoveMisMatches;
settings.SaltConcentration = SaltConcentration;
settings.HybridizationTemperature = HybridizationTemperatureCelsius;
settings.BLASTdna = BLASTdna;
settings.BLASTrna = BLASTrna;
settings.BLASTbatchSize = Parallelization_probeBatchSize;
settings.BLASTsimultaneousParsingOverSequentialParsing = ParsingPreference;
settings.num_parpool_local = feature('numcores');
settings.BlastParameters.reward = reward;
settings.BlastParameters.penalty = penalty;
settings.BlastParameters.wordsize = word_size;
settings.BlastParameters.gapopen = gapopen;
settings.BlastParameters.gapextend = gapextend;
settings.BlastParameters.evalue = evalue;
settings.BlastParameters.num_alignments = num_alignments;
settings.BlastParameters.dust = dust;
settings.TargetBatchSize = Parallelization_targetBatchSize;
settings.MinProbeSize = minProbeSize;
settings.MaxProbeSize = maxProbeSize;
settings.MinHomologySearchTargetSize = MinHomologySearchTargetSize;
settings.N_model = Gibbs_Model;
settings.ExpressionReferenceForDesigningProbes = ExpressionReferenceForProbeDesign;
settings.TMM.logRatioTrim = logRatioTrim;
settings.TMM.sumTrim = sumTrim;
settings.TMM.Acutoff = Acutoff;
settings.TMM.doWeighting = doWeighting;
%% Gene Expression Parameters
settings.DoAllGenesHaveSameExpression = DoAllGenesHaveSameExpression;
settings.HumanSpecific.HumanExpGeneOrTransc = UseGeneOverTranscLevelExpression; % 1 (Gene/EMBL GENEID) , 0 (Transcript/EMBL Transcript ID)
settings.UseRegularDNAExpression = UseRegularDNAExpression;%0 use DNA expression from gene expression track in expression data, 1 set expression to 2 for DNA.
settings.UniformRNAExpression = nullRNAcopynumber;%if assuming no differences in gene's expression sets level.
settings.DNAPloidy = nullDNAcopynumber;


%% Cluster Parameters
settings.clusterStatus = cluster;
settings.isOffline = RunOffline;

%% Selecting Probes Specifications
settings.maxProbes = maxNumberOfProbes;
settings.RemoveProbesWithRibosomalHits = RemoveProbesBindingRibosomalHits;

%% Simulation Settings
settings.SimulationConfiguration.Temperature_Celsius_Model_Vector = Temperature_Celsius_Model_Vector;
settings.SimulationConfiguration.Gibbs_Model_Vector = Gibbs_Model_Vector;
settings.SimulationConfiguration.Dilution_Vector = Dilution_Vector;
settings.SimulationConfiguration.AutoBackground_Mean = AutoBackground_Mean;
settings.SimulationConfiguration.AutoBackground_STD = AutoBackground_STD;
settings.SimulationConfiguration.NumReferenceProbes = NumReferenceProbes;
settings.SimulationConfiguration.SpotIntensity_Mean = SpotIntensity_Mean;
settings.SimulationConfiguration.SpotIntensity_STD = SpotIntensity_STD;
settings.SimulationConfiguration.Tref = Tref;
settings.SimulationConfiguration.Mean_Diameter = Mean_Diameter;
settings.SimulationConfiguration.CellRadius = Rcell;
settings.SimulationConfiguration.SpotRadius = Rspot;
settings.SimulationConfiguration.NumOfZStacks = Nstacks;
settings.SimulationConfiguration.InitialGuessConc = InitialGuessConc;
settings.SimulationConfiguration.ProbeConcentration = ProbeConcentration;
settings.SimulationConfiguration.errThreshold = errThreshold;
settings.SimulationConfiguration.MaxIter = MaxIter;
settings.SimulationConfiguration.Signal_StepSize = Signal_StepSize;
settings.SimulationConfiguration.Signal_MaxValue = Signal_MaxValue;

T_hybrid = HybridizationTemperatureCelsius;
Lmin = minProbeSize; Lmax = maxProbeSize;

%% Generate Folder
if (not(isfolder([saveRoot])))
    mkdir([saveRoot])
end
fprintf('\n')
fprintf("Making probe target design output folder")
fprintf('\n')
if (not(isfolder([saveRoot '(' geneNames ')' '_' settings.rootName])))
    mkdir([saveRoot  '(' geneNames ')' '_' settings.rootName])
end
settings.FolderName = [ '(' geneNames ')' '_' settings.rootName];
settings.FolderRootName = strcat('output',filesep,settings.FolderName);
settings.rootName = strjoin(IncludeAccessionNumbers,'_');
settings.TargetLists = strjoin(IncludeAccessionNumbers,', ');
%% Generate Probes
fprintf('\n') 
fprintf('\n')
fprintf(strcat("Designing probes for"," ",geneNames," ","transcript id:"," ",settings.TargetLists))
fprintf('\n')
fprintf('\n')
try
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_probes' designerName '.mat'],'probes')
    fprintf("Loading probe tile sequences")
catch
    fprintf("Generating probe tile sequences")
    try
        tic
        [init_probes,~,~,~,~] = ...
            BKJH_Probe_Generator([Lmin Lmax],IncludeAccessionNumbers,InclusionSequenceFiles,ExcludeSequenceFiles,ExcludeAccessionNumbers,targetStrand,settings.isOffline,settings.SEQdbRoot);
        probes = init_probes;
        save([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_probes' designerName '.mat'],'probes','-v7.3')
        tEnd = toc;
        fprintf('\n')
        fprintf("Time elapsed to tile probes %g seconds",round(tEnd,3,"significant"))
    catch e
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end
%% BLAST Probes
fprintf('\n')
fprintf('\n')
ProgressBar.deleteAllTimers();
try
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_hits_table' designerName '.mat'],'gene_table')
    fprintf("Loading probe BLAST results table")
catch
    fprintf("BLASTING Probes in batches")
    fprintf('\n')
    try
        tic                                                                                                           %Organism           %(Gene Name)        %(Gene Name)        %ChrNum
        [~,gene_table] = Probe_checker_general_JH10(probes,Organism,geneNames,geneNames,geneChrNum,settings);
        save([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_hits_table' designerName '.mat'],'-mat','gene_table','-v7.3');
        tEnd = toc;fprintf('\n')
        fprintf("Time elapsed to generate probe BLAST results table %g seconds",round(tEnd,3,"significant"))
    catch e
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end

%% Get Gene Expression Information
fprintf('\n')
fprintf('\n')
try
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_ExpressionInfo' designerName '.mat'],'ExpressionMatrix');
    fprintf("Loading BLAST hits gene expression information")
catch
    fprintf("Getting BLAST hits gene expression information")
    fprintf('\n')
    fprintf('\n')
    try
        tic
        [ExpressionMatrix,get_expression_time] = A_JH_GetExpressionInfo_V2(gene_table,settings);
        save([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_ExpressionInfo' designerName '.mat'],'ExpressionMatrix','get_expression_time','settings','-v7.3');
        tEnd = toc;fprintf('\n')
        fprintf("Time elapsed to generate probe BLAST hits gene expression information %g seconds",round(tEnd,3,"significant"))
    catch e
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end
settings.saveRoot = saveRoot;

%% Get Thermodynamic Information (On-Target,Off-Target)
fprintf('\n')
fprintf('\n')
try
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_Tm' num2str(settings.HybridizationTemperature) '_OnOffThermoInfo' designerName '.mat'],'Kon','Kb_Match');
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_dHInfo' designerName '.mat'],'-mat','dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match');
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_dSInfo' designerName '.mat'],'-mat','dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match');
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_TmInfo' designerName '.mat'],'-mat','Tm_on','Tm_Match');
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_dCpInfo' designerName '.mat'],'-mat','dCpon_eq','dCpeq_Match');
    fprintf("Loading probe target thermodynamic information")
catch
    fprintf("Computing probe target thermodynamic information")
    fprintf('\n')
    fprintf('\n')
    try
        tic
        [Kb_Match,Kon,Koff,dHeq_Match,dSeq_Match,dHf_Match,dSf_Match,dHr_Match,dSr_Match,dCpeq_Match,dHon_eq,dSon_eq,dHon_f,dSon_f,dHon_r,dSon_r,dCpon_eq,Tm_on,Tm_Match] = ...
            A_JH_GenerateThermoInfo_V5(probes,gene_table,geneNames,settings);%add Kon Koff
        save([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_Tm' num2str(settings.HybridizationTemperature) '_OnOffThermoInfo' designerName '.mat'],'-mat','Kon','Koff','Kb_Match','-v7.3');
        save([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_dHInfo' designerName '.mat'],'-mat','dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match','-v7.3');
        save([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_dSInfo' designerName '.mat'],'-mat','dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match','-v7.3');
        save([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_TmInfo' designerName '.mat'],'-mat','Tm_on','Tm_Match','-v7.3');
        save([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_dCpInfo' designerName '.mat'],'-mat','dCpon_eq','dCpeq_Match','-v7.3');
        toc
        tEnd = toc;fprintf('\n')
        fprintf("Time elapsed to compute probe target thermodynamic information %g seconds",round(tEnd,3,"significant"))
    catch e
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end

%% Get Binding Site Mapping and Energy
fprintf('\n')
fprintf('\n')
try
    load([settings.FolderRootName filesep settings.GeneName '_binding_hits_map' designerName '.mat'],'DoesProbeBindSite2')
    load([settings.FolderRootName filesep settings.GeneName  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' designerName '.mat'],'Kb_Complement')
    load([settings.FolderRootName filesep settings.GeneName  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' designerName '.mat'],'Kb_mod')
    load([settings.FolderRootName filesep settings.GeneName '_BindingMatrices' designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','dCp_mod')
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
        fprintf("Time elapsed to create probe target binding site maps %g seconds",round(tEnd,3,"significant"))
    catch e
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end

%% Basic Stats for Designing Probes Function
if (ExpressionReferenceForProbeDesign==0)
EKernel = ones(size(ExpressionMatrix,1),1);
else
EKernel = ExpressionMatrix(:,ExpressionReferenceForProbeDesign);
end
Kon = squeeze(Kon(:,Gibbs_Model));
Kb_mod = squeeze(Kb_mod(:,:,:,Gibbs_Model));
DoesProbeBindSite = DoesProbeBindSite2;
FoldName = [];
if (settings.BLASTdna)
    Kb_Complement = squeeze(Kb_Complement(:,:,Gibbs_Model));
else
     dHeq_Complement= [];
    dCp_Complement= [];
    dSeq_Complement= [];
    Kb_Complement = [];
end
fprintf('\n')
fprintf('\n')
try
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_Tm' num2str(T_hybrid) '_BasicDesignerStats' designerName '.mat'],'Tvec_RNA','Svec_RNA','TPvec_RNA','TSvec_RNA','TPvec_logKOFF_RNA','TPvec_logKOFFdivON_RNA','TPvec_logKONdivOFF_RNA',...
        'Nvec_RNAmulti','Tvec_DNA','Svec_DNA','TPvec_DNA','TSvec_DNA','TPvec_logKOFF_DNA','TPvec_logKOFFdivON_DNA','TPvec_logKONdivOFF_DNA','TPvec_logKOFFdivCOMP_DNA','TPvec_logKCOMPdivOFF_DNA',...
        'Off_Score','Specificity_Score','NumRNAOffTargetOptions','Probes_WithNRNAOFF','NumDNAOffTargetOptions','Probes_WithNDNAOFF')
    fprintf("Loading probe target statistics information")
catch
    try
        fprintf("Computing probe target statistics information")
        tic
        [Nvec_RNAmulti,RNAOFF_Score,RNASpecificity_Score,NumRNAOffTargetOptions,Probes_WithNRNAOFF,DNAOFF_Score,DNASpecificity_Score,NumDNAOffTargetOptions,Probes_WithNDNAOFF,Cout] = ...
            A0_BasicDesignerStats([settings.BLASTrna settings.BLASTdna],removeUndesiredIsoformsFromPredictionOffTargets,gene_table,settings,settings.FolderRootName,DoesProbeBindSite,Kon,Kb_mod,Kb_Complement,EKernel);
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
        fprintf("Time elapsed to compute probe target statistics %g seconds",round(tEnd,3,"significant"))
    catch e
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end

%% Selection of TrueProbes Probes
fprintf('\n')
fprintf('\n')
try
    load([saveRoot filesep settings.FolderName filesep settings.FolderName '_chosen.mat'],'chosenProbes')
    fprintf("Loading TrueProbes designed probes")
catch
    try
        fprintf("Designing TrueProbes probe set")
        tic
        chosenProbes = A_ZigZagProbeSelection_V5(probes,gene_table,settings,IncludeSelfHybridizationInProbeSelection,packOptimal_ProbesWithNoOffTargets,Kon,Nvec_RNAmulti,Off_Score,Specificity_Score,Tvec_RNA,Svec_RNA,TPvec_RNA,TSvec_RNA,TPvec_logKOFF_RNA,TPvec_logKOFFdivON_RNA,TPvec_logKONdivOFF_RNA,ExpressionMatrix,DoesProbeBindSite,Kb_mod);
        save([saveRoot filesep settings.FolderName filesep settings.FolderName '_chosen.mat'],'chosenProbes','-v7.3')
        tEnd = toc;fprintf('\n')
        fprintf("Time elapsed to select TrueProbes probes %g seconds",round(tEnd,3,"significant"))
    catch e
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end

%% Print Excel Spreedsheet of Probes
fprintf('\n')
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
filename = [saveRoot filesep settings.FolderName filesep settings.FolderName '_probes_final_' num2str(maxNumberOfProbes) 'max.xlsx'];
if exist(filename,'file')        %delete if already exists
    delete(filename)
end
writetable(T,filename,'Sheet',1,'Range','A1')
size(final_probe_info,1)

%% Get Metric Information (Probes, Final Probe Set)
%Using Concentrations Solve For Equilibrium and Get Distributions and probe set metrics for detection
fprintf('\n')
fprintf('\n')
try
    load([settings.FolderRootName filesep '(' geneNames ')_Tm' num2str(T_hybrid) '_ModelMetrics' designerName '.mat'],'ModelMetrics')
    fprintf("Loading TrueProbes probe set kinetic-model simulation results and RNA-FISH probe design metrics")
catch
    fprintf("Computing TrueProbes probe set kinetic-model simulation results and RNA-FISH probe design metrics")
    try
        tic
        ModelMetrics = ...
            RNAsolver_JH2(chosenProbes,settings,probes,gene_table,ExpressionMatrix,DoesProbeBindSite2,dHeq_mod,dSeq_mod,dCp_mod,dHeq_Complement,dSeq_Complement,dCp_Complement)
        save([settings.FolderRootName filesep '(' geneNames ')_Tm' num2str(T_hybrid) '_ModelMetrics' designerName '.mat'],'ModelMetrics','chosenProbes','settings','-v7.3')
        tEnd = toc;fprintf('\n')
        fprintf("Time elapsed to compute TrueProbes probe metrics %g seconds",round(tEnd,3,"significant"))
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
