function A0_BKJH_ProbeDesign_Wrapper_cluster_V5(id,cluster)
input_file = 'TrueProbes_DesignTargets.csv';
input_parameters = 'TrueProbes_ParameterSettings.xml';
input_databases_file = 'DatabaseLocations.xml';
input_gene_expression_file_locations = 'GeneExpressionDataLocations.xml';
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

%% Folders with files to add to path
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
progbar_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'parwaitbar'));
VarGibbs_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'VarGibbs-4.1'));
hpf_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'HighPrecisionFloat'));
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
addpath(progbar_folders);
addpath(VarGibbs_folders);
addpath(hpf_folders);

%% Parse Input FIles
gene_num = id;
saveRoot = strcat('output',filesep);
input_file_opts = detectImportOptions(input_file);
inputs1 = readmatrix(input_file,input_file_opts);

Organism = inputs1{gene_num,1};
IncludeAccessionNumbers = split(inputs1{gene_num,2},',');
ExcludeAccessionNumbers = split(inputs1{gene_num,3},',');
InclusionSequenceFiles = split(inputs1{gene_num,4},',');
ExcludeSequenceFiles = split(inputs1{gene_num,5},',');
if (sum(cellfun(@isempty,IncludeAccessionNumbers))>0)
    IncludeAccessionNumbers = {};
else
    IncludeAccessionNumbers = reshape(IncludeAccessionNumbers,1,[]);
end
if (sum(cellfun(@isempty,ExcludeAccessionNumbers))>0)
    ExcludeAccessionNumbers = {};
else
    ExcludeAccessionNumbers = reshape(ExcludeAccessionNumbers,1,[]);
end
if (sum(cellfun(@isempty,ExcludeSequenceFiles))>0)
    ExcludeSequenceFiles = {};
else
    ExcludeSequenceFiles = reshape(ExcludeSequenceFiles,1,[]);
end
if (sum(cellfun(@isempty,InclusionSequenceFiles))>0)
    InclusionSequenceFiles= {};
else
    InclusionSequenceFiles = reshape(InclusionSequenceFiles,1,[]);
end
RunOffline = 1;%Is Design Run Offline using databases or looks online to get genbank record
UseFlankingInfo = 0;
%% Parse Input Setting Parameters
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
RemoveProbesBindingOffTargetRibosomalHits = inputsParameterSettings.DesignFiltering_Settings.RemoveProbesBindingOffTargetRibosomalHits;% Filter out probes with targets hits to ribosomal proteins
packOptimal_ProbesWithNoOffTargets = inputsParameterSettings.DesignFiltering_Settings.packOptimal_ProbesWithNoOffTargets;%when designing probes without off-targets do optimal packing to get the most or normal selection not considering packing efficiency
IncludeSelfHybridizationInProbeSelection = inputsParameterSettings.DesignFiltering_Settings.IncludeSelfHybridizationInProbeSelection;%when designing probes consider probe self-hybridization when ranking probes on binding affinity

%% Thermodynamic Settings
Gibbs_Model = inputsParameterSettings.Thermodynamic_Settings.Gibbs_Model; %Which Hybridization model to use for probe design and evaluation
SaltConcentration = inputsParameterSettings.Thermodynamic_Settings.SaltConcentration; %Concentration of Salt in thermodynamic calculations mol/L
HybridizationTemperatureCelsius = inputsParameterSettings.Thermodynamic_Settings.HybridizationTemperature_Celsius; % Temperature for Evaluating Probe Design and Simulations
PrimerConcentration = inputsParameterSettings.Thermodynamic_Settings.PrimerConcentration; % Temperature for Evaluating Probe Design and Simulations
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
parse_seqids =  inputsParameterSettings.MAKEBLASTDB_Settings.parse_seqids; %When making blastdb preserve accession names into blastdb sequence ids (default to true, when getting flanking info if keeping mismatches in binding affinity calcs)
hash_index =  inputsParameterSettings.MAKEBLASTDB_Settings.hash_index; %When making blastdb creates sequence hash values for sequence retrieval (default false, when used can lead to fast exact match retrieval but worse range queries where not exact matches).

%% BLAST Parameters (You Usually will not change)
evalue =  inputsParameterSettings.BLASTN_Settings.evalue; %BLAST expectation value cutoff
dust =  convertStringsToChars(inputsParameterSettings.BLASTN_Settings.dust); % filter BLAST query sequences with DUST
gapextend = inputsParameterSettings.BLASTN_Settings.gapextend; %BLAST cost to extend a gap
gapopen = inputsParameterSettings.BLASTN_Settings.gapopen; %BLAST cost to open a gap
num_alignments = inputsParameterSettings.BLASTN_Settings.numAlignments; %Number of database sequneces to show num_alignments for
penalty = inputsParameterSettings.BLASTN_Settings.penalty; % BLAST penalty for nucleotide mismatch
reward = inputsParameterSettings.BLASTN_Settings.reward; % BLAST reward for a nucleotide match
word_size = inputsParameterSettings.BLASTN_Settings.wordsize; %BLAST word size for wordfinder algorithm
task = inputsParameterSettings.BLASTN_Settings.task; %BLAST task to use for evaluating alignments, megablast, blastn, blastn-short, dc-megablast

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
if (isnumeric(inputsParameterSettings.ModelSimulation_Settings.Dilution_Vector))
    Dilution_Vector =  inputsParameterSettings.ModelSimulation_Settings.Dilution_Vector;
else
    Dilution_Vector =  double(split(inputsParameterSettings.ModelSimulation_Settings.Dilution_Vector,','));
end
if (isnumeric(inputsParameterSettings.ModelSimulation_Settings.Temperature_Celsius_Model_Vector))
    Temperature_Celsius_Model_Vector =  inputsParameterSettings.ModelSimulation_Settings.Temperature_Celsius_Model_Vector;
else
    Temperature_Celsius_Model_Vector =  double(split(inputsParameterSettings.ModelSimulation_Settings.Temperature_Celsius_Model_Vector,','));
end
if (isnumeric(inputsParameterSettings.ModelSimulation_Settings.Gibbs_Model_Vector))
    Gibbs_Model_Vector =  inputsParameterSettings.ModelSimulation_Settings.Gibbs_Model_Vector;
else
    Gibbs_Model_Vector =  double(split(inputsParameterSettings.ModelSimulation_Settings.Gibbs_Model_Vector ,','));
end
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
if (isMATLABReleaseOlderThan("R2022b"))
    msg = 'Error. \n MATLAB must be version 2022b or higher for design software to work.';
    error(msg)
end
if ~(ismcc || isdeployed)
    %#exclude matlab.addons.installedAddons
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
end
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    if ~(ismcc || isdeployed)
        %#exclude matlab.desktop.editor.getActiveFilename
        mfilePath = matlab.desktop.editor.getActiveFilename;
    else
        mFilePath = which(fullfile('A0_BKJH_ProbeDesign_Wrapper_cluster_V5.exe'));
    end
end
if (~strcmp(pwd,extractBefore(mfilePath,strcat(filesep,'A0'))))
    msg = 'Error. The script must be run in the TrueProbes main folder for the code to work properly';
    error(msg)
end

startup;
probes = [];
designerName = '_TrueProbes';
settings.BUILD_STRING = '2025.06.04.00';
settings.VERSION_STRING = 'v1.1.1';
refInfo = IncludeAccessionNumbers{1}(1:2);
if (ismember(refInfo,{'NR','XR','NM','XM'}))
    settings.referenceType = 'RefSeq';
else
    settings.referenceType = 'ENSEMBL';
end
%% Update Location of Databases & Needed Files
settings.EMBLtoNCBI = dictionary(struct2table(readstruct(input_databases_file).EMBL_to_NCBI.row).Organism,struct2table(readstruct(input_databases_file).EMBL_to_NCBI.row).StableIDs);
if (strcmp(settings.referenceType,'RefSeq'))
    DatabaseLocations = struct2table(readstruct(input_databases_file).NCBI.row);
elseif (strcmp(settings.referenceType,'ENSEMBL'))
    DatabaseLocations = struct2table(readstruct(input_databases_file).EMBL.row);
end
for k = 2:size(DatabaseLocations,2)
    settings.(strcat("Loc",DatabaseLocations.Properties.VariableNames{k})) = dictionary(DatabaseLocations.Organism,DatabaseLocations.(DatabaseLocations.Properties.VariableNames{k}));
end
if (isKey(settings.LocRoot_FASTA,Organism))
    settings.SEQdbRoot = char(settings.LocRoot_FASTA(Organism));
else
    settings.SEQdbRoot = char(settings.otherLocRoot);
end

%% Load Annotation File
fprintf("Loading genome and transcriptome annotation files")
fprintf('\n')
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
fprintf("Time elapsed to load annotation files %g seconds",round(tEnd,3,"significant"))
fprintf('\n')
fprintf('\n')

%% Gets List of chromosomes Target Transcripts are on using GFF & GTF
if (strcmp(settings.referenceType,'RefSeq'))
    geneInfo_Table = getGenes(GTFobj,"Transcript",IncludeAccessionNumbers);
    geneNames = char(join(convertCharsToStrings(unique(geneInfo_Table.GeneID)),'_'));
    ids_gff_loc = find(strcmp(GFFobj.Feature,'region').*contains(GFFobj.Attributes,'NC'));
    chr_list = regexp(GFFobj.Attributes(ids_gff_loc),'(?<=Name=)\w*(?=;)','match');
    chrom_list = string(convertCharsToStrings(chr_list));
    ref_list = regexp(GFFobj.Attributes(ids_gff_loc),'(?<=ID=)\w*(?=.)','match');
    reference_list = string(convertCharsToStrings(ref_list));
    ChromosomeToReference = dictionary(reference_list,chrom_list);
    ReferenceToChromosome = dictionary(reference_list,chrom_list);
    geneChrNum = ReferenceToChromosome(extractBefore(string(geneInfo_Table.Reference),'.'));
    geneReference_ID = extractBefore(string(geneInfo_Table.Reference),'.');
elseif (strcmp(settings.referenceType,'ENSEMBL'))
    if (sum(contains(IncludeAccessionNumbers,'ENS'))>0)
    if (sum(ismember(extractBefore(IncludeAccessionNumbers,'.'),getTranscriptNames(GTFobj)))>0)
        geneInfo_Table = getGenes(GTFobj,"Transcript",extractBefore(IncludeAccessionNumbers,'.'));
        ts = 1;
    else
        geneInfo_Table = getGenes(GTFobj,"Transcript",IncludeAccessionNumbers);
        ts = 0;
    end
    else
    try
        geneInfo_Table = getGenes(GTFobj,"Transcript",IncludeAccessionNumbers);
        ts = 0;
    catch
        geneInfo_Table = getGenes(GTFobj,"Transcript",extractBefore(IncludeAccessionNumbers,'.'));
        ts = 1;
    end
    end
    geneNames = char(join(convertCharsToStrings(unique(geneInfo_Table.GeneName)),'_'));
    ids_gff_loc = find(contains(GFFobj.Attributes,'Alias').*contains(GFFobj.Attributes,'NC')); 
    if (isempty(ids_gff_loc))
    ids_gff_loc = find(contains(GFFobj.Attributes,'Alias'));
    end
    chr_list = regexp(GFFobj.Attributes(ids_gff_loc),'(?<=:).*(?=[,;]Alias)','match');
    chrom_list = string(convertCharsToStrings(chr_list));
    geneChrNum = string(geneInfo_Table.Reference);
    try
    NCBI_chrom_ID = string(convertCharsToStrings(regexp(GFFobj.Attributes(ids_gff_loc),'\<N\w*','match')));
    ref_list = regexp(GFFobj.Attributes(ids_gff_loc),'(?<=Alias=).*(?=[,;])','match');
    reference_list = string(convertCharsToStrings(ref_list));
    refparts = split(reference_list,',');
    paired_ID_chrom = sort(refparts,2);
    ReferenceToChromosome = dictionary(paired_ID_chrom(:,1),chrom_list);
    ChromosomeToReference = dictionary(chrom_list,paired_ID_chrom(:,1));
    geneReference_ID = extractBefore(ChromosomeToReference(string(geneInfo_Table.Reference)),'.');
    catch
    end
end
%% Gets List of all isoforms of Target Transcripts using GTF
if (strcmp(settings.referenceType,'RefSeq'))
    transcriptInfo_Table = getTranscripts(GTFobj,"Transcript",IncludeAccessionNumbers);
    all_isoform_transcriptInfo_Table= getTranscripts(GTFobj,"Gene",getGenes(GTFobj,"Transcript",IncludeAccessionNumbers).GeneID);
    Transcript_All_Isoform_IDs = cellfun(@(x) all_isoform_transcriptInfo_Table.Transcript(strcmp(all_isoform_transcriptInfo_Table.GeneID,x)), transcriptInfo_Table.GeneID,'Un',0);
    Transcript_UnDesired_Isoform_IDs = arrayfun(@(x) setdiff(Transcript_All_Isoform_IDs{x},IncludeAccessionNumbers{x}),1:length(IncludeAccessionNumbers),'Un',0);
    Transcript_Joint_UnDesired_Isoforms_IDs = setdiff(all_isoform_transcriptInfo_Table.Transcript,transcriptInfo_Table.Transcript);
    exonInfo_Table = getExons(GTFobj,"Transcript",IncludeAccessionNumbers);
    segmentInfo_Table = getSegments(GTFobj,"Transcript",IncludeAccessionNumbers);
    all_isoform_exonInfo_Table = getExons(GTFobj,"Gene",getGenes(GTFobj,"Transcript",IncludeAccessionNumbers).GeneID);
    all_isoform_segmentInfo_Table = getSegments(GTFobj,"Gene",getGenes(GTFobj,"Transcript",IncludeAccessionNumbers).GeneID);
elseif (strcmp(settings.referenceType,'ENSEMBL'))
    %if (sum(ismember(extractBefore(IncludeAccessionNumbers,'.'),getTranscriptNames(GTFobj)))>0)
    if (ts==1)
    transcriptInfo_Table = getTranscripts(GTFobj,"Transcript",extractBefore(IncludeAccessionNumbers,'.'));
    all_isoform_transcriptInfo_Table= getTranscripts(GTFobj,"Gene",getGenes(GTFobj,"Transcript",extractBefore(IncludeAccessionNumbers,'.')).GeneID);
    Transcript_All_Isoform_IDs = cellfun(@(x) all_isoform_transcriptInfo_Table.Transcript(strcmp(all_isoform_transcriptInfo_Table.GeneID,x)), transcriptInfo_Table.GeneID,'Un',0);
    Transcript_UnDesired_Isoform_IDs = arrayfun(@(x) setdiff(Transcript_All_Isoform_IDs{x},extractBefore(IncludeAccessionNumbers{x},'.')),1:length(IncludeAccessionNumbers),'Un',0);
    Transcript_Joint_UnDesired_Isoforms_IDs = setdiff(all_isoform_transcriptInfo_Table.Transcript,transcriptInfo_Table.Transcript);
    exonInfo_Table = getExons(GTFobj,"Transcript",extractBefore(IncludeAccessionNumbers,'.'));
    segmentInfo_Table = getSegments(GTFobj,"Transcript",extractBefore(IncludeAccessionNumbers,'.'));
    all_isoform_exonInfo_Table = getExons(GTFobj,"Gene",getGenes(GTFobj,"Transcript",extractBefore(IncludeAccessionNumbers,'.')).GeneID);
    all_isoform_segmentInfo_Table = getSegments(GTFobj,"Gene",getGenes(GTFobj,"Transcript",extractBefore(IncludeAccessionNumbers,'.')).GeneID);
    else
    transcriptInfo_Table = getTranscripts(GTFobj,"Transcript",IncludeAccessionNumbers);
    all_isoform_transcriptInfo_Table= getTranscripts(GTFobj,"Gene",getGenes(GTFobj,"Transcript",IncludeAccessionNumbers).GeneID);
    Transcript_All_Isoform_IDs = cellfun(@(x) all_isoform_transcriptInfo_Table.Transcript(strcmp(all_isoform_transcriptInfo_Table.GeneID,x)), transcriptInfo_Table.GeneID,'Un',0);
    Transcript_UnDesired_Isoform_IDs = arrayfun(@(x) setdiff(Transcript_All_Isoform_IDs{x},IncludeAccessionNumbers{x}),1:length(IncludeAccessionNumbers),'Un',0);
    Transcript_Joint_UnDesired_Isoforms_IDs = setdiff(all_isoform_transcriptInfo_Table.Transcript,transcriptInfo_Table.Transcript);     
    exonInfo_Table = getExons(GTFobj,"Transcript",IncludeAccessionNumbers);
    segmentInfo_Table = getSegments(GTFobj,"Transcript",IncludeAccessionNumbers);
    all_isoform_exonInfo_Table = getExons(GTFobj,"Gene",getGenes(GTFobj,"Transcript",IncludeAccessionNumbers).GeneID);
    all_isoform_segmentInfo_Table = getSegments(GTFobj,"Gene",getGenes(GTFobj,"Transcript",IncludeAccessionNumbers).GeneID);
    end
end
%% Gets Ribosomal Transcripts using GTF
ids_gtf = GTFobj.Transcript(strcmp(GTFobj.Feature,'transcript'));
btypes = extractBetween(string(regexp(GTFobj.Attributes(strcmp(GTFobj.Feature,'transcript')),'(?<=transcript_biotype )"\w*"(?=;)','match')),'"','"');
Ribosomal_IDs = unique(ids_gtf(contains(btypes,'rRNA')));
%% Using GTF Gets Transcript IDs and Corresponding Gene Names For Inputing Transcript Level Expression From Gene Level Expression Data
allGTF_Transcript_IDs = getTranscriptNames(GTFobj);
allGTF_Transcript_IDs = allGTF_Transcript_IDs(~cellfun(@isempty,allGTF_Transcript_IDs));

allGTF_Transcript_uniNames = extractBefore(allGTF_Transcript_IDs,'.');
if (sum(ismissing(allGTF_Transcript_uniNames))>0)
    allGTF_Transcript_uniNames(ismissing(allGTF_Transcript_uniNames)) = allGTF_Transcript_IDs(ismissing(allGTF_Transcript_uniNames));
end
allGTF_GeneIDs_IDstring = convertCharsToStrings(getTranscripts(GTFobj,"Transcript",allGTF_Transcript_IDs).GeneID);
allGTF_GeneNames_IDstring = convertCharsToStrings(getTranscripts(GTFobj,"Transcript",allGTF_Transcript_IDs).GeneName);
settings.pairedGTF_TranscriptIDs = convertCharsToStrings(allGTF_Transcript_uniNames(~cellfun(@isempty,allGTF_Transcript_uniNames)));
settings.pairedGTF_GeneIDs = allGTF_GeneIDs_IDstring(~cellfun(@isempty,allGTF_Transcript_uniNames));
settings.pairedGTF_GeneNames_EMBLonly =  allGTF_GeneNames_IDstring(~cellfun(@isempty,allGTF_Transcript_uniNames));
clear GTFobj GFFobj

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
settings.chromosome_IDs = {char(geneInfo_Table.Reference)};
settings.transcript_IDs = IncludeAccessionNumbers;
settings.ribosomal_IDs = Ribosomal_IDs;
settings.transcript_IDs_desired = Transcript_All_Isoform_IDs;
settings.transcript_IDs_undesired = Transcript_UnDesired_Isoform_IDs;
settings.transcript_IDs_joint_undesired = Transcript_Joint_UnDesired_Isoforms_IDs;
settings.ProbeSpacing = minProbeSpacing;
settings.RemoveMisMatches = RemoveMisMatches;
settings.SaltConcentration = SaltConcentration;
settings.HybridizationTemperature = HybridizationTemperatureCelsius;
settings.PrimerConcentration = PrimerConcentration;
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
settings.BlastParameters.task = task;
settings.MakeblastdbParameters.parse_seqids = parse_seqids;
settings.MakeblastdbParameters.hash_index = hash_index;
settings.TargetBatchSize = Parallelization_targetBatchSize;
settings.MinProbeSize = minProbeSize;
settings.MaxProbeSize = maxProbeSize;
settings.MinHomologySearchTargetSize = MinHomologySearchTargetSize;
settings.N_model = Gibbs_Model;
settings.UseFlankingInfo = UseFlankingInfo;
settings.ExpressionReferenceForDesigningProbes = ExpressionReferenceForProbeDesign;


%% Gene Expression Parameters
settings.DoAllGenesHaveSameExpression = DoAllGenesHaveSameExpression;
settings.HumanSpecific.HumanExpGeneOrTransc = UseGeneOverTranscLevelExpression; % 1 (Gene/EMBL GENEID) , 0 (Transcript/EMBL Transcript ID)
settings.UseRegularDNAExpression = UseRegularDNAExpression;%0 use DNA expression from gene expression track in expression data, 1 set expression to 2 for DNA.
settings.UniformRNAExpression = nullRNAcopynumber;%if assuming no differences in gene's expression sets level.
settings.DNAPloidy = nullDNAcopynumber;
settings.TMM.logRatioTrim = logRatioTrim;
settings.TMM.sumTrim = sumTrim;
settings.TMM.Acutoff = Acutoff;
settings.TMM.doWeighting = doWeighting;


%% Cluster Parameters
settings.clusterStatus = cluster;
settings.isOffline = RunOffline;

%% Selecting Probes Specifications
settings.maxProbes = maxNumberOfProbes;
settings.RemoveProbesBindingOffTargetRibosomalHits = RemoveProbesBindingOffTargetRibosomalHits;

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
fprintf("Making probe target design output folder")
fprintf('\n')
fprintf('\n')
if (not(isfolder([saveRoot '(' geneNames ')' '_' settings.rootName])))
    mkdir([saveRoot  '(' geneNames ')' '_' settings.rootName])
end
settings.FolderName = [ '(' geneNames ')' '_' settings.rootName];
settings.FolderRootName = strcat('output',filesep,settings.FolderName);
settings.rootName = strjoin(IncludeAccessionNumbers,'_');
settings.TargetLists = strjoin(IncludeAccessionNumbers,', ');

%% Generate Probes
fprintf(strcat("Designing probes for"," ",strrep(geneNames,'_',', ')," ","transcript id:"," ",settings.TargetLists))
fprintf('\n')
fprintf('\n')
try
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_probes' designerName '.mat'],'probes')
    fprintf("Loading probe tile sequences")
    fprintf('\n')
    fprintf('\n')
catch
    fprintf("Generating probe tile sequences")
    fprintf('\n')
    fprintf('\n')
    tic
    [init_probes,~,~,~,~] = ...
        BKJH_Probe_Generator([Lmin Lmax],IncludeAccessionNumbers,InclusionSequenceFiles,ExcludeSequenceFiles,ExcludeAccessionNumbers,targetStrand,settings.isOffline,settings.SEQdbRoot);
    probes = init_probes;
    save([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_probes' designerName '.mat'],'probes','-v7.3')
    tEnd = toc;
    fprintf("Time elapsed to tile probes %g seconds",round(tEnd,3,"significant"))
    fprintf('\n')
    fprintf('\n')
end

%% BLAST Probes
if (cluster==0)
    if (ismac)
        blastpath = strcat(filesep,'usr',filesep,'local',filesep,'ncbi',filesep,'blast',filesep,'bin');
        addpath(blastpath);
        setenv('PATH',[blastpath ':' getenv('PATH')]);
    else
        curr_dir = pwd;
        root_drive = extractBefore(curr_dir,'Users');
        cd(root_drive);
        if (ispc)
            blast_paths = dir('**/blastn.exe');
        else
            blast_paths = dir('**/blastn');
        end
         valid_path = find(contains(lower({blast_paths.folder}),strcat(filesep,'ncbi',filesep)),1);
         if (isempty(valid_path))
             msg = 'BLAST+ needs to be installed in order for software to work. (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)';
             error(msg)
         else
             blastpath = blast_paths(valid_path).folder;
             addpath(blastpath);
             if (ispc)
             setenv('PATH',[blastpath ';' getenv('PATH')]);
             end
             if (isunix)
             setenv('PATH',[blastpath ':' getenv('PATH')]);
             end
         end
         cd(curr_dir);
    end
end
[status,cmdout] = system('blastn -version');
if status
    msg = 'BLAST+ needs to be installed in order for software to work. (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)';
    error(msg)
else
    fprintf("Version of blastn to be used in probe design:")
    fprintf('\n')
    fprintf(cmdout)
    fprintf('\n')
    fprintf('\n')
end
if (ismac)
    [status,cmdout] = system('which blastn');
    if status
        msg = 'BLAST+ needs to be installed in order for software to work. (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)';
        error(msg)
    else
        fprintf("Location of blastn to be used in probe design:")
        fprintf('\n')
        fprintf(cmdout)
        fprintf('\n')
        fprintf('\n')
    end
end
try
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_hits_table' designerName '.mat'],'gene_table')
    fprintf("Loading probe BLAST results table")
    fprintf('\n')
    fprintf('\n')
catch
    fprintf("BLASTING Probes in batches")
    fprintf('\n')
    fprintf('\n')
    tic
    [~,gene_table] = Probe_checker_general_JH10(probes,Organism,geneNames,geneNames,geneChrNum,settings);
    save([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_hits_table' designerName '.mat'],'-mat','gene_table','-v7.3');
    tEnd = toc;
    fprintf("Time elapsed to generate probe BLAST results table %g seconds",round(tEnd,3,"significant"))
    fprintf('\n')
    fprintf('\n')
end
fprintf("Getting Parser Information from BLAST Databases BLAST")
fprintf('\n')
fprintf('\n')
[RNAdbParser, DNAdbParser] = A3_BlastDBCMD_JH(settings,gene_table);
settings.DNAdbParser = DNAdbParser;
settings.RNAdbParser = RNAdbParser;


%% Get Gene Expression Information
try
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_ExpressionInfo' designerName '.mat'],'ExpressionMatrix');
    fprintf("Loading BLAST hits gene expression information")
    fprintf('\n')
    fprintf('\n')
catch
    fprintf("Getting BLAST hits gene expression information")
    fprintf('\n')
    fprintf('\n')
    tic
    [ExpressionMatrix,get_expression_time] = A_JH_GetExpressionInfo_V3(gene_table,settings,input_gene_expression_file_locations);
    save([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_ExpressionInfo' designerName '.mat'],'ExpressionMatrix','get_expression_time','settings','-v7.3');
    tEnd = toc;fprintf('\n')
    fprintf("Time elapsed to generate probe BLAST hits gene expression information %g seconds",round(tEnd,3,"significant"))
    fprintf('\n')
    fprintf('\n')
end
settings.saveRoot = saveRoot;

%% Get Thermodynamic Information (On-Target,Off-Target)
try
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_Tm' num2str(settings.HybridizationTemperature) '_OnOffThermoInfo' designerName '.mat'],'Kon','Kb_Match');
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_dHInfo' designerName '.mat'],'-mat','dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match');
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_dSInfo' designerName '.mat'],'-mat','dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match');
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_TmInfo' designerName '.mat'],'-mat','Tm_on','Tm_Match');
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_dCpInfo' designerName '.mat'],'-mat','dCpon_eq','dCpeq_Match');
    fprintf("Loading probe target thermodynamic information")
    fprintf('\n')
    fprintf('\n')
catch
    fprintf("Computing probe target thermodynamic information")
    fprintf('\n')
    fprintf('\n')
    tic
    [Kb_Match,Kon,Koff,dHeq_Match,dSeq_Match,dHf_Match,dSf_Match,dHr_Match,dSr_Match,dCpeq_Match,dHon_eq,dSon_eq,dHon_f,dSon_f,dHon_r,dSon_r,dCpon_eq,Tm_on,Tm_Match] = ...
        A_JH_GenerateThermoInfo_V5(probes,gene_table,geneNames,settings);%add Kon Koff
    save([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_Tm' num2str(settings.HybridizationTemperature) '_OnOffThermoInfo' designerName '.mat'],'-mat','Kon','Koff','Kb_Match','-v7.3');
    save([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_dHInfo' designerName '.mat'],'-mat','dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match','-v7.3');
    save([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_dSInfo' designerName '.mat'],'-mat','dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match','-v7.3');
    save([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_TmInfo' designerName '.mat'],'-mat','Tm_on','Tm_Match','-v7.3');
    save([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_dCpInfo' designerName '.mat'],'-mat','dCpon_eq','dCpeq_Match','-v7.3');
    toc
    tEnd = toc;
    fprintf("Time elapsed to compute probe target thermodynamic information %g seconds",round(tEnd,3,"significant"))
    fprintf('\n')
    fprintf('\n')
end

%% Get Binding Site Mapping and Energy
try
    load([settings.FolderRootName filesep '(' geneNames ')_binding_hits_map' designerName '.mat'],'DoesProbeBindSite2')
    load([settings.FolderRootName filesep '(' geneNames ')_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' designerName '.mat'],'Kb_Complement')
    load([settings.FolderRootName filesep '(' geneNames ')_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' designerName '.mat'],'Kb_mod')
    load([settings.FolderRootName filesep '(' geneNames ')_BindingMatrices' designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','dCp_mod')
    if (settings.BLASTdna)
        load([settings.FolderRootName filesep '(' geneNames ')_BindingMatrices2' designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','dCp_Complement')
    end
    fprintf("Loading probe target binding site maps")
    fprintf('\n')
    fprintf('\n')
catch
    fprintf("Generating probe target binding site maps")
    fprintf('\n')
    fprintf('\n')
    tic
    [Kb_mod,Kb_Complement,DoesProbeBindSite2,~,~,~,~,...
        dHeq_mod,dSeq_mod,dHf_mod,dSf_mod,dHr_mod,dSr_mod,~,dCp_mod,...
        dHeq_Complement,dSeq_Complement,dHf_Complement,dSf_Complement,dHr_Complement,dSr_Complement,dCp_Complement] = ...
        A_JH_GetSiteMapping_V6(probes,settings,gene_table,Kb_Match,dHeq_Match,dSeq_Match,dHf_Match,dSf_Match,dHr_Match,dSr_Match,dCpeq_Match,Tm_Match);
    tEnd = toc;
    fprintf("Time elapsed to create probe target binding site maps %g seconds",round(tEnd,3,"significant"))
    fprintf('\n')
    fprintf('\n')
end
clear Kb_Match dCpeq_Match Tm_Match
clear dHeq_Match dHf_Match dHr_Match 
clear dSeq_Match dSf_Match dSr_Match 
%% Basic Stats for Designing Probes Function
if (ExpressionReferenceForProbeDesign==0)
    EKernel = ones(size(ExpressionMatrix,1),1);
else
    EKernel = ExpressionMatrix(:,ExpressionReferenceForProbeDesign);
end
Kon = squeeze(Kon(:,Gibbs_Model));
Kb_mod = squeeze(Kb_mod(:,:,:,Gibbs_Model));
DoesProbeBindSite = DoesProbeBindSite2;
if (settings.BLASTdna)
    Kb_Complement = squeeze(Kb_Complement(:,:,Gibbs_Model));
else
    dHeq_Complement= [];
    dCp_Complement= [];
    dSeq_Complement= [];
    Kb_Complement = [];
end
try
    load([settings.FolderRootName filesep '(' geneNames ')_' settings.rootName '_Tm' num2str(T_hybrid) '_BasicDesignerStats' designerName '.mat'],'Tvec_RNA','Svec_RNA','TPvec_RNA','TSvec_RNA','TPvec_logKOFF_RNA','TPvec_logKOFFdivON_RNA','TPvec_logKONdivOFF_RNA',...
        'Nvec_RNAmulti','Tvec_DNA','Svec_DNA','TPvec_DNA','TSvec_DNA','TPvec_logKOFF_DNA','TPvec_logKOFFdivON_DNA','TPvec_logKONdivOFF_DNA','TPvec_logKOFFdivCOMP_DNA','TPvec_logKCOMPdivOFF_DNA',...
        'Off_Score','Specificity_Score','NumRNAOffTargetOptions','Probes_WithNRNAOFF','NumDNAOffTargetOptions','Probes_WithNDNAOFF')
    fprintf("Loading probe target statistics information")
    fprintf('\n')
    fprintf('\n')
catch
    fprintf("Computing probe target statistics information")
    fprintf('\n')
    fprintf('\n')
    tic
    [Nvec_RNAmulti,RNAOFF_Score,RNASpecificity_Score,NumRNAOffTargetOptions,Probes_WithNRNAOFF,DNAOFF_Score,DNASpecificity_Score,NumDNAOffTargetOptions,Probes_WithNDNAOFF,Cout] = ...
        A0_BasicDesignerStats_V2([settings.BLASTrna settings.BLASTdna],removeUndesiredIsoformsFromPredictionOffTargets,gene_table,settings,settings.FolderRootName,DoesProbeBindSite,Kon,Kb_mod,Kb_Complement,EKernel);
    Tvec_RNA = Cout{1}{1};Svec_RNA = Cout{1}{2};TPvec_RNA = Cout{1}{3};TSvec_RNA = Cout{1}{4};
    TPvec_logKOFF_RNA = Cout{1}{5};TPvec_logKOFFdivON_RNA = Cout{1}{6};TPvec_logKONdivOFF_RNA = Cout{1}{7};
    Tvec_DNA = Cout{2}{1};Svec_DNA = Cout{2}{2};TPvec_DNA = Cout{2}{3};TSvec_DNA = Cout{2}{4};
    TPvec_logKOFF_DNA = Cout{2}{5};TPvec_logKOFFdivON_DNA = Cout{2}{6};TPvec_logKONdivOFF_DNA = Cout{2}{7};
    TPvec_logKOFFdivCOMP_DNA = Cout{2}{8};TPvec_logKCOMPdivOFF_DNA = Cout{2}{9};
    Off_Score = RNAOFF_Score+DNAOFF_Score;
    Specificity_Score = RNASpecificity_Score+DNASpecificity_Score;
    tEnd = toc;
    fprintf("Time elapsed to compute probe target statistics %g seconds",round(tEnd,3,"significant"))
    fprintf('\n')
    fprintf('\n')
end

%% Selection of TrueProbes Probes
try
    load([saveRoot filesep settings.FolderName filesep settings.FolderName '_chosen.mat'],'chosenProbes')
    fprintf("Loading TrueProbes designed probes")
    fprintf('\n')
    fprintf('\n')
catch
    fprintf("Designing TrueProbes probe set")
    fprintf('\n')
    fprintf('\n')
    tic
    chosenProbes = A_ZigZagProbeSelection_V5(probes,gene_table,settings,IncludeSelfHybridizationInProbeSelection,packOptimal_ProbesWithNoOffTargets,Kon,Nvec_RNAmulti,Off_Score,Specificity_Score,Tvec_RNA,Svec_RNA,TPvec_RNA,TSvec_RNA,TPvec_logKOFF_RNA,TPvec_logKOFFdivON_RNA,TPvec_logKONdivOFF_RNA,ExpressionMatrix,DoesProbeBindSite,Kb_mod);
    save([saveRoot filesep settings.FolderName filesep settings.FolderName '_chosen.mat'],'chosenProbes','-v7.3')
    tEnd = toc;
    fprintf("Time elapsed to select TrueProbes probes %g seconds",round(tEnd,3,"significant"))
    fprintf('\n')
    fprintf('\n')
end
gene_table0 = gene_table;
gene_table = sortrows(gene_table,[7 6],'ascend');
gene_table = gene_table(gene_table.Match>=settings.MinHomologySearchTargetSize,:);
MinusStrandedHits = find(contains(gene_table.Strand,'Minus'));
gene_table_NamesZ = convertCharsToStrings(gene_table.Name);
contains_RNA = find(ismember(gene_table_NamesZ,settings.RNAdbParser));
RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);clear contains_RNA
gene_table = gene_table(setdiff(1:size(gene_table,1),RNA_MissedFilteredHits),:);
gene_table.Ax = min(gene_table.SubjectIndices,[],2);
gene_table.Bx = max(gene_table.SubjectIndices,[],2);
gene_table = sortrows(gene_table,[7 13],'ascend');
Names = unique(gene_table.Name);
Names = convertCharsToStrings(Names);
if (and(strcmp(settings.referenceType,"ENSEMBL"),max(double(contains(extractBefore(Names,' '),'ENS')))==0))
    uniNames = extractBefore(Names,' ');
else
    uniNames = extractBefore(Names,'.');
    if (sum(ismissing(uniNames))>0)
        uniNames(ismissing(uniNames)) = extractBefore(Names(ismissing(uniNames)),' ');
    end
end
if (settings.BLASTdna)
DNA_IDs = find(~ismember(Names,settings.DNAdbParser));%IDs
else
DNA_IDs = [];
end
if (settings.BLASTrna)
NonDNA_IDs = find(ismember(Names,settings.RNAdbParser));%IDs
else
NonDNA_IDs =[];
end
ON_IDs_specific = find(ismember(uniNames,extractBefore(settings.transcript_IDs{:},'.')));
ON_IDs_agnostic = find(ismember(uniNames,extractBefore(settings.transcript_IDs_desired{:},'.')));
OFF_IDs = find(~ismember(uniNames,extractBefore(settings.transcript_IDs_desired{:},'.')));

%% Print Excel Spreedsheet of Probes
fprintf("Printing TrueProbes probes to Excel spreadsheet")
fprintf('\n')
fprintf('\n')

Lp_designed_min = min(cell2mat(cellfun(@length,probes(chosenProbes,2),'UniformOutput',false)));
Lp_designed_max = max(cell2mat(cellfun(@length,probes(chosenProbes,2),'UniformOutput',false)));
TL_hits = zeros(length(min(settings.MinHomologySearchTargetSize,Lp_designed_min):Lp_designed_max),size(probes,1));
for L = min(settings.MinHomologySearchTargetSize,Lp_designed_min):Lp_designed_max
    TL_hits(L+1-settings.MinHomologySearchTargetSize,:) = cell2mat(arrayfun(@(z) length(find((gene_table.ProbeNum==z).*(gene_table.Match == L)>0)),1:size(probes,1),'Un',0));
end
final_probe_info = cell(length(chosenProbes),13+length(min(settings.MinHomologySearchTargetSize,Lp_designed_min):Lp_designed_max));
for v = 1:length(chosenProbes)
    Lp = length(probes{chosenProbes(v),2});
    final_probe_info{v,1} = settings.transcript_IDs;
    final_probe_info{v,2} = settings.GeneName;%needs to be different for using multiple genes
    final_probe_info{v,3} = probes{chosenProbes(v),3};
    final_probe_info{v,4} = probes{chosenProbes(v),3}+Lp-1;
    final_probe_info{v,5} = probes{chosenProbes(v),1};
    [rev_comp] = Reverse_complement_probes({probes{chosenProbes(v),2}});  %find reverse complement
    final_probe_info{v,6} = upper(rev_comp{1});
    final_probe_info{v,7} = length(probes{chosenProbes(v),2});
    final_probe_info{v,8} = oligoprop(rev_comp{1},'Salt',SaltConcentration,'Temp',T_hybrid).GC;
    final_probe_info{v,9} = round(oligoprop(rev_comp{1},'Salt',SaltConcentration,'Temp',T_hybrid).Tm(5),2);
    final_probe_info{v,10} = sum(full(squeeze(DoesProbeBindSite2(chosenProbes(v),ON_IDs_agnostic,:))),'all');
    final_probe_info{v,11} = sum(full(squeeze(DoesProbeBindSite2(chosenProbes(v),ON_IDs_specific,:))),'all');
    final_probe_info{v,12} = round(Off_Score(chosenProbes(v)),2);
    final_probe_info{v,13} =  sum(TL_hits(:,chosenProbes(v)))-sum(full(squeeze(DoesProbeBindSite2(chosenProbes(v),ON_IDs_agnostic,:))),'all');
    for Li = 1:length(min(settings.MinHomologySearchTargetSize,Lp_designed_min):Lp_designed_max)
    final_probe_info{v,13+Li} = TL_hits(end-Li+1,chosenProbes(v))-double(Li==1)*sum(full(squeeze(DoesProbeBindSite2(chosenProbes(v),ON_IDs_agnostic,:))),'all');
    end
end
OffTargetMatch_Labels = arrayfun(@(n) char(strcat('Number of '," ",string(n),'nt Off-Target Matches')),Lp_designed_max:-1:min(settings.MinHomologySearchTargetSize,Lp_designed_min),'Un',0);
T = cell2table(final_probe_info,...
    'VariableNames',[{'Probe Target'} {'Transcript Name'} {'Base_Position-start'} {'Base_Position-stop'} {'Probe ID'} {'Probe Sequence'} {'Probe Length'} {'GC_fraction (%)'} {'Probe Tm (°C)'} ...
    {'Number of Isoform-Agnostic On-Target Matches'} {'Number of Isoform-Specific On-Target Matches'} ...
    {'Off-Score'} {'Total Off-Target Matches'} OffTargetMatch_Labels(:)']);
filename = [settings.saveRoot filesep settings.FolderName filesep settings.FolderName '_probes_final_' num2str(settings.maxProbes) 'max.xlsx'];
if exist(filename,'file')        %delete if already exists
    delete(filename)
end
writetable(T,filename,'Sheet',1,'Range','A1')
fprintf("Number of Final Designed Probes: ")
fprintf(num2str(size(final_probe_info,1)))
fprintf('\n')
fprintf('\n')
gene_table = gene_table0;clear gene_table0
%% Get Metric Information (Probes, Final Probe Set)
%Using Concentrations Solve For Equilibrium and Get Distributions and probe set metrics for detection
try
    load([settings.FolderRootName filesep '(' geneNames ')_Tm' num2str(T_hybrid) '_ModelMetrics' designerName '.mat'],'ModelMetrics')
    fprintf("Loading TrueProbes probe set kinetic-model simulation results and RNA-FISH probe design metrics")
    fprintf('\n')
    fprintf('\n')
catch
    fprintf("Computing TrueProbes probe set kinetic-model simulation results and RNA-FISH probe design metrics")
    fprintf('\n')
    fprintf('\n')
    tic
    ModelMetrics = ...
        RNAsolver_JH2(chosenProbes,settings,probes,gene_table,ExpressionMatrix,DoesProbeBindSite2,dHeq_mod,dSeq_mod,dCp_mod,dHeq_Complement,dSeq_Complement,dCp_Complement);
    save([settings.FolderRootName filesep '(' geneNames ')_Tm' num2str(T_hybrid) '_ModelMetrics' designerName '.mat'],'ModelMetrics','chosenProbes','settings','-v7.3')
    tEnd = toc;
    fprintf("Time elapsed to compute TrueProbes probe metrics %g seconds",round(tEnd,3,"significant"))
    fprintf('\n')
    fprintf('\n')
end
%% Get Output File with summary of analysis on designed probes
%% Filter out probes below specific specficity

%% Make/Save Output Figures & Figure Objects
%Output plots and excel or csv table with statistics

end
