function A_IntegrateResults1(id)
input_file = 'TrueProbes_DesignTargets.csv';
input_parameters = 'TrueProbes_ParameterSettings.xml';
input_databases_file = 'DatabaseLocations.xml';
input_gene_expression_file_locations = 'GeneExpressionDataLocations.xml';
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
if ~(ismcc || isdeployed)
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
end

%% Settings Specification
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
HybridizationTemperature = inputsParameterSettings.Thermodynamic_Settings.HybridizationTemperature_Celsius; % Temperature for Evaluating Probe Design and Simulations
if (isempty(Organism))
    msg = 'Error. There must be a organism specified with the design in order to quantify target hits in the reference genome or transcriptome.';
    error(msg)
end
if (isempty(IncludeAccessionNumbers)*isempty(InclusionSequenceFiles)==1)
    msg = 'Error. Must include at least one target Accession Number or target inclusion sequence file.';
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
        if (ismac)
            mfilePath = which(fullfile('A0_BKJH_ProbeDesign_Wrapper_cluster_V5.app'));
        elseif (isunix)
            mfilePath = which(fullfile('A0_BKJH_ProbeDesign_Wrapper_cluster_V5'));
        elseif (ispc)
            mfilePath = which(fullfile('A0_BKJH_ProbeDesign_Wrapper_cluster_V5.exe'));
        end
    end
end
if (~strcmp(pwd,extractBefore(mfilePath,strcat(filesep,'A0'))))
    msg = 'Error. The script must be run in the TrueProbes main folder for the code to work properly';
    error(msg)
end
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
elseif (strcmp(settings.referenceType,'ENSEMBL'))
    if (sum(contains(IncludeAccessionNumbers,'ENS'))>0)
    if (sum(ismember(extractBefore(IncludeAccessionNumbers,'.'),getTranscriptNames(GTFobj)))>0)
        geneInfo_Table = getGenes(GTFobj,"Transcript",extractBefore(IncludeAccessionNumbers,'.'));
    else
        geneInfo_Table = getGenes(GTFobj,"Transcript",IncludeAccessionNumbers);
    end
    else
    try
        geneInfo_Table = getGenes(GTFobj,"Transcript",IncludeAccessionNumbers);
    catch
        geneInfo_Table = getGenes(GTFobj,"Transcript",extractBefore(IncludeAccessionNumbers,'.'));
    end
    end
    geneNames = char(join(convertCharsToStrings(unique(geneInfo_Table.GeneName)),'_'));
end
clear GTFobj
%% Save Settings
Nsoftware = 11;
designerName0 = '_NLPDS';
%% Compile Results from  Supporting Output Files
settings.FolderRootName = strcat(saveRoot,'(',geneNames,')','_',strjoin(IncludeAccessionNumbers,'_'));
settings.rootName = strjoin(IncludeAccessionNumbers,'_');
ProbeDesignResults1.Gene = geneNames;
ProbeDesignResults1.GeneTarget = settings.rootName;
ProbeDesignResults1.Organism = Organism;
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
        gene_table = gene_table(gene_table.Match>=15,:);
        MinusStrandedHits = find(contains(gene_table.Strand,'Minus'));clear gene_table
        RNA_IDs_1 = find(contains(gene_table.Name,'NM_'));
        RNA_IDs_2 = find(contains(gene_table.Name,'NR_'));
        RNA_IDs_3 = find(contains(gene_table.Name,'XM_'));
        RNA_IDs_4 = find(contains(gene_table.Name,'XR_'));
        contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
        RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);clear RNA_IDs_*
        gene_table = gene_table(setdiff(1:size(gene_table,1),RNA_MissedFilteredHits),:);
        gene_table.Ax = min(gene_table.SubjectIndices,[],2);clear MinusStrandedHits
        gene_table.Bx = max(gene_table.SubjectIndices,[],2);clear RNA_MissedFilteredHits
        gene_table3 = sortrows(gene_table,[7 13],'ascend');
        Iz = find(endsWith(gene_table3.Name,', ribosomal RNA'));clear gene_table
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
        gene_table = gene_table(gene_table.Match>=15,:);
        MinusStrandedHits = find(contains(gene_table.Strand,'Minus'));clear gene_table
        RNA_IDs_1 = find(contains(gene_table.Name,'NM_'));
        RNA_IDs_2 = find(contains(gene_table.Name,'NR_'));
        RNA_IDs_3 = find(contains(gene_table.Name,'XM_'));
        RNA_IDs_4 = find(contains(gene_table.Name,'XR_'));
        contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
        RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);clear RNA_IDs_*
        gene_table = gene_table(setdiff(1:size(gene_table,1),RNA_MissedFilteredHits),:);
        gene_table.Ax = min(gene_table.SubjectIndices,[],2);clear MinusStrandedHits
        gene_table.Bx = max(gene_table.SubjectIndices,[],2);clear RNA_MissedFilteredHits
        gene_table3 = sortrows(gene_table,[7 13],'ascend');
        Iz = find(endsWith(gene_table3.Name,', ribosomal RNA'));clear gene_table
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


