function A_IntegrateResults2(id)
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
saveRoot = '/nobackup/p_neuert_lab/Jason/2023-12-08 Probe Designer_JH/';

HybridizationTemperature = inputsParameterSettings.Thermodynamic_Settings.HybridizationTemperature_Celsius; % Temperature for Evaluating Probe Design and Simulations
Nmodel = inputsParameterSettings.Thermodynamic_Settings.Gibbs_Model; %Which Hybridization model to use for probe design and evaluation


%% Compile Results from  Supporting Output Files

Nsoftware = 11;
designerName0 = '_NLPDS';

RoundN = 1;Z = 1;iz = 1;

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

settings.FolderRootName = strcat(saveRoot,'(',geneNames,')','_',strjoin(IncludeAccessionNumbers,'_'));
settings.rootName = strjoin(IncludeAccessionNumbers,'_');

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
        load([settings.FolderRootName '/' inputs1{gene_num,5}  '_Tm' num2str(HybridizationTemperature) '_BindingEnergyMatrix' designerName '.mat'],'Kb_mod')     
        load([settings.FolderRootName '/' inputs1{gene_num,5}  '_Tm' num2str(HybridizationTemperature) '_BindingEnergyMatrix2' designerName '.mat'],'Kb_Complement')      
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
        [K_S,~,~,~,~,~,~,~,...
        K_CD,~,~,~,~,~,~,~,~,~,~,~] = A_JH_GenerateSecondaryStructureInfo_V3(probes,Pset,settings);
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
        load([settings.FolderRootName '/' inputs1{gene_num,5}  '_Tm' num2str(HybridizationTemperature) '_BindingEnergyMatrix' designerName0 '.mat'],'Kb_mod')     
        load([settings.FolderRootName '/' inputs1{gene_num,5}  '_Tm' num2str(HybridizationTemperature) '_BindingEnergyMatrix2' designerName0 '.mat'],'Kb_Complement') 
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
        [K_S,~,~,~,~,~,~,~,...
        K_CD,~,~,~,~,~,~,~,~,~,~,~] = A_JH_GenerateSecondaryStructureInfo_V3(probes,Pset,settings);
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


            