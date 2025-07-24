function A_IntegrateResults3(id)
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
gene_num = id;
HybridizationTemperature = inputsParameterSettings.Thermodynamic_Settings.HybridizationTemperature_Celsius; % Temperature for Evaluating Probe Design and Simulations
Nmodel = inputsParameterSettings.Thermodynamic_Settings.Gibbs_Model; %Which Hybridization model to use for probe design and evaluation

%Save Settings
%%Load Data
Nsoftware = 11;
designerName0 = '_NLPDS';


%Neuert Lab initial Probes,simulated annealing Probes, Old Designer Probes,
% Stellaris, Oligostan, MERFISH, PaintSHOP, Unspecific (design in reverse)

%% Settings Specification
%Save Settings
%% Compile Results from  Supporting Output Files
settings.FolderRootName = strcat(saveRoot,inputs1{gene_num,5},'_',strjoin(inputs1{gene_num,1},'_'));
settings.rootName = strjoin(inputs1{gene_num,1},'_');



load([settings.FolderRootName '/' inputs1{gene_num,5} '_' num2str(HybridizationTemperature) '_DesignResults_CarryOver.mat'],'CarryOver')
for id = 1
ProbeDesignResults3.Gene = inputs1{gene_num,5};
ProbeDesignResults3.GeneTarget = settings.rootName;
ProbeDesignResults3.Organism = inputs1{gene_num,4};
ProbeDesignResults3.FolderName = settings.FolderRootName;
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
        ProbeDesignResults3.Software{id2} = designerName; 
        try
        if (~contains(designerName,'Stellaris'))        
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_probes' designerName '.mat'],'probes');
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_ExpressionInfo' designerName '.mat'],'ExpressionMatrix','settings');
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_Tm' num2str(HybridizationTemperature) '_OnOffThermoInfo' designerName '.mat'],'Kon'); 
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_binding_hits_map' designerName '.mat'],'DoesProbeBindSite2')  
        load([settings.FolderRootName '/' inputs1{gene_num,5}  '_Tm' num2str(HybridizationTemperature) '_BindingEnergyMatrix' designerName '.mat'],'Kb_mod')  
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_Tm' num2str(HybridizationTemperature) '_BasicDesignerStats' designerName '.mat'],'Tvec_RNA','Svec_RNA','TPvec_RNA','TSvec_RNA','TPvec_logKOFF_RNA','TPvec_logKOFFdivON_RNA','TPvec_logKONdivOFF_RNA',...
                'Nvec_RNAmulti','Tvec_DNA','Svec_DNA','TPvec_DNA','TSvec_DNA','TPvec_logKOFF_DNA','TPvec_logKOFFdivON_DNA','TPvec_logKONdivOFF_DNA','TPvec_logKOFFdivCOMP_DNA','TPvec_logKCOMPdivOFF_DNA',...
                'Off_Score','Specificity_Score','NumRNAOffTargetOptions','Probes_WithNRNAOFF','NumDNAOffTargetOptions','Probes_WithNDNAOFF');   
        if (id2==1)
           load([settings.FolderRootName '/' inputs1{gene_num,5} '_' inputs1{gene_num,1}{:}  '_chosen.mat'],'chosenProbes');  
        else
           chosenProbes = 1:size(probes,1);
        end
        AllowableProbes = CarryOver.AllowableProbes{id2};
        ProbesWithRibosomalHits = CarryOver.ProbesWithRibosomalHits{id2};
        OFF_RNAIDs = CarryOver.OFF_RNAIDs{id2};
        Kon1 = squeeze(Kon(:,Nmodel));
        Kb1 = squeeze(Kb_mod(:,:,:,Nmodel));   
        notDone = setdiff(1:length(OFF_RNAIDs),1:size(TPvec_RNA,2));
        Tp = @(x) find(sum(squeeze(DoesProbeBindSite2(:,x,:)),2)>0); 
        Tx =@(y,Z) arrayfun(@(x) find(squeeze(DoesProbeBindSite2(x,y,:))==1),Z,'Un',0);
        for t = notDone
           TPvec_RNA0{t} = Tp(OFF_RNAIDs(t));%does not have multiplicity
           temp_RNAV1b = Tx(OFF_RNAIDs(t),TPvec_RNA0{t});
           temp_RNAV2b = cellfun(@length,temp_RNAV1b);
           TPvec_RNA{t} = cell2mat(arrayfun(@(x) repmat(TPvec_RNA0{t}(x),[1 temp_RNAV2b(x)]),1:length(TPvec_RNA0{t}),'Un',0));
           TSvec_RNA{t} = cell2mat(temp_RNAV1b);%site locations 
           TPvec_logKOFF_RNA{t} = log10(diag(full(squeeze(Kb1(TPvec_RNA{t},OFF_RNAIDs(t),TSvec_RNA{t}))))');
           TPvec_logKOFFdivON_RNA{t} = log10(diag(full(squeeze(Kb1(TPvec_RNA{t},OFF_RNAIDs(t),TSvec_RNA{t}))))'./Kon1(TPvec_RNA{t}));
           TPvec_logKONdivOFF_RNA{t} = log10(Kon1(TPvec_RNA{t})./diag(full(squeeze(Kb1(TPvec_RNA{t},OFF_RNAIDs(t),TSvec_RNA{t}))))');
        end    
        clear TPvec_RNA0 
        Tvec_logKOFF_RNA = arrayfun(@(p) log10(diag(full(squeeze(Kb1(p,Tvec_RNA{p},Svec_RNA{p})))))',1:size(probes,1),'Un',0);
        Tvec_logKOFFdivON_RNA = arrayfun(@(p) log10(diag(full(squeeze(Kb1(p,Tvec_RNA{p},Svec_RNA{p}))))'/Kon1(p)),1:size(probes,1),'Un',0); 
        NumOffTargetOptions = unique(Nvec_RNAmulti(AllowableProbes));
        Probes_WithNOFF_targets = arrayfun(@(i) AllowableProbes(Nvec_RNAmulti(AllowableProbes)==NumOffTargetOptions(i)),1:length(NumOffTargetOptions),'Un',0);
        if (id2>1)
%         AllowMatrix = zeros(length(AllowableProbes),length(AllowableProbes));
%         for u1 = 1:length(AllowableProbes)
%             u = AllowableProbes(u1);
%             for v1 = 1:length(AllowableProbes)
%             v = AllowableProbes(v1);
%             if (u~=v)
%                AllowMatrix(u1,v1) = double(length(intersect([probes{u,3}-spacing:probes{u,3}+length(probes{u,2})-1+spacing],[probes{v,3}-spacing:probes{v,3}+length(probes{v,2})-1+spacing]))>spacing);   
%             end
%             end
%         end
        [~,idx] = sort(Nvec_RNAmulti(chosenProbes),'ascend');
        chosenProbes_Sorted = chosenProbes(idx);
        %order by number of probe off-targets, Kon - Koff;
        %chosenProbes_Sorted = A_ZigZagProbeSelection_V5(AllowableProbes,AllowMatrix,ProbesWithRibosomalHits,OFF_RNAIDs,probes,settings,addSelfProb,packOptimal,Kon,Nvec_RNAmulti,Off_Score,Specificity_Score,Tvec_RNA,Svec_RNA,TPvec_RNA,TSvec_RNA,TPvec_logKOFF_RNA,TPvec_logKOFFdivON_RNA,TPvec_logKONdivOFF_RNA,ExpressionMatrix,DoesProbeBindSite2,Kb_mod);
        else
        chosenProbes_Sorted = chosenProbes;   
        end
        ProbeDesignResultsD.SoftwareResults(id2).Designed_Probes_ZigZagSelectionSorted = chosenProbes_Sorted;
        ProbeDesignResults3.SoftwareResults(id2).settings = settings;
        ProbeDesignResults3.SoftwareResults(id2).probes = probes; clear Probes
        ProbeDesignResults3.SoftwareResults(id2).Designed_Probes = chosenProbes;clear DoesProbeBindSite2 chosenProbes
        ProbeDesignResults3.SoftwareResults(id2).Designed_Probes_ZigZagSelectionSorted = chosenProbes_Sorted;clear chosenProbes_Sorted
        ProbeDesignResults3.SoftwareResults(id2).AllowableProbes = AllowableProbes;
        ProbeDesignResults3.SoftwareResults(id2).ExpressionMatrix = ExpressionMatrix;clear ExpressionMatrix
        ProbeDesignResults3.SoftwareResults(id2).Tvec_RNA = Tvec_RNA;clear Tvec_RNA
        ProbeDesignResults3.SoftwareResults(id2).Svec_RNA = Svec_RNA;clear Svec_RNA
        ProbeDesignResults3.SoftwareResults(id2).TPvec_RNA = TPvec_RNA;clear TPvec_RNA
        ProbeDesignResults3.SoftwareResults(id2).TSvec_RNA = TSvec_RNA;clear TSvec_RNA
        ProbeDesignResults3.SoftwareResults(id2).Tvec_logKOFF_RNA = Tvec_logKOFF_RNA;clear Tvec_logKOFF_RNA
        ProbeDesignResults3.SoftwareResults(id2).Tvec_logKOFFdivON_RNA = Tvec_logKOFFdivON_RNA;clear Tvec_logKOFFdivON_RNA
        ProbeDesignResults3.SoftwareResults(id2).TPvec_logKOFF_RNA = TPvec_logKOFF_RNA;clear TPvec_logKOFF_RNA
        ProbeDesignResults3.SoftwareResults(id2).TPvec_logKOFFdivON_RNA = TPvec_logKOFFdivON_RNA;clear TPvec_logKOFFdivON_RNA
        ProbeDesignResults3.SoftwareResults(id2).TPvec_logKONdivOFF_RNA = TPvec_logKONdivOFF_RNA;clear TPvec_logKONdivOFF_RNA
        ProbeDesignResults3.SoftwareResults(id2).Nvec_RNAmulti = Nvec_RNAmulti;clear Nvec_RNAmulti
        ProbeDesignResults3.SoftwareResults(id2).Off_Score = Off_Score; clear Off_Score
        ProbeDesignResults3.SoftwareResults(id2).Specificity_Score = Specificity_Score; clear Specificity_Score
        ProbeDesignResults3.SoftwareResults(id2).NumRNAOffTargetOptions = NumRNAOffTargetOptions; clear NumRNAOffTargetOptions
        ProbeDesignResults3.SoftwareResults(id2).Probes_WithNRNAOFF = Probes_WithNRNAOFF; clear Probes_WithNRNAOFF
        ProbeDesignResults3.SoftwareResults(id2).NumDNAOffTargetOptions = NumDNAOffTargetOptions; clear umDNAOffTargetOptions
        ProbeDesignResults3.SoftwareResults(id2).Probes_WithNDNAOFF = Probes_WithNDNAOFF; clear Probes_WithNDNAOFF 
        ProbeDesignResults3.SoftwareResults(id2).ProbesWithRibosomalHits = ProbesWithRibosomalHits; clear ProbesWithRibosomalHits
        ProbeDesignResults3.SoftwareResults(id2).Probes_WithNOFF_targets = Probes_WithNOFF_targets; clear Probes_WithNOFF_targets
        ProbeDesignResults3.SoftwareResults(id2).NumOffTargetOptions = NumOffTargetOptions;clear NumOffTargetOptions
        ProbeDesignResults3.SoftwareResults(id2).Tvec_DNA = Tvec_DNA;clear Tvec_DNA
        ProbeDesignResults3.SoftwareResults(id2).Svec_DNA = Svec_DNA;clear Svec_DNA
        ProbeDesignResults3.SoftwareResults(id2).TPvec_DNA = TPvec_DNA;clear TPvec_DNA
        ProbeDesignResults3.SoftwareResults(id2).TSvec_DNA = TSvec_DNA;clear TSvec_DNA
        ProbeDesignResults3.SoftwareResults(id2).TPvec_logKOFF_DNA = TPvec_logKOFF_DNA;clear TPvec_logKOFF_DNA
        ProbeDesignResults3.SoftwareResults(id2).TPvec_logKOFFdivON_DNA = TPvec_logKOFFdivON_DNA;clear TPvec_logKOFFdivON_DNA
        ProbeDesignResults3.SoftwareResults(id2).TPvec_logKONdivOFF_DNA = TPvec_logKONdivOFF_DNA;clear TPvec_logKONdivOFF_DNA
        ProbeDesignResults3.SoftwareResults(id2).TPvec_logKOFFdivCOMP_DNA = TPvec_logKOFFdivCOMP_DNA;clear TPvec_logKOFFdivCOMP_DNA
        ProbeDesignResults3.SoftwareResults(id2).TPvec_logKCOMPdivOFF_DNA = TPvec_logKCOMPdivOFF_DNA;clear TPvec_logKCOMPdivOFF_DNA
        else
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_probes' designerName0 '.mat'],'probes');
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_ExpressionInfo' designerName0 '.mat'],'settings');
        %load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_Tm' num2str(HybridizationTemperature) '_OnOffThermoInfo' designerName0 '.mat'],'Kon'); 
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_Tm' num2str(HybridizationTemperature) '_BasicDesignerStats' designerName0 '.mat'],...
                'Nvec_RNAmulti');      
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_probes' designerName '.mat'],'chosenProbes');
        if (id2>1)
%         AllowMatrix = zeros(length(AllowableProbes),length(AllowableProbes));
%         for u1 = 1:length(AllowableProbes)
%             u = AllowableProbes(u1);
%             for v1 = 1:length(AllowableProbes)
%             v = AllowableProbes(v1);
%             if (u~=v)
%                AllowMatrix(u1,v1) = double(length(intersect([probes{u,3}-spacing:probes{u,3}+length(probes{u,2})-1+spacing],[probes{v,3}-spacing:probes{v,3}+length(probes{v,2})-1+spacing]))>spacing);   
%             end
%             end
%         end
        [~,idx] = sort(Nvec_RNAmulti(chosenProbes),'ascend');
        chosenProbes_Sorted = chosenProbes(idx);
%         chosenProbes_Sorted = A_ZigZagProbeSelection_V5(AllowableProbes,AllowMatrix,ProbesWithRibosomalHits,OFF_RNAIDs,probes,settings,addSelfProb,packOptimal,Kon,Nvec_RNAmulti,Off_Score,Specificity_Score,Tvec_RNA,Svec_RNA,TPvec_RNA,TSvec_RNA,TPvec_logKOFF_RNA,TPvec_logKOFFdivON_RNA,TPvec_logKONdivOFF_RNA,ExpressionMatrix,DoesProbeBindSite2,Kb_mod);
        else
        chosenProbes_Sorted = chosenProbes;    
        end
        ProbeDesignResultsD.SoftwareResults(id2).Designed_Probes_ZigZagSelectionSorted = chosenProbes_Sorted;
        ProbeDesignResults3.SoftwareResults(id2).settings = settings;
        ProbeDesignResults3.SoftwareResults(id2).probes = probes; clear Probes
        ProbeDesignResults3.SoftwareResults(id2).Designed_Probes = chosenProbes;clear chosenProbes
        ProbeDesignResults3.SoftwareResults(id2).Designed_Probes_ZigZagSelectionSorted = chosenProbes_Sorted;clear chosenProbes_Sorted 
        end
        catch
        end
    end 
end
save([settings.FolderRootName '/' inputs1{gene_num,5} '_' num2str(HybridizationTemperature) '_DesignResults3.mat'],'ProbeDesignResults3','-v7.3')
CarryOver4.Designed_Probes_ZigZagSelectionSorted = {ProbeDesignResultsD.SoftwareResults(:).Designed_Probes_ZigZagSelectionSorted};
save([settings.FolderRootName '/' inputs1{gene_num,5} '_' num2str(HybridizationTemperature) '_DesignResults_CarryOver4.mat'],'CarryOver4','-v7.3')





end


