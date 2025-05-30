function A0_BKJH_ProbeComparison_Wrapper_IntegrateResults_V2(id)
inputs1 = ...  
    {... %%RefSeq number and specific transcript name which gives isoform specific
    {'NM_005564.5'},{},{}, 'Human','(LCN2)','(LCN2)','9',{},1 ;... %820bp 
    {'NM_005764.4'},{},{}, 'Human','(PDZK1IP1)','(PDZK1IP1)','1',{},1 ;...  %838bp
    {'NM_006507.4'},{},{}, 'Human','(REG1B)','(REG1B)','2',{},1 ;... % 775bp
    {'NM_001657.4'},{},{}, 'Human','(AREG)','(AREG)','4',{},1 ;... %1234bp
    {'NM_002046.7'},{},{}, 'Human','(GAPDH)','(GAPDH)','12',{},1 ;...%1285bp
    {'NM_000591.4'},{},{}, 'Human','(CD14)','(CD14)','5',{},1;... %1356bp
    {'NM_000733.4'},{},{}, 'Human','(CD3E)','(CD3E)','11',{},1 ;...  %1361bp
    {'NM_000045.4'},{},{}, 'Human','(ARG1)','(ARG1)','6',{},1 ;...%1447bp
    {'NM_020529.3'},{},{}, 'Human','(NFKBIA)','(NFKBIA)','14',{},1;... %1559bp
    {'NM_001251.3'},{},{}, 'Human','(CD68)','(CD68)','17',{},1 ;... % 1705bp
    {'NM_001101.5'},{},{}, 'Human','(ACTB)','(ACTB)','7',{},1 ;... %1812bp
    {'NM_001172.4'},{},{}, 'Human','(ARG2)','(ARG2)','14',{},1 ;...    %1911bp
    {'NM_001770.6'},{},{}, 'Human','(CD19)','(CD19)','16',{},1 ;...   %1918bp
    {'NM_014009.4'},{},{}, 'Human','(FOXP3)','(FOXP3)','X',{},1 ;...  %2264bp
    {'NM_033486.3'},{},{},'Human','(CDK11B)','(CDK11B)','1',{},1;... %2992bp  Has 27 isoforms
    {'NM_001145873.1'},{},{}, 'Human','(CD8A)','(CD8A)','2',{},1 ;... %3177bp
    {'NM_152866.3'},{},{}, 'Human','(MS4A1)','(MS4A1)','11',{},1 ;... %3556bp
    {'NM_203416.4'},{},{}, 'Human','(CD163)','(CD163)','12',{},1 ;... %4154bp
    {'NM_004931.5'},{},{}, 'Human','(CD8B)','(CD8B)','2',{},1 ;...  %4794bp
    {'NM_002838.5'},{},{}, 'Human','(PTPRC)','(PTPRC)','1',{},1 ;... %5357bp 
    };

%Neuert Lab initial Probes,simulated annealing Probes, Old Designer Probes,
% Stellaris, Oligostan, MERFISH, PaintSHOP, Unspecific (design in reverse)

%% Settings Specification
saveRoot = '/gpfs23/scratch/neuertg/Jason/2022-6-24 Probe Designer_JH/';
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
settings.FolderRootName = strcat(saveRoot,inputs1{gene_num,5},'_',strjoin(inputs1{gene_num,1},'_'));
settings.rootName = strjoin(inputs1{gene_num,1},'_');
T_hybrid = HybridizationTemperature;
Lmin = minProbeSize; Lmax = maxProbeSize;
%% Compile Results from  Supporting Output Files
load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_probes.mat'],'probes');
load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_hits.mat'],'gene_hits')
load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_hits_table.mat'],'gene_table') 
load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_ExpressionInfo.mat'],'ExpressionMatrix');
load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_Tm' num2str(HybridizationTemperature) '_RM' num2str(RemoveMisMatches) '_OnOffThermoInfo.mat'],'Kon','Koff'); 
load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_ProbeProperties.mat'],'Probe_Properties','probe_decision_matrix')
load([settings.FolderRootName '/' inputs1{gene_num,5} '_binding_hits_map.mat'],'DoesProbeBindSite2','Num_of_Molecule_Sites')  
load([settings.FolderRootName '/' inputs1{gene_num,5}  '_Tm' num2str(T_hybrid) '_RM' num2str(RemoveMisMatches) '_BindingEnergyMatrix.mat'],'Kb_mod')     
load([settings.FolderRootName '/' inputs1{gene_num,5}  '_Tm' num2str(T_hybrid) '_RM' num2str(RemoveMisMatches) '_BindingEnergyMatrix2.mat'],'Kb_Complement') 
try
    load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_RNASecondaryStructure.mat'],'-mat','Ks_TjLi','IsSiteInLoop')    
catch
    Ks_TjLi = [];
    IsSiteInLoop = [];
end    
%% Get Results from other software
ProbeDesignResults.Gene = inputs1{gene_num,5};
ProbeDesignResults.GeneTarget = settings.rootName;
ProbeDesignResults.Organism = inputs1{gene_num,4};
ProbeDesignResults.FolderName = settings.FolderRootName;
DesignLandscape.probes = probes;
DesignLandscape.gene_hits = gene_hits;
DesignLandscape.gene_hits_table = gene_table;
DesignLandscape.ExpressionMatrix = ExpressionMatrix;
DesignLandscape.Probe_Properties = Probe_Properties;
DesignLandscape.Probe_DecisionMatrix = probe_decision_matrix;
DesignLandscape.BindingSiteMap = DoesProbeBindSite2;
DesignLandscape.NumOfMoleculeSites = Num_of_Molecule_Sites;
DesignLandscape.Kb = Kb_mod;
DesignLandscape.Kb_Complement = Kb_Complement;
DesignLandscape.Kon = Kon;
DesignLandscape.Koff = Koff;
DesignLandscape.Krss_TjLi = Ks_TjLi;
DesignLandscape.rss_IsSiteInLoop = IsSiteInLoop;
ProbeDesignResults.DesignLandscape = DesignLandscape;
for id2 = 1:108
switch id2
    case 1
        designerName = '';
    case 2
        designerName = 'OldNeuertLab';
    case 3
        designerName = 'Stellaris';
    case 4
        designerName = 'OligoStan';
    case 5
        designerName = 'PaintSHOP';
    case 6
        designerName = 'MERFISH';
    case 7
        designerName = 'Unspecific';
    otherwise
        designerName = strcat('Random',num2str(id2-7));  
end
    try  
        if (id2==1)
            ProbeDesignResults.Software{id2} = 'NeuertLab'; 
            load([settings.FolderRootName '/' inputs1{gene_num,5} '_DesignMetrics.mat'],...
            'BindingProbs_Full','BindingProbs_Subset','ProbeSetMetrics','NormalizedMetrics','chosenProbes','settings','Ks','Kd')
            ProbeDesignResults.settings = settings;
        else
           ProbeDesignResults.Software{id2} = designerName; 
           load([settings.FolderRootName '/' inputs1{gene_num,5} '_DesignMetrics' designerName '.mat'],...
            'BindingProbs_Full','BindingProbs_Subset','ProbeSetMetrics','NormalizedMetrics','chosenProbes','Ks','Kd')
        end 
        ProbeDesignResults.SoftwareResults(id2).Designed_Probes = chosenProbes;
        ProbeDesignResults.SoftwareResults(id2).Probe_Ks = Ks;
        ProbeDesignResults.SoftwareResults(id2).Probe_Kd = Kd;
        %what is threshold 5%, 10%, conventional for each gene and values
        %also scaled number of probes above and below threshold
        ProbeDesignResults.SoftwareResults(id2).ProbeSetMetrics = ProbeSetMetrics;
        ProbeDesignResults.SoftwareResults(id2).NormalizedMetrics = NormalizedMetrics;
        ProbeDesignResults.SoftwareResults(id2).BindingProbs_Full = BindingProbs_Full;
        if (~isempty(BindingProbs_Subset{:}))
            for h = 1:size(BindingProbs_Subset,2)
            ProbeDesignResults.SoftwareResults(id2).BindingProbs_Subset{h} = BindingProbs_Subset{h};
            end
            clear BindingProbs_Subset
        else
            try
            %errors?
           for h = 1:size(BindingProbs_Full,2)
                for x1 = 1:size(ProbeSetMetrics{h},1)
                    for x2 = 1:size(ProbeSetMetrics{h},2)
                        ProbeDesignResults.SoftwareResults(id2).ProbeSetMetrics{h}{x1,x2}.Spots_Subset = ProbeSetMetrics{h}{x1,x2}.Spots_Full;
                        ProbeDesignResults.SoftwareResults(id2).ProbeSetMetrics{h}{x1,x2}.Probe_Subset = ProbeSetMetrics{h}{x1,x2}.Probe_Full;
                    end
                end
                ProbeDesignResults.SoftwareResults(id2).BindingProbs_Subset{h} = BindingProbs_Full{h}; 
                ProbeDesignResults.SoftwareResults(id2).BindingProbs_Subset{h}.MoleculeSignal_Subset = BindingProbs_Full{h}.MoleculeSignal_Full; 
                ProbeDesignResults.SoftwareResults(id2).BindingProbs_Subset{h}.MoleculeSignal_Subset_1FLAP = BindingProbs_Full{h}.MoleculeSignal_Full_1FLAP; 
                ProbeDesignResults.SoftwareResults(id2).BindingProbs_Subset{h}.MoleculeSignal_Subset_2FLAP = BindingProbs_Full{h}.MoleculeSignal_Full_2FLAP; 
                ProbeDesignResults.SoftwareResults(id2).BindingProbs_Subset{h}.MoleculeCount_Subset = BindingProbs_Full{h}.MoleculeCount_Full; 
                ProbeDesignResults.SoftwareResults(id2).BindingProbs_Subset{h}.MoleculeCount_Subset_1FLAP = BindingProbs_Full{h}.MoleculeCount_Full_1FLAP; 
                ProbeDesignResults.SoftwareResults(id2).BindingProbs_Subset{h}.MoleculeCount_Subset_2FLAP = BindingProbs_Full{h}.MoleculeCount_Full_2FLAP;  
            end
            catch ME
                disp(ME.message)
            end
        end
        clear chosenProbes Ks Kd ProbeSetMetrics BindingProbs_Full NormalizedMetrics
        %determine special threshold for X%, and conventional
        %ProbeSetMetrics{h}.Spots_Full.P
        %largest drop off. %first derivative largest
    catch
        ProbeDesignResults.Software{id2} = designerName;
        ProbeDesignResults.SoftwareResults(id2).Designed_Probes = [];
        ProbeDesignResults.SoftwareResults(id2).Probe_Ks = [];
        ProbeDesignResults.SoftwareResults(id2).Probe_Kd = [];
        ProbeDesignResults.SoftwareResults(id2).ProbeSetMetrics = [];
        ProbeDesignResults.SoftwareResults(id2).NormalizedMetrics = [];
        ProbeDesignResults.SoftwareResults(id2).BindingProbs_Full = [];
        ProbeDesignResults.SoftwareResults(id2).BindingProbs_Subset = [];
    end
end

save([settings.FolderRootName '/' inputs1{gene_num,5} '_' num2str(HybridizationTemperature) '_' num2str(SaltConcentration) '_' num2str(RemoveMisMatches) '_DesignResults.mat'],'ProbeDesignResults','-v7.3')
end
