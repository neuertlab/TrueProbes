function A_IntegrateResults4(id)
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
HybridizationTemperature = 37;
Nmodel = 4;

%% Compile Results from  Supporting Output Files
settings.FolderRootName = strcat(saveRoot,inputs1{gene_num,5},'_',strjoin(inputs1{gene_num,1},'_'));
settings.rootName = strjoin(inputs1{gene_num,1},'_');


Nsoftware = 11;
designerName0 = '_NLPDS';




%load([settings.FolderRootName '/' inputs1{gene_num,5} '_' num2str(HybridizationTemperature) '_DesignResults_CarryOver.mat'],'CarryOver2')
load([settings.FolderRootName '/' inputs1{gene_num,5} '_' num2str(HybridizationTemperature) '_DesignResults_CarryOver4.mat'],'CarryOver4')
for id = 1
ProbeDesignResults4.Gene = inputs1{gene_num,5};
ProbeDesignResults4.GeneTarget = settings.rootName;
ProbeDesignResults4.Organism = inputs1{gene_num,4};
ProbeDesignResults4.FolderName = settings.FolderRootName;
    for id2 = [1 8 9 2 3 4 5 6 7 10 11]
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
        ProbeDesignResults4.Software{id2} = designerName;  
        if (~contains(designerName,'Stellaris'))    
             try    
        %load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_probes' designerName '.mat'],'probes');
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_ExpressionInfo' designerName '.mat'],'settings');
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_Tm' num2str(HybridizationTemperature) '_ModelMetrics' designerName '.mat'],'ModelMetrics','chosenProbes');  
%         ON_IDs = CarryOver2.ON_IDs{id2};
%         OFF_IDs = CarryOver2.OFF_IDs{id2};
%         UnDesired_Isoforms = CarryOver2.UnDesired_Isoforms{id2};
        chosenProbes_Sorted = CarryOver4.Designed_Probes_ZigZagSelectionSorted{id2};
%         cTargetSites_Bound = ModelMetrics.cTargetSites_Bound;   
        ProbeDesignResults4.SoftwareResults(id2).settings = settings;
        ProbeDesignResults4.SoftwareResults(id2).Designed_Probes = chosenProbes;
        ProbeDesignResults4.SoftwareResults(id2).Designed_Probes_ZigZagSelectionSorted = chosenProbes_Sorted;
        ProbeDesignResults4.SoftwareResults(id2).ModelMetrics = ModelMetrics;clear ModelMetrics      
             catch
             end
%% Get Cumulative Concentrations ordered by Probes from high to low specificity
%         Size_Vec = arrayfun(@(x) size(cTargetSites_Bound,x),4:7);
%         Con_Cumulative = ndSparse.build([length(chosenProbes) Size_Vec],0);    
%         Coff_Cumulative = ndSparse.build([length(chosenProbes) Size_Vec],0);    
%         Con_otherIsos_Cumulative = ndSparse.build([length(chosenProbes) Size_Vec],0);   
%         Con_Fraction_Cumulative = ndSparse.build([length(chosenProbes) Size_Vec],0);    
%         Con_otherIsos_Fraction_Cumulative = ndSparse.build([length(chosenProbes) Size_Vec],0);    
%         Con_Coff_Ratio_Cumulative = ndSparse.build([length(chosenProbes) Size_Vec],0);    
%         Pset_sorted = arrayfun(@(x) find(chosenProbes==x) ,chosenProbes_Sorted);
%         for m=1:length(Pset_sorted) %solve equation first using 1st probe, then solve next given solution,
%             varySet = 1:m;
%             Con_Cumulative(m,:,:,:,:) = sum(sum(sum(cTargetSites_Bound(Pset_sorted(varySet),ON_IDs,:,:,:,:,:),1,'omitnan'),2,'omitnan'),3,'omitnan'); %Desired_Isoforms
%             Coff_Cumulative(m,:,:,:,:) = sum(sum(sum(cTargetSites_Bound(Pset_sorted(varySet),OFF_IDs,:,:,:,:,:),1,'omitnan'),2,'omitnan'),3,'omitnan'); 
%             Con_otherIsos_Cumulative(m,:,:,:,:) = sum(sum(sum(cTargetSites_Bound(Pset_sorted(varySet),UnDesired_Isoforms,:,:,:,:,:),1,'omitnan'),2,'omitnan'),3,'omitnan');
%             Con_Fraction_Cumulative(m,:,:,:,:) = Con_Cumulative(m,:,:,:,:)./(Con_Cumulative(m,:,:,:,:) + Coff_Cumulative(m,:,:,:,:) + Con_otherIsos_Cumulative(m,:,:,:,:));
%             Con_otherIsos_Fraction_Cumulative(m,:,:,:,:) = Con_otherIsos_Cumulative(m,:,:,:,:)./(Con_Cumulative(m,:,:,:,:) + Coff_Cumulative(m,:,:,:,:) + Con_otherIsos_Cumulative(m,:,:,:,:));
%             Con_Coff_Ratio_Cumulative(m,:,:,:,:) = Con_Cumulative(m,:,:,:,:)./Coff_Cumulative(m,:,:,:,:);
%         end 
%         ProbeDesignResults4.SoftwareResults(id2).CumulativeModelMetrics.Sorted_Con_Coff_Ratio = permute(Con_Coff_Ratio_Cumulative,[1 5 2 4 3]);clear Con_Coff_Ratio_Cumulative
%         ProbeDesignResults4.SoftwareResults(id2).CumulativeModelMetrics.Sorted_Con = permute(Con_Cumulative,[1 5 2 4 3]);clear Con_Cumulative
%         ProbeDesignResults4.SoftwareResults(id2).CumulativeModelMetrics.Sorted_Coff = permute(Coff_Cumulative,[1 5 2 4 3]);clear Coff_Cumulative
%         ProbeDesignResults4.SoftwareResults(id2).CumulativeModelMetrics.Sorted_Con_otherIsos = permute(Con_otherIsos_Cumulative,[1 5 2 4 3]);clear Con_otherIsos_Cumulative
%         ProbeDesignResults4.SoftwareResults(id2).CumulativeModelMetrics.Sorted_Con_Fraction = permute(Con_Fraction_Cumulative,[1 5 2 4 3]);clear Con_Fraction_Cumulative
%         ProbeDesignResults4.SoftwareResults(id2).CumulativeModelMetrics.Sorted_Con_otherIsos_Fraction = permute(Con_otherIsos_Fraction_Cumulative,[1 5 2 4 3]);clear Con_otherIsos_Fraction_Cumulative
            else
        %load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_probes' designerName0 '.mat'],'probes');
        try
        s=[settings.FolderRootName '/' inputs1{gene_num,5} '_Tm' num2str(HybridizationTemperature) '_ModelMetrics' designerName '.mat']
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_' settings.rootName '_ExpressionInfo' designerName0 '.mat'],'settings');
        load([settings.FolderRootName '/' inputs1{gene_num,5} '_Tm' num2str(HybridizationTemperature) '_ModelMetrics' designerName '.mat'],'ModelMetrics','chosenProbes');  
        tk = 1;
%         if (curr > 1)
%            tk = length(setdiff(currProbes,chosenProbes))+length(setdiff(chosenProbes,currProbes));
%         end
%         curr = curr + 1;
%         currProbes = chosenProbes;
%         ON_IDs = CarryOver2.ON_IDs{id2};
%         OFF_IDs = CarryOver2.OFF_IDs{id2};
%         UnDesired_Isoforms = CarryOver2.UnDesired_Isoforms{id2};
        chosenProbes_Sorted = CarryOver4.Designed_Probes_ZigZagSelectionSorted{id2};
%         cTargetSites_Bound = ModelMetrics.cTargetSites_Bound;   
        ProbeDesignResults4.SoftwareResults(id2).settings = settings;
        ProbeDesignResults4.SoftwareResults(id2).Designed_Probes = chosenProbes;
        ProbeDesignResults4.SoftwareResults(id2).Designed_Probes_ZigZagSelectionSorted = chosenProbes_Sorted;
        ProbeDesignResults4.SoftwareResults(id2).ModelMetrics = ModelMetrics;clear ModelMetrics   
        catch
        end
        %if (tk>0)
           
%% Get Cumulative Concentrations ordered by Probes from high to low specificity
%         Size_Vec = arrayfun(@(x) size(cTargetSites_Bound,x),4:7);
%         Con_Cumulative = ndSparse.build([length(chosenProbes) Size_Vec],0);    
%         Coff_Cumulative = ndSparse.build([length(chosenProbes) Size_Vec],0);    
%         Con_otherIsos_Cumulative = ndSparse.build([length(chosenProbes) Size_Vec],0);   
%         Con_Fraction_Cumulative = ndSparse.build([length(chosenProbes) Size_Vec],0);    
%         Con_otherIsos_Fraction_Cumulative = ndSparse.build([length(chosenProbes) Size_Vec],0);    
%         Con_Coff_Ratio_Cumulative = ndSparse.build([length(chosenProbes) Size_Vec],0);    
%         Pset_sorted = arrayfun(@(x) find(chosenProbes==x) ,chosenProbes_Sorted);
%         for m=1:length(Pset_sorted) %solve equation first using 1st probe, then solve next given solution,
%             varySet = 1:m;
%             Con_Cumulative(m,:,:,:,:) = sum(sum(sum(cTargetSites_Bound(Pset_sorted(varySet),ON_IDs,:,:,:,:,:),1,'omitnan'),2,'omitnan'),3,'omitnan'); %Desired_Isoforms
%             Coff_Cumulative(m,:,:,:,:) = sum(sum(sum(cTargetSites_Bound(Pset_sorted(varySet),OFF_IDs,:,:,:,:,:),1,'omitnan'),2,'omitnan'),3,'omitnan'); 
%             Con_otherIsos_Cumulative(m,:,:,:,:) = sum(sum(sum(cTargetSites_Bound(Pset_sorted(varySet),UnDesired_Isoforms,:,:,:,:,:),1,'omitnan'),2,'omitnan'),3,'omitnan');
%             Con_Fraction_Cumulative(m,:,:,:,:) = Con_Cumulative(m,:,:,:,:)./(Con_Cumulative(m,:,:,:,:) + Coff_Cumulative(m,:,:,:,:) + Con_otherIsos_Cumulative(m,:,:,:,:));
%             Con_otherIsos_Fraction_Cumulative(m,:,:,:,:) = Con_otherIsos_Cumulative(m,:,:,:,:)./(Con_Cumulative(m,:,:,:,:) + Coff_Cumulative(m,:,:,:,:) + Con_otherIsos_Cumulative(m,:,:,:,:));
%             Con_Coff_Ratio_Cumulative(m,:,:,:,:) = Con_Cumulative(m,:,:,:,:)./Coff_Cumulative(m,:,:,:,:);
%         end 
%         ProbeDesignResults4.SoftwareResults(id2).CumulativeModelMetrics.Sorted_Con_Coff_Ratio = permute(Con_Coff_Ratio_Cumulative,[1 5 2 4 3]);clear Con_Coff_Ratio_Cumulative
%         ProbeDesignResults4.SoftwareResults(id2).CumulativeModelMetrics.Sorted_Con = permute(Con_Cumulative,[1 5 2 4 3]);clear Con_Cumulative
%         ProbeDesignResults4.SoftwareResults(id2).CumulativeModelMetrics.Sorted_Coff = permute(Coff_Cumulative,[1 5 2 4 3]);clear Coff_Cumulative
%         ProbeDesignResults4.SoftwareResults(id2).CumulativeModelMetrics.Sorted_Con_otherIsos = permute(Con_otherIsos_Cumulative,[1 5 2 4 3]);clear Con_otherIsos_Cumulative
%         ProbeDesignResults4.SoftwareResults(id2).CumulativeModelMetrics.Sorted_Con_Fraction = permute(Con_Fraction_Cumulative,[1 5 2 4 3]);clear Con_Fraction_Cumulative
%         ProbeDesignResults4.SoftwareResults(id2).CumulativeModelMetrics.Sorted_Con_otherIsos_Fraction = permute(Con_otherIsos_Fraction_Cumulative,[1 5 2 4 3]);clear Con_otherIsos_Fraction_Cumulative  
        %end
        end
        id2
    end
    save([settings.FolderRootName '/' inputs1{gene_num,5} '_' num2str(HybridizationTemperature) '_DesignResults4.mat'],'ProbeDesignResults4','-v7.3')     
end




end


