function RNA_FISH_ExperimentalAnalysisImageReplicaAggregation_V5(part)
%% load experimental data
no_probe_conditions = {'NoProbes'};
SoftwareNames = {'TrueProbes','Stellaris','OligostanHT','PaintSHOP','MERFISH'};
SoftwareNames_Spelling = {{'trueprobes'},{'stellaris'},{'oligostan','ogliostan'},{'paintshop'},{'merfish'}};
Replica = {'Rep1','Rep2','Rep3','Rep4','Rep5A','Rep5B','Rep6A','Rep6B'};
unique_concentrations = {'_1000','_2000','_4000','MolAdj','NewCalc_MolAdj'};
group_concentrations = {'Cy5_1000','Cy5_2000','Cy5_4000','MolAdj','NewMolAdj'};
replica_num = 1:length(Replica);
thresholdObj = {'Image','Cell'};
thresholdTypes = {'min','q1','mean','median','q3','opt'};
thresholdTitleA = {'Single-Image Threshold:','Single-Cell Threshold'};
thresholdTitleB = {'Min','Quartile 1','Mean','Median','Quartile 3','Optimal'};
n_threshold_options = length(thresholdTypes);

try
    load_pooled_results = 1;
    if (part==1)
        load('pooled_hARF4_RNAFISH_AnalysisQuant.mat','noProbePooledResults_quant','softwarePooledResults_quant');
    end
    if (part==2)
        load('pooled_hARF4_RNAFISH_AnalysisStats.mat','noProbePooledResults_stats','softwarePooledResults_stats');
    end
    if (part==3)
        load('pooled_hARF4_RNAFISH_AnalysisThresholdInfo.mat','noProbePooledResults_thresholdInfo','softwarePooledResults_thresholdInfo');
    end
    if (part==4)
        load('pooled_hARF4_RNAFISH_AnalysisImgDat.mat','noProbePooledResults_imgData','softwarePooledResults_imgData');
        load('pooled_hARF4_RNAFISH_AnalysisCellseg.mat','noProbePooledResults_imgData','softwarePooledResults_imgData');
    end
catch
    load_pooled_results = 0;
end
if (load_pooled_results)
else
    try
        load('pooledIndex_hARF4_RNAFISH_Analysis.mat','noProbeResults','softwareResults');
    catch
        load('hARF4_RNAFISH_Analysis_currentHidden.mat','results');
        resultsA = results;clear results
        load('hARF4_RNAFISH_Analysis_currentMain.mat','results');
        resultsB = results;clear results
        if (~isempty(fieldnames(resultsA)))
            results = [resultsA resultsB];
        else
            results = resultsB;
        end
        clear resultsA resultsB
        load('hARF4_RNAFISH_AnalysisQuantResults_currentHidden.mat','results_quant');
        results_quantA = results_quant;clear results_quant
        load('hARF4_RNAFISH_AnalysisQuantResults_currentMain.mat','results_quant');
        results_quantB = results_quant;clear results_quant 
        if (~isempty(fieldnames(results_quantA)))
            results_quant = [results_quantA results_quantB];
        else
            results_quant = results_quantB;
        end
        clear results_quantA results_quantB
        load('hARF4_RNAFISH_AnalysisCellStats_currentHidden.mat','results_cellStats');
        results_cellStatsA = results_cellStats;clear results_cellStats
        load('hARF4_RNAFISH_AnalysisCellStats_currentMain.mat','results_cellStats');
        results_cellStatsB = results_cellStats;clear results_cellStats
        if (~isempty(fieldnames(results_cellStatsA)))
            results_cellStats = [results_cellStatsA results_cellStatsB];
        else
            results_cellStats = results_cellStatsB;
        end
        clear results_cellStatsA results_cellStatsB
        load('hARF4_RNAFISH_AnalysisSpotStats_currentHidden.mat','results_spotStats');
        results_spotStatsA = results_spotStats;clear results_spotStats
        load('hARF4_RNAFISH_AnalysisSpotStats_currentMain.mat','results_spotStats');
        results_spotStatsB = results_spotStats;clear results_spotStats
        if (~isempty(fieldnames(results_spotStatsA)))
            results_spotStats = [results_spotStatsA results_spotStatsB];
        else
            results_spotStats = results_spotStatsB;
        end
        clear results_spotStatsA results_spotStatsB

        celldata_exists = find(double(~contains({results.file},'70C')).*...
            double(~contains({results.file},'75C')).*...
            double(~contains({results.file},'80C')).*...
            double(~contains({results.file},'85C')).*...
            double(~contains({results.file},'500')));
        existing_files = {results(celldata_exists).file}';
        temp_files = extractAfter(existing_files,'Jurkat_R');
        replica_num_vector = zeros(size(existing_files));
        for f = 1:length(celldata_exists)
            if (str2double(temp_files{f}(1))<5)
                replica_num_vector(f) = str2double(temp_files{f}(1));
            else
                replica_num_vector(f) = find(strcmp(extractAfter(Replica,'Rep'),temp_files{f}(1:2)));
            end
            results(celldata_exists(f)).replica_number = replica_num_vector(f);
        end
        softwareResults = struct();noProbeResults = struct();
        software_conditionQuantStatsExists_matches = cell(length(SoftwareNames_Spelling),length(group_concentrations),length(Replica));
        software_conditionSpotStatsExists_matches = cell(length(SoftwareNames_Spelling),length(group_concentrations),length(Replica));
        software_conditionCellStatsExists_matches = cell(length(SoftwareNames_Spelling),length(group_concentrations),length(Replica));

        %% group all no probe images by replica condition and group all probe images by replica condition, software, and probe concentration
        for ik = 1:length(Replica)
            noProbeResults.(Replica{ik}).ind_results = [results(celldata_exists(double([results_quant(celldata_exists).QuantStatsExist]'==1).*double(replica_num_vector==replica_num(ik)).*contains(existing_files,no_probe_conditions)==1))];
            noProbeResults.(Replica{ik}).ind_results_part_function = @(results_part_file) [results_part_file(celldata_exists(double([results_quant(celldata_exists).QuantStatsExist]'==1).*double(replica_num_vector==replica_num(ik)).*contains(existing_files,no_probe_conditions)==1))];
            noProbeResults.(Replica{ik}).ind_results_part_spots_function = @(results_part_file) [results_part_file(celldata_exists(double([results_spotStats(celldata_exists).SpotStatsExist]'==1).*double(replica_num_vector==replica_num(ik)).*contains(existing_files,no_probe_conditions)==1))];
            noProbeResults.(Replica{ik}).ind_results_part_cells_function = @(results_part_file) [results_part_file(celldata_exists(double([results_cellStats(celldata_exists).CellStatsExist]'==1).*double(replica_num_vector==replica_num(ik)).*contains(existing_files,no_probe_conditions)==1))];
            for ii = 1:length(SoftwareNames_Spelling)
                for ij = 1:length(group_concentrations)
                    if (strcmp(group_concentrations{ij},'MolAdj'))
                        software_conditionQuantStatsExists_matches{ii,ij} = celldata_exists(double([results_quant(celldata_exists).QuantStatsExist]'==1).*double(replica_num_vector==replica_num(ik)).*contains(existing_files,unique_concentrations{ij}).*~contains(existing_files,'NewCalc').*contains(lower(existing_files),SoftwareNames_Spelling{ii})==1);
                        software_conditionSpotStatsExists_matches{ii,ij} = celldata_exists(double([results_spotStats(celldata_exists).SpotStatsExist]'==1).*double(replica_num_vector==replica_num(ik)).*contains(existing_files,unique_concentrations{ij}).*~contains(existing_files,'NewCalc').*contains(lower(existing_files),SoftwareNames_Spelling{ii})==1);
                        software_conditionCellStatsExists_matches{ii,ij} = celldata_exists(double([results_cellStats(celldata_exists).CellStatsExist]'==1).*double(replica_num_vector==replica_num(ik)).*contains(existing_files,unique_concentrations{ij}).*~contains(existing_files,'NewCalc').*contains(lower(existing_files),SoftwareNames_Spelling{ii})==1);
                    else
                        software_conditionQuantStatsExists_matches{ii,ij} = celldata_exists(double([results_quant(celldata_exists).QuantStatsExist]'==1).*double(replica_num_vector==replica_num(ik)).*contains(existing_files,unique_concentrations{ij}).*contains(lower(existing_files),SoftwareNames_Spelling{ii})==1);
                        software_conditionSpotStatsExists_matches{ii,ij} = celldata_exists(double([results_spotStats(celldata_exists).SpotStatsExist]'==1).*double(replica_num_vector==replica_num(ik)).*contains(existing_files,unique_concentrations{ij}).*contains(lower(existing_files),SoftwareNames_Spelling{ii})==1);
                        software_conditionCellStatsExists_matches{ii,ij} = celldata_exists(double([results_cellStats(celldata_exists).CellStatsExist]'==1).*double(replica_num_vector==replica_num(ik)).*contains(existing_files,unique_concentrations{ij}).*contains(lower(existing_files),SoftwareNames_Spelling{ii})==1);
                    end
                    softwareResults.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).ind_results = [results(software_conditionQuantStatsExists_matches{ii,ij})];
                    softwareResults.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).ind_results_part_function = @(results_part_file) [results_part_file(software_conditionQuantStatsExists_matches{ii,ij})];
                    softwareResults.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).ind_results_part_spots_function = @(results_part_file) [results_part_file(software_conditionSpotStatsExists_matches{ii,ij})];
                    softwareResults.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).ind_results_part_cells_function = @(results_part_file) [results_part_file(software_conditionCellStatsExists_matches{ii,ij})];
                end
            end
        end
        save('pooledIndex_hARF4_RNAFISH_Analysis.mat','noProbeResults','softwareResults','-v7.3');
    end  
    if (part==1)
        load('hARF4_RNAFISH_AnalysisQuantResults_currentHidden.mat','results_quant');
        results_quantA = results_quant;clear results_quant
        load('hARF4_RNAFISH_AnalysisQuantResults_currentMain.mat','results_quant');
        results_quantB = results_quant;clear results_quant 
        if (~isempty(fieldnames(results_quantA)))
            results_quant = [results_quantA results_quantB];
        else
            results_quant = results_quantB;
        end
        clear results_quantA results_quantB
    end
    if (part==2)
        load('hARF4_RNAFISH_AnalysisCellStats_currentHidden.mat','results_cellStats');
        results_cellStatsA = results_cellStats;clear results_cellStats
        load('hARF4_RNAFISH_AnalysisCellStats_currentMain.mat','results_cellStats');
        results_cellStatsB = results_cellStats;clear results_cellStats
        if (~isempty(fieldnames(results_cellStatsA)))
            results_cellStats = [results_cellStatsA results_cellStatsB];
        else
            results_cellStats = results_cellStatsB;
        end
        clear results_cellStatsA results_cellStatsB
        load('hARF4_RNAFISH_AnalysisSpotStats_currentHidden.mat','results_spotStats');
        results_spotStatsA = results_spotStats;clear results_spotStats
        load('hARF4_RNAFISH_AnalysisSpotStats_currentMain.mat','results_spotStats');
        results_spotStatsB = results_spotStats;clear results_spotStats
        if (~isempty(fieldnames(results_spotStatsA)))
            results_spotStats = [results_spotStatsA results_spotStatsB];
        else
            results_spotStats = results_spotStatsB;
        end
        clear results_spotStatsA results_spotStatsB
    end
    if (part==3)
        load('hARF4_RNAFISH_imageThresholdsStats_currentHidden.mat','results_ImageThresholdsStats');
        results_ImageThresholdsStatsA = results_ImageThresholdsStats;clear results_ImageThresholdsStats
        load('hARF4_RNAFISH_imageThresholdsStats_currentMain.mat','results_ImageThresholdsStats');
        results_ImageThresholdsStatsB = results_ImageThresholdsStats;clear results_ImageThresholdsStats
        if (~isempty(fieldnames(results_ImageThresholdsStatsA)))
            results_ImageThresholdsStats = [results_ImageThresholdsStatsA results_ImageThresholdsStatsB];
        else
            results_ImageThresholdsStats = results_ImageThresholdsStatsB;
        end
        clear results_ImageThresholdsStatsA results_ImageThresholdsStatsB
        load('hARF4_RNAFISH_CellThresholdStats_currentHidden.mat','results_CellThresholdsStats');
        results_CellThresholdsStatsA = results_CellThresholdsStats;clear results_CellThresholdsStats
        load('hARF4_RNAFISH_CellThresholdStats_currentMain.mat','results_CellThresholdsStats');
        results_CellThresholdsStatsB = results_CellThresholdsStats;clear results_CellThresholdsStats
        if (~isempty(fieldnames(results_CellThresholdsStatsA)))
            results_CellThresholdsStats = [results_CellThresholdsStatsA results_CellThresholdsStatsB];
        else
            results_CellThresholdsStats = results_CellThresholdsStatsB;
        end
        clear results_CellThresholdsStatsA results_CellThresholdsStatsB
        load('hARF4_RNAFISH_AnalysisCallTables_currentHidden.mat','results_callTable');
        results_callTableA = results_callTable;clear results_callTable
        load('hARF4_RNAFISH_AnalysisCallTables_currentMain.mat','results_callTable');
        results_callTableB = results_callTable;clear results_callTable
        if (~isempty(fieldnames(results_callTableA)))
            results_callTable = [results_callTableA results_callTableB];
        else
            results_callTable = results_callTableB;
        end
        clear results_callTableA results_callTableB
    end
    if (part==4)
        load('hARF4_RNAFISH_CellSegStats_currentHidden.mat','results_cellseg');
        results_cellsegA = results_cellseg;clear results_cellseg
        load('hARF4_RNAFISH_CellSegStats_currentMain.mat','results_cellseg');
        results_cellsegB = results_cellseg;clear results_cellseg
        if (~isempty(fieldnames(results_cellsegA)))
            results_cellseg = [results_cellsegA results_cellsegB];
        else
            results_cellseg = results_cellsegB;
        end
        clear results_cellsegA results_cellsegB
        load('hARF4_RNAFISH_AnalysisCellStats_currentHidden.mat','results_cellStats');
        results_cellStatsA = results_cellStats;clear results_cellStats
        load('hARF4_RNAFISH_AnalysisCellStats_currentMain.mat','results_cellStats');
        results_cellStatsB = results_cellStats;clear results_cellStats
        if (~isempty(fieldnames(results_cellStatsA)))
            results_cellStats = [results_cellStatsA results_cellStatsB];
        else
            results_cellStats = results_cellStatsB;
        end
        clear results_cellStatsA results_cellStatsB
        load('hARF4_RNAFISH_AnalysisImgData_currentHidden.mat','results_imgData');
        results_imgDataA = results_imgData;clear results_imgData
        load('hARF4_RNAFISH_AnalysisImgData_currentMain.mat','results_imgData');
        results_imgDataB = results_imgData;clear results_imgData
        if (~isempty(fieldnames(results_imgDataA)))
            results_imgData = [results_imgDataA results_imgDataB];
        else
            results_imgData = results_imgDataB;
        end
        clear results_imgDataA results_imgDataB
    end
    noProbePooledResults_quant = struct();
    noProbePooledResults_stats = struct();
    noProbePooledResults_thresholdInfo = struct();
    noProbePooledResults_callTable = struct();
    noProbePooledResults_cellseg = struct();
    noProbePooledResults_imgData = struct();
    %% combine individual no_probe replica result structs together into single structures
    for ik = 1:length(Replica)
        n_reps = length(noProbeResults.(Replica{ik}).ind_results);
        if (n_reps>0)
            if (part==1)
                temp_quant_results = noProbeResults.(Replica{ik}).ind_results_part_function(results_quant);
                existing_quant = find([temp_quant_results.QuantStatsExist]==1);
                num_cells_withoutspots = [temp_quant_results(existing_quant).num_cells_withoutspots];
                struct_quant_results_spotTable = vertcat(temp_quant_results(existing_quant).quant_results_spotTable);
                if (~isempty(struct_quant_results_spotTable))
                    image_number_spotQuantTable = cell2mat(arrayfun(@(x) x*ones(1,size(temp_quant_results(x).quant_results_spotTable,1)),existing_quant,'Un',0))';
                    struct_quant_results_spotTable(:,'image') = array2table(image_number_spotQuantTable);
                end
                noProbePooledResults_quant.(Replica{ik}).num_cells_withoutspots = num_cells_withoutspots;
                noProbePooledResults_quant.(Replica{ik}).struct_quant_results_spotTable = struct_quant_results_spotTable;
                spotQuantTable_in_cell_X_in_image_Y = @(x,y) find(double(table2array(struct_quant_results_spotTable(:,'cell'))==x).*double(table2array(struct_quant_results_spotTable(:,'image'))==y)==1);
                spotQuantTable_nas_nuc_cyt_type_in_cell_X_in_image_Y = @(x,y,ai,bi,ci) find(double(table2array(struct_quant_results_spotTable(:,'cell'))==x).*...
                    double(table2array(struct_quant_results_spotTable(:,'image'))==y).*...
                    double(table2array(struct_quant_results_spotTable(:,'nascent_flag'))==ai).*...
                    double(table2array(struct_quant_results_spotTable(:,'nucRNA'))==bi).*...
                    double(table2array(struct_quant_results_spotTable(:,'cytoRNA'))==ci)==1);
                check_spotQuantTable_if_spot_nas_nuc_cyt_type = @(W,x,y,z) find(double(table2array(struct_quant_results_spotTable(W,'nascent_flag'))==x).*...
                    double(table2array(struct_quant_results_spotTable(W,'nucRNA'))==y).*...
                    double(table2array(struct_quant_results_spotTable(W,'cytoRNA'))==z)==1);
                check_spotQuantTable_meetThreshold = @(W,Thr) find(double(table2array(struct_quant_results_spotTable(W,'dropout_thresh'))>Thr)==1);
                get_spotQuantTable_meetThreshold = @(W,Thr) W(double(table2array(struct_quant_results_spotTable(W,'dropout_thresh'))>Thr)==1);
                getSpotQuantTableMeetingThresholds = @(imageList,cellList,ThrList,ai,bi,ci) arrayfun(@(x) get_spotQuantTable_meetThreshold(spotQuantTable_nas_nuc_cyt_type_in_cell_X_in_image_Y(cellList(x),imageList(x),ai,bi,ci),ThrList(x))',1:length(ThrList),'Un',0);
                numSpotQuantTable_at_Thr_in_cell_image = @(imageList,cellList,ThrList) CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(ci) CATnWrapper(arrayfun(@(bi) arrayfun(@(ai) length(check_spotQuantTable_meetThreshold(check_spotQuantTable_if_spot_nas_nuc_cyt_type(spotQuantTable_in_cell_X_in_image_Y(cellList(x),imageList(x)),ai,bi,ci),ThrList(x))),0:1)',0:1,'Un',0),2),0:1,'Un',0),3),1:length(imageList),'Un',0),4);
                numSpotQuantTable_max_at_cell_image = @(imageList,cellList) CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(ci) CATnWrapper(arrayfun(@(bi) arrayfun(@(ai) length(check_spotQuantTable_if_spot_nas_nuc_cyt_type(spotQuantTable_in_cell_X_in_image_Y(cellList(x),imageList(x)),ai,bi,ci)),0:1),0:1,'Un',0),1),0:1,'Un',0),3),1:length(imageList),'Un',0),4);
                noProbePooledResults_quant.(Replica{ik}).numSpotQuantTable_at_Thr_in_cell_image = numSpotQuantTable_at_Thr_in_cell_image;
                noProbePooledResults_quant.(Replica{ik}).numSpotQuantTable_max_at_cell_image = numSpotQuantTable_max_at_cell_image;
                noProbePooledResults_quant.(Replica{ik}).getSpotQuantTableMeetingThresholds = getSpotQuantTableMeetingThresholds;
            end
            if (part==2)
                temp_ImageCellStats = noProbeResults.(Replica{ik}).ind_results_part_cells_function(results_cellStats);
                existing_cellStats = find([temp_ImageCellStats.CellStatsExist]==1);
                struct_ImageCellStats = vertcat(temp_ImageCellStats(existing_cellStats).ImageCellStats);
                if (~isempty(struct_ImageCellStats))
                    image_number_cellStatTable = cell2mat(arrayfun(@(x) x*ones(1,size(temp_ImageCellStats(x).ImageCellStats,1)),existing_cellStats,'Un',0))';
                    struct_ImageCellStats(:,'image') = array2table(image_number_cellStatTable);
                end
                noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats = struct_ImageCellStats;
                existing_cellStats_IntCorrected = find([temp_ImageCellStats.CellStatsExist_IntCorrected]==1);
                struct_ImageCellStats_IntCorrected = vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).ImageCellStats_IntCorrected);
                if (~isempty(struct_ImageCellStats_IntCorrected))
                    image_number_cellStatTable_IntCorrected = cell2mat(arrayfun(@(x) x*ones(1,size(temp_ImageCellStats(x).ImageCellStats_IntCorrected,1)),existing_cellStats,'Un',0))';
                    struct_ImageCellStats_IntCorrected(:,'image') = array2table(image_number_cellStatTable_IntCorrected);
                end
                noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats_IntCorrected = struct_ImageCellStats_IntCorrected;
                % if (length(existing_cellStats)>1)
                % struct_ImageCamBkgFullStats = struct2table(vertcat(temp_ImageCellStats(existing_cellStats).image_cambkg_Full));
                % else
                % struct_ImageCamBkgFullStats = struct2table(vertcat(temp_ImageCellStats(existing_cellStats).image_cambkg_Full),'AsArray',1);    
                % end
                % if (~isempty(struct_ImageCamBkgFullStats))
                %     image_number_CamBkgFullStats = cell2mat(arrayfun(@(x) x,existing_cellStats,'Un',0))';
                %     struct_ImageCamBkgFullStats(:,'image') = array2table(image_number_CamBkgFullStats);
                % end
                % noProbePooledResults_stats.(Replica{ik}).struct_ImageCamBkgFullStats = struct_ImageCamBkgFullStats;
                % if (length(existing_cellStats)>1)
                % struct_ImageCamBkgTrimStats = struct2table(vertcat(temp_ImageCellStats(existing_cellStats).image_cambkg_Trim));
                % else
                % struct_ImageCamBkgTrimStats = struct2table(vertcat(temp_ImageCellStats(existing_cellStats).image_cambkg_Trim),'AsArray',1);    
                % end
                % if (~isempty(struct_ImageCamBkgTrimStats))
                %     image_number_CamBkgTrimStats = cell2mat(arrayfun(@(x) x,existing_cellStats,'Un',0))';
                %     struct_ImageCamBkgTrimStats(:,'image') = array2table(image_number_CamBkgTrimStats);
                % end
                % noProbePooledResults_stats.(Replica{ik}).struct_ImageCamBkgTrimStats = struct_ImageCamBkgTrimStats;
                % if (length(existing_cellStats)>1)
                % struct_ImageCellBkgFullStats = struct2table(vertcat(temp_ImageCellStats(existing_cellStats).image_cellbkg_Full));
                % else
                % struct_ImageCellBkgFullStats = struct2table(vertcat(temp_ImageCellStats(existing_cellStats).image_cellbkg_Full),'AsArray',1);    
                % end
                % if (~isempty(struct_ImageCellBkgFullStats))
                %     image_number_CellBkgFullStats = cell2mat(arrayfun(@(x) x,existing_cellStats,'Un',0))';
                %     struct_ImageCellBkgFullStats(:,'image') = array2table(image_number_CellBkgFullStats);
                % end
                % noProbePooledResults_stats.(Replica{ik}).struct_ImageCellBkgFullStats = struct_ImageCellBkgFullStats;
                % if (length(existing_cellStats)>1)
                % struct_ImageCellBkgTrimStats = struct2table(vertcat(temp_ImageCellStats(existing_cellStats).image_cellbkg_Trim));
                % else
                % struct_ImageCellBkgTrimStats = struct2table(vertcat(temp_ImageCellStats(existing_cellStats).image_cellbkg_Trim),'AsArray',1);    
                % end
                % if (~isempty(struct_ImageCellBkgTrimStats))
                %     image_number_CellBkgTrimStats = cell2mat(arrayfun(@(x) x,existing_cellStats,'Un',0))';
                %     struct_ImageCellBkgTrimStats(:,'image') = array2table(image_number_CellBkgTrimStats);
                % end
                % noProbePooledResults_stats.(Replica{ik}).struct_ImageCellBkgTrimStats = struct_ImageCellBkgTrimStats; 
                % if (length(existing_cellStats_IntCorrected)>1)
                % struct_ImageCamBkgFullStats_IntCorrected = struct2table(vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).image_cambkg_Full_IntCorrected));
                % else
                % struct_ImageCamBkgFullStats_IntCorrected = struct2table(vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).image_cambkg_Full_IntCorrected),'AsArray',1);    
                % end
                % if (~isempty(struct_ImageCamBkgFullStats_IntCorrected))
                %     image_number_CamBkgFullStats_IntCorrected = cell2mat(arrayfun(@(x) x,existing_cellStats_IntCorrected,'Un',0))';
                %     struct_ImageCamBkgFullStats_IntCorrected(:,'image') = array2table(image_number_CamBkgFullStats_IntCorrected);
                % end
                % noProbePooledResults_stats.(Replica{ik}).struct_ImageCamBkgFullStats_IntCorrected = struct_ImageCamBkgFullStats_IntCorrected;
                % if (length(existing_cellStats_IntCorrected)>1)
                % struct_ImageCamBkgTrimStats_IntCorrected = struct2table(vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).image_cambkg_Trim_IntCorrected));
                % else
                % struct_ImageCamBkgTrimStats_IntCorrected = struct2table(vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).image_cambkg_Trim_IntCorrected),'AsArray',1);    
                % end
                % if (~isempty(struct_ImageCamBkgTrimStats_IntCorrected))
                %     image_number_CamBkgTrimStats_IntCorrected = cell2mat(arrayfun(@(x) x,existing_cellStats_IntCorrected,'Un',0))';
                %     struct_ImageCamBkgTrimStats_IntCorrected(:,'image') = array2table(image_number_CamBkgTrimStats_IntCorrected);
                % end
                % noProbePooledResults_stats.(Replica{ik}).struct_ImageCamBkgTrimStats_IntCorrected = struct_ImageCamBkgTrimStats_IntCorrected;
                % if (length(existing_cellStats_IntCorrected)>1)
                % struct_ImageCellBkgFullStats_IntCorrected = struct2table(vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).image_cellbkg_Full_IntCorrected));
                % else
                % struct_ImageCellBkgFullStats_IntCorrected = struct2table(vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).image_cellbkg_Full_IntCorrected),'AsArray',1);    
                % end
                % if (~isempty(struct_ImageCellBkgFullStats_IntCorrected))
                %     image_number_CellBkgFullStats_IntCorrected = cell2mat(arrayfun(@(x) x,existing_cellStats_IntCorrected,'Un',0))';
                %     struct_ImageCellBkgFullStats_IntCorrected(:,'image') = array2table(image_number_CellBkgFullStats_IntCorrected);
                % end
                % noProbePooledResults_stats.(Replica{ik}).struct_ImageCellBkgFullStats_IntCorrected = struct_ImageCellBkgFullStats_IntCorrected;
                % if (length(existing_cellStats_IntCorrected)>1)
                % struct_ImageCellBkgTrimStats_IntCorrected = struct2table(vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).image_cellbkg_Trim_IntCorrected));
                % else
                % struct_ImageCellBkgTrimStats_IntCorrected = struct2table(vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).image_cellbkg_Trim_IntCorrected),'AsArray',1);    
                % end
                % if (~isempty(struct_ImageCellBkgTrimStats_IntCorrected))
                %     image_number_CellBkgTrimStats_IntCorrected = cell2mat(arrayfun(@(x) x,existing_cellStats_IntCorrected,'Un',0))';
                %     struct_ImageCellBkgTrimStats_IntCorrected(:,'image') = array2table(image_number_CellBkgTrimStats_IntCorrected);
                % end
                % noProbePooledResults_stats.(Replica{ik}).struct_ImageCellBkgTrimStats_IntCorrected = struct_ImageCellBkgTrimStats_IntCorrected;
                temp_ImageSpotStats = noProbeResults.(Replica{ik}).ind_results_part_spots_function(results_spotStats);
                existing_spotStats = find([temp_ImageSpotStats.SpotStatsExist]==1);
                struct_ImageSpotStats = vertcat(temp_ImageSpotStats(existing_spotStats).ImageSpotStats);
                if (~isempty(struct_ImageSpotStats))
                    image_number_spotStatTable = cell2mat(arrayfun(@(x) x*ones(1,size(temp_ImageSpotStats(x).ImageSpotStats,1)),existing_spotStats,'Un',0))';
                    struct_ImageSpotStats(:,'image') = array2table(image_number_spotStatTable);
                end
                noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats = struct_ImageSpotStats;
                existing_spotStats_IntCorrected = find([temp_ImageSpotStats.SpotStatsExist_IntCorrected]==1);
                struct_ImageSpotStats_IntCorrected = vertcat(temp_ImageSpotStats(existing_spotStats_IntCorrected).ImageSpotStats_IntCorrected);
                if (~isempty(struct_ImageSpotStats_IntCorrected))
                    image_number_spotStatTable_IntCorrected = cell2mat(arrayfun(@(x) x*ones(1,size(temp_ImageSpotStats(x).ImageSpotStats_IntCorrected,1)),existing_spotStats_IntCorrected,'Un',0))';
                    struct_ImageSpotStats_IntCorrected(:,'image') = array2table(image_number_spotStatTable_IntCorrected);
                end
                noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats_IntCorrected = struct_ImageSpotStats_IntCorrected;
                spotStatTable_nas_nuc_cyt_type_in_cell_X_in_image_Y = @(x,y,ai,bi,ci) find(double(table2array(struct_ImageSpotStats(:,'cell'))==x).*...
                    double(table2array(struct_ImageSpotStats(:,'image'))==y).*...
                    double(table2array(struct_ImageSpotStats(:,'nascent_flag'))==ai).*...
                    double(table2array(struct_ImageSpotStats(:,'nuc_flag'))==bi).*...
                    double(table2array(struct_ImageSpotStats(:,'nuc_flag'))==1-ci)==1);
                check_spotStatTable_if_spot_nas_nuc_cyt_type = @(W,x,y,z) find(double(table2array(struct_ImageSpotStats(W,'nascent_flag'))==x).*...
                    double(table2array(struct_ImageSpotStats(W,'nuc_flag'))==y).*...
                    double(table2array(struct_ImageSpotStats(W,'nuc_flag'))==1-z)==1);
                check_spotStatTable_meetThreshold = @(W,Thr) find(double(table2array(struct_ImageSpotStats(W,'dropout_thresh'))>Thr)==1);
                get_spotStatTable_meetThreshold = @(W,Thr) W(double(table2array(struct_ImageSpotStats(W,'dropout_thresh'))>Thr)==1);
                getSpotStatTableMeetingThresholds = @(imageList,cellList,ThrList,ai,bi,ci) arrayfun(@(x) get_spotStatTable_meetThreshold(spotStatTable_nas_nuc_cyt_type_in_cell_X_in_image_Y(cellList(x),imageList(x),ai,bi,ci),ThrList(x))',1:length(ThrList),'Un',0);
                numSpotStatTable_at_Thr_in_cell_image = @(imageList,cellList,ThrList) CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(ci) CATnWrapper(arrayfun(@(bi) arrayfun(@(ai) length(check_spotStatTable_meetThreshold(check_spotStatTable_if_spot_nas_nuc_cyt_type(spotStatTable_in_cell_X_in_image_Y(cellList(x),imageList(x)),ai,bi,ci),ThrList(x))),0:1)',0:1,'Un',0),2),0:1,'Un',0),3),1:length(imageList),'Un',0),4);
                numSpotStatTable_max_at_cell_image = @(imageList,cellList) CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(ci) CATnWrapper(arrayfun(@(bi) arrayfun(@(ai) length(check_spotStatTable_if_spot_nas_nuc_cyt_type(spotStatTable_in_cell_X_in_image_Y(cellList(x),imageList(x)),ai,bi,ci)),0:1),0:1,'Un',0),1),0:1,'Un',0),3),1:length(imageList),'Un',0),4);
                noProbePooledResults_stats.(Replica{ik}).numSpotStatTable_at_Thr_in_cell_image = numSpotStatTable_at_Thr_in_cell_image;
                noProbePooledResults_stats.(Replica{ik}).numSpotStatTable_max_at_cell_image = numSpotStatTable_max_at_cell_image;
                noProbePooledResults_stats.(Replica{ik}).getSpotStatTableMeetingThresholds = getSpotStatTableMeetingThresholds;
                if (~isempty(noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats))
                noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats(:,{'max_spot_int'}) = ...
                    array2table(cellfun(@(x) double(x.max),noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats.noBkgSubStats));
                noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats(:,{'mean_bkg_int'}) = ...
                    array2table(cellfun(@(x) double(x.mean),noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats.localBkgStats));
                noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats(:,{'Sig_minus_bkg'}) = ...
                    array2table(cellfun(@(x) double(x.max),noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats.noBkgSubStats)) - array2table(cellfun(@(x) double(x.mean),noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats.localBkgStats));
                end
                if (~isempty(noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats))
                noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats(:,{'mean_cell_bkg'}) = ...
                    array2table(cellfun(@(x) double(x.mean),noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats.cell_bkg));
                noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats(:,{'mean_cam_bkg'}) = ...
                    array2table(cellfun(@(x) double(x.mean),noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats.local_img_bkg));
                noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats(:,{'mean_cell_minus_cam_bkg'}) = ...
                    array2table(cellfun(@(x) double(x.mean),noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats.cell_bkg)) - array2table(cellfun(@(x) double(x.mean),noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats.local_img_bkg));
                end
                if (~isempty(noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats_IntCorrected))
                noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats_IntCorrected(:,{'max_spot_int'}) = ...
                    array2table(cellfun(@(x) double(x.max),noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats_IntCorrected.noBkgSubStats));
                noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats_IntCorrected(:,{'mean_bkg_int'}) = ...
                    array2table(cellfun(@(x) double(x.mean),noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats_IntCorrected.localBkgStats));
                noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats_IntCorrected(:,{'Sig_minus_bkg'}) = ...
                    array2table(cellfun(@(x) double(x.max),noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats_IntCorrected.noBkgSubStats)) - array2table(cellfun(@(x) double(x.mean),noProbePooledResults_stats.(Replica{ik}).struct_ImageSpotStats_IntCorrected.localBkgStats));
                end
                if (~isempty(noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats_IntCorrected))
                noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats_IntCorrected(:,{'mean_cell_bkg'}) = ...
                    array2table(cellfun(@(x) double(x.mean),noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats_IntCorrected.cell_bkg));
                noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats_IntCorrected(:,{'mean_cam_bkg'}) = ...
                    array2table(cellfun(@(x) double(x.mean),noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats_IntCorrected.local_img_bkg));
                noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats_IntCorrected(:,{'mean_cell_minus_cam_bkg'}) = ...
                    array2table(cellfun(@(x) double(x.mean),noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats_IntCorrected.cell_bkg)) - array2table(cellfun(@(x) double(x.mean),noProbePooledResults_stats.(Replica{ik}).struct_ImageCellStats_IntCorrected.local_img_bkg));
                spotStatTable_in_cell_X_in_image_Y = @(x,y) find(double(table2array(struct_ImageSpotStats(:,'cell'))==x).*double(table2array(struct_ImageSpotStats(:,'image'))==y)==1);
                end
            end
            if (part==3)
                %% threshold candidate info for each image
                n_reps = size(noProbeResults.(Replica{ik}).ind_results_part_cells_function(results_ImageThresholdsStats),2);
                temp_ImageThresholdsStats = noProbeResults.(Replica{ik}).ind_results_part_cells_function(results_ImageThresholdsStats);
                image_sugg_m_table = cell2mat(arrayfun(@(x) temp_ImageThresholdsStats(x).sugg_m',1:n_reps,'Un',0))';
                image_threshold_mean = [temp_ImageThresholdsStats.mean_overall];
                image_threshold_med = [temp_ImageThresholdsStats.med_overall];
                image_threshold_std = [temp_ImageThresholdsStats.std_overall];
                image_threshold_min = [temp_ImageThresholdsStats.min_overall];
                image_threshold_max = [temp_ImageThresholdsStats.max_overall];
                image_threshold_opt =  [temp_ImageThresholdsStats.threshold];
                image_threshold_spot_count_curve_x = cell2mat(arrayfun(@(x) temp_ImageThresholdsStats(x).threshold_spot_count_curve.x',1:n_reps,'Un',0))';
                image_threshold_spot_count_curve_y = cell2mat(arrayfun(@(x) temp_ImageThresholdsStats(x).threshold_spot_count_curve.y',1:n_reps,'Un',0))';
                image_threshold_score_curve_x = cell2mat(arrayfun(@(x) temp_ImageThresholdsStats(x).threshold_score_curve.x',1:n_reps,'Un',0))';
                image_threshold_score_curve_score = cell2mat(arrayfun(@(x) temp_ImageThresholdsStats(x).threshold_score_curve.score',1:n_reps,'Un',0))';
                image_threshold_score_curve_score_m = cell2mat(arrayfun(@(x) temp_ImageThresholdsStats(x).threshold_score_curve.score_m',1:n_reps,'Un',0))';
                image_threshold_score_curve_score_f = cell2mat(arrayfun(@(x) temp_ImageThresholdsStats(x).threshold_score_curve.score_f',1:n_reps,'Un',0))';
                image_threshold_score_curve_score_fri = cell2mat(arrayfun(@(x) temp_ImageThresholdsStats(x).threshold_score_curve.score_fri',1:n_reps,'Un',0))';
                noProbePooledResults_thresholdInfo.(Replica{ik}).image_sugg_m_table = image_sugg_m_table;
                noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_mean = image_threshold_mean;
                noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_med = image_threshold_med;
                noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_std = image_threshold_std;
                noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_min = image_threshold_min;
                noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_max = image_threshold_max;
                noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_opt = image_threshold_opt;
                noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_spot_count_curve_x = image_threshold_spot_count_curve_x;
                noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_spot_count_curve_y = image_threshold_spot_count_curve_y;
                noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_score_curve_x = image_threshold_score_curve_x;
                noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_score_curve_score = image_threshold_score_curve_score;
                noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_score_curve_score_m = image_threshold_score_curve_score_m;
                noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_score_curve_score_f = image_threshold_score_curve_score_f;
                noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_score_curve_score_fri = image_threshold_score_curve_score_fri;
                %% threshold info for single-cells
                temp_CellThresholdsStats = noProbeResults.(Replica{ik}).ind_results_part_cells_function(results_CellThresholdsStats);
                existing_cellThresholdsStats = find([temp_CellThresholdsStats.CellThresholdsExist]==1);
                cell_sugg_m_table = cell2mat(arrayfun(@(x) temp_CellThresholdsStats(x).scThresholdSuggestions.sugg_m',existing_cellThresholdsStats,'Un',0))';
                cell_threshold_mean = cell2mat(arrayfun(@(x) temp_CellThresholdsStats(x).scThresholdSuggestions.mean_overall',existing_cellThresholdsStats,'Un',0))';
                cell_threshold_med = cell2mat(arrayfun(@(x) temp_CellThresholdsStats(x).scThresholdSuggestions.med_overall',existing_cellThresholdsStats,'Un',0))';
                cell_threshold_std = cell2mat(arrayfun(@(x) temp_CellThresholdsStats(x).scThresholdSuggestions.std_overall',existing_cellThresholdsStats,'Un',0))';
                cell_threshold_min = cell2mat(arrayfun(@(x) temp_CellThresholdsStats(x).scThresholdSuggestions.min_overall',existing_cellThresholdsStats,'Un',0))';
                cell_threshold_max = cell2mat(arrayfun(@(x) temp_CellThresholdsStats(x).scThresholdSuggestions.max_overall',existing_cellThresholdsStats,'Un',0))';
                cell_threshold_opt = cell2mat(arrayfun(@(x) temp_CellThresholdsStats(x).scThresholdSuggestions.threshold',existing_cellThresholdsStats,'Un',0))';
                cell_threshold_spot_count_curve_x = cell2mat(arrayfun(@(v) [temp_CellThresholdsStats(v).scThresholdSuggestions.x{:}],existing_cellThresholdsStats,'Un',0))';
                cell_threshold_spot_count_curve_y = cell2mat(arrayfun(@(x) [temp_CellThresholdsStats(x).scThresholdSuggestions.spot_counts{:}],existing_cellThresholdsStats,'Un',0))';
                cell_threshold_score_curve_x = cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(w) temp_CellThresholdsStats(x).scThresholdSuggestions.value_table{w}.threshold_value, ...
                    1:size(temp_CellThresholdsStats(x).scThresholdSuggestions.value_table,1),'Un',0)),existing_cellThresholdsStats,'Un',0))';
                cell_threshold_score_curve_score = cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(w) temp_CellThresholdsStats(x).scThresholdSuggestions.value_table{w}.score, ...
                    1:size(temp_CellThresholdsStats(x).scThresholdSuggestions.value_table,1),'Un',0)),existing_cellThresholdsStats,'Un',0))';
                cell_threshold_score_curve_score_m = cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(w) temp_CellThresholdsStats(x).scThresholdSuggestions.value_table{w}.score_m, ...
                    1:size(temp_CellThresholdsStats(x).scThresholdSuggestions.value_table,1),'Un',0)),existing_cellThresholdsStats,'Un',0))';
                cell_threshold_score_curve_score_f = cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(w) temp_CellThresholdsStats(x).scThresholdSuggestions.value_table{w}.score_f, ...
                    1:size(temp_CellThresholdsStats(x).scThresholdSuggestions.value_table,1),'Un',0)),existing_cellThresholdsStats,'Un',0))';
                cell_threshold_score_curve_score_fri = cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(w) temp_CellThresholdsStats(x).scThresholdSuggestions.value_table{w}.score_fri, ...
                    1:size(temp_CellThresholdsStats(x).scThresholdSuggestions.value_table,1),'Un',0)),existing_cellThresholdsStats,'Un',0))';
                noProbePooledResults_thresholdInfo.(Replica{ik}).cell_sugg_m_table = cell_sugg_m_table;
                noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_mean = cell_threshold_mean;
                noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_med = cell_threshold_med;
                noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_std = cell_threshold_std;
                noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_min = cell_threshold_min;
                noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_max = cell_threshold_max;
                noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_opt = cell_threshold_opt;
                noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_spot_count_curve_x = cell_threshold_spot_count_curve_x;
                noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_spot_count_curve_y = cell_threshold_spot_count_curve_y;
                noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_score_curve_x = cell_threshold_score_curve_x;
                noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_score_curve_score = cell_threshold_score_curve_score;
                noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_score_curve_score_m = cell_threshold_score_curve_score_m;
                noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_score_curve_score_f = cell_threshold_score_curve_score_f;
                noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_score_curve_score_fri = cell_threshold_score_curve_score_fri;
                temp_callTable = noProbeResults.(Replica{ik}).ind_results_part_cells_function(results_callTable);
                struct_callTable = vertcat(temp_callTable(:).call_table);
                if (~isempty(struct_callTable))
                    image_number_callTable = cell2mat(arrayfun(@(x) x*ones(1,size(temp_callTable(x).call_table,1)),1:n_reps,'Un',0))';
                    struct_callTable(:,'image') = array2table(image_number_callTable);
                end
                in_image_Y = @(y) find(double(table2array(struct_callTable(:,'image'))==y)==1);
                no_cell_in_image_Y = @(y) find(double(table2array(struct_callTable(:,'cell'))==0).*...
                    double(table2array(struct_callTable(:,'image'))==y)==1);
                in_cells_in_image_Y = @(y) find(double(table2array(struct_callTable(:,'cell'))>0).*...
                    double(table2array(struct_callTable(:,'image'))==y)==1);
                in_cell_X_in_image_Y = @(x,y) find(double(table2array(struct_callTable(:,'cell'))==x).*...
                    double(table2array(struct_callTable(:,'image'))==y)==1);
                extended_callTable_spot_count_curve_in_image_Y_y = cell2mat(arrayfun(@(Y) ...
                    numel(in_image_Y(Y))-histcounts(table2array(struct_callTable(in_image_Y(Y),'dropout_thresh')),'Normalization','cumcount','BinEdges',[-Inf 1:size(image_threshold_spot_count_curve_x,2)])',1:n_reps,'Un',0))';
                unique_image_and_cell_pair = unique(double(table2array(struct_callTable(find(double(table2array(struct_callTable(:,'cell'))>0)),{'image','cell'}))),'rows');
                extended_callTable_spot_count_curve_no_cell_in_image_Y_y = cell2mat(arrayfun(@(Y) ...
                    numel(no_cell_in_image_Y(Y))-histcounts(table2array(struct_callTable(no_cell_in_image_Y(Y),'dropout_thresh')),'Normalization','cumcount','BinEdges',[-Inf 1:size(image_threshold_spot_count_curve_x,2)])',1:n_reps,'Un',0))';
                extended_callTable_spot_count_curve_in_cells_in_image_Y_y = cell2mat(arrayfun(@(Y) ...
                    numel(in_cells_in_image_Y(Y))-histcounts(table2array(struct_callTable(in_cells_in_image_Y(Y),'dropout_thresh')),'Normalization','cumcount','BinEdges',[-Inf 1:size(image_threshold_spot_count_curve_x,2)])',1:n_reps,'Un',0))';
                extended_callTable_spot_count_curve_in_cell_X_in_image_Y_y = cell2mat(arrayfun(@(Y) ...
                    numel(in_cell_X_in_image_Y(unique_image_and_cell_pair(Y,2),unique_image_and_cell_pair(Y,1)))-histcounts(table2array(struct_callTable(in_cell_X_in_image_Y(unique_image_and_cell_pair(Y,2),unique_image_and_cell_pair(Y,1)),'dropout_thresh')),'Normalization','cumcount','BinEdges',[-Inf 1:size(image_threshold_spot_count_curve_x,2)])',1:size(unique_image_and_cell_pair,1),'Un',0))';
                cell_spot_count_at_thrList_for_unique_image_and_cell_pair = @(cellList,ThrList) extended_callTable_spot_count_curve_in_cell_X_in_image_Y_y(sub2ind(...
                    size(extended_callTable_spot_count_curve_in_cell_X_in_image_Y_y),cellList,ceil(ThrList)));
                callTable_in_cell_X_in_image_Y = @(x,y) find(double(table2array(struct_callTable(:,'cell'))==x).*double(table2array(struct_callTable(:,'image'))==y)==1);
                check_callTable_meetThreshold = @(W,Thr) find(double(table2array(struct_callTable(W,'dropout_thresh'))>=Thr)==1);
                get_callTable_meetThreshold = @(W,Thr) W(double(table2array(struct_callTable(W,'dropout_thresh'))>=Thr)==1);
                getCallTableMeetingThresholds = @(imageList,cellList,ThrList) arrayfun(@(x) get_callTable_meetThreshold(callTable_in_cell_X_in_image_Y(cellList(x),imageList(x)),ThrList(x))',1:length(ThrList),'Un',0);
                numCallTable_at_Thr_in_cell_image = @(imageList,cellList,ThrList) arrayfun(@(x) length(check_callTable_meetThreshold(callTable_in_cell_X_in_image_Y(cellList(x),imageList(x)),ThrList(x))),1:length(imageList));
                is_main_spot = A_BHJH_ApplySpotMergingGeneric(struct_callTable,'isnap_x','isnap_y','isnap_z','intensity');
                struct_callTable(:,'is_main_spot') = array2table(is_main_spot);
                in_cell_X_in_image_Y_merged = @(x,y) find(double(table2array(struct_callTable(:,'cell'))==x).*...
                    double(table2array(struct_callTable(:,'image'))==y).*...
                    double(table2array(struct_callTable(:,'is_main_spot'))==1)==1);
                callTable_merged_spot_count_curve_in_cell_X_in_image_Y_y = cell2mat(arrayfun(@(Y) ...
                    numel(in_cell_X_in_image_Y_merged(unique_image_and_cell_pair(Y,2),unique_image_and_cell_pair(Y,1)))-histcounts(table2array(struct_callTable(in_cell_X_in_image_Y_merged(unique_image_and_cell_pair(Y,2),unique_image_and_cell_pair(Y,1)),'dropout_thresh')),'Normalization','cumcount','BinEdges',[-Inf 1:size(image_threshold_spot_count_curve_x,2)])',1:size(unique_image_and_cell_pair,1),'Un',0))';
                spotmerged_callTable_in_cell_X_in_image_Y= @(x,y) find(double(table2array(struct_callTable(:,'cell'))==x).*double(table2array(struct_callTable(:,'image'))==y).*double(table2array(struct_callTable(:,'is_main_spot'))==1)==1);
                spotmerged_check_callTable_meetThreshold = @(W,Thr) find(double(table2array(struct_callTable(W,'dropout_thresh'))>=Thr).*double(table2array(struct_callTable(W,'is_main_spot'))==1)==1);
                spotmerged_get_callTable_meetThresholds = @(W,Thr) W(double(table2array(struct_callTable(W,'dropout_thresh'))>=Thr).*double(table2array(struct_callTable(W,'is_main_spot'))==1)==1);
                spotmerged_numCallTable_at_Thr_in_cell_image = @(imageList,cellList,ThrList) arrayfun(@(x) length(spotmerged_check_callTable_meetThreshold(spotmerged_callTable_in_cell_X_in_image_Y(cellList(x),imageList(x)),ThrList(x))),1:length(imageList));
                noProbePooledResults_thresholdInfo.(Replica{ik}).unique_image_and_cell_pair = unique_image_and_cell_pair;
                noProbePooledResults_thresholdInfo.(Replica{ik}).extended_callTable_spot_count_curve_in_cell_X_in_image_Y_y = extended_callTable_spot_count_curve_in_cell_X_in_image_Y_y;
                noProbePooledResults_thresholdInfo.(Replica{ik}).numCallTable_at_Thr_in_cell_image = numCallTable_at_Thr_in_cell_image;
                noProbePooledResults_thresholdInfo.(Replica{ik}).callTable_merged_spot_count_curve_in_cell_X_in_image_Y_y = callTable_merged_spot_count_curve_in_cell_X_in_image_Y_y;
                noProbePooledResults_thresholdInfo.(Replica{ik}).spotmerged_numCallTable_at_Thr_in_cell_image = spotmerged_numCallTable_at_Thr_in_cell_image;
                noProbePooledResults_callTable.(Replica{ik}).struct_callTable = struct_callTable;
                noProbePooledResults_callTable.(Replica{ik}).unique_image_and_cell_pair = unique_image_and_cell_pair;
                noProbePooledResults_callTable.(Replica{ik}).extended_callTable_spot_count_curve_in_image_Y_y = extended_callTable_spot_count_curve_in_image_Y_y;
                noProbePooledResults_callTable.(Replica{ik}).extended_callTable_spot_count_curve_no_cell_in_image_Y_y = extended_callTable_spot_count_curve_no_cell_in_image_Y_y;
                noProbePooledResults_callTable.(Replica{ik}).extended_callTable_spot_count_curve_in_cells_in_image_Y_y = extended_callTable_spot_count_curve_in_cells_in_image_Y_y;
                noProbePooledResults_callTable.(Replica{ik}).extended_callTable_spot_count_curve_in_cell_X_in_image_Y_y = extended_callTable_spot_count_curve_in_cell_X_in_image_Y_y;
                noProbePooledResults_callTable.(Replica{ik}).cell_spot_count_at_thrList_for_unique_image_and_cell_pair = cell_spot_count_at_thrList_for_unique_image_and_cell_pair;
                noProbePooledResults_callTable.(Replica{ik}).numCallTable_at_Thr_in_cell_image = numCallTable_at_Thr_in_cell_image;
                noProbePooledResults_callTable.(Replica{ik}).getCallTableMeetingThresholds = getCallTableMeetingThresholds;
                noProbePooledResults_callTable.(Replica{ik}).callTable_merged_spot_count_curve_in_cell_X_in_image_Y_y = callTable_merged_spot_count_curve_in_cell_X_in_image_Y_y;
                noProbePooledResults_callTable.(Replica{ik}).spotmerged_numCallTable_at_Thr_in_cell_image = spotmerged_numCallTable_at_Thr_in_cell_image;
                noProbePooledResults_callTable.(Replica{ik}).spotmerged_getCallTableMeetingThresholds = spotmerged_get_callTable_meetThresholds;
            end
            if (part==4)
                subfields = {'cell_mask',...
                    'nucmask_max_z','nucmask_opt_z','nucmask_mid_z',...
                    'TRANS_max_z','TRANS_opt_z','TRANS_mid_z',...
                    'DAPI_max_z','DAPI_opt_z','DAPI_mid_z',...
                    'PROBE_max_z','PROBE_opt_z','PROBE_mid_z'};
                temp_ImageCellStats = noProbeResults.(Replica{ik}).ind_results_part_cells_function(results_cellStats);
                existing_cellStats = find([temp_ImageCellStats.CellStatsExist]==1);
                struct_ImageCellStats = vertcat(temp_ImageCellStats(existing_cellStats).ImageCellStats);
                try
                    temp_cellseg = noProbeResults.(Replica{ik}).ind_results_part_cells_function(results_cellseg);
                    temp_cellseg = vertcat(temp_cellseg(:).cellSeg);
                    if (size(temp_cellseg,1)>1)
                        struct_cellSeg = struct2table(temp_cellseg);
                    else
                        struct_cellSeg = struct2table(temp_cellseg,'AsArray',1);
                    end
                    if (~isempty(struct_cellSeg))
                        struct_cellSeg(:,'image') = array2table(existing_cellStats');
                    end
                    noProbePooledResults_cellseg.(Replica{ik}).struct_cellSeg = struct_cellSeg;
                catch e
                    disp(e)
                    for v = 1:length(e.stack)
                        e.stack(v)
                    end
                end
                temp_ImgData = noProbeResults.(Replica{ik}).ind_results_part_cells_function(results_imgData);
                if (size(temp_ImgData,2)>1)
                    struct_ImgData = struct2table(temp_ImgData);
                else
                    struct_ImgData = struct2table(temp_ImgData,'AsArray',1);
                end
                if (~isempty(struct_ImgData))
                    image_number_imgData = 1:size(struct_ImgData,1);
                    struct_ImgData(:,'image') = array2table(image_number_imgData');
                end
                noProbePooledResults_imgData.(Replica{ik}).struct_ImgData = struct_ImgData;
                if (~isempty(struct_ImageCellStats))
                    image_number_cellStatTable = cell2mat(arrayfun(@(x) x*ones(1,size(temp_ImageCellStats(x).ImageCellStats,1)),existing_cellStats,'Un',0))';
                    struct_ImageCellStats(:,'image') = array2table(image_number_cellStatTable);
                    length_Y = arrayfun(@(w) length([struct_ImageCellStats.box{w}.top:struct_ImageCellStats.box{w}.bottom]),1:size(struct_ImageCellStats,1));
                    length_X = arrayfun(@(w) length([struct_ImageCellStats.box{w}.left:struct_ImageCellStats.box{w}.right]),1:size(struct_ImageCellStats,1));
                    for f = 1:length(subfields)
                        try
                            struct_ImageCellStats(:,subfields{f}) = cell2table(arrayfun(@(w) struct_ImgData.(subfields{f}){...
                                struct_ImageCellStats.image(w)}(struct_ImageCellStats.box{w}.top:struct_ImageCellStats.box{w}.bottom,...
                                struct_ImageCellStats.box{w}.left:struct_ImageCellStats.box{w}.right),1:size(struct_ImageCellStats,1),'Un',0)');
                        catch
                        end
                    end
                end
                %have load's image (then stores intensity values in box
                %around each spot, so see raw image of each spot)
                noProbePooledResults_imgData.(Replica{ik}).struct_SingleCellImgData = struct_ImageCellStats;
            end
            %% function gives threshold curves for subsets
            %% function gives list cell counts at given threshold for each cell
            %gets position of spots meeting thresholds  3edc
            %number_at_thr(#nascent,#nuclear,$cytoplasmic)
            %number_max(#nascent,#nuclear,$cytoplasmic,cell)
            %% storing no probe main information
            %count and spot position information
            noProbePooledResults_quant.(Replica{ik}).n_reps = n_reps;
            noProbePooledResults_stats.(Replica{ik}).n_reps = n_reps;
            noProbePooledResults_thresholdInfo.(Replica{ik}).n_reps = n_reps;
            noProbePooledResults_callTable.(Replica{ik}).n_reps = n_reps;
            noProbePooledResults_imgData.(Replica{ik}).n_reps = n_reps;
        else
            noProbePooledResults_quant.(Replica{ik}).n_reps = 0;
            noProbePooledResults_stats.(Replica{ik}).n_reps = 0;
            noProbePooledResults_thresholdInfo.(Replica{ik}).n_reps = 0;
            noProbePooledResults_callTable.(Replica{ik}).n_reps = 0;
            noProbePooledResults_imgData.(Replica{ik}).n_reps = 0;
        end
    end
    %% combine individual probe replica results structs together into single structures
    softwarePooledResults_quant = struct();
    softwarePooledResults_stats = struct();
    softwarePooledResults_thresholdInfo = struct();
    softwarePooledResults_callTable = struct();
    softwarePooledResults_cellseg = struct();
    softwarePooledResults_imgData = struct();
    for ii = 1:length(SoftwareNames_Spelling)
        for ij = 1:length(unique_concentrations)
            for ik = 1:length(Replica)
                n_reps = length(softwareResults.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).ind_results);
                if (n_reps>0)
                    if (part==1)
                        temp_quant_results = softwareResults.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).ind_results_part_function(results_quant);
                        existing_quant = find([temp_quant_results.QuantStatsExist]==1);
                        num_cells_withoutspots = [temp_quant_results(existing_quant).num_cells_withoutspots];
                        struct_quant_results_spotTable = vertcat(temp_quant_results(existing_quant).quant_results_spotTable);
                        if (~isempty(struct_quant_results_spotTable))
                            image_number_spotQuantTable = cell2mat(arrayfun(@(x) x*ones(1,size(temp_quant_results(x).quant_results_spotTable,1)),existing_quant,'Un',0))';
                            struct_quant_results_spotTable(:,'image') = array2table(image_number_spotQuantTable);
                        end
                        softwarePooledResults_quant.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).num_cells_withoutspots = num_cells_withoutspots;
                        softwarePooledResults_quant.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_quant_results_spotTable = struct_quant_results_spotTable;
                        spotQuantTable_in_cell_X_in_image_Y = @(x,y) find(double(table2array(struct_quant_results_spotTable(:,'cell'))==x).*double(table2array(struct_quant_results_spotTable(:,'image'))==y)==1);
                        spotQuantTable_nas_nuc_cyt_type_in_cell_X_in_image_Y = @(x,y,ai,bi,ci) find(double(table2array(struct_quant_results_spotTable(:,'cell'))==x).*...
                            double(table2array(struct_quant_results_spotTable(:,'image'))==y).*...
                            double(table2array(struct_quant_results_spotTable(:,'nascent_flag'))==ai).*...
                            double(table2array(struct_quant_results_spotTable(:,'nucRNA'))==bi).*...
                            double(table2array(struct_quant_results_spotTable(:,'cytoRNA'))==ci)==1);
                        check_spotQuantTable_if_spot_nas_nuc_cyt_type = @(W,x,y,z) find(double(table2array(struct_quant_results_spotTable(W,'nascent_flag'))==x).*...
                            double(table2array(struct_quant_results_spotTable(W,'nucRNA'))==y).*...
                            double(table2array(struct_quant_results_spotTable(W,'cytoRNA'))==z)==1);
                        check_spotQuantTable_meetThreshold = @(W,Thr) find(double(table2array(struct_quant_results_spotTable(W,'dropout_thresh'))>Thr)==1);
                        get_spotQuantTable_meetThreshold = @(W,Thr) W(double(table2array(struct_quant_results_spotTable(W,'dropout_thresh'))>Thr)==1);
                        getSpotQuantTableMeetingThresholds = @(imageList,cellList,ThrList,ai,bi,ci) arrayfun(@(x) get_spotQuantTable_meetThreshold(spotQuantTable_nas_nuc_cyt_type_in_cell_X_in_image_Y(cellList(x),imageList(x),ai,bi,ci),ThrList(x))',1:length(ThrList),'Un',0);
                        numSpotQuantTable_at_Thr_in_cell_image = @(imageList,cellList,ThrList) CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(ci) CATnWrapper(arrayfun(@(bi) arrayfun(@(ai) length(check_spotQuantTable_meetThreshold(check_spotQuantTable_if_spot_nas_nuc_cyt_type(spotQuantTable_in_cell_X_in_image_Y(cellList(x),imageList(x)),ai,bi,ci),ThrList(x))),0:1)',0:1,'Un',0),2),0:1,'Un',0),3),1:length(imageList),'Un',0),4);
                        numSpotQuantTable_max_at_cell_image = @(imageList,cellList) CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(ci) CATnWrapper(arrayfun(@(bi) arrayfun(@(ai) length(check_spotQuantTable_if_spot_nas_nuc_cyt_type(spotQuantTable_in_cell_X_in_image_Y(cellList(x),imageList(x)),ai,bi,ci)),0:1),0:1,'Un',0),1),0:1,'Un',0),3),1:length(imageList),'Un',0),4);
                        softwarePooledResults_quant.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).numSpotQuantTable_at_Thr_in_cell_image = numSpotQuantTable_at_Thr_in_cell_image;
                        softwarePooledResults_quant.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).numSpotQuantTable_max_at_cell_image = numSpotQuantTable_max_at_cell_image;
                        softwarePooledResults_quant.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).getSpotQuantTableMeetingThresholds = getSpotQuantTableMeetingThresholds;
                    end
                    if (part==2)
                        temp_ImageCellStats = softwareResults.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).ind_results_part_cells_function(results_cellStats);
                        existing_cellStats = find([temp_ImageCellStats.CellStatsExist]==1);
                        struct_ImageCellStats = vertcat(temp_ImageCellStats(existing_cellStats).ImageCellStats);
                        if (~isempty(struct_ImageCellStats))
                            image_number_cellStatTable = cell2mat(arrayfun(@(x) x*ones(1,size(temp_ImageCellStats(x).ImageCellStats,1)),existing_cellStats,'Un',0))';
                            struct_ImageCellStats(:,'image') = array2table(image_number_cellStatTable);
                        end
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats = struct_ImageCellStats;
                        existing_cellStats_IntCorrected = find([temp_ImageCellStats.CellStatsExist_IntCorrected]==1);
                        struct_ImageCellStats_IntCorrected = vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).ImageCellStats_IntCorrected);
                        if (~isempty(struct_ImageCellStats_IntCorrected))
                            image_number_cellStatTable_IntCorrected = cell2mat(arrayfun(@(x) x*ones(1,size(temp_ImageCellStats(x).ImageCellStats_IntCorrected,1)),existing_cellStats_IntCorrected,'Un',0))';
                            struct_ImageCellStats_IntCorrected(:,'image') = array2table(image_number_cellStatTable_IntCorrected);
                        end
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats_IntCorrected = struct_ImageCellStats_IntCorrected;
                        % if (length(existing_cellStats)>1)
                        % struct_ImageCamBkgFullStats = struct2table(vertcat(temp_ImageCellStats(existing_cellStats).image_cambkg_Full));
                        % else
                        % struct_ImageCamBkgFullStats = struct2table(vertcat(temp_ImageCellStats(existing_cellStats).image_cambkg_Full),'AsArray',1);    
                        % end
                        % if (~isempty(struct_ImageCamBkgFullStats))
                        %     image_number_CamBkgFullStats = cell2mat(arrayfun(@(x) x,existing_cellStats,'Un',0))';
                        %     struct_ImageCamBkgFullStats(:,'image') = array2table(image_number_CamBkgFullStats);
                        % end
                        % softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCamBkgFullStats = struct_ImageCamBkgFullStats;
                        % if (length(existing_cellStats)>1)
                        % struct_ImageCamBkgTrimStats = struct2table(vertcat(temp_ImageCellStats(existing_cellStats).image_cambkg_Trim));
                        % else
                        % struct_ImageCamBkgTrimStats = struct2table(vertcat(temp_ImageCellStats(existing_cellStats).image_cambkg_Trim),'AsArray',1);
                        % end
                        % if (~isempty(struct_ImageCamBkgTrimStats))
                        %     image_number_CamBkgTrimStats = cell2mat(arrayfun(@(x) x,existing_cellStats,'Un',0))';
                        %     struct_ImageCamBkgTrimStats(:,'image') = array2table(image_number_CamBkgTrimStats);
                        % end
                        % softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCamBkgTrimStats = struct_ImageCamBkgTrimStats;
                        % if (length(existing_cellStats)>1)
                        % struct_ImageCellBkgFullStats = struct2table(vertcat(temp_ImageCellStats(existing_cellStats).image_cellbkg_Full));
                        % else
                        % struct_ImageCellBkgFullStats = struct2table(vertcat(temp_ImageCellStats(existing_cellStats).image_cellbkg_Full),'AsArray',1);    
                        % end
                        % if (~isempty(struct_ImageCellBkgFullStats))
                        %     image_number_CellBkgFullStats = cell2mat(arrayfun(@(x) x,existing_cellStats,'Un',0))';
                        %     struct_ImageCellBkgFullStats(:,'image') = array2table(image_number_CellBkgFullStats);
                        % end
                        % softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellBkgFullStats = struct_ImageCellBkgFullStats;
                        % if (length(existing_cellStats)>1)
                        % struct_ImageCellBkgTrimStats = struct2table(vertcat(temp_ImageCellStats(existing_cellStats).image_cellbkg_Trim));
                        % else
                        % struct_ImageCellBkgTrimStats = struct2table(vertcat(temp_ImageCellStats(existing_cellStats).image_cellbkg_Trim),'AsArray',1);    
                        % end
                        % if (~isempty(struct_ImageCellBkgTrimStats))
                        %     image_number_CellBkgTrimStats = cell2mat(arrayfun(@(x) x,existing_cellStats,'Un',0))';
                        %     struct_ImageCellBkgTrimStats(:,'image') = array2table(image_number_CellBkgTrimStats);
                        % end
                        % softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellBkgTrimStats = struct_ImageCellBkgTrimStats;
                        % if (length(existing_cellStats_IntCorrected)>1)
                        % struct_ImageCamBkgFullStats_IntCorrected = struct2table(vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).image_cambkg_Full_IntCorrected));
                        % else
                        % struct_ImageCamBkgFullStats_IntCorrected = struct2table(vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).image_cambkg_Full_IntCorrected),'AsArray',1);    
                        % end
                        % if (~isempty(struct_ImageCamBkgFullStats_IntCorrected))
                        %     image_number_CamBkgFullStats_IntCorrected = cell2mat(arrayfun(@(x) x,existing_cellStats_IntCorrected,'Un',0))';
                        %     struct_ImageCamBkgFullStats_IntCorrected(:,'image') = array2table(image_number_CamBkgFullStats_IntCorrected);
                        % end
                        % softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCamBkgFullStats_IntCorrected = struct_ImageCamBkgFullStats_IntCorrected;
                        % if (length(existing_cellStats_IntCorrected)>1)
                        % struct_ImageCamBkgTrimStats_IntCorrected = struct2table(vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).image_cambkg_Trim_IntCorrected));
                        % else
                        % struct_ImageCamBkgTrimStats_IntCorrected = struct2table(vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).image_cambkg_Trim_IntCorrected),'AsArray',1);
                        % end
                        % if (~isempty(struct_ImageCamBkgTrimStats_IntCorrected))
                        %     image_number_CamBkgTrimStats_IntCorrected = cell2mat(arrayfun(@(x) x,existing_cellStats_IntCorrected,'Un',0))';
                        %     struct_ImageCamBkgTrimStats_IntCorrected(:,'image') = array2table(image_number_CamBkgTrimStats_IntCorrected);
                        % end
                        % softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCamBkgTrimStats_IntCorrected = struct_ImageCamBkgTrimStats_IntCorrected;
                        % if (length(existing_cellStats_IntCorrected)>1)
                        % struct_ImageCellBkgFullStats_IntCorrected = struct2table(vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).image_cellbkg_Full_IntCorrected));
                        % else
                        % struct_ImageCellBkgFullStats_IntCorrected = struct2table(vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).image_cellbkg_Full_IntCorrected),'AsArray',1);    
                        % end
                        % if (~isempty(struct_ImageCellBkgFullStats_IntCorrected))
                        %     image_number_CellBkgFullStats_IntCorrected = cell2mat(arrayfun(@(x) x,existing_cellStats_IntCorrected,'Un',0))';
                        %     struct_ImageCellBkgFullStats_IntCorrected(:,'image') = array2table(image_number_CellBkgFullStats_IntCorrected);
                        % end
                        % softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellBkgFullStats_IntCorrected = struct_ImageCellBkgFullStats_IntCorrected;
                        % if (length(existing_cellStats_IntCorrected)>1)
                        % struct_ImageCellBkgTrimStats_IntCorrected = struct2table(vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).image_cellbkg_Trim_IntCorrected));
                        % else
                        % struct_ImageCellBkgTrimStats_IntCorrected = struct2table(vertcat(temp_ImageCellStats(existing_cellStats_IntCorrected).image_cellbkg_Trim_IntCorrected),'AsArray',1);    
                        % end
                        % if (~isempty(struct_ImageCellBkgTrimStats_IntCorrected))
                        %     image_number_CellBkgTrimStats_IntCorrected = cell2mat(arrayfun(@(x) x,existing_cellStats_IntCorrected,'Un',0))';
                        %     struct_ImageCellBkgTrimStats_IntCorrected(:,'image') = array2table(image_number_CellBkgTrimStats_IntCorrected);
                        % end
                        %softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellBkgTrimStats_IntCorrected = struct_ImageCellBkgTrimStats_IntCorrected;
                        temp_ImageSpotStats = softwareResults.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).ind_results_part_spots_function(results_spotStats);
                        existing_spotStats = find([temp_ImageSpotStats.SpotStatsExist]==1);
                        struct_ImageSpotStats = vertcat(temp_ImageSpotStats(existing_spotStats).ImageSpotStats);
                        if (~isempty(struct_ImageSpotStats))
                            image_number_spotStatTable = cell2mat(arrayfun(@(x) x*ones(1,size(temp_ImageSpotStats(x).ImageSpotStats,1)),existing_spotStats,'Un',0))';
                            struct_ImageSpotStats(:,'image') = array2table(image_number_spotStatTable);
                        end
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats = struct_ImageSpotStats;
                        existing_spotStats_IntCorrected = find([temp_ImageSpotStats.SpotStatsExist_IntCorrected]==1);
                        struct_ImageSpotStats_IntCorrected = vertcat(temp_ImageSpotStats(existing_spotStats_IntCorrected).ImageSpotStats_IntCorrected);
                        if (~isempty(struct_ImageSpotStats_IntCorrected))
                            image_number_spotStatTable_IntCorrected = cell2mat(arrayfun(@(x) x*ones(1,size(temp_ImageSpotStats(x).ImageSpotStats_IntCorrected,1)),existing_spotStats_IntCorrected,'Un',0))';
                            struct_ImageSpotStats_IntCorrected(:,'image') = array2table(image_number_spotStatTable_IntCorrected);
                        end
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats_IntCorrected = struct_ImageSpotStats_IntCorrected;
                        spotStatTable_in_cell_X_in_image_Y = @(x,y) find(double(table2array(struct_ImageSpotStats(:,'cell'))==x).*double(table2array(struct_ImageSpotStats(:,'image'))==y)==1);
                        spotStatTable_nas_nuc_cyt_type_in_cell_X_in_image_Y = @(x,y,ai,bi,ci) find(double(table2array(struct_ImageSpotStats(:,'cell'))==x).*...
                            double(table2array(struct_ImageSpotStats(:,'image'))==y).*...
                            double(table2array(struct_ImageSpotStats(:,'nascent_flag'))==ai).*...
                            double(table2array(struct_ImageSpotStats(:,'nuc_flag'))==bi).*...
                            double(table2array(struct_ImageSpotStats(:,'nuc_flag'))==1-ci)==1);
                        check_spotStatTable_if_spot_nas_nuc_cyt_type = @(W,x,y,z) find(double(table2array(struct_ImageSpotStats(W,'nascent_flag'))==x).*...
                            double(table2array(struct_ImageSpotStats(W,'nuc_flag'))==y).*...
                            double(table2array(struct_ImageSpotStats(W,'nuc_flag'))==1-z)==1);
                        check_spotStatTable_meetThreshold = @(W,Thr) find(double(table2array(struct_ImageSpotStats(W,'dropout_thresh'))>Thr)==1);
                        get_spotStatTable_meetThreshold = @(W,Thr) W(double(table2array(struct_ImageSpotStats(W,'dropout_thresh'))>Thr)==1);
                        getSpotStatTableMeetingThresholds = @(imageList,cellList,ThrList,ai,bi,ci) arrayfun(@(x) get_spotStatTable_meetThreshold(spotStatTable_nas_nuc_cyt_type_in_cell_X_in_image_Y(cellList(x),imageList(x),ai,bi,ci),ThrList(x))',1:length(ThrList),'Un',0);
                        numSpotStatTable_at_Thr_in_cell_image = @(imageList,cellList,ThrList) CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(ci) CATnWrapper(arrayfun(@(bi) arrayfun(@(ai) length(check_spotStatTable_meetThreshold(check_spotStatTable_if_spot_nas_nuc_cyt_type(spotStatTable_in_cell_X_in_image_Y(cellList(x),imageList(x)),ai,bi,ci),ThrList(x))),0:1)',0:1,'Un',0),2),0:1,'Un',0),3),1:length(imageList),'Un',0),4);
                        numSpotStatTable_max_at_cell_image = @(imageList,cellList) CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(ci) CATnWrapper(arrayfun(@(bi) arrayfun(@(ai) length(check_spotStatTable_if_spot_nas_nuc_cyt_type(spotStatTable_in_cell_X_in_image_Y(cellList(x),imageList(x)),ai,bi,ci)),0:1),0:1,'Un',0),1),0:1,'Un',0),3),1:length(imageList),'Un',0),4);
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).numSpotStatTable_at_Thr_in_cell_image = numSpotStatTable_at_Thr_in_cell_image;
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).numSpotStatTable_max_at_cell_image = numSpotStatTable_max_at_cell_image;
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).getSpotStatTableMeetingThresholds = getSpotStatTableMeetingThresholds;
                        if (~isempty(softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats))
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats(:,{'max_spot_int'}) = ...
                            array2table(cellfun(@(x) double(x.max),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats.noBkgSubStats));
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats(:,{'mean_bkg_int'}) = ...
                            array2table(cellfun(@(x) double(x.mean),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats.localBkgStats));
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats(:,{'Sig_minus_bkg'}) = ...
                            array2table(cellfun(@(x) double(x.max),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats.noBkgSubStats)) - array2table(cellfun(@(x) double(x.mean),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats.localBkgStats));
                        end
                        if (~isempty(softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats))
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats(:,{'mean_cell_bkg'}) = ...
                            array2table(cellfun(@(x) double(x.mean),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats.cell_bkg));
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats(:,{'mean_cam_bkg'}) = ...
                            array2table(cellfun(@(x) double(x.mean),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats.local_img_bkg));
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats(:,{'mean_cell_minus_cam_bkg'}) = ...
                            array2table(cellfun(@(x) double(x.mean),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats.cell_bkg)) - array2table(cellfun(@(x) double(x.mean),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats.local_img_bkg));
                         end
                        if (~isempty(softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats_IntCorrected))
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats_IntCorrected(:,{'max_spot_int'}) = ...
                            array2table(cellfun(@(x) double(x.max),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats_IntCorrected.noBkgSubStats));
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats_IntCorrected(:,{'mean_bkg_int'}) = ...
                            array2table(cellfun(@(x) double(x.mean),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats_IntCorrected.localBkgStats));
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats_IntCorrected(:,{'Sig_minus_bkg'}) = ...
                            array2table(cellfun(@(x) double(x.max),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats_IntCorrected.noBkgSubStats)) - array2table(cellfun(@(x) double(x.mean),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageSpotStats_IntCorrected.localBkgStats));
                        end
                        if (~isempty(softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats_IntCorrected))
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats_IntCorrected(:,{'mean_cell_bkg'}) = ...
                            array2table(cellfun(@(x) double(x.mean),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats_IntCorrected.cell_bkg));
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats_IntCorrected(:,{'mean_cam_bkg'}) = ...
                            array2table(cellfun(@(x) double(x.mean),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats_IntCorrected.local_img_bkg));
                        softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats_IntCorrected(:,{'mean_cell_minus_cam_bkg'}) = ...
                            array2table(cellfun(@(x) double(x.mean),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats_IntCorrected.cell_bkg)) - array2table(cellfun(@(x) double(x.mean),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImageCellStats_IntCorrected.local_img_bkg));
                        end
                    end
                    if (part==3)
                        n_reps = size(softwareResults.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).ind_results_part_cells_function(results_ImageThresholdsStats),2);
                        %% threshold candidate info for each image
                        temp_ImageThresholdsStats = softwareResults.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).ind_results_part_cells_function(results_ImageThresholdsStats);
                        image_sugg_m_table = cell2mat(arrayfun(@(x) temp_ImageThresholdsStats(x).sugg_m',1:n_reps,'Un',0))';
                        image_threshold_mean = [temp_ImageThresholdsStats.mean_overall];
                        image_threshold_med = [temp_ImageThresholdsStats.med_overall];
                        image_threshold_std = [temp_ImageThresholdsStats.std_overall];
                        image_threshold_min = [temp_ImageThresholdsStats.min_overall];
                        image_threshold_max = [temp_ImageThresholdsStats.max_overall];
                        image_threshold_opt =  [temp_ImageThresholdsStats.threshold];
                        image_threshold_spot_count_curve_x = cell2mat(arrayfun(@(x) temp_ImageThresholdsStats(x).threshold_spot_count_curve.x',1:n_reps,'Un',0))';
                        image_threshold_spot_count_curve_y = cell2mat(arrayfun(@(x) temp_ImageThresholdsStats(x).threshold_spot_count_curve.y',1:n_reps,'Un',0))';
                        image_threshold_score_curve_x = cell2mat(arrayfun(@(x) temp_ImageThresholdsStats(x).threshold_score_curve.x',1:n_reps,'Un',0))';
                        image_threshold_score_curve_score = cell2mat(arrayfun(@(x) temp_ImageThresholdsStats(x).threshold_score_curve.score',1:n_reps,'Un',0))';
                        image_threshold_score_curve_score_m = cell2mat(arrayfun(@(x) temp_ImageThresholdsStats(x).threshold_score_curve.score_m',1:n_reps,'Un',0))';
                        image_threshold_score_curve_score_f = cell2mat(arrayfun(@(x) temp_ImageThresholdsStats(x).threshold_score_curve.score_f',1:n_reps,'Un',0))';
                        image_threshold_score_curve_score_fri = cell2mat(arrayfun(@(x) temp_ImageThresholdsStats(x).threshold_score_curve.score_fri',1:n_reps,'Un',0))';
                        %% storing image threshold information
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_sugg_m_table = image_sugg_m_table;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_mean = image_threshold_mean;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_med = image_threshold_med;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_std = image_threshold_std;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_min = image_threshold_min;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_max = image_threshold_max;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_opt = image_threshold_opt;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_spot_count_curve_x = image_threshold_spot_count_curve_x;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_spot_count_curve_y = image_threshold_spot_count_curve_y;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_score_curve_x = image_threshold_score_curve_x;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_score_curve_score = image_threshold_score_curve_score;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_score_curve_score_m = image_threshold_score_curve_score_m;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_score_curve_score_f = image_threshold_score_curve_score_f;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_score_curve_score_fri = image_threshold_score_curve_score_fri;
                        %% threshold info for single-cells
                        temp_CellThresholdsStats = softwareResults.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).ind_results_part_cells_function(results_CellThresholdsStats);
                        existing_cellThresholdsStats = find([temp_CellThresholdsStats.CellThresholdsExist]==1);
                        cell_sugg_m_table = cell2mat(arrayfun(@(x) temp_CellThresholdsStats(x).scThresholdSuggestions.sugg_m',existing_cellThresholdsStats,'Un',0))';
                        cell_threshold_mean = cell2mat(arrayfun(@(x) temp_CellThresholdsStats(x).scThresholdSuggestions.mean_overall',existing_cellThresholdsStats,'Un',0))';
                        cell_threshold_med = cell2mat(arrayfun(@(x) temp_CellThresholdsStats(x).scThresholdSuggestions.med_overall',existing_cellThresholdsStats,'Un',0))';
                        cell_threshold_std = cell2mat(arrayfun(@(x) temp_CellThresholdsStats(x).scThresholdSuggestions.std_overall',existing_cellThresholdsStats,'Un',0))';
                        cell_threshold_min = cell2mat(arrayfun(@(x) temp_CellThresholdsStats(x).scThresholdSuggestions.min_overall',existing_cellThresholdsStats,'Un',0))';
                        cell_threshold_max = cell2mat(arrayfun(@(x) temp_CellThresholdsStats(x).scThresholdSuggestions.max_overall',existing_cellThresholdsStats,'Un',0))';
                        cell_threshold_opt = cell2mat(arrayfun(@(x) temp_CellThresholdsStats(x).scThresholdSuggestions.threshold',existing_cellThresholdsStats,'Un',0))';
                        cell_threshold_spot_count_curve_x = cell2mat(arrayfun(@(v) [temp_CellThresholdsStats(v).scThresholdSuggestions.x{:}],existing_cellThresholdsStats,'Un',0))';
                        cell_threshold_spot_count_curve_y = cell2mat(arrayfun(@(x) [temp_CellThresholdsStats(x).scThresholdSuggestions.spot_counts{:}],existing_cellThresholdsStats,'Un',0))';
                        cell_threshold_score_curve_x = cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(w) temp_CellThresholdsStats(x).scThresholdSuggestions.value_table{w}.threshold_value, ...
                            1:size(temp_CellThresholdsStats(x).scThresholdSuggestions.value_table,1),'Un',0)),existing_cellThresholdsStats,'Un',0))';
                        cell_threshold_score_curve_score = cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(w) temp_CellThresholdsStats(x).scThresholdSuggestions.value_table{w}.score, ...
                            1:size(temp_CellThresholdsStats(x).scThresholdSuggestions.value_table,1),'Un',0)),existing_cellThresholdsStats,'Un',0))';
                        cell_threshold_score_curve_score_m = cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(w) temp_CellThresholdsStats(x).scThresholdSuggestions.value_table{w}.score_m, ...
                            1:size(temp_CellThresholdsStats(x).scThresholdSuggestions.value_table,1),'Un',0)),existing_cellThresholdsStats,'Un',0))';
                        cell_threshold_score_curve_score_f = cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(w) temp_CellThresholdsStats(x).scThresholdSuggestions.value_table{w}.score_f, ...
                            1:size(temp_CellThresholdsStats(x).scThresholdSuggestions.value_table,1),'Un',0)),existing_cellThresholdsStats,'Un',0))';
                        cell_threshold_score_curve_score_fri = cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(w) temp_CellThresholdsStats(x).scThresholdSuggestions.value_table{w}.score_fri, ...
                            1:size(temp_CellThresholdsStats(x).scThresholdSuggestions.value_table,1),'Un',0)),existing_cellThresholdsStats,'Un',0))';
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_sugg_m_table = cell_sugg_m_table;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_mean = cell_threshold_mean;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_med = cell_threshold_med;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_std = cell_threshold_std;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_min = cell_threshold_min;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_max = cell_threshold_max;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_opt = cell_threshold_opt;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_spot_count_curve_x = cell_threshold_spot_count_curve_x;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_spot_count_curve_y = cell_threshold_spot_count_curve_y;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_score_curve_x = cell_threshold_score_curve_x;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_score_curve_score = cell_threshold_score_curve_score;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_score_curve_score_m = cell_threshold_score_curve_score_m;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_score_curve_score_f = cell_threshold_score_curve_score_f;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_score_curve_score_fri = cell_threshold_score_curve_score_fri;
                        temp_callTable = softwareResults.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).ind_results_part_cells_function(results_callTable);
                        struct_callTable = vertcat(temp_callTable(:).call_table);
                        if (~isempty(struct_callTable))
                            image_number_callTable = cell2mat(arrayfun(@(x) x*ones(1,size(temp_callTable(x).call_table,1)),1:n_reps,'Un',0))';
                            struct_callTable(:,'image') = array2table(image_number_callTable);
                        end
                        in_image_Y = @(y) find(double(table2array(struct_callTable(:,'image'))==y)==1);
                        no_cell_in_image_Y = @(y) find(double(table2array(struct_callTable(:,'cell'))==0).*...
                            double(table2array(struct_callTable(:,'image'))==y)==1);
                        in_cells_in_image_Y = @(y) find(double(table2array(struct_callTable(:,'cell'))>0).*...
                            double(table2array(struct_callTable(:,'image'))==y)==1);
                        in_cell_X_in_image_Y = @(x,y) find(double(table2array(struct_callTable(:,'cell'))==x).*...
                            double(table2array(struct_callTable(:,'image'))==y)==1);
                        unique_image_and_cell_pair = unique(double(table2array(struct_callTable(find(double(table2array(struct_callTable(:,'cell'))>0)),{'image','cell'}))),'rows');
                        extended_callTable_spot_count_curve_in_image_Y_y = cell2mat(arrayfun(@(Y) ...
                            numel(in_image_Y(Y))-histcounts(table2array(struct_callTable(in_image_Y(Y),'dropout_thresh')),'Normalization','cumcount','BinEdges',[-Inf 1:size(image_threshold_spot_count_curve_x,2)])',1:n_reps,'Un',0))';
                        extended_callTable_spot_count_curve_no_cell_in_image_Y_y = cell2mat(arrayfun(@(Y) ...
                            numel(no_cell_in_image_Y(Y))-histcounts(table2array(struct_callTable(no_cell_in_image_Y(Y),'dropout_thresh')),'Normalization','cumcount','BinEdges',[-Inf 1:size(image_threshold_spot_count_curve_x,2)])',1:n_reps,'Un',0))';
                        extended_callTable_spot_count_curve_in_cells_in_image_Y_y = cell2mat(arrayfun(@(Y) ...
                            numel(in_cells_in_image_Y(Y))-histcounts(table2array(struct_callTable(in_cells_in_image_Y(Y),'dropout_thresh')),'Normalization','cumcount','BinEdges',[-Inf 1:size(image_threshold_spot_count_curve_x,2)])',1:n_reps,'Un',0))';
                        extended_callTable_spot_count_curve_in_cell_X_in_image_Y_y = cell2mat(arrayfun(@(Y) ...
                            numel(in_cell_X_in_image_Y(unique_image_and_cell_pair(Y,2),unique_image_and_cell_pair(Y,1)))-histcounts(table2array(struct_callTable(in_cell_X_in_image_Y(unique_image_and_cell_pair(Y,2),unique_image_and_cell_pair(Y,1)),'dropout_thresh')),'Normalization','cumcount','BinEdges',[-Inf 1:size(image_threshold_spot_count_curve_x,2)])',1:size(unique_image_and_cell_pair,1),'Un',0))';
                        cell_spot_count_at_thrList_for_unique_image_and_cell_pair = @(cellList,ThrList) extended_callTable_spot_count_curve_in_cell_X_in_image_Y_y(sub2ind(...
                            size(extended_callTable_spot_count_curve_in_cell_X_in_image_Y_y),cellList,ceil(ThrList)));
                        callTable_in_cell_X_in_image_Y = @(x,y) find(double(table2array(struct_callTable(:,'cell'))==x).*double(table2array(struct_callTable(:,'image'))==y)==1);
                        check_callTable_meetThreshold = @(W,Thr) find(double(table2array(struct_callTable(W,'dropout_thresh'))>=Thr)==1);
                        get_callTable_meetThreshold = @(W,Thr) W(double(table2array(struct_callTable(W,'dropout_thresh'))>=Thr)==1);
                        getCallTableMeetingThresholds = @(imageList,cellList,ThrList) arrayfun(@(x) get_callTable_meetThreshold(callTable_in_cell_X_in_image_Y(cellList(x),imageList(x)),ThrList(x))',1:length(ThrList),'Un',0);
                        numCallTable_at_Thr_in_cell_image = @(imageList,cellList,ThrList) arrayfun(@(x) length(check_callTable_meetThreshold(callTable_in_cell_X_in_image_Y(cellList(x),imageList(x)),ThrList(x))),1:length(imageList));
                        is_main_spot = A_BHJH_ApplySpotMergingGeneric(struct_callTable,'isnap_x','isnap_y','isnap_z','intensity');
                        struct_callTable(:,'is_main_spot') = array2table(is_main_spot);
                        in_cell_X_in_image_Y_merged = @(x,y) find(double(table2array(struct_callTable(:,'cell'))==x).*...
                            double(table2array(struct_callTable(:,'image'))==y).*...
                            double(table2array(struct_callTable(:,'is_main_spot'))==1)==1);                        
                        callTable_merged_spot_count_curve_in_cell_X_in_image_Y_y = cell2mat(arrayfun(@(Y) ...
                            numel(in_cell_X_in_image_Y_merged(unique_image_and_cell_pair(Y,2),unique_image_and_cell_pair(Y,1)))-histcounts(table2array(struct_callTable(in_cell_X_in_image_Y_merged(unique_image_and_cell_pair(Y,2),unique_image_and_cell_pair(Y,1)),'dropout_thresh')),'Normalization','cumcount','BinEdges',[-Inf 1:size(image_threshold_spot_count_curve_x,2)])',1:size(unique_image_and_cell_pair,1),'Un',0))';
                        spotmerged_callTable_in_cell_X_in_image_Y= @(x,y) find(double(table2array(struct_callTable(:,'cell'))==x).*double(table2array(struct_callTable(:,'image'))==y).*double(table2array(struct_callTable(:,'is_main_spot'))==1)==1);
                        spotmerged_check_callTable_meetThreshold = @(W,Thr) find(double(table2array(struct_callTable(W,'dropout_thresh'))>=Thr).*double(table2array(struct_callTable(W,'is_main_spot'))==1)==1);
                        spotmerged_get_callTable_meetThresholds = @(W,Thr) W(double(table2array(struct_callTable(W,'dropout_thresh'))>=Thr).*double(table2array(struct_callTable(W,'is_main_spot'))==1)==1);
                        spotmerged_numCallTable_at_Thr_in_cell_image = @(imageList,cellList,ThrList) arrayfun(@(x) length(spotmerged_check_callTable_meetThreshold(spotmerged_callTable_in_cell_X_in_image_Y(cellList(x),imageList(x)),ThrList(x))),1:length(imageList));
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).unique_image_and_cell_pair = unique_image_and_cell_pair;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).callTable_spot_count_curve_in_cell_X_in_image_Y_y = extended_callTable_spot_count_curve_in_cell_X_in_image_Y_y;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).numCallTable_at_Thr_in_cell_image = numCallTable_at_Thr_in_cell_image;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).callTable_merged_spot_count_curve_in_cell_X_in_image_Y_y = callTable_merged_spot_count_curve_in_cell_X_in_image_Y_y;
                        softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).spotmerged_numCallTable_at_Thr_in_cell_image = spotmerged_numCallTable_at_Thr_in_cell_image;
                        softwarePooledResults_callTable.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_callTable = struct_callTable;                        
                        softwarePooledResults_callTable.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).unique_image_and_cell_pair = unique_image_and_cell_pair;
                        softwarePooledResults_callTable.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).extended_callTable_spot_count_curve_in_image_Y_y = extended_callTable_spot_count_curve_in_image_Y_y;
                        softwarePooledResults_callTable.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).extended_callTable_spot_count_curve_no_cell_in_image_Y_y = extended_callTable_spot_count_curve_no_cell_in_image_Y_y;
                        softwarePooledResults_callTable.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).extended_callTable_spot_count_curve_in_cells_in_image_Y_y = extended_callTable_spot_count_curve_in_cells_in_image_Y_y;
                        softwarePooledResults_callTable.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).extended_callTable_spot_count_curve_in_cell_X_in_image_Y_y = extended_callTable_spot_count_curve_in_cell_X_in_image_Y_y;
                        softwarePooledResults_callTable.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_spot_count_at_thrList_for_unique_image_and_cell_pair = cell_spot_count_at_thrList_for_unique_image_and_cell_pair;
                        softwarePooledResults_callTable.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).numCallTable_at_Thr_in_cell_image = numCallTable_at_Thr_in_cell_image;
                        softwarePooledResults_callTable.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).getCallTableMeetingThresholds = getCallTableMeetingThresholds;
                        softwarePooledResults_callTable.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).callTable_merged_spot_count_curve_in_cell_X_in_image_Y_y = callTable_merged_spot_count_curve_in_cell_X_in_image_Y_y;
                        softwarePooledResults_callTable.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).spotmerged_numCallTable_at_Thr_in_cell_image = spotmerged_numCallTable_at_Thr_in_cell_image;
                        softwarePooledResults_callTable.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).spotmerged_getCallTableMeetingThresholds = spotmerged_get_callTable_meetThresholds;
                    end
                    if (part==4)
                        subfields = {'cell_mask',...
                            'nucmask_max_z','nucmask_opt_z','nucmask_mid_z',...
                            'TRANS_max_z','TRANS_opt_z','TRANS_mid_z',...
                            'DAPI_max_z','DAPI_opt_z','DAPI_mid_z',...
                            'PROBE_max_z','PROBE_opt_z','PROBE_mid_z'};
                        temp_ImageCellStats = softwareResults.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).ind_results_part_cells_function(results_cellStats);
                        existing_cellStats = find([temp_ImageCellStats.CellStatsExist]==1);
                        struct_ImageCellStats = vertcat(temp_ImageCellStats(existing_cellStats).ImageCellStats);
                        try
                            temp_cellseg = softwareResults.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).ind_results_part_cells_function(results_cellseg);
                            temp_cellseg = vertcat(temp_cellseg(:).cellSeg);
                            if (size(temp_cellseg,1)>1)
                                struct_cellSeg = struct2table(temp_cellseg);
                            else
                                struct_cellSeg = struct2table(temp_cellseg,'AsArray',1);
                            end
                            if (~isempty(struct_cellSeg))
                                struct_cellSeg(:,'image') = array2table(existing_cellStats');
                            end
                            softwarePooledResults_cellseg.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_cellSeg = struct_cellSeg;
                        catch e
                            disp(e)
                            for v = 1:length(e.stack)
                                e.stack(v)
                            end
                        end
                        temp_ImgData = softwareResults.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).ind_results_part_cells_function(results_imgData);
                        if (size(temp_ImgData,2)>1)
                            struct_ImgData = struct2table(temp_ImgData);
                        else
                            struct_ImgData = struct2table(temp_ImgData,'AsArray',1);
                        end
                        if (~isempty(struct_ImgData))
                            image_number_imgData = 1:size(struct_ImgData,1);
                            struct_ImgData(:,'image') = array2table(image_number_imgData');
                        end
                        softwarePooledResults_imgData.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_ImgData = struct_ImgData;
                        if (~isempty(struct_ImageCellStats))
                            image_number_cellStatTable = cell2mat(arrayfun(@(x) x*ones(1,size(temp_ImageCellStats(x).ImageCellStats,1)),existing_cellStats,'Un',0))';
                            struct_ImageCellStats(:,'image') = array2table(image_number_cellStatTable);
                            length_Y = arrayfun(@(w) length([struct_ImageCellStats.box{w}.top:struct_ImageCellStats.box{w}.bottom]),1:size(struct_ImageCellStats,1));
                            length_X = arrayfun(@(w) length([struct_ImageCellStats.box{w}.left:struct_ImageCellStats.box{w}.right]),1:size(struct_ImageCellStats,1));
                            for f = 1:length(subfields)
                                try
                                    struct_ImageCellStats(:,subfields{f}) = cell2table(arrayfun(@(w) struct_ImgData.(subfields{f}){...
                                        struct_ImageCellStats.image(w)}(struct_ImageCellStats.box{w}.top:struct_ImageCellStats.box{w}.bottom,...
                                        struct_ImageCellStats.box{w}.left:struct_ImageCellStats.box{w}.right),1:size(struct_ImageCellStats,1),'Un',0)');
                                catch
                                end
                            end
                        end
                        softwarePooledResults_imgData.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).struct_SingleCellImgData = struct_ImageCellStats;
                    end
                    %% function gives threshold curves for subsets
                    %% function gives list cell counts at given threshold for each cell
                    % %each cell (list number of nuclear,cytoplasmic, and nascent spots)
                    %gets position of spots meeting thresholds [has to be altered]
                    %number_at_thr(#nascent,#nuclear,$cytoplasmic)
                    %number_max(#nascent,#nuclear,$cytoplasmic,cell)
                    %% storing main info
                    softwarePooledResults_quant.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).n_reps = n_reps;
                    softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).n_reps = n_reps;
                    softwarePooledResults_callTable.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).n_reps = n_reps;
                    softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).n_reps = n_reps;
                    softwarePooledResults_cellseg.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).n_reps = n_reps;
                    softwarePooledResults_imgData.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).n_reps = n_reps;
                else
                    softwarePooledResults_quant.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).n_reps = 0;
                    softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).n_reps = 0;
                    softwarePooledResults_callTable.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).n_reps = 0;
                    softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).n_reps = 0;
                    softwarePooledResults_cellseg.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).n_reps = 0;
                    softwarePooledResults_imgData.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).n_reps = 0;
                end
            end
        end
    end
    if (part==1)
        save('pooled_hARF4_RNAFISH_AnalysisQuant.mat','noProbePooledResults_quant','softwarePooledResults_quant','-v7.3');
    end
    if (part==2)
        save('pooled_hARF4_RNAFISH_AnalysisStats.mat','noProbePooledResults_stats','softwarePooledResults_stats','-v7.3');
    end
    if (part==3)
        save('pooled_hARF4_RNAFISH_AnalysisThresholdInfo.mat','noProbePooledResults_thresholdInfo','softwarePooledResults_thresholdInfo','-v7.3');
        save('pooled_hARF4_RNAFISH_AnalysisCallTable.mat','noProbePooledResults_callTable','softwarePooledResults_callTable','-v7.3');
    end
    if (part==4)
        save('pooled_hARF4_RNAFISH_AnalysisImgData.mat','noProbePooledResults_imgData','softwarePooledResults_imgData','-v7.3');
        save('pooled_hARF4_RNAFISH_AnalysisCellseg.mat','noProbePooledResults_cellseg','softwarePooledResults_cellseg','-v7.3');
    end
end
clear softwareResults noProbeResults
%% handle functions for getting software image thresholds using min, max, mean,mean+std,med, optimal
if (part==3)
    try
        load('pooled_hARF4_RNAFISH_AnalysisThresholdInfo_ThresholdTypes.mat','noProbePooledResults_thresholdInfo','softwarePooledResults_thresholdInfo');
    catch
        quartile_function = @(A,Q) prctile(A,Q);
        min_imageThreshold = @(ii,ij,ik) cell2mat(arrayfun(@(k) min(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_sugg_m_table(k,:),[],'omitnan'),1:size(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_sugg_m_table,1),'Un',0));
        min_imageThreshold_image = @(ii,ij,ik,ki) cell2mat(arrayfun(@(k) min(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_sugg_m_table(k,:),[],'omitnan'), ki,'Un',0));
        mean_std_imageThreshold = @(ii,ij,ik,up_or_down) cell2mat(arrayfun(@(k) mean(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_sugg_m_table(k,:),'omitnan')+...
            up_or_down*std(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_sugg_m_table(k,:),'omitnan'),1:size(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_sugg_m_table,1),'Un',0));
        mean_std_imageThreshold_image = @(ii,ij,ik,up_or_down,ki) cell2mat(arrayfun(@(k) mean(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_sugg_m_table(k,:),'omitnan')+...
            up_or_down*std(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_sugg_m_table(k,:),'omitnan'),ki,'Un',0));
        quartile_imageThreshold = @(ii,ij,ik,Q) cell2mat(arrayfun(@(k) quartile_function(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_sugg_m_table(k,:),Q),1:size(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_sugg_m_table,1),'Un',0));
        quartile_imageThreshold_image = @(ii,ij,ik,Q,ki) cell2mat(arrayfun(@(k) quartile_function(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_sugg_m_table(k,:),Q),ki,'Un',0));
        med_imageThreshold = @(ii,ij,ik) softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_med;
        med_imageThreshold_image = @(ii,ij,ik,ki) softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_med(ki);
        opt_imageThreshold = @(ii,ij,ik) softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_opt;
        opt_imageThreshold_image = @(ii,ij,ik,ki) softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).image_threshold_opt(ki);
        %% handle functions for getting no-probe image thresholds using min, max, mean,mean+std,med, optimal
        min_imageThreshold_noprobe = @(ik) cell2mat(arrayfun(@(k) min(noProbePooledResults_thresholdInfo.(Replica{ik}).image_sugg_m_table(k,:),[],'omitnan'),1:size(noProbePooledResults_thresholdInfo.(Replica{ik}).image_sugg_m_table,1),'Un',0));
        min_imageThreshold_image_noprobe = @(ik,ki) cell2mat(arrayfun(@(k) min(noProbePooledResults_thresholdInfo.(Replica{ik}).image_sugg_m_table(k,:),[],'omitnan'),ki,'Un',0));
        mean_std_imageThreshold_noprobe = @(ik,up_or_down) cell2mat(arrayfun(@(k) mean(noProbePooledResults_thresholdInfo.(Replica{ik}).image_sugg_m_table(k,:),'omitnan')+...
            up_or_down*std(noProbePooledResults_thresholdInfo.(Replica{ik}).image_sugg_m_table(k,:),'omitnan'),1:size(noProbePooledResults_thresholdInfo.(Replica{ik}).image_sugg_m_table,1),'Un',0));
        mean_std_imageThreshold_image_noprobe = @(ik,up_or_down,ki) cell2mat(arrayfun(@(k) mean(noProbePooledResults_thresholdInfo.(Replica{ik}).image_sugg_m_table(k,:),'omitnan')+...
            up_or_down*std(noProbePooledResults_thresholdInfo.(Replica{ik}).image_sugg_m_table(k,:),'omitnan'),ki,'Un',0));
        quartile_imageThreshold_noprobe = @(ik,Q) cell2mat(arrayfun(@(k) quartile_function(noProbePooledResults_thresholdInfo.(Replica{ik}).image_sugg_m_table(k,:),Q),1:size(noProbePooledResults_thresholdInfo.(Replica{ik}).image_sugg_m_table,1),'Un',0));
        quartile_imageThreshold_image_noprobe = @(ik,Q,ki) cell2mat(arrayfun(@(k) quartile_function(noProbePooledResults_thresholdInfo.(Replica{ik}).image_sugg_m_table(k,:),Q),ki,'Un',0));
        med_imageThreshold_noprobe = @(ik) noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_med;
        med_imageThreshold_image_noprobe = @(ik,ki) noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_med(ki);
        opt_imageThreshold_noprobe = @(ik) noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_opt;
        opt_imageThreshold_image_noprobe = @(ik,ki) noProbePooledResults_thresholdInfo.(Replica{ik}).image_threshold_opt(ki);
        %% handle functions for getting software cell thresholds using min, max, mean,mean+std, med, optimal
        min_cellThreshold = @(ii,ij,ik) cell2mat(arrayfun(@(k) min(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_sugg_m_table(k,:),[],'omitnan'),1:size(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_sugg_m_table,1),'Un',0));
        min_cellThreshold_cell = @(ii,ij,ik,ki) cell2mat(arrayfun(@(k) min(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_sugg_m_table(k,:),[],'omitnan'), ki,'Un',0));
        mean_std_cellThreshold = @(ii,ij,ik,up_or_down) cell2mat(arrayfun(@(k) mean(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_sugg_m_table(k,:),'omitnan')+...
            up_or_down*std(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_sugg_m_table(k,:),'omitnan'),1:size(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_sugg_m_table,1),'Un',0));
        mean_std_cellThreshold_cell = @(ii,ij,ik,up_or_down,ki) cell2mat(arrayfun(@(k) mean(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_sugg_m_table(k,:),'omitnan')+...
            up_or_down*std(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_sugg_m_table(k,:),'omitnan'),ki,'Un',0));
        quartile_cellThreshold = @(ii,ij,ik,Q) cell2mat(arrayfun(@(k) quartile_function(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_sugg_m_table(k,:),Q),1:size(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_sugg_m_table,1),'Un',0));
        quartile_cellThreshold_cell = @(ii,ij,ik,Q,ki) cell2mat(arrayfun(@(k) quartile_function(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_sugg_m_table(k,:),Q),ki,'Un',0));
        med_cellThreshold = @(ii,ij,ik) softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_med;
        med_cellThreshold_cell = @(ii,ij,ik,ki) softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_med(ki);
        opt_cellThreshold = @(ii,ij,ik) softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_opt;
        opt_cellThreshold_cell = @(ii,ij,ik,ki) softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).cell_threshold_opt(ki);
        %% handle functions for getting no-probe cell thresholds using min, max, mean,mean+std, med, optimal
        min_cellThreshold_noprobe = @(ik) cell2mat(arrayfun(@(k) min(noProbePooledResults_thresholdInfo.(Replica{ik}).cell_sugg_m_table(k,:),[],'omitnan'),1:size(noProbePooledResults_thresholdInfo.(Replica{ik}).cell_sugg_m_table,1),'Un',0));
        min_cellThreshold_cell_noprobe = @(ik,ki) cell2mat(arrayfun(@(k) min(noProbePooledResults_thresholdInfo.(Replica{ik}).cell_sugg_m_table(k,:),[],'omitnan'),ki,'Un',0));
        mean_std_cellThreshold_noprobe = @(ik,up_or_down) cell2mat(arrayfun(@(k) mean(noProbePooledResults_thresholdInfo.(Replica{ik}).cell_sugg_m_table(k,:),'omitnan')+...
            up_or_down*std(noProbePooledResults_thresholdInfo.(Replica{ik}).cell_sugg_m_table(k,:),'omitnan'),1:size(noProbePooledResults_thresholdInfo.(Replica{ik}).cell_sugg_m_table,1),'Un',0));
        mean_std_cellThreshold_cell_noprobe = @(ik,up_or_down,ki) cell2mat(arrayfun(@(k) mean(noProbePooledResults_thresholdInfo.(Replica{ik}).cell_sugg_m_table(k,:),'omitnan')+...
            up_or_down*std(noProbePooledResults_thresholdInfo.(Replica{ik}).cell_sugg_m_table(k,:),'omitnan'),ki,'Un',0));
        quartile_cellThreshold_noprobe = @(ik,Q) cell2mat(arrayfun(@(k) quartile_function(noProbePooledResults_thresholdInfo.(Replica{ik}).cell_sugg_m_table(k,:),Q),1:size(noProbePooledResults_thresholdInfo.(Replica{ik}).cell_sugg_m_table,1),'Un',0));
        quartile_cellThreshold_cell_noprobe = @(ik,Q,ki) cell2mat(arrayfun(@(k) quartile_function(noProbePooledResults_thresholdInfo.(Replica{ik}).cell_sugg_m_table(k,:),Q),ki,'Un',0));
        med_cellThreshold_noprobe = @(ik) noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_med;
        med_cellThreshold_cell_noprobe = @(ik,ki) noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_med(ki);
        opt_cellThreshold_noprobe = @(ik) noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_opt;
        opt_cellThreshold_cell_noprobe = @(ik,ki) noProbePooledResults_thresholdInfo.(Replica{ik}).cell_threshold_opt(ki);
        for vi = 1:2*n_threshold_options
            %create handle function for threshold type to keep minimum line usage
            switch vi
                case 1 %min threshold per image
                    single_threshold_handle = @(Sof,Con,Rep,k) min_imageThreshold_image(Sof,Con,Rep,k);
                    threshold_handle = @(Sof,Con,Rep) min_imageThreshold(Sof,Con,Rep);
                    %%goes from image_and_cell_pair list find threshold for each cell in list of pairss
                    single_threshold_handle_noprobe = @(Rep,k) min_imageThreshold_image_noprobe(Rep,k);
                    threshold_handle_noprobe = @(Rep) min_imageThreshold_noprobe(Rep);
                case 2 %Q25 threshold per image
                    Percentile = 25;
                    single_threshold_handle = @(Sof,Con,Rep,k) quartile_imageThreshold_image(Sof,Con,Rep,Percentile,k);
                    threshold_handle = @(Sof,Con,Rep) quartile_imageThreshold(Sof,Con,Rep,Percentile);
                    %%goes from image_and_cell_pair list find threshold for each cell in list of pairss
                    single_threshold_handle_noprobe = @(Rep,k) quartile_imageThreshold_image_noprobe(Rep,Percentile,k);
                    threshold_handle_noprobe = @(Rep) quartile_imageThreshold_noprobe(Rep,Percentile);
                case 3 %mean threshold per image
                    spin = 0;
                    single_threshold_handle = @(Sof,Con,Rep,k) mean_std_imageThreshold_image(Sof,Con,Rep,spin,k);
                    threshold_handle = @(Sof,Con,Rep) mean_std_imageThreshold(Sof,Con,Rep,spin);
                    %%goes from image_and_cell_pair list find threshold for each cell in list of pairss
                    single_threshold_handle_noprobe = @(Rep,k) mean_std_imageThreshold_image_noprobe(Rep,spin,k);
                    threshold_handle_noprobe = @(Rep) mean_std_imageThreshold_noprobe(Rep,spin);
                case 4 %med threshold per image
                    single_threshold_handle = @(Sof,Con,Rep,k) med_imageThreshold_image(Sof,Con,Rep,k);
                    threshold_handle = @(Sof,Con,Rep) med_imageThreshold(Sof,Con,Rep);
                    %%goes from image_and_cell_pair list find threshold for each cell in list of pairss
                    single_threshold_handle_noprobe = @(Rep,k) med_imageThreshold_image_noprobe(Rep,k);
                    threshold_handle_noprobe = @(Rep) med_imageThreshold_noprobe(Rep);
                case 5 %Q75 threshold per image
                    Percentile = 75;
                    single_threshold_handle = @(Sof,Con,Rep,k) quartile_imageThreshold_image(Sof,Con,Rep,Percentile,k);
                    threshold_handle = @(Sof,Con,Rep) quartile_imageThreshold(Sof,Con,Rep,Percentile);
                    %%goes from image_and_cell_pair list find threshold for each cell in list of pairss
                    single_threshold_handle_noprobe = @(Rep,k) quartile_imageThreshold_image_noprobe(Rep,Percentile,k);
                    threshold_handle_noprobe = @(Rep) quartile_imageThreshold_noprobe(Rep,Percentile);
                case 6 %opt threshold per image
                    single_threshold_handle = @(Sof,Con,Rep,k) opt_imageThreshold_image(Sof,Con,Rep,k);
                    threshold_handle = @(Sof,Con,Rep) opt_imageThreshold(Sof,Con,Rep);
                    %%goes from image_and_cell_pair list find threshold for each cell in list of pairss
                    single_threshold_handle_noprobe = @(Rep,k) opt_imageThreshold_image_noprobe(Rep,k);
                    threshold_handle_noprobe = @(Rep) opt_imageThreshold_noprobe(Rep);
                case 7 %min threshold per cell
                    single_threshold_handle = @(Sof,Con,Rep,k) min_cellThreshold_cell(Sof,Con,Rep,k);
                    threshold_handle = @(Sof,Con,Rep) min_cellThreshold(Sof,Con,Rep);
                    %%goes from image_and_cell_pair list find threshold for each cell in list of pairss
                    single_threshold_handle_noprobe = @(Rep,k) min_cellThreshold_cell_noprobe(Rep,k);
                    threshold_handle_noprobe = @(Rep) min_cellThreshold_noprobe(Rep);
                case 8 %Q25 threshold per cell
                    Percentile = 25;
                    single_threshold_handle = @(Sof,Con,Rep,k) quartile_cellThreshold_cell(Sof,Con,Rep,Percentile,k);
                    threshold_handle = @(Sof,Con,Rep) quartile_cellThreshold(Sof,Con,Rep,Percentile);
                    %%goes from image_and_cell_pair list find threshold for each cell in list of pairs
                    single_threshold_handle_noprobe = @(Rep,k) quartile_cellThreshold_cell_noprobe(Rep,Percentile,k);
                    threshold_handle_noprobe = @(Rep) quartile_cellThreshold_noprobe(Rep,Percentile);
                case 9 %mean threshold per cell
                    spin = 0;
                    single_threshold_handle = @(Sof,Con,Rep,k) mean_std_cellThreshold_cell(Sof,Con,Rep,spin,k);
                    threshold_handle = @(Sof,Con,Rep) mean_std_cellThreshold(Sof,Con,Rep,spin);
                    %%goes from image_and_cell_pair list find threshold for each cell in list of pairs
                    single_threshold_handle_noprobe = @(Rep,k) mean_std_cellThreshold_cell_noprobe(Rep,spin,k);
                    threshold_handle_noprobe = @(Rep) mean_std_cellThreshold_noprobe(Rep,spin);
                case 10 %med threshold per cell
                    single_threshold_handle = @(Sof,Con,Rep,k) med_cellThreshold_cell(Sof,Con,Rep,k);
                    threshold_handle = @(Sof,Con,Rep) med_cellThreshold(Sof,Con,Rep);
                    single_threshold_handle_noprobe = @(Rep,k) med_cellThreshold_cell_noprobe(Rep,k);
                    threshold_handle_noprobe = @(Rep) med_cellThreshold_noprobe(Rep);
                case 11 %Q75 threshold per cell
                    Percentile = 75;
                    single_threshold_handle = @(Sof,Con,Rep,k) quartile_cellThreshold_cell(Sof,Con,Rep,Percentile,k);
                    threshold_handle = @(Sof,Con,Rep) quartile_cellThreshold(Sof,Con,Rep,Percentile);
                    %%goes from image_and_cell_pair list find threshold for each cell in list of pairs
                    single_threshold_handle_noprobe = @(Rep,k) quartile_cellThreshold_cell_noprobe(Rep,Percentile,k);
                    threshold_handle_noprobe = @(Rep) quartile_cellThreshold_noprobe(Rep,Percentile);
                case 12 %opt threshold per cell
                    single_threshold_handle = @(Sof,Con,Rep,k) opt_cellThreshold_cell(Sof,Con,Rep,k);
                    threshold_handle = @(Sof,Con,Rep) opt_cellThreshold(Sof,Con,Rep);
                    single_threshold_handle_noprobe = @(Rep,k) opt_cellThreshold_cell_noprobe(Rep,k);
                    threshold_handle_noprobe = @(Rep) opt_cellThreshold_noprobe(Rep);
            end
            if (vi<=n_threshold_options)%image threshold passed to single cells
                %goes from image_and_cell_pair list find threshold for each cell in list of pairs
                cell_threshold_handle = @(Sof,Con,Rep,unique_IC_pair_image_num)...
                    floor(cell2mat(arrayfun(@(k) single_threshold_handle(Sof,Con,Rep,k)*ones(1,sum(double(unique_IC_pair_image_num==k))),1:length(threshold_handle(Sof,Con,Rep)),'Un',0)));
                cell_threshold_handle_noprobe = @(Rep,unique_IC_pair_image_num)...
                    floor(cell2mat(arrayfun(@(k) single_threshold_handle_noprobe(Rep,k)*ones(1,sum(double(unique_IC_pair_image_num==k))),1:length(threshold_handle_noprobe(Rep)),'Un',0)));
                %gets thresholds for each spot in for each cell and image from list
                ob_num = 1;
                ob_type = vi;
            else
                cell_threshold_handle = @(Sof,Con,Rep,unique_IC_pair_image_num)...
                    floor(cell2mat(arrayfun(@(k) single_threshold_handle(Sof,Con,Rep,k),1:length(threshold_handle(Sof,Con,Rep)),'Un',0)));
                cell_threshold_handle_noprobe = @(Rep,unique_IC_pair_image_num)...
                    floor(cell2mat(arrayfun(@(k) single_threshold_handle_noprobe(Rep,k),1:length(threshold_handle_noprobe(Rep)),'Un',0)));
                ob_num = 2;
                ob_type = vi-n_threshold_options;
            end
            for ii = 1:length(SoftwareNames_Spelling)
                for ij = 1:length(unique_concentrations)
                    for ik = 1:length(Replica)
                        try
                            n_reps = softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).n_reps;
                            if (n_reps>0)
                                softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).(thresholdObj{ob_num}).(thresholdTypes{ob_type}).cell_threshold_handle  = cell_threshold_handle;
                            end
                        catch
                        end
                    end
                end
            end
            for ik = 1:length(Replica)
                try
                    n_reps = noProbePooledResults_thresholdInfo.(Replica{ik}).n_reps;
                    if (n_reps>0)
                        noProbePooledResults_thresholdInfo.(Replica{ik}).(thresholdObj{ob_num}).(thresholdTypes{ob_type}).cell_threshold_handle  = cell_threshold_handle_noprobe;
                    end
                catch
                end
            end
        end
        %% odds that no spots detected at low thresholds in a cell think near none.
        for ik = 1:length(Replica)
            n_reps = noProbePooledResults_thresholdInfo.(Replica{ik}).n_reps;
            if (n_reps>0)
                try
                    unique_image_and_cell_pair = noProbePooledResults_thresholdInfo.(Replica{ik}).unique_image_and_cell_pair;
                    for vi = 1:2*n_threshold_options
                        if (vi<=n_threshold_options)
                            ob_num = 1;
                            ob_type = vi;
                        else
                            ob_num = 2;
                            ob_type = vi-n_threshold_options;
                        end
                        cell_threshold_handle = ...
                            noProbePooledResults_thresholdInfo.(Replica{ik}).(thresholdObj{ob_num}).(thresholdTypes{ob_type}).cell_threshold_handle;
                        curr_cell_image_thresholds = cell_threshold_handle(ik,unique_image_and_cell_pair(:,1));
                        curr_cell_image_thresholds(curr_cell_image_thresholds<0) = 1;
                        curr_cell_image_thresholds(isnan(curr_cell_image_thresholds)) = 1;
                        noProbePooledResults_thresholdInfo.(Replica{ik}).(thresholdObj{ob_num}).(thresholdTypes{ob_type}).thresholds = curr_cell_image_thresholds;
                    end
                catch
                end
            end
        end
        for ii = 1:length(SoftwareNames_Spelling)
            for ij = 1:length(unique_concentrations)
                for ik = 1:length(Replica)
                    try
                        n_reps = softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).n_reps;
                        if (n_reps>0)
                            unique_image_and_cell_pair = softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).unique_image_and_cell_pair;
                            %% get id of cell for cell level stats
                            for vi = 1:2*n_threshold_options
                                if (vi<=n_threshold_options)
                                    ob_num = 1;
                                    ob_type = vi;
                                else
                                    ob_num = 2;
                                    ob_type = vi-n_threshold_options;
                                end
                                cell_threshold_handle = ...
                                    softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).(thresholdObj{ob_num}).(thresholdTypes{ob_type}).cell_threshold_handle;
                                curr_cell_image_thresholds = cell_threshold_handle(ii,ij,ik,unique_image_and_cell_pair(:,1));
                                curr_cell_image_thresholds(curr_cell_image_thresholds<0) = 1;
                                curr_cell_image_thresholds(isnan(curr_cell_image_thresholds)) = 1;
                                softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{ik}).(thresholdObj{ob_num}).(thresholdTypes{ob_type}).thresholds = curr_cell_image_thresholds;
                            end
                        end
                    catch
                    end
                end
            end
        end
        save('pooled_hARF4_RNAFISH_AnalysisThresholdInfo_ThresholdTypes.mat','noProbePooledResults_thresholdInfo','softwarePooledResults_thresholdInfo','-v7.3');
    end
end
end