function RNA_FISH_ExperimentalAnalysisStructure_V5(mmm,part)
%low use pool
%RNAThreshold has plot functions for plotting curves involved in spot
%detection
%parse results for each summary file
%check what happens when no probe image comes up
%problem with using lbl_mid for 2d cell is it is just from dapi so might
%miss cytoplasmic RNA, as gates being outside cell
%update to combine folders for ProbeTest_20240827 and ProbeTest_20240906
%For ToolThreshold need to combine multiple tools one.
%why not use spot_mask tif for the actual spot mask?
root2 = '../../hospelb/RNAFISH/ProbeTest_ThLow/JAJH_ProbeTest/IntensityDistro/nodpc/';
root4 = '../../hospelb/RNAFISH/ProbeTest_ThLow/JAJH_ProbeTest/IntensityDistro/';

switch mmm
    case 1
        root1 = '../../hospelb/RNAFISH/ProbeTest_ThLow/_hidey_hole/';
        root3 = '../../hospelb/RNAFISH/Analysis/_hidey_hole/';
        root1a = '../../hospelb/RNAFISH/ProbeTest_ThLow/JAJH_ProbeTest/';
        root3a = '../../hospelb/RNAFISH/Analysis/JAJH_ProbeTest/';
        foldName = 'currentHidden';
    case 2
        root1 = '../../hospelb/RNAFISH/ProbeTest_ThLow/JAJH_ProbeTest/';
        root3 = '../../hospelb/RNAFISH/Analysis/JAJH_ProbeTest/';
        root1a = '../../hospelb/RNAFISH/ProbeTest_ThLow/_hidey_hole/';
        root3a = '../../hospelb/RNAFISH/Analysis/_hidey_hole/';
        foldName = 'currentMain';
end
summarys0 = dir(root1);
all_summary_files = {summarys0([summarys0(:).isdir]==1).name}';
all_summary_files(ismember(all_summary_files,{'.','..','IntensityDistro'})) = [];
all_summary_files(contains(all_summary_files,'70C')) = [];
all_summary_files(contains(all_summary_files,'75C')) = [];
all_summary_files(contains(all_summary_files,'80C')) = [];
all_summary_files(contains(all_summary_files,'85C')) = [];
all_summary_files(contains(all_summary_files,'500')) = [];
all_summary_files(contains(all_summary_files,'R1')) = [];
all_summary_files(contains(all_summary_files,'R2')) = [];
all_summary_files(contains(all_summary_files,'R3')) = [];
all_summary_files(contains(all_summary_files,'R4')) = [];
results = struct();
results_cellseg = struct();
results_callTable = struct();
results_quant = struct();
results_cellStats = struct();
results_spotStats = struct();
results_ImageThresholdsStats = struct();
results_CellThresholdsStats = struct();
results_imgData = struct();
xdim = 2048; ydim = 2048; zdim = 67;
if (part==1)
    save(strcat('hARF4_RNAFISH_Analysis_',foldName,'.mat'),'results','-v7.3');
    save(strcat('hARF4_RNAFISH_imageThresholdsStats_',foldName,'.mat'),'results_ImageThresholdsStats','-v7.3');
end
if (part==2)
    save(strcat('hARF4_RNAFISH_CellSegStats_',foldName,'.mat'),'results_cellseg','-v7.3');
end
if (part==3)
    save(strcat('hARF4_RNAFISH_AnalysisCallTables_',foldName,'.mat'),'results_callTable','-v7.3');
    save(strcat('hARF4_RNAFISH_CellThresholdStats_',foldName,'.mat'),'results_CellThresholdsStats','-v7.3');
end
if (part==4)
    save(strcat('hARF4_RNAFISH_AnalysisQuantResults_',foldName,'.mat'),'results_quant','-v7.3');
end
if (part==5)
    save(strcat('hARF4_RNAFISH_AnalysisSpotStats_',foldName,'.mat'),'results_spotStats','-v7.3');
end
if (part==6)
    save(strcat('hARF4_RNAFISH_AnalysisCellStats_',foldName,'.mat'),'results_cellStats','-v7.3');
end
if (part==7)
    save(strcat('hARF4_RNAFISH_AnalysisImgData_',foldName,'.mat'),'results_imgData','-v7.3');
end
for f = 1:length(all_summary_files)
    try
        results(f).file = all_summary_files{f};
        results_callTable(f).file = all_summary_files{f};
        results_quant(f).file = all_summary_files{f};
        results_cellStats(f).file = all_summary_files{f};
        results_spotStats(f).file = all_summary_files{f};
        results_ImageThresholdsStats(f).file = all_summary_files{f};
        results_CellThresholdsStats(f).file = all_summary_files{f};
        results_cellseg(f).file = all_summary_files{f};
        results_imgData(f).file = all_summary_files{f};
        cellSegFile = strcat('CellSeg_',all_summary_files{f},'.mat');
        load([root1 all_summary_files{f} '/CH3/' all_summary_files{f} '_CH3_spotCall_rnaspotsrun.mat'],'run_info');
        results(f).run_info = run_info;
        results(f).type_target = run_info.meta.type_target;
        results(f).type_targetmol = run_info.meta.type_targetmol;
        results(f).is_no_probe = run_info.meta.noProbe_flag;
        results_ImageThresholdsStats(f).sugg_m = run_info.threshold_results.pool.sugg_m;
        results_ImageThresholdsStats(f).sugg_f = run_info.threshold_results.pool.sugg_f;
        results_ImageThresholdsStats(f).sugg_fri = run_info.threshold_results.pool.sugg_fri;
        results_ImageThresholdsStats(f).mean_overall = run_info.threshold_results.thstats.mean_overall;
        results_ImageThresholdsStats(f).med_overall = run_info.threshold_results.thstats.med_overall;
        results_ImageThresholdsStats(f).std_overall = run_info.threshold_results.thstats.std_overall;
        results_ImageThresholdsStats(f).min_overall = run_info.threshold_results.thstats.min_overall;
        results_ImageThresholdsStats(f).max_overall = run_info.threshold_results.thstats.max_overall;
        results_ImageThresholdsStats(f).threshold = run_info.intensity_threshold;
        results_ImageThresholdsStats(f).window_scores = run_info.threshold_results.window_scores';
        results_ImageThresholdsStats(f).threshold_curve_spline_info = run_info.threshold_results.test_winsc;
        results_ImageThresholdsStats(f).threshold_score_curve.x = [];
        results_ImageThresholdsStats(f).threshold_score_curve.score = [];
        results_ImageThresholdsStats(f).threshold_score_curve.score_m = [];
        results_ImageThresholdsStats(f).threshold_score_curve.score_f = [];
        results_ImageThresholdsStats(f).threshold_score_curve.score_fri = [];
        results_ImageThresholdsStats(f).threshold_spot_count_curve.x = [];
        results_ImageThresholdsStats(f).threshold_spot_count_curve.y = [];
        results_ImageThresholdsStats(f).spot_count_curve.x = [];
        results_ImageThresholdsStats(f).spot_count_curve.y = [];
        results_ImageThresholdsStats(f).bkg_masked_spot_count_curve.x = [];
        results_ImageThresholdsStats(f).bkg_masked_spot_count_curve.y = [];
        try
            results_ImageThresholdsStats(f).threshold_score_curve.x = run_info.threshold_results.thstats.value_table.threshold_value';
            results_ImageThresholdsStats(f).threshold_score_curve.score = run_info.threshold_results.thstats.value_table.score';
            results_ImageThresholdsStats(f).threshold_score_curve.score_m = run_info.threshold_results.thstats.value_table.score_m';
            results_ImageThresholdsStats(f).threshold_score_curve.score_f = run_info.threshold_results.thstats.value_table.score_f';
            results_ImageThresholdsStats(f).threshold_score_curve.score_fri = run_info.threshold_results.thstats.value_table.score_fri';
        catch
        end
        try
            results_ImageThresholdsStats(f).threshold_spot_count_curve.x = run_info.threshold_results.x';
            results_ImageThresholdsStats(f).threshold_spot_count_curve.y = run_info.threshold_results.spot_counts';
        catch
        end
        try
            load([root3 all_summary_files{f} '/CH3/' all_summary_files{f} '_CH3_spotCall_bkgmasked_spotTable.mat'],'spot_table')
            bkgmasked_spot_table = spot_table; clear spot_table
            results_ImageThresholdsStats(f).bkg_masked_spot_count_curve.x = bkgmasked_spot_table(:,1)';
            results_ImageThresholdsStats(f).bkg_masked_spot_count_curve.y = bkgmasked_spot_table(:,2)';
        catch
            try
                load([root3 all_summary_files{f} '/CH3/' all_summary_files{f} '_CH3_spotCall_bkgmasked_spotTable.mat'],'spot_table')
                bkgmasked_spot_table = spot_table; clear spot_table
                results_ImageThresholdsStats(f).bkg_masked_spot_count_curve.x = bkgmasked_spot_table(:,1)';
                results_ImageThresholdsStats(f).bkg_masked_spot_count_curve.y = bkgmasked_spot_table(:,2)';
            catch
            end
        end
        try
            load([root3 all_summary_files{f} '/CH3/' all_summary_files{f} '_CH3_spotCall_spotTable.mat'],'spot_table')
            results_ImageThresholdsStats(f).spot_count_curve.x = spot_table(:,1)';
            results_ImageThresholdsStats(f).spot_count_curve.y = spot_table(:,2)';
        catch
            try
                load([root3a all_summary_files{f} '/CH3/' all_summary_files{f} '_CH3_spotCall_spotTable.mat'],'spot_table')
                results_ImageThresholdsStats(f).spot_count_curve.x = spot_table(:,1)';
                results_ImageThresholdsStats(f).spot_count_curve.y = spot_table(:,2)';
            catch
            end
        end
        results_ImageThresholdsStats(f).thresholdresultParameters = run_info.th_params;
    catch e
        disp(e)
        for v = 1:length(e.stack)
            e.stack(v)
        end
        all_summary_files{f}
    end
    if (part==1)
        save(strcat('hARF4_RNAFISH_Analysis_',foldName,'.mat'),'results','-v7.3');
        save(strcat('hARF4_RNAFISH_imageThresholdsStats_',foldName,'.mat'),'results_ImageThresholdsStats','-v7.3');
    end
    if (part==2)
        results_cellseg(f).cellSegParameters = [];
        results_cellseg(f).nucleiSegParameters = [];
        results_cellseg(f).cellSeg = [];
        results_cellseg(f).x_trim = [];
        results_cellseg(f).y_trim = [];
        results_cellseg(f).cell_z_min = [];
        results_cellseg(f).cell_z_max = [];
        results_cellseg(f).nuclei_z_min = [];
        results_cellseg(f).nuclei_z_max = [];
        results_cellseg(f).min_cell_size = [];
        results_cellseg(f).max_cell_size = [];
        results_cellseg(f).min_nuclei_size = [];
        results_cellseg(f).max_nuclei_size = [];
        try  %just uses Image Specific folder
            try
                load([root3 all_summary_files{f} '/' cellSegFile],'nucleiSeg','cellSeg')
            catch
                load([root3a all_summary_files{f} '/' cellSegFile],'nucleiSeg','cellSeg')
            end
            cellsegfields = fieldnames(cellSeg);
            for ii = 1:13
                cellSegParameters.(cellsegfields{ii}) = cellSeg.(cellsegfields{ii});
            end
            results_cellseg(f).cellSegParameters = cellSegParameters; clear cellSegParameters
            results_cellseg(f).nucleiSegParameters = nucleiSeg.params;
            results_cellseg(f).cellSeg = cellSeg;
            results_cellseg(f).x_trim = cellSeg.x_trim;
            results_cellseg(f).y_trim = cellSeg.y_trim;
            results_cellseg(f).cell_z_min = cellSeg.z_min;
            results_cellseg(f).cell_z_max = cellSeg.z_max;
            results_cellseg(f).nuclei_z_min = nucleiSeg.params.z_min;
            results_cellseg(f).nuclei_z_max = nucleiSeg.params.z_max;
            results_cellseg(f).min_cell_size = cellSeg.min_cell_size;
            results_cellseg(f).max_cell_size = cellSeg.max_cell_size;
            results_cellseg(f).min_nuclei_size = nucleiSeg.params.min_nucleus_size;
            results_cellseg(f).max_nuclei_size = nucleiSeg.params.max_nucleus_size;
        catch e
            disp(e)
            for v = 1:length(e.stack)
                e.stack(v)
            end
        end
        clear nucleiSeg cellSeg
        save(strcat('hARF4_RNAFISH_CellSegStats_',foldName,'.mat'),'results_cellseg','-v7.3');
    end
    if (part==3)
        results_callTable(f).call_table = [];
        results_CellThresholdsStats(f).extended_callTable_spot_count_curve_in_cell_X_y = [];
        results_CellThresholdsStats(f).scThresholdSuggestions = [];
        %if no cellSeg file then
        try
            load([root3 all_summary_files{f} '/' cellSegFile],'cellSeg')
        catch
            load([root3a all_summary_files{f} '/' cellSegFile],'cellSeg')
        end
        try
            try
                load([root3 all_summary_files{f} '/CH3/' all_summary_files{f} '_CH3_spotCall_callTable.mat'],'call_table')
                results_callTable(f).call_table = call_table;
                results_callTable(f).call_table = RNACoords.applyCellSegMask(call_table, cellSeg.cell_mask);
            catch
                try
                    load([root3a all_summary_files{f} '/CH3/' all_summary_files{f} '_CH3_spotCall_callTable.mat'],'call_table')
                    results_callTable(f).call_table = call_table;
                    results_callTable(f).call_table = RNACoords.applyCellSegMask(call_table, cellSeg.cell_mask);
                catch  e
                    disp(e)
                    for v = 1:length(e.stack)
                        e.stack(v)
                    end
                end
            end
            current_callTable = results_callTable(f).call_table;
            if (~isempty(current_callTable))
                in_cell_X = @(x) find(double(table2array(current_callTable(:,'cell'))==x));
                threshold_spot_count_curve_x = run_info.threshold_results.x';
                n_cells = length(setdiff(unique(table2array(current_callTable(:,'cell'))),0));
                results_CellThresholdsStats(f).CellThresholdsExist = 0;
                if (n_cells>0)
                    results_CellThresholdsStats(f).CellThresholdsExist = 1;
                    extended_callTable_spot_count_curve_in_cell_X_y = cell2mat(arrayfun(@(Y) ...
                        numel(in_cell_X(Y))-histcounts(table2array(current_callTable(in_cell_X(Y),'dropout_thresh')),'Normalization','cumcount','BinEdges',[-Inf 1:length(threshold_spot_count_curve_x)])',1:n_cells,'Un',0))';
                    param_struct_vector = arrayfun(@(x) RNAThreshold.genEmptyThresholdParamStruct(),1:n_cells,'Un',0);
                    for nn = 1:n_cells
                        param_struct_vector{nn} = run_info.th_params;
                        param_struct_vector{nn}.sample_spot_table = [threshold_spot_count_curve_x' extended_callTable_spot_count_curve_in_cell_X_y(nn,:)'];
                    end
                    scThresholdSuggestions = arrayfun(@(nn) RNAThreshold.scoreThresholdSuggestions(RNAThreshold.estimateThreshold(param_struct_vector{nn})),1:n_cells,'Un',0);
                    scThresholdSuggestions =  [scThresholdSuggestions{:}];
                    subfield_groups = {'pool','thstats'};
                    for v = 1:length(subfield_groups)
                        subfields = fieldnames(scThresholdSuggestions(1).(subfield_groups{v}));
                        for subf = 1:length(subfields)
                            subf_vals = arrayfun(@(nn) scThresholdSuggestions(nn).(subfield_groups{v}).(subfields{subf}),1:size(scThresholdSuggestions,2),'UniformOutput',0);
                            [scThresholdSuggestions.(subfields{subf})] = subf_vals{:};
                        end
                    end
                    scThresholdSuggestions = rmfield(scThresholdSuggestions,subfield_groups);
                    results_CellThresholdsStats(f).extended_callTable_spot_count_curve_in_cell_X_y =extended_callTable_spot_count_curve_in_cell_X_y;
                    if (size(scThresholdSuggestions,2)>1)
                        results_CellThresholdsStats(f).scThresholdSuggestions = struct2table(scThresholdSuggestions);
                    else
                        results_CellThresholdsStats(f).scThresholdSuggestions = struct2table(scThresholdSuggestions,'AsArray',1);
                    end
                end
            end
        catch e
            disp(e)
            for v = 1:length(e.stack)
                e.stack(v)
            end
            all_summary_files{f}
        end
        clear bkgmasked_call_table call_table_cell call_table cellSeg
        save(strcat('hARF4_RNAFISH_AnalysisCallTables_',foldName,'.mat'),'results_callTable','-v7.3');
        clear param_struct_vector extended_callTable_spot_count_curve_in_cell_X_y scThresholdSuggestions
        save(strcat('hARF4_RNAFISH_CellThresholdStats_',foldName,'.mat'),'results_CellThresholdsStats','-v7.3');
    end
    if (part==4)
        results_quant(f).gaussian_radius = [];
        results_quant(f).quantresultParameters = [];
        results_quant(f).quant_results_spotTable = [];
        results_quant(f).spotcount_nuc = [];
        results_quant(f).spotcount_cyt = [];
        results_quant(f).spotcount_total = [];
        results_quant(f).nucNascentCount = [];
        results_quant(f).nucCount = [];
        results_quant(f).cytoCount = [];
        results_quant(f).nucNascentCloud = [];
        results_quant(f).nucCloud = [];
        results_quant(f).cytoCloud = [];
        results_quant(f).num_cells_withoutspots =  [];
        try
            results_quant(f).QuantStatsExist = 0;
            results_quant(f).NoSpots = 1;
            try
                load([root1 all_summary_files{f} '/CH3/' all_summary_files{f} '_CH3_quantData.mat'],'quant_results')
            catch
                load([root1a all_summary_files{f} '/CH3/' all_summary_files{f} '_CH3_quantData.mat'],'quant_results')
            end
            results_quant(f).gaussian_radius = quant_results.gaussian_radius;
            quantresultsfields = fieldnames(quant_results);
            for ii = [1:16 18 19]
                quantresultsParameters.(quantresultsfields{ii}) = quant_results.(quantresultsfields{ii});
            end
            results_quant(f).quantresultParameters = quantresultsParameters;
            if (size(quant_results.cellData,2)>0)
                results_quant(f).spotcount_nuc = [quant_results.cellData.spotcount_nuc];
                results_quant(f).spotcount_cyt = [quant_results.cellData.spotcount_total]-[quant_results.cellData.spotcount_nuc];
                results_quant(f).spotcount_total = [quant_results.cellData.spotcount_total];
                results_quant(f).nucNascentCount = [quant_results.cellData.nucNascentCount];
                results_quant(f).nucCount = [quant_results.cellData.nucCount];
                results_quant(f).cytoCount = [quant_results.cellData.cytoCount];
                results_quant(f).nucNascentCloud = [quant_results.cellData.cytoCloud];
                results_quant(f).nucCloud = [quant_results.cellData.nucCloud];
                results_quant(f).cytoCloud = [quant_results.cellData.nucNascentCloud];
                celldata = quant_results.cellData;
                results_quant(f).num_cells_withoutspots = sum(arrayfun(@(x) isempty(celldata(x).spotTable),1:size(celldata,2))==1);
                cell_spotDataExists = find(cell2mat(cellfun(@isempty, {celldata.spotTable},'Un',0))==0);
                celldata_spotZFits = [celldata(cell_spotDataExists).spotZFits];
                celldata_spotTable = vertcat(celldata(cell_spotDataExists).spotTable);
                cell_num = cell2mat(arrayfun(@(n) double(celldata(n).cell_number)*ones(1,celldata(n).spotcount_total),cell_spotDataExists,'Un',0))';
                celldata_cell_loc = arrayfun(@(vv) celldata(vv).cell_loc,cell_num,'Un',0);
                if (~isempty(cell_spotDataExists))
                celldata_spotTable(:,'cell') = array2table(cell_num);
                celldata_spotTable(:,'spotZFits') = cell2table(celldata_spotZFits');
                celldata_spotTable(:,'cell_loc') = cell2table(celldata_cell_loc);
                results_quant(f).NoSpots = 0;
                end
                results_quant(f).quant_results_spotTable = celldata_spotTable;
                results_quant(f).QuantStatsExist = 1;
            end
        catch e
            disp(e)
            for v = 1:length(e.stack)
                e.stack(v)
            end
            all_summary_files{f}
        end
        clear quant_results celldata celldata_spotTable cell_num celldata_cell_loc cell_spotDataExists
        save(strcat('hARF4_RNAFISH_AnalysisQuantResults_',foldName,'.mat'),'results_quant','-v7.3');
    end
    if (part==5)
        results_spotStats(f).intensityStatThresholdResults = [];
        results_spotStats(f).ImageSpotStats = [];
        try
            load([root2 all_summary_files{f} '_CH3_istats.mat'],'intensityStats');
        catch
        end
        try
            imThresholdStat = intensityStats.threshold_results;
            subfield_groups = {'pool','thstats'};
            for v = 1:length(subfield_groups)
                subfields = fieldnames(imThresholdStat(1).(subfield_groups{v}));
                for subf = 1:length(subfields)
                    subf_vals = arrayfun(@(nn) imThresholdStat(nn).(subfield_groups{v}).(subfields{subf}),1:size(imThresholdStat,2),'UniformOutput',0);
                    [imThresholdStat.(subfields{subf})] = subf_vals{:};
                end
            end
            imThresholdStat = rmfield(imThresholdStat,subfield_groups);
            results_spotStats(f).intensityStatThresholdResults = imThresholdStat;
        catch
        end
        try
            results_spotStats(f).SpotStatsExist = 0;
            if (size(intensityStats.spotProfile.spotData,1)>0)
                ImageSpotStats = intensityStats.spotProfile.spotData;
                xi = double(ImageSpotStats.x)-xdim/2;
                yi = double(ImageSpotStats.y)-ydim/2;
                zi = double(ImageSpotStats.z)-zdim/2;
                ri = sqrt(xi.^2+yi.^2+zi.^2);
                thetai = acosd(zi./ri);
                phii = sign(yi).*acosd(xi./sqrt(xi.^2+yi.^2));
                rp = sqrt(xi.^2+yi.^2);
                thetap = atand(yi./xi);
                ImageSpotStats.spot_polar_r = rp;
                ImageSpotStats.spot_sphere_r = ri;
                ImageSpotStats.spot_polar_theta = thetap;
                ImageSpotStats.spot_sphere_theta = thetai;
                ImageSpotStats.spot_sphere_phi = phii;
                results_spotStats(f).ImageSpotStats = ImageSpotStats;
                results_spotStats(f).SpotStatsExist = 1;
            end
        catch e
            disp(e)
            for v = 1:length(e.stack)
                e.stack(v)
            end
            all_summary_files{f}
        end
        clear ImageSpotStats intensityStats imThresholdStat
        try
            load([root4 all_summary_files{f} '_CH3_istats.mat'],'intensityStats');
        catch
        end
        try
            imThresholdStat = intensityStats.threshold_results;
            subfield_groups = {'pool','thstats'};
            for v = 1:length(subfield_groups)
                subfields = fieldnames(imThresholdStat(1).(subfield_groups{v}));
                for subf = 1:length(subfields)
                    subf_vals = arrayfun(@(nn) imThresholdStat(nn).(subfield_groups{v}).(subfields{subf}),1:size(imThresholdStat,2),'UniformOutput',0);
                    [imThresholdStat.(subfields{subf})] = subf_vals{:};
                end
            end
            imThresholdStat = rmfield(imThresholdStat,subfield_groups);
            results_spotStats(f).intensityStatThresholdResults = imThresholdStat;
        catch
        end
        try
            results_spotStats(f).SpotStatsExist_IntCorrected = 0;
            if (size(intensityStats.spotProfile.spotData,1)>0)
                ImageSpotStats = intensityStats.spotProfile.spotData;
                xi = double(ImageSpotStats.x)-xdim/2;
                yi = double(ImageSpotStats.y)-ydim/2;
                zi = double(ImageSpotStats.z)-zdim/2;
                ri = sqrt(xi.^2+yi.^2+zi.^2);
                thetai = acosd(zi./ri);
                phii = sign(yi).*acosd(xi./sqrt(xi.^2+yi.^2));
                rp = sqrt(xi.^2+yi.^2);
                thetap = atand(yi./xi);
                ImageSpotStats.spot_polar_r = rp;
                ImageSpotStats.spot_sphere_r = ri;
                ImageSpotStats.spot_polar_theta = thetap;
                ImageSpotStats.spot_sphere_theta = thetai;
                ImageSpotStats.spot_sphere_phi = phii;
                results_spotStats(f).ImageSpotStats_IntCorrected = ImageSpotStats;
                results_spotStats(f).SpotStatsExist_IntCorrected = 1;
            end
        catch e
            disp(e)
            for v = 1:length(e.stack)
                e.stack(v)
            end
            all_summary_files{f}
        end
        clear ImageSpotStats intensityStats imThresholdStat
        save(strcat('hARF4_RNAFISH_AnalysisSpotStats_',foldName,'.mat'),'results_spotStats','-v7.3');
    end
    if (part==6)
        results_cellStats(f).ImageCellStats = [];
        try
            load([root2 all_summary_files{f} '_CH3_istats.mat'],'intensityStats');
        catch
        end
        try
            results_cellStats(f).image_cambkg_Full = intensityStats.statsRawFull.ibkg;
            results_cellStats(f).image_cambkg_Trim = intensityStats.statsRawZTrim.ibkg;
            results_cellStats(f).image_cellbkg_Full = intensityStats.statsRawFull.cbkg;
            results_cellStats(f).image_cellbkg_Trim = intensityStats.statsRawZTrim.cbkg;
        catch e
            disp(e)
            for v = 1:length(e.stack)
                e.stack(v)
            end
            all_summary_files{f}
        end
        try
            results_cellStats(f).CellStatsExist = 0;
            if (size(intensityStats.cellProfile,1)>0)
                ImageCellStats = intensityStats.cellProfile;
                cell_centroid_x = intensityStats.cellProfile.centroid_x;
                cell_centroid_y = intensityStats.cellProfile.centroid_y;
                cxi = cell_centroid_x-xdim/2;
                cyi = cell_centroid_y-ydim/2;
                crp = sqrt(cxi.^2+cyi.^2);
                cthetap = atand(cyi./cxi);
                ImageCellStats(:,'cell_polar_r') = array2table(crp);
                ImageCellStats(:,'cell_polar_theta') = array2table(cthetap);
                results_cellStats(f).ImageCellStats = ImageCellStats;
                results_cellStats(f).CellStatsExist = 1;
            end
        catch e
            disp(e)
            for v = 1:length(e.stack)
                e.stack(v)
            end
            all_summary_files{f}
        end
        clear intensityStats ImageCellStats
        try
            load([root4 all_summary_files{f} '_CH3_istats.mat'],'intensityStats');
        catch
        end
        try
            results_cellStats(f).image_cambkg_Full_IntCorrected = intensityStats.statsRawFull.ibkg;
            results_cellStats(f).image_cambkg_Trim_IntCorrected = intensityStats.statsRawZTrim.ibkg;
            results_cellStats(f).image_cellbkg_Full_IntCorrected = intensityStats.statsRawFull.cbkg;
            results_cellStats(f).image_cellbkg_Trim_IntCorrected = intensityStats.statsRawZTrim.cbkg;
        catch e
            disp(e)
            for v = 1:length(e.stack)
                e.stack(v)
            end
            all_summary_files{f}
        end
        try
            results_cellStats(f).CellStatsExist_IntCorrected = 0;
            if (size(intensityStats.cellProfile,1)>0)
                ImageCellStats = intensityStats.cellProfile;
                cell_centroid_x = intensityStats.cellProfile.centroid_x;
                cell_centroid_y = intensityStats.cellProfile.centroid_y;
                cxi = cell_centroid_x-xdim/2;
                cyi = cell_centroid_y-ydim/2;
                crp = sqrt(cxi.^2+cyi.^2);
                cthetap = atand(cyi./cxi);
                ImageCellStats(:,'cell_polar_r') = array2table(crp);
                ImageCellStats(:,'cell_polar_theta') = array2table(cthetap);
                results_cellStats(f).ImageCellStats_IntCorrected = ImageCellStats;
                results_cellStats(f).CellStatsExist_IntCorrected = 1;
            end
        catch e
            disp(e)
            for v = 1:length(e.stack)
                e.stack(v)
            end
            all_summary_files{f}
        end
        clear intensityStats ImageCellStats
        save(strcat('hARF4_RNAFISH_AnalysisCellStats_',foldName,'.mat'),'results_cellStats','-v7.3');
    end
    if (part==7)
        try
            run_info.channels.dna_ch = 2;
            imgPath = run_info.paths.img_path;
            csPath = run_info.paths.cellseg_path;
            results_imgData(f).nucmask_max_z = [];
            results_imgData(f).nucmask_opt_z = [];
            results_imgData(f).nucmask_mid_z= [];
            results_imgData(f).TRANS_max_z = [];
            results_imgData(f).TRANS_opt_z = [];
            results_imgData(f).TRANS_mid_z = [];
            results_imgData(f).DAPI_max_z = [];
            results_imgData(f).DAPI_opt_z = [];
            results_imgData(f).DAPI_mid_z = [];
            results_imgData(f).PROBE_max_z = [];
            results_imgData(f).PROBE_opt_z = [];
            results_imgData(f).PROBE_mid_z = [];
            try
                if ~isempty(csPath)
                    if isfile(csPath)
                        cell_mask = CellSeg.openCellMask(csPath);
                        cell_mask = (cell_mask > 0);
                        nuc_mask = CellSeg.openNucMask(csPath);
                    end
                end
                [chn,~] =  LoadTif(imgPath,run_info.channels.total_ch,run_info.channels.light_ch,1);
                [y,x] = find(cell_mask);
                contrast_z = cell2mat(arrayfun(@(zi) range(chn{run_info.channels.light_ch,1}(sub2ind(size(chn{run_info.channels.light_ch,1}),y,x,zi*ones(length(x),1))),'all'),1:run_info.dims.z_max,'Un',0));
                contrast_zs = smooth(contrast_z,3);
                contrast_trim1 = 10:length(contrast_zs)-10;
                contrast_trim2 = contrast_trim1(islocalmax(contrast_zs(contrast_trim1)));
                contrast_trim = min(contrast_trim2):max(contrast_trim2);
                zlayer_middle = round(run_info.dims.idims_sample.z ./ 2);
                zlayer_opt_contrast = contrast_trim(find(contrast_zs(contrast_trim)==min(contrast_zs(contrast_trim)),1));
                results_imgData(f).contrast_z = contrast_z;
                results_imgData(f).contrast_zsmooth = contrast_zs;
                results_imgData(f).zopt = zlayer_opt_contrast;
                results_imgData(f).zmid = zlayer_middle;
                results_imgData(f).cell_mask = cell_mask;
                results_imgData(f).nucmask_max_z = max(nuc_mask,[],3);
                results_imgData(f).TRANS_max_z = max(chn{run_info.channels.light_ch,1},[],3);
                results_imgData(f).nucmask_mid_z= squeeze(nuc_mask(:,:,zlayer_middle));
                results_imgData(f).TRANS_mid_z = squeeze(chn{run_info.channels.light_ch,1}(:,:,zlayer_middle));
                if (~isempty(zlayer_opt_contrast))
                    results_imgData(f).nucmask_opt_z = squeeze(nuc_mask(:,:,zlayer_opt_contrast));
                    results_imgData(f).TRANS_opt_z = squeeze(chn{run_info.channels.light_ch,1}(:,:,zlayer_opt_contrast));
                end
                clear chn
                [chn,~] =  LoadTif(imgPath,run_info.channels.total_ch,run_info.channels.dna_ch,1);
                results_imgData(f).DAPI_max_z = max(chn{run_info.channels.dna_ch,1},[],3);
                results_imgData(f).DAPI_mid_z  = squeeze(chn{run_info.channels.dna_ch,1}(:,:,zlayer_middle));
                if (~isempty(zlayer_opt_contrast))
                    results_imgData(f).DAPI_opt_z  = squeeze(chn{run_info.channels.dna_ch,1}(:,:,zlayer_opt_contrast));
                end
                clear chn
                [chn,~] =  LoadTif(imgPath,run_info.channels.total_ch,run_info.channels.rna_ch,1);
                try
                    try
                        load([root1 all_summary_files{f} '/CH3/' all_summary_files{f} '_CH3_quantData.mat'],'quant_results')
                    catch
                        load([root1a all_summary_files{f} '/CH3/' all_summary_files{f} '_CH3_quantData.mat'],'quant_results')
                    end
                    results_imgData(f).PROBE_max_z = max(chn{run_info.channels.rna_ch,1},[],3);
                    results_imgData(f).PROBE_mid_z = squeeze(chn{run_info.channels.rna_ch,1}(:,:,zlayer_middle));
                    if (~isempty(zlayer_opt_contrast))
                        results_imgData(f).PROBE_opt_z = squeeze(chn{run_info.channels.rna_ch,1}(:,:,zlayer_opt_contrast));
                    end
                    %RNAQuant requires RNACloud, RNASpot RNAUtils,
                    %RNA_Clouds, RNA_Threshold_Common SingleCell
                    % [img_clean, ~, ~, ~, ~] = RNAQuant.RNAProcess3Dim(chn{run_info.channels.rna_ch,1}, quant_results.threshold, quant_results.small_obj_size,quant_results.connect_size, quant_results.gaussian_radius);
                    % results_imgData(f).PROBE_max_z = max(img_clean,[],3);
                    % results_imgData(f).PROBE_mid_z = squeeze(img_clean(:,:,zlayer_middle));
                    % if (~isempty(zlayer_opt_contrast))
                    %     results_imgData(f).PROBE_opt_z = squeeze(img_clean(:,:,zlayer_opt_contrast));
                    % end
                    % results_imgData(f).PROBE_dpcleaned = 1; clear img_clean
                catch
                    results_imgData(f).PROBE_max_z = max(chn{run_info.channels.rna_ch,1},[],3);
                    results_imgData(f).PROBE_mid_z = squeeze(chn{run_info.channels.rna_ch,1}(:,:,zlayer_middle));
                    if (~isempty(zlayer_opt_contrast))
                        results_imgData(f).PROBE_opt_z = squeeze(chn{run_info.channels.rna_ch,1}(:,:,zlayer_opt_contrast));
                    end
                    results_imgData(f).PROBE_dpcleaned = 0;
                end
                clear chn
            catch e
                disp(e)
                for v = 1:length(e.stack)
                    e.stack(v)
                end
                all_summary_files{f}
            end
        catch e
            disp(e)
            for v = 1:length(e.stack)
                e.stack(v)
            end
            all_summary_files{f}
        end
        clear cell_mask nuc_mask zlayer_middle zlayer_opt_contrast quant_results
        clear contrast_zc contrast_z contrast_zs contrast_trim1contrast_trim2 contrast_trim
        save(strcat('hARF4_RNAFISH_AnalysisImgData_',foldName,'.mat'),'results_imgData','-v7.3');
    end
    clear run_info
end


end





