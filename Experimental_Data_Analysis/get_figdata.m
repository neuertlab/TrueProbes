%https://www.mathworks.com/matlabcentral/answers/1985039-size-plot-and-legend-position-programmatically-for-publication     
%outlier cell vol is done per replica condition not of all cells at once as
%a group.
%[p,tbl,stats] = kruskalwallis(x,[],'off')
%re-check if badImages are bad, and results with them in it?
% 1 counts number of spots, 2, counts per cell, 3 per on cell cells
max_bkg_cutoff = 1000;
%look at ImageCellStats max background is super high.
%TrueProbes 1 maybe one cell, has some overlap with debris.


%images where 1) bad cells signal, but not air bubble.
% throw out cells within specific region, or detecting (background)
%variant of cell 
%if filter out abnormally small spots xgw and ygw do by software where
%number of probes could affect result
%error where spots in bad segmented cell, 
zlow = 0;
zhigh = 70;
xlow = 512-1;
xhigh = 1024+512+1;
ylow = 512-1;
yhigh = 1024+512+1;
xlow = 0;
ylow = 0;
xhigh = 2049;
yhigh = 2049;
xgwlow = 0.6;
xgwhigh = 40;
ygwlow = 0.6;
ygwhigh = 40;
%are cases where cells are at top stacks.
pAlpha_threshold = 0.05;

%xgw and ygw >1
no_probe_conditions = {'NoProbes'};
SoftwareNames = {'TrueProbes','Stellaris','OligostanHT','PaintSHOP','MERFISH'};
SoftwareNames_Spelling = {{'trueprobes'},{'stellaris'},{'oligostan','ogliostan'},{'paintshop'},{'merfish'}};
Conditions = {SoftwareNames{:},no_probe_conditions{:}};
int_level_types = {'max','mean','median'};
thresholdObj = {'Image','Cell'};
intType = {'MaxInt','MeanInt'};
fixType = {'UnfixedThr','fixedThrReplica','fixThrAll'};
illumType = {'UnfixedIllum','fixedllum'};
thresholdTypes = {'min','q1','mean','median','q3','opt'};
SubRNA_Categories = {'All','Nascent_Nuclear','Mature_Cytoplasmic','Mature Nuclear','Nuclear','Cytoplasmic','Mature'};
bkgType = {'Cell','Local'};
Replica = {'Rep1','Rep2','Rep3','Rep4','Rep5A','Rep5B','Rep6A','Rep6B'};
replica_num = 1:length(Replica);
concentrations_label = {'1:1000','1:2000','1:4000','Mol Adjusted 1:285','MolAdj 1:2000'};
unique_concentrations = {'_1000','_2000','_4000','MolAdj','NewCalc_MolAdj'};
group_concentrations = {'Cy5_1000','Cy5_2000','Cy5_4000','MolAdj','NewMolAdj'};
nascent_status = {[0 1],[1],[0],[0],[0 1],[0 1],[0]};
nuc_status =        {[0 1],[1],[0],[1],[1],[0],[0 1]};

y_vector = @(ik,var) cell2mat(arrayfun(@(x)var{x,ik}',1:length(SoftwareNames),'Un',0));
g1_vector = @(ik,var) string(SoftwareNames(cell2mat(arrayfun(@(x)x*ones(size(var{x,ik}))',1:length(SoftwareNames),'Un',0))));
g2_vector = @(ik,var) string(Replica(Replicas_To_Show(ik)*cell2mat(arrayfun(@(x) ones(size(var{x,ik}))',1:length(SoftwareNames),'Un',0))));
y_multi_vector = @(var) cell2mat(arrayfun(@(ik) cell2mat(arrayfun(@(x)var{x,ik}',1:length(SoftwareNames),'Un',0)),1:length(Replicas_To_Show),'Un',0));
g1_multi_vector = @(var) string(SoftwareNames(cell2mat(arrayfun(@(ik) cell2mat(arrayfun(@(x)x*ones(size(var{x,ik}))',1:length(SoftwareNames),'Un',0)),1:length(Replicas_To_Show),'Un',0))));
g2_multi_vector = @(var) string(Replica(Replicas_To_Show(cell2mat(arrayfun(@(ik) ik*cell2mat(arrayfun(@(x) ones(size(var{x,ik}))',1:length(SoftwareNames),'Un',0)),1:length(Replicas_To_Show),'Un',0)))));
y_vector2 = @(ik,var) cell2mat(arrayfun(@(x)var{x,ik},1:length(SoftwareNames),'Un',0));
g1_vector2 = @(ik,var) string(SoftwareNames(cell2mat(arrayfun(@(x)x*ones(size(var{x,ik})),1:length(SoftwareNames),'Un',0))));
g2_vector2 = @(ik,var) string(Replica(Replicas_To_Show(ik)*cell2mat(arrayfun(@(x) ones(size(var{x,ik})),1:length(SoftwareNames),'Un',0))));
y_multi_vector2 = @(var) cell2mat(arrayfun(@(ik) cell2mat(arrayfun(@(x)var{x,ik},1:length(SoftwareNames),'Un',0)),1:length(Replicas_To_Show),'Un',0));
g1_multi_vector2 = @(var) string(SoftwareNames(cell2mat(arrayfun(@(ik) cell2mat(arrayfun(@(x)x*ones(size(var{x,ik})),1:length(SoftwareNames),'Un',0)),1:length(Replicas_To_Show),'Un',0))));
g2_multi_vector2 = @(var) string(Replica(Replicas_To_Show(cell2mat(arrayfun(@(ik) ik*cell2mat(arrayfun(@(x) ones(size(var{x,ik})),1:length(SoftwareNames),'Un',0)),1:length(Replicas_To_Show),'Un',0)))));


threshold_labelC = cell(1,length(Replicas_To_Show));     
numcell_labelC = cell(1,length(Replicas_To_Show));  
numcell_label = cell(1,length(length(Replicas_To_Show)));
numcells_rm_label = cell(1,length(length(Replicas_To_Show)));
numcell_withspots_label = cell(1,length(length(Replicas_To_Show)));      
numcell_nospots_label = cell(1,length(length(Replicas_To_Show)));     
threshold_label = cell(1,length(length(Replicas_To_Show)));     
numspot_label = cell(1,length(length(Replicas_To_Show)));
numspots_rm_label = cell(1,length(length(Replicas_To_Show)));
meancell_spotcount_withspots = cell(1,length(length(Replicas_To_Show)));      
meancell_spotcount = cell(1,length(length(Replicas_To_Show)));    
noprobe_numcell_label = cell(1,length(length(Replicas_To_Show)));
noprobe_numcells_rm_label = cell(1,length(length(Replicas_To_Show)));
noprobe_numcell_withspots_label = cell(1,length(length(Replicas_To_Show)));      
noprobe_numcell_nospots_label = cell(1,length(length(Replicas_To_Show)));     
noprobe_threshold_label = cell(1,length(length(Replicas_To_Show)));     
noprobe_numspot_label = cell(1,length(length(Replicas_To_Show)));
noprobe_numspots_rm_label = cell(1,length(length(Replicas_To_Show)));
noprobe_meancell_spotcount_withspots = cell(1,length(length(Replicas_To_Show)));      
noprobe_meancell_spotcount = cell(1,length(length(Replicas_To_Show)));     
ij_vector = zeros(1,length(Replicas_To_Show));
if (vi<=length(thresholdTypes))
    ob_num = 1;
    ob_type = vi;
else
    ob_num = 2;
    ob_type = vi-length(thresholdTypes);
end
for ik = 1:length(Replicas_To_Show)
    if (Replicas_To_Show(ik)==6)
        ij_vector(ik) = 2;
    else
        ij_vector(ik) = 5;
    end
end
software_and_replica_pairs = allcomb(1:length(Replicas_To_Show),1:length(SoftwareNames));
signal_values = cell(length(SoftwareNames_Spelling),length(Replicas_To_Show));
backgd_values = cell(length(SoftwareNames_Spelling),length(Replicas_To_Show));
signal_minus_backgd_values = cell(length(SoftwareNames_Spelling),length(Replicas_To_Show));
SNR_values = cell(length(SoftwareNames_Spelling),length(Replicas_To_Show));
spot_count_each_type_values = cell(length(SoftwareNames_Spelling),length(Replicas_To_Show));
spot_count_values = cell(length(SoftwareNames_Spelling),length(Replicas_To_Show));
spot_count_above_noprobe_callTable_merged_values = cell(length(SoftwareNames_Spelling),length(Replicas_To_Show));
spot_count_above_noprobe_callTable_unmerged_values = cell(length(SoftwareNames_Spelling),length(Replicas_To_Show));
cell_above_noprobe_spot_detection_threshold_values = cell(length(SoftwareNames_Spelling),length(Replicas_To_Show));
single_cell_spot_count_curves_merged_values = cell(length(Conditions),length(Replicas_To_Show));
single_cell_spot_count_curves_unmerged_values = cell(length(Conditions),length(Replicas_To_Show));
cell_spot_detection_threshold_values = cell(length(Conditions),length(Replicas_To_Show));
spot_count_callTable_merged_values = cell(length(Conditions),length(Replicas_To_Show));
spot_count_callTable_unmerged_values = cell(length(Conditions),length(Replicas_To_Show));
cell_keep_in_values = cell(length(Conditions),length(Replicas_To_Show));
cellbkg_mean_values = cell(length(Conditions),length(Replicas_To_Show));
cellbkg_fano_values = cell(length(Conditions),length(Replicas_To_Show));
noprobe_signal_values = cell(1,length(Replicas_To_Show));
noprobe_backgd_values = cell(1,length(Replicas_To_Show));
noprobe_signal_minus_backgd_values = cell(1,length(Replicas_To_Show));
noprobe_SNR_values = cell(1,length(Replicas_To_Show));
noprobe_spot_count_values = cell(1,length(Replicas_To_Show));
signal_k2statMatrix = cell(1,length(Replicas_To_Show));
backgd_k2statMatrix = cell(1,length(Replicas_To_Show));
signal_minus_backgd_k2statMatrix = cell(1,length(Replicas_To_Show));
SNR_k2statMatrix = cell(1,length(Replicas_To_Show));
spotcount_k2statMatrix = cell(1,length(Replicas_To_Show));
signal_pvalMatrix = cell(1,length(Replicas_To_Show));
backgd_pvalMatrix = cell(1,length(Replicas_To_Show));
signal_minus_backgd_pvalMatrix = cell(1,length(Replicas_To_Show));
SNR_pvalMatrix = cell(1,length(Replicas_To_Show));
spotcount_pvalMatrix = cell(1,length(Replicas_To_Show));
signal_hMatrix = cell(1,length(Replicas_To_Show));
backgd_hMatrix = cell(1,length(Replicas_To_Show));
signal_minus_backgd_hMatrix = cell(1,length(Replicas_To_Show));
SNR_hMatrix = cell(1,length(Replicas_To_Show));
spotcount_hMatrix = cell(1,length(Replicas_To_Show));
cellbkg_mean_k2statMatrix = cell(1,length(Replicas_To_Show));
cellbkg_fano_k2statMatrix = cell(1,length(Replicas_To_Show));
cellbkg_mean_pvalMatrix = cell(1,length(Replicas_To_Show));
cellbkg_fano_pvalMatrix = cell(1,length(Replicas_To_Show));
cellbkg_mean_hMatrix = cell(1,length(Replicas_To_Show));
cellbkg_fano_hMatrix = cell(1,length(Replicas_To_Show));
dup_fraction = zeros(length(SoftwareNames_Spelling),length(Replicas_To_Show));
signal_hMatrix_multi1 = zeros(length(SoftwareNames));
signal_pvalMatrix_multi1 = zeros(length(SoftwareNames));
signal_k2statMatrix_multi1 = zeros(length(SoftwareNames));
backgd_hMatrix_multi1 = zeros(length(SoftwareNames));
backgd_pvalMatrix_multi1 = zeros(length(SoftwareNames));
backgd_k2statMatrix_multi1 = zeros(length(SoftwareNames));
signal_minus_backgd_hMatrix_multi1 = zeros(length(SoftwareNames));
signal_minus_backgd_pvalMatrix_multi1 = zeros(length(SoftwareNames));
signal_minus_backgd_k2statMatrix_multi1 = zeros(length(SoftwareNames));
SNR_hMatrix_multi1 = zeros(length(SoftwareNames));
SNR_pvalMatrix_multi1 = zeros(length(SoftwareNames));
SNR_k2statMatrix_multi1 = zeros(length(SoftwareNames));
spotcount_hMatrix_multi1 = zeros(length(SoftwareNames));
spotcount_pvalMatrix_multi1 = zeros(length(SoftwareNames));
spotcount_k2statMatrix_multi1 = zeros(length(SoftwareNames));
cellbkg_mean_hMatrix_multi1 = zeros(length(Conditions));
cellbkg_mean_pvalMatrix_multi1 = zeros(length(Conditions));
cellbkg_mean_k2statMatrix_multi1 = zeros(length(Conditions));
cellbkg_fano_hMatrix_multi1 = zeros(length(Conditions));
cellbkg_fano_pvalMatrix_multi1 = zeros(length(Conditions));
cellbkg_fano_k2statMatrix_multi1 = zeros(length(Conditions));
signal_hMatrix_multi2 = zeros(length(Replicas_To_Show)*length(SoftwareNames));
signal_pvalMatrix_multi2 = zeros(length(Replicas_To_Show)*length(SoftwareNames));
signal_k2statMatrix_multi2 = zeros(length(Replicas_To_Show)*length(SoftwareNames));
backgd_hMatrix_multi2 = zeros(length(Replicas_To_Show)*length(SoftwareNames));
backgd_pvalMatrix_multi2 = zeros(length(Replicas_To_Show)*length(SoftwareNames));
backgd_k2statMatrix_multi2 = zeros(length(Replicas_To_Show)*length(SoftwareNames));
signal_minus_backgd_hMatrix_multi2 = zeros(length(Replicas_To_Show)*length(SoftwareNames));
signal_minus_backgd_pvalMatrix_multi2 = zeros(length(Replicas_To_Show)*length(SoftwareNames));
signal_minus_backgd_k2statMatrix_multi2 = zeros(length(Replicas_To_Show)*length(SoftwareNames));
SNR_hMatrix_multi2 = zeros(length(Replicas_To_Show)*length(SoftwareNames));
SNR_pvalMatrix_multi2 = zeros(length(Replicas_To_Show)*length(SoftwareNames));
SNR_k2statMatrix_multi2 = zeros(length(Replicas_To_Show)*length(SoftwareNames));
spotcount_hMatrix_multi2 = zeros(length(Replicas_To_Show)*length(SoftwareNames));
spotcount_pvalMatrix_multi2 = zeros(length(Replicas_To_Show)*length(SoftwareNames));
spotcount_k2statMatrix_multi2 = zeros(length(Replicas_To_Show)*length(SoftwareNames));
cellbkg_mean_hMatrix_multi2 = zeros(length(Replicas_To_Show)*length(Conditions));
cellbkg_mean_pvalMatrix_multi2 = zeros(length(Replicas_To_Show)*length(Conditions));
cellbkg_mean_k2statMatrix_multi2 = zeros(length(Replicas_To_Show)*length(Conditions));
cellbkg_fano_hMatrix_multi2 = zeros(length(Replicas_To_Show)*length(Conditions));
cellbkg_fano_pvalMatrix_multi2 = zeros(length(Replicas_To_Show)*length(Conditions));
cellbkg_fano_k2statMatrix_multi2 = zeros(length(Replicas_To_Show)*length(Conditions));

allThrs = cell2mat(arrayfun(@(ii) cell2mat(arrayfun(@(x) softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij_vector(x)}).(Replica{Replicas_To_Show(x)}).(thresholdObj{ob_num}).(thresholdTypes{ob_type}).thresholds,1:length(Replicas_To_Show),'Un',0)),1:length(SoftwareNames),'Un',0));
for ik = 1:length(Replicas_To_Show)
    ij = ij_vector(ik);
    allThrs_Replica = cell2mat(arrayfun(@(ii) cell2mat(arrayfun(@(x) softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij_vector(x)}).(Replica{Replicas_To_Show(x)}).(thresholdObj{ob_num}).(thresholdTypes{ob_type}).thresholds,ik,'Un',0)),1:length(SoftwareNames),'Un',0));
     if (noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).n_reps>0)
        cell_From_Image = noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats.image;
        unique_image_and_cell_pair = noProbePooledResults_thresholdInfo.(Replica{Replicas_To_Show(ik)}).unique_image_and_cell_pair;
        cells_planeVol = cellfun(@(x) sum(double(x.histo_y)),noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats.cell_bkg);
        cells_with_outlier_Vol = find(isoutlier(cells_planeVol));
        cells_without_outlier_Vol = find(~isoutlier(cells_planeVol));
        if (illumination_corrected)
            ImageCellStat_image_and_cell_pair = table2array(noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats_IntCorrected(:,{'image','cellNo'}));
            max_bkg = cellfun(@(x) double(x.max),noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats_IntCorrected.cell_bkg);
        else
            ImageCellStat_image_and_cell_pair = table2array(noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats(:,{'image','cellNo'}));
            max_bkg = cellfun(@(x) double(x.max),noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats.cell_bkg);
        end
        rm_cells = find(max_bkg>=max_bkg_cutoff);
        BadImage_cells = find(ismember(unique_image_and_cell_pair,ImageCellStat_image_and_cell_pair(rm_cells,:),'rows'));
        cells_to_keep_in = setdiff(1:length(cell_From_Image),BadImage_cells);
        noprobe_numcell_label{ik}= num2str(length(cells_to_keep_in));
        noprobe_numcells_rm_label{ik}= num2str(length(union(cells_with_outlier_Vol',BadImage_cells)));
        cell_keep_in_values{6,ik} = cells_to_keep_in;
        cellThr_noprobe = noProbePooledResults_thresholdInfo.(Replica{Replicas_To_Show(ik)}).(thresholdObj{ob_num}).(thresholdTypes{ob_type}).thresholds;
        if (vi<=length(thresholdTypes))
            thrList = reshape(noProbePooledResults_thresholdInfo.(Replica{Replicas_To_Show(ik)}).image_sugg_m_table,[],1);
            noprobe_threshold_label{ik} = strcat(num2str(mean(thrList,'omitnan'),'%1.0f'),char(177),num2str(std(thrList,'omitnan'),'%1.0f'));
        else
            thrList = reshape(noProbePooledResults_thresholdInfo.(Replica{Replicas_To_Show(ik)}).cell_sugg_m_table,[],1);
            noprobe_threshold_label{ik}= strcat(num2str(mean(cellThr_noprobe(cells_to_keep_in),'omitnan'),'%1.0f'),char(177),num2str(std(cellThr_noprobe(cells_to_keep_in),'omitnan'),'%1.0f'));
        end
        cellThr_noprobe_kept_cells = cellThr_noprobe(cells_to_keep_in);
        cell_spot_detection_threshold_values{length(Conditions),ik} = cellThr_noprobe_kept_cells;
        call_table_spot_count_no_merging = noProbePooledResults_thresholdInfo.(Replica{Replicas_To_Show(ik)}).numCallTable_at_Thr_in_cell_image(unique_image_and_cell_pair(cells_to_keep_in,1),unique_image_and_cell_pair(cells_to_keep_in,2),cellThr_noprobe(cells_to_keep_in));
        spot_count_callTable_unmerged_values{length(Conditions),ik} = call_table_spot_count_no_merging;
        single_cell_spot_count_curves_no_merging = noProbePooledResults_thresholdInfo.(Replica{Replicas_To_Show(ik)}).extended_callTable_spot_count_curve_in_cell_X_in_image_Y_y(cells_to_keep_in,:);
        single_cell_spot_count_curves_unmerged_values{length(Conditions),ik} = single_cell_spot_count_curves_no_merging;
        call_table_spot_count_merged = noProbePooledResults_thresholdInfo.(Replica{Replicas_To_Show(ik)}).spotmerged_numCallTable_at_Thr_in_cell_image(unique_image_and_cell_pair(cells_to_keep_in,1),unique_image_and_cell_pair(cells_to_keep_in,2),cellThr_noprobe(cells_to_keep_in));
        spot_count_callTable_merged_values{length(Conditions),ik} = call_table_spot_count_merged;
        single_cell_spot_count_curves_merged = noProbePooledResults_thresholdInfo.(Replica{Replicas_To_Show(ik)}).callTable_merged_spot_count_curve_in_cell_X_in_image_Y_y(cells_to_keep_in,:);
        single_cell_spot_count_curves_merged_values{length(Conditions),ik} = single_cell_spot_count_curves_merged;
        if (illumination_corrected)
            cells_cellbkg_mean = cellfun(@(x) x.mean,noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats_IntCorrected.cell_bkg(cells_to_keep_in));
            cell_bkg_std = cellfun(@(x) x.stdev,noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats_IntCorrected.cell_bkg(cells_to_keep_in));
        else
            cells_cellbkg_mean = cellfun(@(x) x.mean,noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats.cell_bkg(cells_to_keep_in));
            cell_bkg_std = cellfun(@(x) x.stdev,noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats.cell_bkg(cells_to_keep_in));
        end
        cells_cellbkg_var = cell_bkg_std.^2;
        cells_cellbkg_fano = cells_cellbkg_var./cells_cellbkg_mean;
        cellbkg_mean_values{length(Conditions),ik} = cells_cellbkg_mean;
        cellbkg_fano_values{length(Conditions),ik} = cells_cellbkg_fano;
        if (~isempty(noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats))
        spot_image_cell_num_pairs = [double(noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats.image) ...
            double(noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats.cell)];
        [~,spot_image_cell_num_list] = ismember(spot_image_cell_num_pairs,unique_image_and_cell_pair,'rows');
        %spot background all pixels in bounding box.
        if (illumination_corrected)
            spot_int = cellfun(@(x) double(x.(int_level_types{int_level})),noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats_IntCorrected.noBkgSubStats);
            if (local_over_cell_bkg)
                local_spot_bkg_int_mean = cellfun(@(x) double(x.mean),noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats_IntCorrected.localBkgStats);
                local_spot_bkg_int_std = cellfun(@(x) double(x.stdev),noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats_IntCorrected.localBkgStats);
                spot_bkg_int_mean = local_spot_bkg_int_mean;
                spot_bkg_int_std = local_spot_bkg_int_std;
            else
                cell_bkg_mean = cellfun(@(x) x.mean,noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats_IntCorrected.cell_bkg);
                cell_bkg_std = cellfun(@(x) x.stdev,noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats_IntCorrected.cell_bkg);
                cell_spot_bkg_int_mean = cell_bkg_mean(spot_image_cell_num_list);
                cell_spot_bkg_int_std = cell_bkg_std(spot_image_cell_num_list);
                spot_bkg_int_mean = cell_spot_bkg_int_mean;
                spot_bkg_int_std = cell_spot_bkg_int_std;
            end
        else
            spot_int = cellfun(@(x) double(x.(int_level_types{int_level})),noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats.noBkgSubStats);
            spot_int_tot = cellfun(@(x) sum(double(x.histo_x).*double(x.histo_y)),noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats.noBkgSubStats);
            if (local_over_cell_bkg)
                local_spot_bkg_int_mean = cellfun(@(x) double(x.mean),noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats.localBkgStats);
                local_spot_bkg_int_std = cellfun(@(x) double(x.stdev),noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats.localBkgStats);
                spot_bkg_int_mean = local_spot_bkg_int_mean;
                spot_bkg_int_std = local_spot_bkg_int_std;
            else
                cell_bkg_mean = cellfun(@(x) x.mean,noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats.cell_bkg);
                cell_bkg_std = cellfun(@(x) x.stdev,noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats.cell_bkg);
                cell_spot_bkg_int_mean = cell_bkg_mean(spot_image_cell_num_list);
                cell_spot_bkg_int_std = cell_bkg_std(spot_image_cell_num_list);
                spot_bkg_int_mean = cell_spot_bkg_int_mean;
                spot_bkg_int_std = cell_spot_bkg_int_std;
            end
        end
        spot_bkg_corr_int_mean = spot_int-spot_bkg_int_mean;
        spot_int_snr = (spot_int-spot_bkg_int_mean)./spot_bkg_int_std;
        obj.spotTable = noProbePooledResults_quant.(Replica{Replicas_To_Show(ik)}).struct_quant_results_spotTable;
        obj = A_BH_ApplySpotMerging(obj);
        ImageSpotStats_Table = noProbePooledResults_stats.(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats;
        ImageSpotStats_Table(:,'is_dup') = obj.spotTable(:,'is_a_duplicate');clear obj
        unique_location_cell_rows_base_func = @(cell_rows) double(table2array(unique(ImageSpotStats_Table(cell_rows,{'x','y','z'}),'rows')));
        unique_location_cell_rows_func = @(cell_rows,C) cell_rows(arrayfun(@(ind) find(ismember(double(table2array(ImageSpotStats_Table(cell_rows,{'x','y','z'}))),C(ind,:),'rows'),1),1:size(C,1)));
        %Spot Stats
        cell_rows_base_func = @(nn) find(ismember(table2array(ImageSpotStats_Table(:,{'image','cell'})),unique_image_and_cell_pair(cells_to_keep_in(nn),:),'rows'));
        cell_rows_func = @(nn)  unique_location_cell_rows_func(cell_rows_base_func(nn),unique_location_cell_rows_base_func(cell_rows_base_func(nn)));
        inner_cell_rows_func = @(cell_rows,nn) cell_rows(double(table2array(ImageSpotStats_Table(cell_rows,'dropout_thresh'))>=cellThr_noprobe(cells_to_keep_in(nn))).* ...
            double(ismember(table2array(ImageSpotStats_Table(cell_rows,'nascent_flag')),nascent_status{ui})).*...
            double(ismember(table2array(ImageSpotStats_Table(cell_rows,'nuc_flag')),nuc_status{ui})).*...
            double(ismember(table2array(ImageSpotStats_Table(cell_rows,'is_dup')),0)).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'xgw'))>xgwlow).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'xgw'))<xgwhigh).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'ygw'))>ygwlow).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'ygw'))<ygwhigh).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'x'))>xlow).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'x'))<xhigh).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'y'))>ylow).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'y'))<yhigh).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'z'))>zlow).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'z'))<zhigh)==1).';
        %cellThr_noprobe_kept_cells
        keep_spots = cell2mat(arrayfun(@(nn) inner_cell_rows_func(cell_rows_func(nn),nn),1:length(cells_to_keep_in),'Un',0));
        % Cell Stats
        cells_keep_rows_base_func = @(nn) find(ismember(table2array(ImageSpotStats_Table(:,{'image','cell'})),unique_image_and_cell_pair(cells_to_keep_in(nn),:),'rows'));
        cells_keep_rows_func = @(nn)  unique_location_cell_rows_func(cells_keep_rows_base_func(nn),unique_location_cell_rows_base_func(cells_keep_rows_base_func(nn)));
        inner_cells_keep_rows_func = @(cell_rows,nn) sum(double(table2array(ImageSpotStats_Table(cell_rows,'dropout_thresh'))>=cellThr_noprobe(cells_to_keep_in(nn))).* ...
            double(ismember(table2array(ImageSpotStats_Table(cell_rows,'nascent_flag')),nascent_status{ui})).*...
            double(ismember(table2array(ImageSpotStats_Table(cell_rows,'nuc_flag')),nuc_status{ui})).*...
            double(ismember(table2array(ImageSpotStats_Table(cell_rows,'is_dup')),0)).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'xgw'))>xgwlow).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'xgw'))<xgwhigh).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'ygw'))>ygwlow).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'ygw'))<ygwhigh).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'x'))>xlow).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'x'))<xhigh).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'y'))>ylow).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'y'))<yhigh).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'z'))>zlow).*...
            double(table2array(ImageSpotStats_Table(cell_rows,'z'))<zhigh));
        cell_spot_count = cell2mat(arrayfun(@(nn) inner_cells_keep_rows_func(cells_keep_rows_func(nn),nn),1:length(cells_to_keep_in),'Un',0));
        noprobe_numspot_label{ik} = num2str(length(keep_spots));
        noprobe_numspots_rm_label{ik} = num2str(size(ImageSpotStats_Table,1)-length(keep_spots));
        noprobe_numcell_withspots_label{ik}= num2str(sum(cell_spot_count>0));
        noprobe_numcell_nospots_label{ik}= num2str(length(sum(cell_spot_count==0)));
        noprobe_meancell_spotcount_withspots{ik} = mean(cell_spot_count(cell_spot_count>0));
        noprobe_meancell_spotcount{ik} = mean(cell_spot_count);
        clear ImageSpotStats_Table
        if (~isempty(keep_spots))
            curr_vector = spot_int(keep_spots);
            curr_vector(isnan(curr_vector))=[];
            noprobe_signal_values{ik} = curr_vector;
            curr_vector = spot_bkg_int_mean(keep_spots);
            curr_vector(isnan(curr_vector))=[];
            noprobe_backgd_values{ik} = curr_vector;
            curr_vector = spot_bkg_corr_int_mean(keep_spots);
            curr_vector(isnan(curr_vector))=[];
            noprobe_signal_minus_backgd_values{ik} = curr_vector;
            curr_vector = spot_int_snr(keep_spots);
            curr_vector(isnan(curr_vector))=[];
            noprobe_SNR_values{ik} = curr_vector;
            curr_vector = cell_spot_count;
            noprobe_spot_count_values{ik} = curr_vector';
        end
        end
    end
    for ii = 1:length(SoftwareNames_Spelling)
        if (softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).n_reps>0)
            cell_From_Image = softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats.image;
            unique_image_and_cell_pair = softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).unique_image_and_cell_pair;
            spot_image_cell_num_pairs = [double(softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats.image) ...
                double(softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats.cell)];
            [~,spot_image_cell_num_list] = ismember(spot_image_cell_num_pairs,unique_image_and_cell_pair,'rows');
            cells_planeVol = cellfun(@(x) sum(double(x.histo_y)),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats.cell_bkg);
            cells_with_outlier_Vol = find(isoutlier(cells_planeVol));
            cells_without_outlier_Vol = find(~isoutlier(cells_planeVol));
            %spot background all pixels in bounding box.
            if (illumination_corrected)
                spot_int = cellfun(@(x) double(x.(int_level_types{int_level})),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats_IntCorrected.noBkgSubStats);
                if (local_over_cell_bkg)
                    local_spot_bkg_int_mean = cellfun(@(x) double(x.mean),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats_IntCorrected.localBkgStats);
                    local_spot_bkg_int_std = cellfun(@(x) double(x.stdev),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats_IntCorrected.localBkgStats);
                    spot_bkg_int_mean = local_spot_bkg_int_mean;
                    spot_bkg_int_std = local_spot_bkg_int_std;
                else
                    cell_bkg_mean = cellfun(@(x) x.mean,softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats_IntCorrected.cell_bkg);
                    cell_bkg_std = cellfun(@(x) x.stdev,softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats_IntCorrected.cell_bkg);
                    cell_spot_bkg_int_mean = cell_bkg_mean(spot_image_cell_num_list);
                    cell_spot_bkg_int_std = cell_bkg_std(spot_image_cell_num_list);
                    spot_bkg_int_mean = cell_spot_bkg_int_mean;
                    spot_bkg_int_std = cell_spot_bkg_int_std;
                end
            else
                spot_int = cellfun(@(x) double(x.(int_level_types{int_level})),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats.noBkgSubStats);
                spot_int_tot = cellfun(@(x) sum(double(x.histo_x).*double(x.histo_y)),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats.noBkgSubStats);
                if (local_over_cell_bkg)
                    local_spot_bkg_int_mean = cellfun(@(x) double(x.mean),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats.localBkgStats);
                    local_spot_bkg_int_std = cellfun(@(x) double(x.stdev),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats.localBkgStats);
                    spot_bkg_int_mean = local_spot_bkg_int_mean;
                    spot_bkg_int_std = local_spot_bkg_int_std;
                else
                    cell_bkg_mean = cellfun(@(x) x.mean,softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats.cell_bkg);
                    cell_bkg_std = cellfun(@(x) x.stdev,softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats.cell_bkg);
                    cell_spot_bkg_int_mean = cell_bkg_mean(spot_image_cell_num_list);
                    cell_spot_bkg_int_std = cell_bkg_std(spot_image_cell_num_list);
                    spot_bkg_int_mean = cell_spot_bkg_int_mean;
                    spot_bkg_int_std = cell_spot_bkg_int_std;
                end
            end
            spot_bkg_corr_int_mean = spot_int-spot_bkg_int_mean;
            spot_int_snr = (spot_int-spot_bkg_int_mean)./spot_bkg_int_std;
            obj.spotTable = softwarePooledResults_quant.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_quant_results_spotTable;
            obj = A_BH_ApplySpotMerging(obj);
           ImageSpotStats_Table = softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageSpotStats;
           ImageSpotStats_Table(:,'is_dup') = obj.spotTable(:,'is_a_duplicate');clear obj
            dup_fraction(ii,ik) = sum(ImageSpotStats_Table.is_dup)/size(ImageSpotStats_Table,1);
            unique_location_cell_rows_base_func = @(cell_rows) double(table2array(unique(ImageSpotStats_Table(cell_rows,{'x','y','z'}),'rows')));
            unique_location_cell_rows_func = @(cell_rows,C) cell_rows(arrayfun(@(ind) find(ismember(double(table2array(ImageSpotStats_Table(cell_rows,{'x','y','z'}))),C(ind,:),'rows'),1),1:size(C,1)));
            %Spot Stats
            if (illumination_corrected)
                ImageCellStat_image_and_cell_pair = table2array(softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats_IntCorrected(:,{'image','cellNo'}));
                max_bkg = cellfun(@(x) double(x.max),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats_IntCorrected.cell_bkg);
            else
                ImageCellStat_image_and_cell_pair = table2array(softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats(:,{'image','cellNo'}));
                max_bkg = cellfun(@(x) double(x.max),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats.cell_bkg);
            end
            rm_cells = find(max_bkg>=max_bkg_cutoff);
            BadImage_cells = find(ismember(unique_image_and_cell_pair,ImageCellStat_image_and_cell_pair(rm_cells,:),'rows'));
            cells_to_keep_in = setdiff(1:length(cell_From_Image),BadImage_cells);
            if (illumination_corrected)
                    cells_cellbkg_mean = cellfun(@(x) x.mean,softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats_IntCorrected.cell_bkg(cells_to_keep_in));
                    cell_bkg_std = cellfun(@(x) x.stdev,softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats_IntCorrected.cell_bkg(cells_to_keep_in));
             else
                    cells_cellbkg_mean = cellfun(@(x) x.mean,softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats.cell_bkg(cells_to_keep_in));
                    cell_bkg_std = cellfun(@(x) x.stdev,softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats.cell_bkg(cells_to_keep_in));
             end
             cells_cellbkg_var = cell_bkg_std.^2;
             cells_cellbkg_fano = cells_cellbkg_var./cells_cellbkg_mean;
             cellbkg_mean_values{ii,ik} = cells_cellbkg_mean;
             cellbkg_fano_values{ii,ik} = cells_cellbkg_fano;
            switch useFixedThreshold
                case 0
                    cellThr = softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).(thresholdObj{ob_num}).(thresholdTypes{ob_type}).thresholds;
                    if (vi<=length(thresholdTypes))
                        thrList = reshape(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).image_sugg_m_table,[],1);
                        threshold_label{ik}{ii} = strcat(num2str(mean(thrList,'omitnan'),'%1.0f'),char(177),num2str(std(thrList,'omitnan'),'%1.0f'));
                    else
                        thrList = reshape(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).cell_sugg_m_table,[],1);
                        threshold_label{ik}{ii} = strcat(num2str(mean(cellThr(cells_to_keep_in),'omitnan'),'%1.0f'),char(177),num2str(std(cellThr(cells_to_keep_in),'omitnan'),'%1.0f'));
                    end
                case 1
                    cellThr = mean(allThrs_Replica,'omitnan')*ones(size(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).(thresholdObj{ob_num}).(thresholdTypes{ob_type}).thresholds));
                    threshold_label{ik}{ii} = num2str(mean(allThrs_Replica,'omitnan'),'%1.0f');
                case 2
                    cellThr = mean(allThrs,'omitnan')*ones(size(softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).(thresholdObj{ob_num}).(thresholdTypes{ob_type}).thresholds));
                    threshold_label{ik}{ii} = num2str(mean(allThrs,'omitnan'),'%1.0f');
            end
            cell_rows_base_func = @(nn) find(ismember(table2array(ImageSpotStats_Table(:,{'image','cell'})),unique_image_and_cell_pair(cells_to_keep_in(nn),:),'rows'));
            cell_rows_func = @(nn)  unique_location_cell_rows_func(cell_rows_base_func(nn),unique_location_cell_rows_base_func(cell_rows_base_func(nn)));
          % double(table2array(ImageSpotStats_Table(cell_rows,'dropout_thresh'))>=mean(cellThr_noprobe_kept_cells)).*
            inner_cell_rows_func = @(cell_rows,nn) cell_rows(double(table2array(ImageSpotStats_Table(cell_rows,'dropout_thresh'))>=cellThr(cells_to_keep_in(nn))).* ...
                double(ismember(table2array(ImageSpotStats_Table(cell_rows,'nascent_flag')),nascent_status{ui})).*...
                double(ismember(table2array(ImageSpotStats_Table(cell_rows,'nuc_flag')),nuc_status{ui})).*...
                double(ismember(table2array(ImageSpotStats_Table(cell_rows,'is_dup')),0)).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'xgw'))>xgwlow).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'xgw'))<xgwhigh).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'ygw'))>ygwlow).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'ygw'))<ygwhigh).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'x'))>xlow).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'x'))<xhigh).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'y'))>ylow).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'y'))<yhigh).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'z'))>zlow).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'z'))<zhigh)==1).';
            keep_spots = cell2mat(arrayfun(@(nn) inner_cell_rows_func(cell_rows_func(nn),nn),1:length(cells_to_keep_in),'Un',0));
            % Cell Stats
            if (illumination_corrected)
                ImageCellStat_image_and_cell_pair = table2array(softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats_IntCorrected(:,{'image','cellNo'}));
                max_bkg = cellfun(@(x) double(x.max),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats_IntCorrected.cell_bkg);
            else
                ImageCellStat_image_and_cell_pair = table2array(softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats(:,{'image','cellNo'}));
                max_bkg = cellfun(@(x) double(x.max),softwarePooledResults_stats.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).struct_ImageCellStats.cell_bkg);
            end
            rm_cells = find(max_bkg>=max_bkg_cutoff);
            BadImage_cells = find(ismember(unique_image_and_cell_pair,ImageCellStat_image_and_cell_pair(rm_cells,:),'rows'));
            cells_to_keep_in = setdiff(cells_without_outlier_Vol',BadImage_cells);
            cell_spot_detection_threshold_values{ii,ik} = cellThr(cells_to_keep_in);
            cellThr1 = cellThr;
            cellThr1(cellThr<=ceil(mean(cellThr_noprobe_kept_cells)))=ceil(mean(cellThr_noprobe_kept_cells))+1;
            cellThr_above_noprobe = cellThr1(cells_to_keep_in);clear cellThr1
            cell_above_noprobe_spot_detection_threshold_values{ii,ik} = cellThr_above_noprobe;
            call_table_spot_count_no_merging = softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).numCallTable_at_Thr_in_cell_image(unique_image_and_cell_pair(cells_to_keep_in,1),unique_image_and_cell_pair(cells_to_keep_in,2),cellThr(cells_to_keep_in));
            spot_count_callTable_unmerged_values{ii,ik} = call_table_spot_count_no_merging;
            single_cell_spot_count_curves_no_merging = softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).callTable_spot_count_curve_in_cell_X_in_image_Y_y(cells_to_keep_in,:);
            single_cell_spot_count_curves_unmerged_values{ii,ik} = single_cell_spot_count_curves_no_merging;
            call_table_above_noprobe_spot_count_no_merging = softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).numCallTable_at_Thr_in_cell_image(unique_image_and_cell_pair(cells_to_keep_in,1),unique_image_and_cell_pair(cells_to_keep_in,2),cellThr_above_noprobe);
            spot_count_above_noprobe_callTable_unmerged_values{ii,ik} = call_table_above_noprobe_spot_count_no_merging;
            call_table_spot_count_merged = softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).spotmerged_numCallTable_at_Thr_in_cell_image(unique_image_and_cell_pair(cells_to_keep_in,1),unique_image_and_cell_pair(cells_to_keep_in,2),cellThr(cells_to_keep_in));
            spot_count_callTable_merged_values{ii,ik} = call_table_spot_count_merged;
            call_table_above_noprobe_spot_count_merged = softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).spotmerged_numCallTable_at_Thr_in_cell_image(unique_image_and_cell_pair(cells_to_keep_in,1),unique_image_and_cell_pair(cells_to_keep_in,2),cellThr_above_noprobe);
            spot_count_above_noprobe_callTable_merged_values{ii,ik} = call_table_above_noprobe_spot_count_merged;
            single_cell_spot_count_curves_merged = softwarePooledResults_thresholdInfo.(SoftwareNames{ii}).(group_concentrations{ij}).(Replica{Replicas_To_Show(ik)}).callTable_merged_spot_count_curve_in_cell_X_in_image_Y_y(cells_to_keep_in,:);
            single_cell_spot_count_curves_merged_values{ii,ik} = single_cell_spot_count_curves_merged;
            cells_keep_rows_base_func = @(nn) find(ismember(table2array(ImageSpotStats_Table(:,{'image','cell'})),unique_image_and_cell_pair(cells_to_keep_in(nn),:),'rows'));
            cells_keep_rows_func = @(nn)  unique_location_cell_rows_func(cells_keep_rows_base_func(nn),unique_location_cell_rows_base_func(cells_keep_rows_base_func(nn)));
            inner_cells_keep_rows_func = @(cell_rows,nn) sum(double(table2array(ImageSpotStats_Table(cell_rows,'dropout_thresh'))>=cellThr(cells_to_keep_in(nn))).* ...
                double(ismember(table2array(ImageSpotStats_Table(cell_rows,'nascent_flag')),nascent_status{ui})).*...
                double(ismember(table2array(ImageSpotStats_Table(cell_rows,'nuc_flag')),nuc_status{ui})).*...
                double(ismember(table2array(ImageSpotStats_Table(cell_rows,'is_dup')),0)).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'xgw'))>xgwlow).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'xgw'))<xgwhigh).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'ygw'))>ygwlow).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'ygw'))<ygwhigh).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'x'))>xlow).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'x'))<xhigh).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'y'))>ylow).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'y'))<yhigh).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'z'))>zlow).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'z'))<zhigh));
            cell_spot_count = cell2mat(arrayfun(@(nn) inner_cells_keep_rows_func(cells_keep_rows_func(nn),nn),1:length(cells_to_keep_in),'Un',0));
            inner_cells_keep_rows_type_func = @(cell_rows,nn,nuc_state,nascent_state) sum(double(table2array(ImageSpotStats_Table(cell_rows,'dropout_thresh'))>=cellThr(cells_to_keep_in(nn))).* ...
                double(ismember(table2array(ImageSpotStats_Table(cell_rows,'nascent_flag')),nascent_state)).*...
                double(ismember(table2array(ImageSpotStats_Table(cell_rows,'nuc_flag')),nuc_state)).*...
                double(ismember(table2array(ImageSpotStats_Table(cell_rows,'is_dup')),0)).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'xgw'))>xgwlow).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'xgw'))<xgwhigh).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'ygw'))>ygwlow).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'ygw'))<ygwhigh).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'x'))>xlow).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'x'))<xhigh).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'y'))>ylow).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'y'))<yhigh).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'z'))>zlow).*...
                double(table2array(ImageSpotStats_Table(cell_rows,'z'))<zhigh));
              cell_state_spot_count_func = @(nuc_state,nas_state) cell2mat(arrayfun(@(nn) inner_cells_keep_rows_type_func(cells_keep_rows_func(nn),nn,nuc_state,nas_state),1:length(cells_to_keep_in),'Un',0));
            state_vector = {[0 0],[1 0],[0 1],[1 1]};%cytomatu, nucmat, cytonas, nucnasc
            spot_count_each_type_values{ii,ik} = CATnWrapper(cellfun(@(x) cell_state_spot_count_func(x(1),x(2))',state_vector,'Un',0),2);
            cell_keep_in_values{ii,ik} = cells_to_keep_in;
            numspot_label{ik}{ii} = num2str(length(keep_spots));
            numspots_rm_label{ik}{ii} = num2str(size(ImageSpotStats_Table,1)-length(keep_spots));
            numcell_label{ik}{ii} = num2str(length(cells_to_keep_in));
            numcells_rm_label{ik}{ii} = num2str(length(union(cells_with_outlier_Vol',BadImage_cells)));
            numcell_withspots_label{ik}{ii} = num2str(sum(cell_spot_count>0));
            numcell_nospots_label{ik}{ii} = num2str(length(sum(cell_spot_count==0)));
            meancell_spotcount_withspots{ik}(ii) = mean(cell_spot_count(cell_spot_count>0));
            meancell_spotcount{ik}(ii) = mean(cell_spot_count);
            clear ImageSpotStats_Table
            if (length(keep_spots)>1)
                curr_vector = spot_int(keep_spots);
                curr_vector(isnan(curr_vector))=[];
                signal_values{ii,ik} = curr_vector;
                curr_vector = spot_bkg_int_mean(keep_spots);
                curr_vector(isnan(curr_vector))=[];
                backgd_values{ii,ik} = curr_vector;
                curr_vector = spot_bkg_corr_int_mean(keep_spots);
                curr_vector(isnan(curr_vector))=[];
                signal_minus_backgd_values{ii,ik} = curr_vector;
                curr_vector = spot_int_snr(keep_spots);
                curr_vector(isnan(curr_vector))=[];
                SNR_values{ii,ik} = curr_vector;
                curr_vector = cell_spot_count;
                spot_count_values{ii,ik} = curr_vector';
            end
        end
    end   
end
thresholdType = {'Single-Image','Single-Cell'};
thresholdLevels = {'Min','Quartile 1','Mean','Median','Quartile 3','Optimal'};
backgroundTypes = {'Cell','Spot'};
SpotType = {'All','Nascent Nuclear','Mature Cytoplasmic','Mature Nuclear','Nuclear','Cytoplasmic','Mature'};
ThresholdGroup = {'Software/Replica-Specific','Replica-Specific/Software-Agnostic','Replica/Software-Agnostic'};
figData = struct();
figData.configuration.SpotType = SpotType{ui};
figData.configuration.illumination_corrected = logical(illumination_corrected);
figData.configuration.intensityLevel = int_level_types{int_level};
figData.configuration.backgroundType = backgroundTypes{local_over_cell_bkg+1};
figData.configuration.ThresholdType = thresholdType{ob_num};
figData.configuration.ThresholdLevel = thresholdLevels{ob_type};
figData.configuration.ThresholdGrouping = ThresholdGroup{useFixedThreshold+1};
figData.signal_values = signal_values;
figData.backgd_values = backgd_values;
figData.signal_minus_backgd_values = signal_minus_backgd_values;
figData.SNR_values = SNR_values;
figData.spot_count_values = spot_count_values;
figData.spot_count_each_type_values = spot_count_each_type_values;
figData.cell_keep_in_values = cell_keep_in_values;
figData.spot_count_callTable_merged_values = spot_count_callTable_merged_values;
figData.spot_count_callTable_unmerged_values = spot_count_callTable_unmerged_values;
figData.spot_count_above_noprobe_callTable_merged_values = spot_count_above_noprobe_callTable_merged_values;
figData.spot_count_above_noprobe_callTable_unmerged_values = spot_count_above_noprobe_callTable_unmerged_values;
figData.single_cell_spot_count_curves_unmerged_values = single_cell_spot_count_curves_unmerged_values;
figData.single_cell_spot_count_curves_merged_values = single_cell_spot_count_curves_merged_values;
figData.cell_spot_detection_threshold_values = cell_spot_detection_threshold_values;
figData.cell_above_noprobe_spot_detection_threshold_values = cell_above_noprobe_spot_detection_threshold_values;
figData.cell_backgd_mean_values = cellbkg_mean_values;
figData.cell_backgd_fano_values = cellbkg_fano_values;
figData.noprobe_signal_values = noprobe_signal_values;
figData.noprobe_backgd_values = noprobe_backgd_values;
figData.noprobe_signal_minus_backgd_values = noprobe_signal_minus_backgd_values;
figData.noprobe_SNR_values = noprobe_SNR_values;
figData.noprobe_spot_count_values = noprobe_spot_count_values;
figData.labels.numspot_label = numspot_label;
figData.labels.numspots_rm_label = numspots_rm_label;
figData.labels.numcell_label = numcell_label;
figData.labels.numcells_rm_label = numcells_rm_label;
figData.labels.numcell_withspots_label = numcell_withspots_label;
figData.labels.numcell_nospots_label = numcell_nospots_label;
figData.labels.threshold_label = threshold_label;
figData.labels.meancell_spotcount_withspots = meancell_spotcount_withspots;
figData.labels.meancell_spotcount_nospots = meancell_spotcount;
figData.labels.noprobe_numspot_label = noprobe_numspot_label;
figData.labels.noprobe_numspots_rm_label = noprobe_numspots_rm_label;
figData.labels.noprobe_numcell_label = noprobe_numcell_label;
figData.labels.noprobe_numcells_rm_label = noprobe_numcells_rm_label;
figData.labels.noprobe_numcell_withspots_label = noprobe_numcell_withspots_label;
figData.labels.noprobe_numcell_nospots_label = noprobe_numcell_nospots_label;
figData.labels.noprobe_threshold_label = noprobe_threshold_label;
figData.labels.noprobe_meancell_spotcount_withspots = noprobe_meancell_spotcount_withspots;
figData.labels.noprobe_meancell_spotcount_nospots = noprobe_meancell_spotcount;


if (~skip_statistical_tests)
compareType = "bonferroni";
anova1_results = struct();
anova1_results.signal_pval = cell(1,length(Replicas_To_Show));
anova1_results.signal_stats = cell(1,length(Replicas_To_Show));
anova1_results.signal_AnalysisOfVariance = cell(1,length(Replicas_To_Show));
anova1_results.signal_CI = cell(1,length(Replicas_To_Show));
anova1_results.signal_multiple_comparisons = cell(1,length(Replicas_To_Show));
anova1_results.backgd_pval = cell(1,length(Replicas_To_Show));
anova1_results.backgd_stats = cell(1,length(Replicas_To_Show));
anova1_results.backgd_AnalysisOfVariance = cell(1,length(Replicas_To_Show));
anova1_results.backgd_CI = cell(1,length(Replicas_To_Show));
anova1_results.backgd_multiple_comparisons = cell(1,length(Replicas_To_Show));
anova1_results.signal_minus_backgd_pval = cell(1,length(Replicas_To_Show));
anova1_results.signal_minus_backgd_stats = cell(1,length(Replicas_To_Show));
anova1_results.signal_minus_backgd_AnalysisOfVariance = cell(1,length(Replicas_To_Show));
anova1_results.signal_minus_backgd_CI = cell(1,length(Replicas_To_Show));
anova1_results.signal_minus_backgd_multiple_comparisons = cell(1,length(Replicas_To_Show));
anova1_results.SNR_pval = cell(1,length(Replicas_To_Show));
anova1_results.SNR_stats = cell(1,length(Replicas_To_Show));
anova1_results.SNR_AnalysisOfVariance = cell(1,length(Replicas_To_Show));
anova1_results.SNR_CI = cell(1,length(Replicas_To_Show));
anova1_results.SNR_multiple_comparisons = cell(1,length(Replicas_To_Show));
anova1_results.spotcount_pval = cell(1,length(Replicas_To_Show));
anova1_results.spotcount_stats = cell(1,length(Replicas_To_Show));
anova1_results.spotcount_AnalysisOfVariance = cell(1,length(Replicas_To_Show));
anova1_results.spotcount_CI = cell(1,length(Replicas_To_Show));
anova1_results.spotcount_multiple_comparisons = cell(1,length(Replicas_To_Show));
for ik = 1:length(Replicas_To_Show)
    signal_k2statMatrix{ik} = NaN(length(SoftwareNames_Spelling));
    backgd_k2statMatrix{ik} = NaN(length(SoftwareNames_Spelling));
    signal_minus_backgd_k2statMatrix{ik} = NaN(length(SoftwareNames_Spelling));
    SNR_k2statMatrix{ik} = NaN(length(SoftwareNames_Spelling));
    spotcount_k2statMatrix{ik} = NaN(length(SoftwareNames_Spelling));
    signal_pvalMatrix{ik} = NaN(length(SoftwareNames_Spelling));
    backgd_pvalMatrix{ik} = NaN(length(SoftwareNames_Spelling));
    signal_minus_backgd_pvalMatrix{ik} = NaN(length(SoftwareNames_Spelling));
    SNR_pvalMatrix{ik} = NaN(length(SoftwareNames_Spelling));
    spotcount_pvalMatrix{ik} = NaN(length(SoftwareNames_Spelling));
    signal_hMatrix{ik} = NaN(length(SoftwareNames_Spelling));
    backgd_hMatrix{ik} = NaN(length(SoftwareNames_Spelling));
    signal_minus_backgd_hMatrix{ik} = NaN(length(SoftwareNames_Spelling));
    SNR_hMatrix{ik} = NaN(length(SoftwareNames_Spelling));
    spotcount_hMatrix{ik} = NaN(length(SoftwareNames_Spelling));
    for ii = 1:length(SoftwareNames_Spelling)
        for il = 1:length(SoftwareNames_Spelling)
            [h,p,ks2stat] = kstest2(signal_values{ii,ik},signal_values{il,ik},'Alpha',pAlpha_threshold);
            signal_hMatrix{ik}(ii,il) = h;signal_pvalMatrix{ik}(ii,il) = p;signal_k2statMatrix{ik}(ii,il) = ks2stat;
            [h,p,ks2stat] = kstest2(backgd_values{ii,ik},backgd_values{il,ik},'Alpha',pAlpha_threshold);
            backgd_hMatrix{ik}(ii,il) = h;backgd_pvalMatrix{ik}(ii,il) = p;backgd_k2statMatrix{ik}(ii,il) = ks2stat;
            [h,p,ks2stat] = kstest2(signal_minus_backgd_values{ii,ik},signal_minus_backgd_values{il,ik},'Alpha',pAlpha_threshold);
            signal_minus_backgd_hMatrix{ik}(ii,il) = h;signal_minus_backgd_pvalMatrix{ik}(ii,il) = p;signal_minus_backgd_k2statMatrix{ik}(ii,il) = ks2stat;
            [h,p,ks2stat] = kstest2(SNR_values{ii,ik},SNR_values{il,ik},'Alpha',pAlpha_threshold);
            SNR_hMatrix{ik}(ii,il) = h;SNR_pvalMatrix{ik}(ii,il) = p;SNR_k2statMatrix{ik}(ii,il) = ks2stat;
            [h,p,ks2stat] = kstest2(spot_count_values{ii,ik},spot_count_values{il,ik},'Alpha',pAlpha_threshold);
            spotcount_hMatrix{ik}(ii,il) = h;spotcount_pvalMatrix{ik}(ii,il) = p;spotcount_k2statMatrix{ik}(ii,il) = ks2stat;
        end
    end
    [pval,tbl0,stats] = anova1(y_vector(ik,signal_values),g1_vector(ik,signal_values),'off');
    [compare,means,~,gnames] = multcompare(stats,"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
    tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
    tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
    anova1_results.signal_pval{ik} = pval;anova1_results.signal_stats{ik} = stats;
    anova1_results.signal_AnalysisOfVariance{ik} = tbl0;anova1_results.signal_CI{ik} = tbl1;anova1_results.signal_multiple_comparisons{ik} = tbl2;
    anova1_results.signal_multiple_comparisons{ik} = sortrows(anova1_results.signal_multiple_comparisons{ik},"P-value","ascend");
    [pval,tbl0,stats] = anova1(y_vector(ik,backgd_values),g1_vector(ik,backgd_values),'off');
    [compare,means,~,gnames] = multcompare(stats,"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
    tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
    tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
    anova1_results.backgd_pval{ik} = pval;anova1_results.backgd_stats{ik} = stats;
    anova1_results.backgd_AnalysisOfVariance{ik} = tbl0;anova1_results.backgd_CI{ik} = tbl1;anova1_results.backgd_multiple_comparisons{ik} = tbl2;
    anova1_results.backgd_multiple_comparisons{ik} = sortrows(anova1_results.backgd_multiple_comparisons{ik},"P-value","ascend");
    [pval,tbl0,stats] = anova1(y_vector(ik,signal_minus_backgd_values),g1_vector(ik,signal_minus_backgd_values),'off');
    [compare,means,~,gnames] = multcompare(stats,"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
    tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
    tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
    anova1_results.signal_minus_backgd_pval{ik} = pval;anova1_results.signal_minus_backgd_stats{ik} = stats;
    anova1_results.signal_minus_backgd_AnalysisOfVariance{ik} = tbl0;anova1_results.signal_minus_backgd_CI{ik} = tbl1;anova1_results.signal_minus_backgd_multiple_comparisons{ik} = tbl2;
    anova1_results.signal_minus_backgd_multiple_comparisons{ik} = sortrows(anova1_results.signal_minus_backgd_multiple_comparisons{ik},"P-value","ascend");
    [pval,tbl0,stats] = anova1(y_vector(ik,SNR_values),g1_vector(ik,SNR_values),'off');
    [compare,means,~,gnames] = multcompare(stats,"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
    tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
    tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
    anova1_results.SNR_pval{ik} = pval;anova1_results.SNR_stats{ik} = stats;
    anova1_results.SNR_AnalysisOfVariance{ik} = tbl0;anova1_results.SNR_CI{ik} = tbl1;anova1_results.SNR_multiple_comparisons{ik} = tbl2;
    anova1_results.SNR_multiple_comparisons{ik} = sortrows(anova1_results.SNR_multiple_comparisons{ik},"P-value","ascend");
    [pval,tbl0,stats] = anova1(y_vector(ik,spot_count_values),g1_vector(ik,spot_count_values),'off');
    [compare,means,~,gnames] = multcompare(stats,"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
    tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
    tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
    anova1_results.spotcount_pval{ik} = pval;anova1_results.spotcount_stats{ik} = stats;
    anova1_results.spotcount_AnalysisOfVariance{ik} = tbl0;anova1_results.spotcount_CI{ik} = tbl1;anova1_results.spotcount_multiple_comparisons{ik} = tbl2;
    anova1_results.spotcount_multiple_comparisons{ik} = sortrows(anova1_results.spotcount_multiple_comparisons{ik},"P-value","ascend");
    [pval,tbl0,stats] = anova1(y_vector2(ik,spot_count_callTable_merged_values),g1_vector2(ik,spot_count_callTable_merged_values),'off');
    [compare,means,~,gnames] = multcompare(stats,"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
    tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
    tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
    anova1_results.callTableSpotCount_pval{ik} = pval;anova1_results.callTableSpotCount_stats{ik} = stats;
    anova1_results.callTableSpotCount_AnalysisOfVariance{ik} = tbl0;anova1_results.callTableSpotCount_CI{ik} = tbl1;anova1_results.callTableSpotCount_multiple_comparisons{ik} = tbl2;
    anova1_results.callTableSpotCount_multiple_comparisons{ik} = sortrows(anova1_results.callTableSpotCount_multiple_comparisons{ik},"P-value","ascend");
    [pval,tbl0,stats] = anova1(y_vector2(ik,spot_count_above_noprobe_callTable_merged_values),g1_vector2(ik,spot_count_above_noprobe_callTable_merged_values),'off');
    [compare,means,~,gnames] = multcompare(stats,"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
    tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
    tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
    anova1_results.callTableAboveNoProbeSpotCount_pval{ik} = pval;anova1_results.callTableAboveNoProbeSpotCount_stats{ik} = stats;
    anova1_results.callTableAboveNoProbeSpotCount_AnalysisOfVariance{ik} = tbl0;anova1_results.callTableAboveNoProbeSpotCount_CI{ik} = tbl1;anova1_results.callTableAboveNoProbeSpotCount_multiple_comparisons{ik} = tbl2;
    anova1_results.callTableAboveNoProbeSpotCount_multiple_comparisons{ik} = sortrows(anova1_results.callTableAboveNoProbeSpotCount_multiple_comparisons{ik},"P-value","ascend");
    [pval,tbl0,stats] = anova1(y_vector2(ik,cell_spot_detection_threshold_values),g1_vector2(ik,cell_spot_detection_threshold_values),'off');
    [compare,means,~,gnames] = multcompare(stats,"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
    tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
    tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
    anova1_results.DetectionThreshold_pval{ik} = pval;anova1_results.DetectionThreshold_stats{ik} = stats;
    anova1_results.DetectionThreshold_AnalysisOfVariance{ik} = tbl0;anova1_results.DetectionThreshold_CI{ik} = tbl1;anova1_results.DetectionThreshold_multiple_comparisons{ik} = tbl2;
    anova1_results.DetectionThreshold_multiple_comparisons{ik} = sortrows(anova1_results.DetectionThreshold_multiple_comparisons{ik},"P-value","ascend");
    [pval,tbl0,stats] = anova1(y_vector2(ik,cell_above_noprobe_spot_detection_threshold_values),g1_vector2(ik,cell_above_noprobe_spot_detection_threshold_values),'off');
    [compare,means,~,gnames] = multcompare(stats,"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
    tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
    tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
    anova1_results.AboveNoProbeDetectionThreshold_pval{ik} = pval;anova1_results.AboveNoProbeDetectionThreshold_stats{ik} = stats;
    anova1_results.AboveNoProbeDetectionThreshold_AnalysisOfVariance{ik} = tbl0;anova1_results.AboveNoProbeDetectionThreshold_CI{ik} = tbl1;anova1_results.AboveNoProbeDetectionThreshold_multiple_comparisons{ik} = tbl2;
    anova1_results.AboveNoProbeDetectionThreshold_multiple_comparisons{ik} = sortrows(anova1_results.AboveNoProbeDetectionThreshold_multiple_comparisons{ik},"P-value","ascend");
end
numDims = {'dim1','dim12'};
anova2_results = struct();
for dims = 1:2
[pval,tbl0,stats,~] = anovan(y_multi_vector(signal_values),{g1_multi_vector(signal_values),g2_multi_vector(signal_values)},'model','interaction','varnames',{'Software','Replica'},"Display","off");
[compare,means,~,gnames] = multcompare(stats,"Dimension",[1:dims],"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
anova2_results.(numDims{dims}).signal_pval = pval;anova2_results.(numDims{dims}).signal_stats = stats;
anova2_results.(numDims{dims}).signal_AnalysisOfVariance = tbl0;anova2_results.(numDims{dims}).signal_CI = tbl1;
anova2_results.(numDims{dims}).signal_multiple_comparisons = tbl2;
anova2_results.(numDims{dims}).signal_multiple_comparisons = sortrows(anova2_results.(numDims{dims}).signal_multiple_comparisons,"P-value","ascend");
[pval,tbl0,stats,~] = anovan(y_multi_vector(backgd_values),{g1_multi_vector(backgd_values),g2_multi_vector(backgd_values)},'model','interaction','varnames',{'Software','Replica'},"Display","off");
[compare,means,~,gnames] = multcompare(stats,"Dimension",[1:dims],"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
anova2_results.(numDims{dims}).backgd_pval = pval;anova2_results.(numDims{dims}).backgd_stats = stats;
anova2_results.(numDims{dims}).backgd_AnalysisOfVariance = tbl0;anova2_results.(numDims{dims}).backgd_CI = tbl1;
anova2_results.(numDims{dims}).backgd_multiple_comparisons = tbl2;
anova2_results.(numDims{dims}).backgd_multiple_comparisons = sortrows(anova2_results.(numDims{dims}).backgd_multiple_comparisons,"P-value","ascend");
[pval,tbl0,stats,~]= anovan(y_multi_vector(signal_minus_backgd_values),{g1_multi_vector(signal_minus_backgd_values),g2_multi_vector(signal_minus_backgd_values)},'model','interaction','varnames',{'Software','Replica'},"Display","off");
[compare,means,~,gnames] = multcompare(stats,"Dimension",[1:dims],"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
anova2_results.(numDims{dims}).signal_minus_backgd_pval = pval;anova2_results.(numDims{dims}).signal_minus_backgd_stats = stats;
anova2_results.(numDims{dims}).signal_minus_backgd_AnalysisOfVariance = tbl0;anova2_results.(numDims{dims}).signal_minus_backgd_CI = tbl1;
anova2_results.(numDims{dims}).signal_minus_backgd_multiple_comparisons = tbl2;
anova2_results.(numDims{dims}).signal_minus_backgd_multiple_comparisons = sortrows(anova2_results.(numDims{dims}).signal_minus_backgd_multiple_comparisons,"P-value","ascend");
[pval,tbl0,stats,~] = anovan(y_multi_vector(SNR_values),{g1_multi_vector(SNR_values),g2_multi_vector(SNR_values)},'model','interaction','varnames',{'Software','Replica'},"Display","off");
[compare,means,~,gnames] = multcompare(stats,"Dimension",[1:dims],"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
anova2_results.(numDims{dims}).SNR_pval = pval;anova2_results.(numDims{dims}).SNR_stats = stats;
anova2_results.(numDims{dims}).SNR_AnalysisOfVariance = tbl0;anova2_results.(numDims{dims}).SNR_CI = tbl1;
anova2_results.(numDims{dims}).SNR_multiple_comparisons = tbl2;
anova2_results.(numDims{dims}).SNR_multiple_comparisons = sortrows(anova2_results.(numDims{dims}).SNR_multiple_comparisons,"P-value","ascend");
[pval,tbl0,stats,~]= anovan(y_multi_vector(spot_count_values),{g1_multi_vector(spot_count_values),g2_multi_vector(spot_count_values)},'model','interaction','varnames',{'Software','Replica'},"Display","off");
[compare,means,~,gnames] = multcompare(stats,"Dimension",[1:dims],"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
anova2_results.(numDims{dims}).spotcount_pval = pval;anova2_results.(numDims{dims}).spotcount_stats = stats;
anova2_results.(numDims{dims}).spotcount_AnalysisOfVariance = tbl0;anova2_results.(numDims{dims}).spotcount_CI = tbl1;
anova2_results.(numDims{dims}).spotcount_multiple_comparisons = tbl2;
anova2_results.(numDims{dims}).spotcount_multiple_comparisons = sortrows(anova2_results.(numDims{dims}).spotcount_multiple_comparisons,"P-value","ascend");
[pval,tbl0,stats,~]= anovan(y_multi_vector2(spot_count_callTable_merged_values),{g1_multi_vector2(spot_count_callTable_merged_values),g2_multi_vector2(spot_count_callTable_merged_values)},'model','interaction','varnames',{'Software','Replica'},"Display","off");
[compare,means,~,gnames] = multcompare(stats,"Dimension",[1:dims],"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
anova2_results.(numDims{dims}).callTableSpotCount_pval = pval;anova2_results.(numDims{dims}).callTableSpotCount_stats = stats;
anova2_results.(numDims{dims}).callTableSpotCount_AnalysisOfVariance = tbl0;anova2_results.(numDims{dims}).callTableSpotCount_CI = tbl1;
anova2_results.(numDims{dims}).callTableSpotCount_multiple_comparisons = tbl2;
anova2_results.(numDims{dims}).callTableSpotCount_multiple_comparisons = sortrows(anova2_results.(numDims{dims}).callTableSpotCount_multiple_comparisons,"P-value","ascend");
[pval,tbl0,stats,~]= anovan(y_multi_vector2(spot_count_above_noprobe_callTable_merged_values),{g1_multi_vector2(spot_count_above_noprobe_callTable_merged_values),g2_multi_vector2(spot_count_above_noprobe_callTable_merged_values)},'model','interaction','varnames',{'Software','Replica'},"Display","off");
[compare,means,~,gnames] = multcompare(stats,"Dimension",[1:dims],"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
anova2_results.(numDims{dims}).callTableAboveNoProbeSpotCount_pval = pval;anova2_results.(numDims{dims}).callTableAboveNoProbeSpotCount_stats = stats;
anova2_results.(numDims{dims}).callTableAboveNoProbeSpotCount_AnalysisOfVariance = tbl0;anova2_results.(numDims{dims}).callTableAboveNoProbeSpotCount_CI = tbl1;
anova2_results.(numDims{dims}).callTableAboveNoProbeSpotCount_multiple_comparisons = tbl2;
anova2_results.(numDims{dims}).callTableAboveNoProbeSpotCount_multiple_comparisons = sortrows(anova2_results.(numDims{dims}).callTableAboveNoProbeSpotCount_multiple_comparisons,"P-value","ascend");
[pval,tbl0,stats,~]= anovan(y_multi_vector2(cell_spot_detection_threshold_values),{g1_multi_vector2(cell_spot_detection_threshold_values),g2_multi_vector2(cell_spot_detection_threshold_values)},'model','interaction','varnames',{'Software','Replica'},"Display","off");
[compare,means,~,gnames] = multcompare(stats,"Dimension",[1:dims],"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
anova2_results.(numDims{dims}).DetectionThreshold_pval = pval;anova2_results.(numDims{dims}).DetectionThreshold_stats = stats;
anova2_results.(numDims{dims}).DetectionThreshold_AnalysisOfVariance = tbl0;anova2_results.(numDims{dims}).DetectionThreshold_CI = tbl1;
anova2_results.(numDims{dims}).DetectionThreshold_multiple_comparisons = tbl2;
anova2_results.(numDims{dims}).DetectionThreshold_multiple_comparisons = sortrows(anova2_results.(numDims{dims}).DetectionThreshold_multiple_comparisons,"P-value","ascend");
[pval,tbl0,stats,~]= anovan(y_multi_vector2(cell_above_noprobe_spot_detection_threshold_values),{g1_multi_vector2(cell_above_noprobe_spot_detection_threshold_values),g2_multi_vector2(cell_above_noprobe_spot_detection_threshold_values)},'model','interaction','varnames',{'Software','Replica'},"Display","off");
[compare,means,~,gnames] = multcompare(stats,"Dimension",[1:dims],"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
anova2_results.(numDims{dims}).AboveNoProbeDetectionThreshold_pval = pval;anova2_results.(numDims{dims}).AboveNoProbeDetectionThreshold_stats = stats;
anova2_results.(numDims{dims}).AboveNoProbeDetectionThreshold_AnalysisOfVariance = tbl0;anova2_results.(numDims{dims}).AboveNoProbeDetectionThreshold_CI = tbl1;
anova2_results.(numDims{dims}).AboveNoProbeDetectionThreshold_multiple_comparisons = tbl2;
anova2_results.(numDims{dims}).AboveNoProbeDetectionThreshold_multiple_comparisons = sortrows(anova2_results.(numDims{dims}).AboveNoProbeDetectionThreshold_multiple_comparisons,"P-value","ascend");
end
for ii = 1:length(SoftwareNames_Spelling)
        for il = 1:length(SoftwareNames_Spelling)
            [h,p,ks2stat] = kstest2(cell2mat(cellfun(@(x) x', signal_values(ii,:),'Un',0))',cell2mat(cellfun(@(x) x', signal_values(il,:),'Un',0))','Alpha',pAlpha_threshold);
            signal_hMatrix_multi1(ii,il) = h;
            signal_pvalMatrix_multi1(ii,il) = p;
            signal_k2statMatrix_multi1(ii,il) = ks2stat;
            [h,p,ks2stat] = kstest2(cell2mat(cellfun(@(x) x', backgd_values(ii,:),'Un',0))',cell2mat(cellfun(@(x) x', backgd_values(il,:),'Un',0))','Alpha',pAlpha_threshold);
            backgd_hMatrix_multi1(ii,il) = h;
            backgd_pvalMatrix_multi1(ii,il) = p;
            backgd_k2statMatrix_multi1(ii,il) = ks2stat;
            [h,p,ks2stat] = kstest2(cell2mat(cellfun(@(x) x', signal_minus_backgd_values(ii,:),'Un',0))',cell2mat(cellfun(@(x) x', signal_minus_backgd_values(il,:),'Un',0))','Alpha',pAlpha_threshold);
            signal_minus_backgd_hMatrix_multi1(ii,il) = h;
            signal_minus_backgd_pvalMatrix_multi1(ii,il) = p;
            signal_minus_backgd_k2statMatrix_multi1(ii,il) = ks2stat;
            [h,p,ks2stat] = kstest2(cell2mat(cellfun(@(x) x', SNR_values(ii,:),'Un',0))',cell2mat(cellfun(@(x) x', SNR_values(il,:),'Un',0))','Alpha',pAlpha_threshold);
            SNR_hMatrix_multi1(ii,il) = h;
            SNR_pvalMatrix_multi1(ii,il) = p;
            SNR_k2statMatrix_multi1(ii,il) = ks2stat;
            [h,p,ks2stat] = kstest2(cell2mat(cellfun(@(x) x', spot_count_values(ii,:),'Un',0))',cell2mat(cellfun(@(x) x', spot_count_values(il,:),'Un',0))','Alpha',pAlpha_threshold);
            spotcount_hMatrix_multi1(ii,il) = h;
            spotcount_pvalMatrix_multi1(ii,il) = p;
            spotcount_k2statMatrix_multi1(ii,il) = ks2stat;
        end
 end
for ii = 1:size(software_and_replica_pairs,1)
    for il = 1:size(software_and_replica_pairs,1)
        [h,p,ks2stat] = kstest2(signal_values{software_and_replica_pairs(ii,2),software_and_replica_pairs(ii,1)},signal_values{software_and_replica_pairs(il,2),software_and_replica_pairs(il,1)},'Alpha',pAlpha_threshold);
        signal_hMatrix_multi2(ii,il) = h;
        signal_pvalMatrix_multi2(ii,il) = p;
        signal_k2statMatrix_multi2(ii,il) = ks2stat;
        [h,p,ks2stat] = kstest2(backgd_values{software_and_replica_pairs(ii,2),software_and_replica_pairs(ii,1)},backgd_values{software_and_replica_pairs(il,2),software_and_replica_pairs(il,1)},'Alpha',pAlpha_threshold);
        backgd_hMatrix_multi2(ii,il) = h;
        backgd_pvalMatrix_multi2(ii,il) = p;
        backgd_k2statMatrix_multi2(ii,il) = ks2stat;
        [h,p,ks2stat] = kstest2(signal_minus_backgd_values{software_and_replica_pairs(ii,2),software_and_replica_pairs(ii,1)},signal_minus_backgd_values{software_and_replica_pairs(il,2),software_and_replica_pairs(il,1)},'Alpha',pAlpha_threshold);
        signal_minus_backgd_hMatrix_multi2(ii,il) = h;
        signal_minus_backgd_pvalMatrix_multi2(ii,il) = p;
        signal_minus_backgd_k2statMatrix_multi2(ii,il) = ks2stat;
        [h,p,ks2stat] = kstest2(SNR_values{software_and_replica_pairs(ii,2),software_and_replica_pairs(ii,1)},SNR_values{software_and_replica_pairs(il,2),software_and_replica_pairs(il,1)},'Alpha',pAlpha_threshold);
        SNR_hMatrix_multi2(ii,il) = h;
        SNR_pvalMatrix_multi2(ii,il) = p;
        SNR_k2statMatrix_multi2(ii,il) = ks2stat;
        [h,p,ks2stat] = kstest2(spot_count_values{software_and_replica_pairs(ii,2),software_and_replica_pairs(ii,1)},spot_count_values{software_and_replica_pairs(il,2),software_and_replica_pairs(il,1)},'Alpha',pAlpha_threshold);
        spotcount_hMatrix_multi2(ii,il) = h;
        spotcount_pvalMatrix_multi2(ii,il) = p;
        spotcount_k2statMatrix_multi2(ii,il) = ks2stat;
    end
end
within_replica_matrix = triu(ones(length(SoftwareNames_Spelling)))-diag(ones(1,length(SoftwareNames_Spelling)));
[first_software,second_software] = find(within_replica_matrix);
within_replica_loc = sub2ind(size(within_replica_matrix),first_software,second_software);
first_software_list = cellfun(@(x) strcat('Software=',x),SoftwareNames(first_software),'Un',0);
second_software_list = cellfun(@(x) strcat('Software=',x),SoftwareNames(second_software),'Un',0);
between_replica_matrix = triu(ones(length(SoftwareNames_Spelling)*length(Replicas_To_Show)))-diag(ones(1,length(SoftwareNames_Spelling)*length(Replicas_To_Show)));
[first_pair_value,second_pair_value] = find(between_replica_matrix);
between_replica_loc = sub2ind(size(between_replica_matrix),first_pair_value,second_pair_value);
first_software_and_replica_list = arrayfun(@(x) ...
    strcat('Software=',SoftwareNames{software_and_replica_pairs(first_pair_value(x),2)},...
    ',Replica=',Replica{Replicas_To_Show(software_and_replica_pairs(first_pair_value(x),1))})...
    ,1:length(first_pair_value),'Un',0);
second_software_and_replica_list = arrayfun(@(x) ...
    strcat('Software=',SoftwareNames{software_and_replica_pairs(second_pair_value(x),2)},...
    ',Replica=',Replica{Replicas_To_Show(software_and_replica_pairs(second_pair_value(x),1))})...
    ,1:length(second_pair_value),'Un',0);
ks2stat_signal_comparison_table = cell(1,length(Replicas_To_Show));
ks2stat_backgd_comparison_table = cell(1,length(Replicas_To_Show));
ks2stat_signal_minus_backgd_comparison_table = cell(1,length(Replicas_To_Show));
ks2stat_SNR_comparison_table = cell(1,length(Replicas_To_Show));
ks2stat_spotcount_comparison_table = cell(1,length(Replicas_To_Show));
for ik = 1:length(Replicas_To_Show)
ks2stat_signal_comparison_table{ik}= table(first_software_list',second_software_list',signal_k2statMatrix{ik}(within_replica_loc),signal_pvalMatrix{ik}(within_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_backgd_comparison_table{ik}= table(first_software_list',second_software_list',backgd_k2statMatrix{ik}(within_replica_loc),backgd_pvalMatrix{ik}(within_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_signal_minus_backgd_comparison_table{ik}= table(first_software_list',second_software_list',signal_minus_backgd_k2statMatrix{ik}(within_replica_loc),signal_minus_backgd_pvalMatrix{ik}(within_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_SNR_comparison_table{ik}= table(first_software_list',second_software_list',SNR_k2statMatrix{ik}(within_replica_loc),SNR_pvalMatrix{ik}(within_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_spotcount_comparison_table{ik}= table(first_software_list',second_software_list',spotcount_k2statMatrix{ik}(within_replica_loc),spotcount_pvalMatrix{ik}(within_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_signal_comparison_table{ik} = sortrows(ks2stat_signal_comparison_table{ik},"P-value","ascend");
ks2stat_backgd_comparison_table{ik} = sortrows(ks2stat_backgd_comparison_table{ik},"P-value","ascend");
ks2stat_signal_minus_backgd_comparison_table{ik} = sortrows(ks2stat_signal_minus_backgd_comparison_table{ik},"P-value","ascend");
ks2stat_SNR_comparison_table{ik} = sortrows(ks2stat_SNR_comparison_table{ik},"P-value","ascend");
ks2stat_spotcount_comparison_table{ik} = sortrows(ks2stat_spotcount_comparison_table{ik},"P-value","ascend");
end
ks2stat_multi1_signal_comparison_table = table(first_software_list',second_software_list',signal_k2statMatrix_multi1(within_replica_loc),signal_pvalMatrix_multi1(within_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_multi2_signal_comparison_table = table(first_software_and_replica_list',second_software_and_replica_list',signal_k2statMatrix_multi2(between_replica_loc),signal_pvalMatrix_multi2(between_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_multi1_backgd_comparison_table = table(first_software_list',second_software_list',backgd_k2statMatrix_multi1(within_replica_loc),backgd_pvalMatrix_multi1(within_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_multi2_backgd_comparison_table = table(first_software_and_replica_list',second_software_and_replica_list',backgd_k2statMatrix_multi2(between_replica_loc),backgd_pvalMatrix_multi2(between_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_multi1_signal_minus_backgd_comparison_table = table(first_software_list',second_software_list',signal_minus_backgd_k2statMatrix_multi1(within_replica_loc),signal_minus_backgd_pvalMatrix_multi1(within_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_multi2_signal_minus_backgd_comparison_table = table(first_software_and_replica_list',second_software_and_replica_list',signal_minus_backgd_k2statMatrix_multi2(between_replica_loc),signal_minus_backgd_pvalMatrix_multi2(between_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_multi1_SNR_comparison_table = table(first_software_list',second_software_list',SNR_k2statMatrix_multi1(within_replica_loc),SNR_pvalMatrix_multi1(within_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_multi2_SNR_comparison_table = table(first_software_and_replica_list',second_software_and_replica_list',SNR_k2statMatrix_multi2(between_replica_loc),SNR_pvalMatrix_multi2(between_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_multi1_spotcount_comparison_table = table(first_software_list',second_software_list',spotcount_k2statMatrix_multi1(within_replica_loc),spotcount_pvalMatrix_multi1(within_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_multi2_spotcount_comparison_table = table(first_software_and_replica_list',second_software_and_replica_list',spotcount_k2statMatrix_multi2(between_replica_loc),spotcount_pvalMatrix_multi2(between_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_multi1_signal_comparison_table = sortrows(ks2stat_multi1_signal_comparison_table,"P-value","ascend");
ks2stat_multi2_signal_comparison_table = sortrows(ks2stat_multi2_signal_comparison_table,"P-value","ascend");
ks2stat_multi1_backgd_comparison_table = sortrows(ks2stat_multi1_backgd_comparison_table,"P-value","ascend");
ks2stat_multi2_backgd_comparison_table = sortrows(ks2stat_multi2_backgd_comparison_table,"P-value","ascend");
ks2stat_multi1_signal_minus_backgd_comparison_table = sortrows(ks2stat_multi1_signal_minus_backgd_comparison_table,"P-value","ascend");
ks2stat_multi2_signal_minus_backgd_comparison_table = sortrows(ks2stat_multi2_signal_minus_backgd_comparison_table,"P-value","ascend");
ks2stat_multi1_SNR_comparison_table = sortrows(ks2stat_multi1_SNR_comparison_table,"P-value","ascend");
ks2stat_multi2_SNR_comparison_table = sortrows(ks2stat_multi2_SNR_comparison_table,"P-value","ascend");
ks2stat_multi1_spotcount_comparison_table = sortrows(ks2stat_multi1_spotcount_comparison_table,"P-value","ascend");
ks2stat_multi2_spotcount_comparison_table = sortrows(ks2stat_multi2_spotcount_comparison_table,"P-value","ascend");
software_and_replica_pairs = allcomb(1:length(Replicas_To_Show),1:length(Conditions));
y_vector = @(ik,var) cell2mat(arrayfun(@(x)var{x,ik}',1:length(Conditions),'Un',0));
g1_vector = @(ik,var) string(Conditions(cell2mat(arrayfun(@(x)x*ones(size(var{x,ik}))',1:length(Conditions),'Un',0))));
g2_vector = @(ik,var) string(Replica(Replicas_To_Show(ik)*cell2mat(arrayfun(@(x) ones(size(var{x,ik}))',1:length(Conditions),'Un',0))));
y_multi_vector = @(var) cell2mat(arrayfun(@(ik) cell2mat(arrayfun(@(x)var{x,ik}',1:length(Conditions),'Un',0)),1:length(Replicas_To_Show),'Un',0));
g1_multi_vector = @(var) string(Conditions(cell2mat(arrayfun(@(ik) cell2mat(arrayfun(@(x)x*ones(size(var{x,ik}))',1:length(Conditions),'Un',0)),1:length(Replicas_To_Show),'Un',0))));
g2_multi_vector = @(var) string(Replica(Replicas_To_Show(cell2mat(arrayfun(@(ik) ik*cell2mat(arrayfun(@(x) ones(size(var{x,ik}))',1:length(Conditions),'Un',0)),1:length(Replicas_To_Show),'Un',0)))));
for ik = 1:length(Replicas_To_Show)
    cellbkg_mean_k2statMatrix{ik} = NaN(length(Conditions));
    cellbkg_fano_k2statMatrix{ik} = NaN(length(Conditions));
    cellbkg_mean_pvalMatrix{ik} = NaN(length(Conditions));
    cellbkg_fano_pvalMatrix{ik} = NaN(length(Conditions));
    cellbkg_mean_hMatrix{ik} = NaN(length(Conditions));
    cellbkg_fano_hMatrix{ik} = NaN(length(Conditions));
    for ii = 1:length(Conditions)
        for il = 1:length(Conditions)
            [h,p,ks2stat] = kstest2(cellbkg_mean_values{ii,ik},cellbkg_mean_values{il,ik},'Alpha',pAlpha_threshold);
            cellbkg_mean_hMatrix{ik}(ii,il) = h;cellbkg_mean_pvalMatrix{ik}(ii,il) = p;cellbkg_mean_k2statMatrix{ik}(ii,il) = ks2stat;
            [h,p,ks2stat] = kstest2(cellbkg_fano_values{ii,ik},cellbkg_fano_values{il,ik},'Alpha',pAlpha_threshold);
            cellbkg_fano_hMatrix{ik}(ii,il) = h;cellbkg_fano_pvalMatrix{ik}(ii,il) = p;cellbkg_fano_k2statMatrix{ik}(ii,il) = ks2stat;
        end
    end
    [pval,tbl0,stats] = anova1(y_vector(ik,cellbkg_mean_values),g1_vector(ik,cellbkg_mean_values),'off');
    [compare,means,~,gnames] = multcompare(stats,"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
    tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
    tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
    anova1_results.cellbkg_mean_pval{ik} = pval;anova1_results.cellbkg_mean_stats{ik} = stats;
    anova1_results.cellbkg_mean_AnalysisOfVariance{ik} = tbl0;anova1_results.cellbkg_mean_CI{ik} = tbl1;anova1_results.cellbkg_mean_multiple_comparisons{ik} = tbl2;
    anova1_results.cellbkg_mean_multiple_comparisons{ik} = sortrows(anova1_results.cellbkg_mean_multiple_comparisons{ik},"P-value","ascend");
    [pval,tbl0,stats] = anova1(y_vector(ik,cellbkg_fano_values),g1_vector(ik,cellbkg_fano_values),'off');
    [compare,means,~,gnames] = multcompare(stats,"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
    tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
    tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
    anova1_results.cellbkg_fano_pval{ik} = pval;anova1_results.cellbkg_fano_stats{ik} = stats;
    anova1_results.cellbkg_fano_AnalysisOfVariance{ik} = tbl0;anova1_results.cellbkg_fano_CI{ik} = tbl1;anova1_results.cellbkg_fano_multiple_comparisons{ik} = tbl2;
    anova1_results.cellbkg_fano_multiple_comparisons{ik} = sortrows(anova1_results.cellbkg_fano_multiple_comparisons{ik},"P-value","ascend");
end
numDims = {'dim1','dim12'};
for dims = 1:2
[pval,tbl0,stats,~] = anovan(y_multi_vector(cellbkg_mean_values),{g1_multi_vector(cellbkg_mean_values),g2_multi_vector(cellbkg_mean_values)},'model','interaction','varnames',{'Software','Replica'},"Display","off");
[compare,means,~,gnames] = multcompare(stats,"Dimension",[1:dims],"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
anova2_results.(numDims{dims}).cellbkg_mean_pval = pval;anova2_results.(numDims{dims}).cellbkg_mean_stats = stats;
anova2_results.(numDims{dims}).cellbkg_mean_AnalysisOfVariance = tbl0;anova2_results.(numDims{dims}).cellbkg_mean_CI = tbl1;
anova2_results.(numDims{dims}).cellbkg_mean_multiple_comparisons = tbl2;
anova2_results.(numDims{dims}).cellbkg_mean_multiple_comparisons = sortrows(anova2_results.(numDims{dims}).cellbkg_mean_multiple_comparisons,"P-value","ascend");
[pval,tbl0,stats,~] = anovan(y_multi_vector(cellbkg_fano_values),{g1_multi_vector(cellbkg_fano_values),g2_multi_vector(cellbkg_fano_values)},'model','interaction','varnames',{'Software','Replica'},"Display","off");
[compare,means,~,gnames] = multcompare(stats,"Dimension",[1:dims],"CriticalValueType",compareType,"Alpha",pAlpha_threshold,"Display","off");
tbl1 = array2table(means,"RowNames",gnames,"VariableNames",["Mean","Standard Error"]);
tbl2 = array2table(compare,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A") = gnames(tbl2.("Group A"));tbl2.("Group B") = gnames(tbl2.("Group B"));
anova2_results.(numDims{dims}).cellbkg_fano_pval = pval;anova2_results.(numDims{dims}).cellbkg_fano_stats = stats;
anova2_results.(numDims{dims}).cellbkg_fano_AnalysisOfVariance = tbl0;anova2_results.(numDims{dims}).cellbkg_fano_CI = tbl1;
anova2_results.(numDims{dims}).cellbkg_fano_multiple_comparisons = tbl2;
anova2_results.(numDims{dims}).cellbkg_fano_multiple_comparisons = sortrows(anova2_results.(numDims{dims}).cellbkg_fano_multiple_comparisons,"P-value","ascend");
end
for ii = 1:length(Conditions)
        for il = 1:length(Conditions)
            [h,p,ks2stat] = kstest2(cell2mat(cellfun(@(x) x', cellbkg_mean_values(ii,:),'Un',0))',cell2mat(cellfun(@(x) x', cellbkg_mean_values(il,:),'Un',0))','Alpha',pAlpha_threshold);
            cellbkg_mean_hMatrix_multi1(ii,il) = h;
            cellbkg_mean_pvalMatrix_multi1(ii,il) = p;
            cellbkg_mean_k2statMatrix_multi1(ii,il) = ks2stat;
            [h,p,ks2stat] = kstest2(cell2mat(cellfun(@(x) x', cellbkg_fano_values(ii,:),'Un',0))',cell2mat(cellfun(@(x) x', cellbkg_fano_values(il,:),'Un',0))','Alpha',pAlpha_threshold);
            cellbkg_fano_hMatrix_multi1(ii,il) = h;
            cellbkg_fano_pvalMatrix_multi1(ii,il) = p;
            cellbkg_fano_k2statMatrix_multi1(ii,il) = ks2stat;
        end
 end
for ii = 1:size(software_and_replica_pairs,1)
    for il = 1:size(software_and_replica_pairs,1)
        [h,p,ks2stat] = kstest2(cellbkg_mean_values{software_and_replica_pairs(ii,2),software_and_replica_pairs(ii,1)},cellbkg_mean_values{software_and_replica_pairs(il,2),software_and_replica_pairs(il,1)},'Alpha',pAlpha_threshold);
        cellbkg_mean_hMatrix_multi2(ii,il) = h;
        cellbkg_mean_pvalMatrix_multi2(ii,il) = p;
        cellbkg_mean_k2statMatrix_multi2(ii,il) = ks2stat;
        [h,p,ks2stat] = kstest2(cellbkg_fano_values{software_and_replica_pairs(ii,2),software_and_replica_pairs(ii,1)},cellbkg_fano_values{software_and_replica_pairs(il,2),software_and_replica_pairs(il,1)},'Alpha',pAlpha_threshold);
        cellbkg_fano_hMatrix_multi2(ii,il) = h;
        cellbkg_fano_pvalMatrix_multi2(ii,il) = p;
        cellbkg_fano_k2statMatrix_multi2(ii,il) = ks2stat;
    end
end
within_replica_matrix = triu(ones(length(Conditions)))-diag(ones(1,length(Conditions)));
[first_software,second_software] = find(within_replica_matrix);
within_replica_loc = sub2ind(size(within_replica_matrix),first_software,second_software);
first_software_list = cellfun(@(x) strcat('Software=',x),Conditions(first_software),'Un',0);
second_software_list = cellfun(@(x) strcat('Software=',x),Conditions(second_software),'Un',0);
between_replica_matrix = triu(ones(length(Conditions)*length(Replicas_To_Show)))-diag(ones(1,length(Conditions)*length(Replicas_To_Show)));
[first_pair_value,second_pair_value] = find(between_replica_matrix);
between_replica_loc = sub2ind(size(between_replica_matrix),first_pair_value,second_pair_value);
first_software_and_replica_list = arrayfun(@(x) ...
    strcat('Software=',Conditions{software_and_replica_pairs(first_pair_value(x),2)},...
    ',Replica=',Replica{Replicas_To_Show(software_and_replica_pairs(first_pair_value(x),1))})...
    ,1:length(first_pair_value),'Un',0);
second_software_and_replica_list = arrayfun(@(x) ...
    strcat('Software=',Conditions{software_and_replica_pairs(second_pair_value(x),2)},...
    ',Replica=',Replica{Replicas_To_Show(software_and_replica_pairs(second_pair_value(x),1))})...
    ,1:length(second_pair_value),'Un',0);
ks2stat_cellbkg_mean_comparison_table = cell(1,length(Replicas_To_Show));
ks2stat_cellbkg_fano_comparison_table = cell(1,length(Replicas_To_Show));
for ik = 1:length(Replicas_To_Show)
ks2stat_cellbkg_mean_comparison_table{ik}= table(first_software_list',second_software_list',cellbkg_mean_k2statMatrix{ik}(within_replica_loc),cellbkg_mean_pvalMatrix{ik}(within_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_cellbkg_fano_comparison_table{ik}= table(first_software_list',second_software_list',cellbkg_fano_k2statMatrix{ik}(within_replica_loc),cellbkg_fano_pvalMatrix{ik}(within_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_cellbkg_mean_comparison_table{ik} = sortrows(ks2stat_cellbkg_mean_comparison_table{ik},"P-value","ascend");
ks2stat_cellbkg_fano_comparison_table{ik} = sortrows(ks2stat_cellbkg_fano_comparison_table{ik},"P-value","ascend");
end
ks2stat_multi1_cellbkg_mean_comparison_table = table(first_software_list',second_software_list',cellbkg_mean_k2statMatrix_multi1(within_replica_loc),cellbkg_mean_pvalMatrix_multi1(within_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_multi2_cellbkg_mean_comparison_table = table(first_software_and_replica_list',second_software_and_replica_list',cellbkg_mean_k2statMatrix_multi2(between_replica_loc),cellbkg_mean_pvalMatrix_multi2(between_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_multi1_cellbkg_fano_comparison_table = table(first_software_list',second_software_list',cellbkg_fano_k2statMatrix_multi1(within_replica_loc),cellbkg_fano_pvalMatrix_multi1(within_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_multi2_cellbkg_fano_comparison_table = table(first_software_and_replica_list',second_software_and_replica_list',cellbkg_fano_k2statMatrix_multi2(between_replica_loc),cellbkg_fano_pvalMatrix_multi2(between_replica_loc),'VariableNames',{'Group A','Group B','KS2-statistic','P-value'});
ks2stat_multi1_cellbkg_mean_comparison_table = sortrows(ks2stat_multi1_cellbkg_mean_comparison_table,"P-value","ascend");
ks2stat_multi2_cellbkg_mean_comparison_table = sortrows(ks2stat_multi2_cellbkg_mean_comparison_table,"P-value","ascend");
ks2stat_multi1_cellbkg_fano_comparison_table = sortrows(ks2stat_multi1_cellbkg_fano_comparison_table,"P-value","ascend");
ks2stat_multi2_cellbkg_fano_comparison_table = sortrows(ks2stat_multi2_cellbkg_fano_comparison_table,"P-value","ascend");
ks2stat_results = struct();
ks2stat_results.signal_hMatrix = signal_hMatrix;            
ks2stat_results.signal_pvalMatrix = signal_pvalMatrix;
ks2stat_results.signal_k2statMatrix = signal_k2statMatrix;
ks2stat_results.backgd_hMatrix = backgd_hMatrix;
ks2stat_results.backgd_pvalMatrix = backgd_pvalMatrix;
ks2stat_results.backgd_k2statMatrix = backgd_k2statMatrix;
ks2stat_results.signal_minus_backgd_hMatrix = signal_minus_backgd_hMatrix;
ks2stat_results.signal_minus_backgd_pvalMatrix = signal_minus_backgd_pvalMatrix;
ks2stat_results.signal_minus_backgd_k2statMatrix = signal_minus_backgd_k2statMatrix;
ks2stat_results.SNR_hMatrix = SNR_hMatrix; 
ks2stat_results.SNR_pvalMatrix = SNR_pvalMatrix;
ks2stat_results.SNR_k2statMatrix = SNR_k2statMatrix;
ks2stat_results.spotcount_hMatrix = spotcount_hMatrix;
ks2stat_results.spotcount_pvalMatrix = spotcount_pvalMatrix;
ks2stat_results.spotcount_k2statMatrix = spotcount_k2statMatrix;
ks2stat_results.ks2stat_signal_comparison_table = ks2stat_signal_comparison_table;
ks2stat_results.ks2stat_backgd_comparison_table = ks2stat_backgd_comparison_table;
ks2stat_results.ks2stat_signal_minus_backgd_comparison_table = ks2stat_signal_minus_backgd_comparison_table;
ks2stat_results.ks2stat_SNR_comparison_table = ks2stat_SNR_comparison_table;
ks2stat_results.ks2stat_spotcount_comparison_table = ks2stat_spotcount_comparison_table;
ks2stat_results.ks2stat_multi1_signal_comparison_table = ks2stat_multi1_signal_comparison_table;
ks2stat_results.ks2stat_multi2_signal_comparison_table = ks2stat_multi2_signal_comparison_table;
ks2stat_results.ks2stat_multi1_backgd_comparison_table = ks2stat_multi1_backgd_comparison_table;
ks2stat_results.ks2stat_multi2_backgd_comparison_table = ks2stat_multi2_backgd_comparison_table;
ks2stat_results.ks2stat_multi1_signal_minus_backgd_comparison_table = ks2stat_multi1_signal_minus_backgd_comparison_table;
ks2stat_results.ks2stat_multi2_signal_minus_backgd_comparison_table = ks2stat_multi2_signal_minus_backgd_comparison_table;
ks2stat_results.ks2stat_multi1_SNR_comparison_table = ks2stat_multi1_SNR_comparison_table;
ks2stat_results.ks2stat_multi2_SNR_comparison_table = ks2stat_multi2_SNR_comparison_table;
ks2stat_results.ks2stat_multi1_spotcount_comparison_table = ks2stat_multi1_spotcount_comparison_table;
ks2stat_results.ks2stat_multi2_spotcount_comparison_table = ks2stat_multi2_spotcount_comparison_table;
ks2stat_results.cellbkg_mean_hMatrix = cellbkg_mean_hMatrix;            
ks2stat_results.cellbkg_mean_pvalMatrix = cellbkg_mean_pvalMatrix;
ks2stat_results.cellbkg_mean_k2statMatrix = cellbkg_mean_k2statMatrix;
ks2stat_results.cellbkg_fano_hMatrix = cellbkg_fano_hMatrix;
ks2stat_results.cellbkg_fano_pvalMatrix = cellbkg_fano_pvalMatrix;
ks2stat_results.cellbkg_fano_k2statMatrix = cellbkg_fano_k2statMatrix;
ks2stat_results.ks2stat_cellbkg_mean_comparison_table = ks2stat_cellbkg_mean_comparison_table;
ks2stat_results.ks2stat_cellbkg_fano_comparison_table = ks2stat_cellbkg_fano_comparison_table;
ks2stat_results.ks2stat_multi1_cellbkg_mean_comparison_table = ks2stat_multi1_cellbkg_mean_comparison_table;
ks2stat_results.ks2stat_multi2_cellbkg_mean_comparison_table = ks2stat_multi2_cellbkg_mean_comparison_table;
ks2stat_results.ks2stat_multi1_cellbkg_fano_comparison_table = ks2stat_multi1_cellbkg_fano_comparison_table;
ks2stat_results.ks2stat_multi2_cellbkg_fano_comparison_table = ks2stat_multi2_cellbkg_fano_comparison_table;
figData.ks2stat_results = ks2stat_results;clear ks2stat_results
figData.anova1_results = anova1_results;clear anova1_results
figData.anova2_results = anova2_results;clear anova2_results


for ii = 1:length(SoftwareNames)
    for ik = 1:length(Replicas_To_Show)
          P_joint = zeros(max(figData.spot_count_each_type_values{ii,ik}(:,2))+1,...
              max(figData.spot_count_each_type_values{ii,ik}(:,1))+1);
          joint_expr_options = unique(figData.spot_count_each_type_values{ii,ik}(:,[2 1]),'rows');
          for il = 1:size(joint_expr_options,1)
          P_joint(joint_expr_options(il,1)+1,joint_expr_options(il,2)+1) = sum(ismember(figData.spot_count_each_type_values{ii,ik}(:,[2 1]),joint_expr_options(il,:),'rows'))/size(figData.spot_count_each_type_values{ii,ik},1);
          end
    end
end

end





%save('figData.mat','figData','-v7.3')







