function [Nvec_RNAmulti,RNAOFF_Score,RNASpecificity_Score,NumRNAOffTargetOptions,Probes_WithNRNAOFF,DNAOFF_Score,DNASpecificity_Score,NumDNAOffTargetOptions,Probes_WithNDNAOFF,Cout] = A0_BasicDesignerStats_V2(targetTypes,removeUndesiredIsos,gene_table,settings,FoldName,DoesProbeBindSite,Kon,Kb,Kb_Complement,EKernel)
% This function determine the metrics and statistics used for selecting RNA-FISH probes.
most_recent_num_local = settings.num_parpool_local;
designerName = settings.designerName;
FolderRootName = settings.FolderRootName;
TranscriptName = settings.GeneName;
T_hybrid = settings.HybridizationTemperature;
Cout = cell(1,2);
Cout{1} = cell(1,7);
Cout{2} = cell(1,9);
if (settings.clusterStatus)
    most_recent_num = double(string(getenv('SLURM_JOB_CPUS_PER_NODE')));
else
    most_recent_num = most_recent_num_local;
end
probeBatchSize = settings.BLASTbatchSize;
targetBatchSize = settings.TargetBatchSize;
RNASpecificity_Score = zeros(1,size(DoesProbeBindSite,1));RNAOFF_Score = zeros(1,size(DoesProbeBindSite,1));
DNASpecificity_Score = zeros(1,size(DoesProbeBindSite,1));DNAOFF_Score = zeros(1,size(DoesProbeBindSite,1));
NumRNAOffTargetOptions = [];Probes_WithNRNAOFF = [];
NumDNAOffTargetOptions = [];Probes_WithNDNAOFF = [];
isRNA = targetTypes(1);
isDNA = targetTypes(2);
AllowableProbes = 1:size(DoesProbeBindSite,1);
gene_table = sortrows(gene_table,[7 6],'ascend');
gene_table = gene_table(gene_table.Match>=settings.MinHomologySearchTargetSize,:);
MinusStrandedHits = find(contains(gene_table.Strand,'Minus'));
gene_table_NamesZ = convertCharsToStrings(gene_table.Name);
contains_RNA = find(ismember(gene_table_NamesZ,settings.RNAdbParser));clear gene_table_NamesZ
RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);clear contains_RNA
gene_table = gene_table(setdiff(1:size(gene_table,1),RNA_MissedFilteredHits),:);clear RNA_MissedFilteredHits
gene_table.Ax = min(gene_table.SubjectIndices,[],2);
gene_table.Bx = max(gene_table.SubjectIndices,[],2);
gene_table = sortrows(gene_table,[7 13],'ascend');
Names = unique(gene_table.Name);clear gene_table
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
ON_RNAIDs = find(ismember(uniNames,extractBefore(settings.transcript_IDs_desired{:},'.')));
OFF_RNAIDs = setdiff(NonDNA_IDs,ON_RNAIDs);
ON_RNAIDs_Isos =  find(ismember(uniNames,extractBefore(settings.transcript_IDs_desired{:},'.')));
Desired_Isoforms =  find(ismember(uniNames,extractBefore(settings.transcript_IDs{:},'.')));
UnDesired_Isoforms = setdiff(ON_RNAIDs_Isos,Desired_Isoforms);
OFF_RNAIDs_minusIsos = setdiff(OFF_RNAIDs,UnDesired_Isoforms);
if (removeUndesiredIsos)
    OFF_RNAIDs = OFF_RNAIDs_minusIsos;
end
%Finds Off-targets and off-target binding sites
% finds list of each off-target both number, site and location, as well as Koff and Kon equilibrium constants
Tvec_RNA = cell(1,size(DoesProbeBindSite,1));
TSvec_RNA = cell(1,size(DoesProbeBindSite,1));
TPvec_RNA = cell(1,size(DoesProbeBindSite,1));
Svec_RNA = cell(1,size(DoesProbeBindSite,1));
Nvec_RNAsingle = zeros(1,size(DoesProbeBindSite,1));
Nvec_RNAmulti = zeros(1,size(DoesProbeBindSite,1));
NTPvec_RNAsingle = zeros(1,size(DoesProbeBindSite,1));
NTPvec_RNAmulti = zeros(1,size(DoesProbeBindSite,1));
Tvec_logKOFF_RNA = cell(1,size(DoesProbeBindSite,1));
Tvec_logKOFFdivON_RNA = cell(1,size(DoesProbeBindSite,1));
Tvec_logKONdivOFF_RNA = cell(1,size(DoesProbeBindSite,1));
TPvec_logKOFF_RNA = cell(1,length(OFF_RNAIDs));
TPvec_logKOFFdivON_RNA = cell(1,length(OFF_RNAIDs));
TPvec_logKONdivOFF_RNA = cell(1,length(OFF_RNAIDs));
Tvec_DNA = cell(1,size(DoesProbeBindSite,1));
TSvec_DNA = cell(1,size(DoesProbeBindSite,1));
TPvec_DNA = cell(1,size(DoesProbeBindSite,1));
Svec_DNA = cell(1,size(DoesProbeBindSite,1));
Nvec_DNAsingle = zeros(1,size(DoesProbeBindSite,1));
Nvec_DNAmulti = zeros(1,size(DoesProbeBindSite,1));
NTPvec_DNAsingle = zeros(1,size(DoesProbeBindSite,1));
NTPvec_DNAmulti = zeros(1,size(DoesProbeBindSite,1));
Tvec_logKOFF_DNA = cell(1,size(DoesProbeBindSite,1));
Tvec_logKOFFdivON_DNA = cell(1,size(DoesProbeBindSite,1));
TPvec_logKOFF_DNA = cell(1,length(DNA_IDs));
TPvec_logKOFFdivON_DNA = cell(1,length(DNA_IDs));
TPvec_logKONdivOFF_DNA = cell(1,length(DNA_IDs));
TPvec_logKOFFdivCOMP_DNA = cell(1,length(DNA_IDs));
TPvec_logKCOMPdivOFF_DNA = cell(1,length(DNA_IDs));
if (isRNA)
    if (~isfile([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_BasicStats_RNA_ByProbe_Info' settings.designerName '.mat']))%check if temp file exists
        N_Probes = size(DoesProbeBindSite,1);
        N_ProbeBatches = ceil(N_Probes/probeBatchSize);
        R = mod(N_Probes,probeBatchSize);
        probeBatch = cell(1,N_ProbeBatches);
        if (R==0)
            for k = 1:N_ProbeBatches
                probeBatch{k} = probeBatchSize*(k-1)+1:probeBatchSize*k;
            end
        else
            for k = 1:N_ProbeBatches-1
                probeBatch{k} = probeBatchSize*(k-1)+1:probeBatchSize*k;
            end
            probeBatch{N_ProbeBatches} = probeBatchSize*(N_ProbeBatches-1)+1:probeBatchSize*(N_ProbeBatches-1)+R;
        end
        ResultsExist = zeros(1,N_ProbeBatches);
        ResultsDate = cell(1,N_ProbeBatches);
        fprintf("Check if probe RNA off-target statistic files exist")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_ProbeBatches,'WaitMessage','Checking');
        parfor i = 1:N_ProbeBatches
            if (isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsRNA_Info_probebatch' num2str(i) '.mat']))%check if temp file exists
                d = dir([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsRNA_Info_probebatch' num2str(i) '.mat']);
                if (d.bytes>0)%check size greater than zero
                    ResultsExist(i) = 1;
                end
                ResultsDate{i} = datetime(d.date);
            end
            progress(wb);
        end
        wb.delete();
        fprintf('\n')
        fprintf('\n')
        Results_NotMade = find(ResultsExist==0);
        Results_Made = find(ResultsExist==1);
        %Sort get most recent ResultsMade GeneHitsMade and GeneHitsTable Made and add to probe_check_list
        if (length(Results_Made)<=most_recent_num)
            results_check1 = Results_Made;
        else
            Results_RecentMade_Dates(:,1) = ResultsDate(Results_Made);
            Results_RecentMade_Dates(:,2) = num2cell(Results_Made);
            Results_RecentMade_Dates = table2timetable(cell2table(Results_RecentMade_Dates));
            Results_RecentMade_Dates = sortrows(Results_RecentMade_Dates,1,'descend');
            Results_RecentMade_Dates.Properties.VariableNames = {'ID'};
            results_check1 = Results_RecentMade_Dates.ID(1:most_recent_num).';
            clear Results_RecentMade_Dates
        end
        batch_nums_to_check = union(Results_NotMade,results_check1);
        Kb_List = parallel.pool.Constant(Kb);
        Kon_List = parallel.pool.Constant(Kon);
        probeBatch_List = parallel.pool.Constant(probeBatch);
        OFF_RNAIDs_List = parallel.pool.Constant(OFF_RNAIDs);
        DoesProbeBindSite_List = parallel.pool.Constant(DoesProbeBindSite);
        fprintf(strcat("Number of Probe Batches: ",string(N_ProbeBatches)))
        fprintf('\n')
        fprintf('\n')
        fprintf("Batches to Compute Probe RNA Target Statistics: ")
        fprintf(num2str(batch_nums_to_check1))
        fprintf('\n')
        fprintf('\n')
        fprintf("Computing probe batch RNA off-target statistics")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(length(batch_nums_to_check),'WaitMessage','Computing');
        parfor v = 1:length(batch_nums_to_check)
            Js =  @(x) find(sum(squeeze(sum(DoesProbeBindSite_List.Value(x,:,:),1)),2)>0);
            Sx = @(x,Z) arrayfun(@(y) find(squeeze(DoesProbeBindSite_List.Value(x,y,:))==1),Z,'Un',0);
            Js_OFFRNA = @(x) OFF_RNAIDs_List.Value(ismember(OFF_RNAIDs_List.Value,Js(x)));
            Js_OFFRNAi = @(x,y) OFF_RNAIDs_List.Value(find(cumsum(ismember(OFF_RNAIDs_List.Value,Js(x)))==y,1));
            temp_designer_stats_rna_p = cell(1,length(probeBatch_List.Value{batch_nums_to_check(v)}));
            for w = 1:length(probeBatch_List.Value{batch_nums_to_check(v)})
                p = probeBatch_List.Value{batch_nums_to_check(v)}(w);
                temp_designer_stats_rna_p{w}{1} = length(Js_OFFRNA(p));
                temp_designer_stats_rna_p{w}{2} = sum(cellfun(@length,Sx(p,Js_OFFRNA(p))));
                temp_designer_stats_rna_p{w}{3} = cell2mat(Sx(p,Js_OFFRNA(p)));
                temp_designer_stats_rna_p{w}{4} = cell2mat(arrayfun(@(x) repmat(Js_OFFRNAi(p,x),[1 cellfun(@length,Sx(p,Js_OFFRNAi(p,x)))]),1:length(Js_OFFRNA(p)),'Un',0));
                temp_designer_stats_rna_p{w}{5} = log10(diag(full(squeeze(Kb_List.Value(p,temp_designer_stats_rna_p{w}{4},temp_designer_stats_rna_p{w}{3}))))');
                temp_designer_stats_rna_p{w}{6} = log10(diag(full(squeeze(Kb_List.Value(p,temp_designer_stats_rna_p{w}{4},temp_designer_stats_rna_p{w}{3}))))'/Kon_List.Value(p));
                temp_designer_stats_rna_p{w}{7} = log10(Kon_List.Value(p)./diag(full(squeeze(Kb_List.Value(p,temp_designer_stats_rna_p{w}{4},temp_designer_stats_rna_p{w}{3}))))');
            end
            parsave_partial_designer_stats([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsRNA_Info_probebatch' num2str(batch_nums_to_check(v)) '.mat'],temp_designer_stats_rna_p);
            progress(wb);
        end
        wb.delete();
        fprintf('\n')
        fprintf('\n')
        Nvec_RNAsingle = struct('Nvec_RNAsingle',cell(1,N_ProbeBatches));
        Nvec_RNAmulti = struct('Nvec_RNAmulti',cell(1,N_ProbeBatches));
        Svec_RNA = struct('Svec_RNA',cell(1,N_ProbeBatches));
        Tvec_RNA = struct('Tvec_RNA',cell(1,N_ProbeBatches));
        Tvec_logKOFF_RNA = struct('Tvec_logKOFF_RNA',cell(1,N_ProbeBatches));
        Tvec_logKOFFdivON_RNA =struct('Tvec_logKOFFdivON_RNA',cell(1,N_ProbeBatches));
        Tvec_logKONdivOFF_RNA = struct('Tvec_logKONdivOFF_RNA',cell(1,N_ProbeBatches));
        fprintf("Aggregating probe batch RNA off-target statistics")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_ProbeBatches+7,'WaitMessage','Aggregating');
        parfor v = 1:N_ProbeBatches
            if isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsRNA_Info_probebatch' num2str(v) '.mat'])
                partial_designer_stats_rna_p_tmp = load([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsRNA_Info_probebatch' num2str(v) '.mat']).partial_designer_stats_tmp;
                Nvec_RNAsingle(v).Nvec_RNAsingle= cell2mat(cellfun(@(x) x{1},partial_designer_stats_rna_p_tmp,'Un',0));
                Nvec_RNAmulti(v).Nvec_RNAmulti = cell2mat(cellfun(@(x) x{2},partial_designer_stats_rna_p_tmp,'Un',0));
                Svec_RNA(v).Svec_RNA = cellfun(@(x) x{3},partial_designer_stats_rna_p_tmp,'Un',0);
                Tvec_RNA(v).Tvec_RNA = cellfun(@(x) x{4},partial_designer_stats_rna_p_tmp,'Un',0);
                Tvec_logKOFF_RNA(v).Tvec_logKOFF_RNA = cellfun(@(x) x{5},partial_designer_stats_rna_p_tmp,'Un',0);
                Tvec_logKOFF_RNA(v).Tvec_logKOFF_RNA = cellfun(@(x) x{6},partial_designer_stats_rna_p_tmp,'Un',0);
                Tvec_logKONdivOFF_RNA(v).Tvec_logKONdivOFF_RNA= cellfun(@(x) x{7},partial_designer_stats_rna_p_tmp,'Un',0);
            end
            progress(wb);
        end
        Nvec_RNAsingle = vertcat(Nvec_RNAsingle(:).Nvec_RNAsingle);progress(wb);
        Nvec_RNAmulti = vertcat(Nvec_RNAmulti(:).Nvec_RNAmulti);progress(wb);
        Svec_RNA = vertcat(Svec_RNA(:).Svec_RNA);progress(wb);
        Tvec_RNA = vertcat(Tvec_RNA(:).Tvec_RNA);progress(wb);
        Tvec_logKOFF_RNA = vertcat(Tvec_logKOFF_RNA(:).Tvec_logKOFF_RNA);progress(wb);
        Tvec_logKOFFdivON_RNA = vertcat(Tvec_logKOFFdivON_RNA(:).Tvec_logKOFFdivON_RNA);progress(wb);
        Tvec_logKONdivOFF_RNA = vertcat(Tvec_logKONdivOFF_RNA(:).Tvec_logKONdivOFF_RNA);progress(wb);
        wb.delete();
        fprintf('\n')
        fprintf('\n')
        save([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_BasicStats_RNA_ByProbe_Info' settings.designerName '.mat'],...
            'Nvec_RNAsingle','Nvec_RNAmulti','Svec_RNA','Tvec_RNA','Tvec_logKOFF_RNA','Tvec_logKOFFdivON_RNA','Tvec_logKONdivOFF_RNA','-v7.3')
        fprintf('\n')
        fprintf('\n')
        fprintf("Deleting temporary probe batch RNA off-target statistics files")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_ProbeBatches,'WaitMessage', 'Deleting');
        parfor v = 1:N_ProbeBatches
            if exist([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsRNA_Info_probebatch' num2str(v) '.mat'],'file')        %delete temp mat file if already exists
                delete([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsRNA_Info_probebatch' num2str(v) '.mat']);
            end
            progress(wb);
        end
        wb.delete();
        fprintf('\n')
        fprintf('\n')
    else
        load([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_BasicStats_RNA_ByProbe_Info' settings.designerName '.mat'],...
            'Nvec_RNAsingle','Nvec_RNAmulti','Svec_RNA','Tvec_RNA','Tvec_logKOFF_RNA','Tvec_logKOFFdivON_RNA','Tvec_logKONdivOFF_RNA')
    end
    if (~isfile([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_BasicStats_RNA_ByTarget_Info' settings.designerName '.mat']))%check if temp file exists
        N_RNATargetBatches = ceil(length(OFF_RNAIDs)/targetBatchSize);
        R = mod(length(OFF_RNAIDs),targetBatchSize);
        rnaTargetBatches = cell(1,N_RNATargetBatches);
        if (R==0)
            for k = 1:N_RNATargetBatches
                rnaTargetBatches{k} = targetBatchSize*(k-1)+1:targetBatchSize*k;
            end
        else
            for k = 1:N_RNATargetBatches-1
                rnaTargetBatches{k} = targetBatchSize*(k-1)+1:targetBatchSize*k;
            end
            rnaTargetBatches{N_RNATargetBatches} = targetBatchSize*(N_RNATargetBatches-1)+1:targetBatchSize*(N_RNATargetBatches-1)+R;
        end
        ResultsExist2 = zeros(1,N_RNATargetBatches);
        ResultsDate2 = cell(1,N_RNATargetBatches);
        fprintf("Check if RNA off-target by probe statistic files exist")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_RNATargetBatches,'WaitMessage','Checking');
        parfor i = 1:N_RNATargetBatches
            if (isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsRNA_Info_targetbatch' num2str(i) '.mat']))%check if temp file exists
                d = dir([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsRNA_Info_targetbatch' num2str(i) '.mat']);
                if (d.bytes>0)%check size greater than zero
                    ResultsExist2(i) = 1;
                end
                ResultsDate2{i} = datetime(d.date);
            end
            progress(wb);
        end
        wb.delete();
        fprintf('\n')
        fprintf('\n')
        Results_NotMade2 = find(ResultsExist2==0);
        Results_Made2 = find(ResultsExist2==1);
        %Sort get most recent ResultsMade GeneHitsMade and GeneHitsTable Made and add to probe_check_list
        if (length(Results_Made2)<=most_recent_num)
            results_check2 = Results_Made2;
        else
            Results_RecentMade_Dates2(:,1) = ResultsDate2(Results_Made2);
            Results_RecentMade_Dates2(:,2) = num2cell(Results_Made2);
            Results_RecentMade_Dates2 = table2timetable(cell2table(Results_RecentMade_Dates2));
            Results_RecentMade_Dates2 = sortrows(Results_RecentMade_Dates2,1,'descend');
            Results_RecentMade_Dates2.Properties.VariableNames = {'ID'};
            results_check2 = Results_RecentMade_Dates2.ID(1:most_recent_num).';
            clear Results_RecentMade_Dates
        end
        batch_nums_to_check2 = union(Results_NotMade2,results_check2);
        Kb_List = parallel.pool.Constant(Kb);
        Kon_List = parallel.pool.Constant(Kon);
        rnaTargetBatches_List = parallel.pool.Constant(rnaTargetBatches);
        OFF_RNAIDs_List = parallel.pool.Constant(OFF_RNAIDs);
        DoesProbeBindSite_List = parallel.pool.Constant(DoesProbeBindSite);
        fprintf(strcat("Number of RNA Target Batches: ",string(N_RNATargetBatches)))
        fprintf('\n')
        fprintf('\n')
        fprintf("Batches to Compute RNA Target Statistics: ")
        fprintf(num2str(batch_nums_to_check2))
        fprintf('\n')
        fprintf('\n')
        fprintf("Computing batch RNA off-target by probe statistics")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(length(batch_nums_to_check2),'WaitMessage','Computing');
        parfor v = 1:length(batch_nums_to_check2)
            Tp = @(x) find(sum(squeeze(DoesProbeBindSite_List.Value(:,x,:)),2)>0);
            Tx =@(y,Z) arrayfun(@(x) find(squeeze(DoesProbeBindSite_List.Value(x,y,:))==1),Z,'Un',0);
            temp_designer_stats_rna_t = cell(1,length(rnaTargetBatches_List.Value{batch_nums_to_check2(v)}));
            for w = 1:length(rnaTargetBatches_List.Value{batch_nums_to_check2(v)})
                temp_RNAV0b = Tp(OFF_RNAIDs_List.Value(rnaTargetBatches_List.Value{batch_nums_to_check2(v)}(w)));%does not have multiplicity
                temp_RNAV1b = Tx(OFF_RNAIDs_List.Value(rnaTargetBatches_List.Value{batch_nums_to_check2(v)}(w)),temp_RNAV0b);
                temp_RNAV2b = cellfun(@length,temp_RNAV1b);
                temp_designer_stats_rna_t{w}{1} = length(temp_RNAV0b);
                temp_designer_stats_rna_t{w}{2} = sum(temp_RNAV2b);
                temp_designer_stats_rna_t{w}{3} = cell2mat(arrayfun(@(x) repmat(temp_RNAV0b(x),[1 temp_RNAV2b(x)]),1:length(temp_RNAV0b),'Un',0));
                temp_designer_stats_rna_t{w}{4} = cell2mat(temp_RNAV1b);%site locations
                temp_designer_stats_rna_t{w}{5} = log10(diag(full(squeeze(Kb_List.Value(temp_designer_stats_rna_t{w}{3},OFF_RNAIDs_List.Value(rnaTargetBatches_List.Value{batch_nums_to_check2(v)}(w)),temp_designer_stats_rna_t{w}{4}))))');
                temp_designer_stats_rna_t{w}{6} = log10(diag(full(squeeze(Kb_List.Value(temp_designer_stats_rna_t{w}{3},OFF_RNAIDs_List.Value(rnaTargetBatches_List.Value{batch_nums_to_check2(v)}(w)),temp_designer_stats_rna_t{w}{4}))))'./Kon_List.Value(temp_designer_stats_rna_t{w}{3}));
                temp_designer_stats_rna_t{w}{7} = log10(Kon_List.Value(temp_designer_stats_rna_t{w}{3})./diag(full(squeeze(Kb_List.Value(temp_designer_stats_rna_t{w}{3},OFF_RNAIDs_List.Value(rnaTargetBatches_List.Value{batch_nums_to_check2(v)}(w)),temp_designer_stats_rna_t{w}{4}))))');
            end
            parsave_partial_designer_stats([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsRNA_Info_targetbatch' num2str(batch_nums_to_check2(v)) '.mat'],temp_designer_stats_rna_t);
            progress(wb);
        end
        wb.delete();
        NTPvec_RNAsingle = struct('NTPvec_RNAsingle',cell(1,N_RNATargetBatches));
        NTPvec_RNAmulti = struct('NTPvec_RNAmulti',cell(1,N_RNATargetBatches));
        TPvec_RNA = struct('TPvec_RNA',cell(1,N_RNATargetBatches));
        TSvec_RNA = struct('TSvec_RNA',cell(1,N_RNATargetBatches));
        TPvec_logKOFF_RNA = struct('TPvec_logKOFF_RNA',cell(1,N_RNATargetBatches));
        TPvec_logKOFFdivON_RNA =struct('TPvec_logKOFFdivON_RNA',cell(1,N_RNATargetBatches));
        TPvec_logKONdivOFF_RNA = struct('TPvec_logKONdivOFF_RNA',cell(1,N_RNATargetBatches));
        fprintf("Aggregating batch RNA off-target by target statistics")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_RNATargetBatches+7,'WaitMessage','Aggregating');
        parfor v = 1:N_RNATargetBatches
            if isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsRNA_Info_targetbatch' num2str(v) '.mat'])
                partial_designer_stats_rna_t_tmp = load([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsRNA_Info_targetbatch' num2str(v) '.mat']).partial_designer_stats_tmp;
                NTPvec_RNAsingle(v).NTPvec_RNAsingle = cell2mat(cellfun(@(x) x{1},partial_designer_stats_rna_t_tmp,'Un',0));
                NTPvec_RNAmulti(v).NTPvec_RNAmulti = cell2mat(cellfun(@(x) x{2},partial_designer_stats_rna_t_tmp,'Un',0));
                TPvec_RNA(v).TPvec_RNA = cellfun(@(x) x{3},partial_designer_stats_rna_t_tmp,'Un',0);
                TSvec_RNA(v).TSvec_RNA = cellfun(@(x) x{4},partial_designer_stats_rna_t_tmp,'Un',0);
                TPvec_logKOFF_RNA(v).TPvec_logKOFF_RNA = cellfun(@(x) x{5},partial_designer_stats_rna_t_tmp,'Un',0);
                TPvec_logKOFFdivON_RNA(v).TPvec_logKOFFdivON_RNA = cellfun(@(x) x{6},partial_designer_stats_rna_t_tmp,'Un',0);
                TPvec_logKONdivOFF_RNA(v).TPvec_logKONdivOFF_RNA= cellfun(@(x) x{7},partial_designer_stats_rna_t_tmp,'Un',0);
            end
            progress(wb);
        end
        NTPvec_RNAsingle = vertcat(NTPvec_RNAsingle(:).NTPvec_RNAsingle);progress(wb);
        NTPvec_RNAmulti = vertcat(NTPvec_RNAmulti(:).NTPvec_RNAmulti);progress(wb);
        TPvec_RNA = vertcat(TPvec_RNA(:).TPvec_RNA);progress(wb);
        TSvec_RNA = vertcat(TSvec_RNA(:).TSvec_RNA);progress(wb);
        TPvec_logKOFF_RNA = vertcat(TPvec_logKOFF_RNA(:).TPvec_logKOFF_RNA);progress(wb);
        TPvec_logKOFFdivON_RNA = vertcat(TPvec_logKOFFdivON_RNA(:).TPvec_logKOFFdivON_RNA);progress(wb);
        TPvec_logKONdivOFF_RNA = vertcat(TPvec_logKONdivOFF_RNA(:).TPvec_logKONdivOFF_RNA);progress(wb);
        wb.delete();
        fprintf('\n')
        fprintf('\n')
        save([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_BasicStats_RNA_ByTarget_Info' settings.designerName '.mat'],...
            'NTPvec_RNAsingle','NTPvec_RNAmulti','TSvec_RNA','TPvec_RNA','TPvec_logKOFF_RNA','TPvec_logKOFFdivON_RNA','TPvec_logKONdivOFF_RNA','-v7.3')
        fprintf('\n')
        fprintf('\n')
        fprintf("Deleting temporary target batch RNA off-target statistics files")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_RNATargetBatches,'WaitMessage', 'Deleting');
        parfor v = 1:N_RNATargetBatches
            if exist([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsRNA_Info_targetbatch' num2str(v) '.mat'],'file')        %delete temp mat file if already exists
                delete([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsRNA_Info_targetbatch' num2str(v) '.mat']);
            end
            progress(wb);
        end
        wb.delete();
        fprintf('\n')
        fprintf('\n')
    else
        load([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_BasicStats_RNA_ByTarget_Info' settings.designerName '.mat'],...
            'NTPvec_RNAsingle','NTPvec_RNAmulti','TSvec_RNA','TPvec_RNA','TPvec_logKOFF_RNA','TPvec_logKOFFdivON_RNA','TPvec_logKONdivOFF_RNA')
    end
    Cout{1}{1} = Tvec_RNA;
    Cout{1}{2} = Svec_RNA;
    Cout{1}{3} = TPvec_RNA;
    Cout{1}{4} = TSvec_RNA;
    Cout{1}{5} = TPvec_logKOFF_RNA;
    Cout{1}{6} = TPvec_logKOFFdivON_RNA;
    Cout{1}{7} = TPvec_logKONdivOFF_RNA;
    Probes_With_RNAOFF = AllowableProbes(Nvec_RNAmulti(AllowableProbes)>0);
    NumRNAOffTargetOptions = unique(Nvec_RNAmulti(AllowableProbes));
    Probes_WithNRNAOFF = arrayfun(@(i) AllowableProbes(Nvec_RNAmulti(AllowableProbes)==NumRNAOffTargetOptions(i)),1:length(NumRNAOffTargetOptions),'Un',0);
    fprintf("Converting RNA Statistics into Probe Specificity and OFF-target Scores")
    fprintf('\n')
    fprintf('\n')
    wb = parwaitbar(length(Probes_With_RNAOFF),'WaitMessage','Converting');
    for v = 1:length(Probes_With_RNAOFF)
        temp_T = Tvec_RNA{Probes_With_RNAOFF(v)};
        temp_KOFF = Tvec_logKOFF_RNA{Probes_With_RNAOFF(v)};
        temp_KOFFdivON = Tvec_logKOFFdivON_RNA{Probes_With_RNAOFF(v)};
        RNASpecificity_Score(Probes_With_RNAOFF(v)) = dot(EKernel(temp_T)',temp_KOFFdivON);
        RNAOFF_Score(Probes_With_RNAOFF(v)) = dot(EKernel(temp_T)',temp_KOFF);
        progress(wb);
    end
    wb.delete();
    fprintf('\n')
    fprintf('\n')
end
if (isDNA)
    if (~isfile([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_BasicStats_DNA_ByProbe_Info' settings.designerName '.mat']))%check if temp file exists
        N_Probes = size(DoesProbeBindSite,1);
        N_ProbeBatches = ceil(N_Probes/probeBatchSize);
        R = mod(N_Probes,probeBatchSize);
        probeBatch = cell(1,N_ProbeBatches);
        if (R==0)
            for k = 1:N_ProbeBatches
                probeBatch{k} = probeBatchSize*(k-1)+1:probeBatchSize*k;
            end
        else
            for k = 1:N_ProbeBatches-1
                probeBatch{k} = probeBatchSize*(k-1)+1:probeBatchSize*k;
            end
            probeBatch{N_ProbeBatches} = probeBatchSize*(N_ProbeBatches-1)+1:probeBatchSize*(N_ProbeBatches-1)+R;
        end
        ResultsExist3 = zeros(1,N_ProbeBatches);
        ResultsDate3 = cell(1,N_ProbeBatches);
        fprintf("Check if probe DNA off-target statistic files exist")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_ProbeBatches,'WaitMessage','Checking');
        parfor i = 1:N_ProbeBatches
            if (isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsDNA_Info_probebatch' num2str(i) '.mat']))%check if temp file exists
                d = dir([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsDNA_Info_probebatch' num2str(i) '.mat']);
                if (d.bytes>0)%check size greater than zero
                    ResultsExist3(i) = 1;
                end
                ResultsDate3{i} = datetime(d.date);
            end
            progress(wb);
        end
        wb.delete();
        fprintf('\n')
        fprintf('\n')
        Results_NotMade3 = find(ResultsExist3==0);
        Results_Made3 = find(ResultsExist3==1);
        %Sort get most recent ResultsMade GeneHitsMade and GeneHitsTable Made and add to probe_check_list
        if (length(Results_Made3)<=most_recent_num)
            results_check3 = Results_Made3;
        else
            Results_RecentMade_Dates3(:,1) = ResultsDate3(Results_Made3);
            Results_RecentMade_Dates3(:,2) = num2cell(Results_Made3);
            Results_RecentMade_Dates3 = table2timetable(cell2table(Results_RecentMade_Dates3));
            Results_RecentMade_Dates3 = sortrows(Results_RecentMade_Dates3,1,'descend');
            Results_RecentMade_Dates3.Properties.VariableNames = {'ID'};
            results_check3 = Results_RecentMade_Dates3.ID(1:most_recent_num).';
            clear Results_RecentMade_Dates3
        end
        batch_nums_to_check3 = union(Results_NotMade3,results_check3);
        Kb_List = parallel.pool.Constant(Kb);
        Kon_List = parallel.pool.Constant(Kon);
        Kb_Complement_List = parallel.pool.Constant(Kb_Complement);
        probeBatch_List = parallel.pool.Constant(probeBatch);
        DNA_IDs_List = parallel.pool.Constant(DNA_IDs);
        DoesProbeBindSite_List = parallel.pool.Constant(DoesProbeBindSite);
        fprintf(strcat("Number of Probe Batches: ",string(N_ProbeBatches)))
        fprintf('\n')
        fprintf('\n')
        fprintf("Batches to Compute Probe DNA Target Statistics: ")
        fprintf(num2str(batch_nums_to_check3))
        fprintf('\n')
        fprintf('\n')
        fprintf("Computing probe batch DNA off-target statistics")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(length(batch_nums_to_check3),'WaitMessage','Computing');
        parfor v = 1:length(batch_nums_to_check3)
            Js = @(x) find(sum(squeeze(sum(DoesProbeBindSite_List.Value(x,:,:),1)),2)>0);
            Sx =@(x,Z) arrayfun(@(y) find(squeeze(DoesProbeBindSite_List.Value(x,y,:))==1),Z,'Un',0);
            Js_OFFDNA = @(x)DNA_IDs_List.Value(ismember(DNA_IDs_List.Value,Js(x)));
            Js_OFFDNAi = @(x,y)DNA_IDs_List.Value(find(cumsum(ismember(DNA_IDs_List.Value,Js(x)))==y,1));
            temp_designer_stats_dna_p = cell(1,length(probeBatch_List.Value{batch_nums_to_check3(v)}));
            for w = 1:length(probeBatch_List.Value{batch_nums_to_check3(v)})
                p = probeBatch_List.Value{batch_nums_to_check3(v)}(w);
                temp_designer_stats_dna_p{w}{1} = length(Js_OFFDNA(p));
                temp_designer_stats_dna_p{w}{2} = sum(cellfun(@length,Sx(p,Js_OFFDNA(p))));
                temp_designer_stats_dna_p{w}{3} = cell2mat(Sx(p,Js_OFFDNA(p)));
                temp_designer_stats_dna_p{w}{4} = cell2mat(arrayfun(@(x) repmat(Js_OFFDNAi(p,x),[1 cellfun(@length,Sx(p,Js_OFFDNAi(p,x)))]),1:length(Js_OFFDNA(p)),'Un',0));
                temp_designer_stats_dna_p{w}{5} = log10(diag(full(squeeze(Kb_List.Value(p,temp_designer_stats_dna_p{4},temp_designer_stats_dna_p{3}))))');
                temp_designer_stats_dna_p{w}{6} = log10(diag(full(squeeze(Kb_List.Value(p,temp_designer_stats_dna_p{4},temp_designer_stats_dna_p{3}))))'/Kon_List.Value(p));
                temp_designer_stats_dna_p{w}{7} = log10(Kon_List.Value(p)./diag(full(squeeze(Kb_List.Value(p,temp_designer_stats_dna_p{4},temp_designer_stats_dna_p{w}{3}))))');
                temp_designer_stats_dna_p{w}{8} = log10(diag(full(squeeze(Kb_List.Value(p,temp_designer_stats_dna_p{w}{4},temp_designer_stats_dna_p{w}{3}))))'./diag(full(squeeze(Kb_Complement_List.Value(temp_designer_stats_dna_p{w}{4},temp_designer_stats_dna_p{w}{3}))))');
                temp_designer_stats_dna_p{w}{9} = log10(diag(full(squeeze(Kb_Complement_List.Value(temp_designer_stats_dna_p{w}{4},temp_designer_stats_dna_p{w}{3}))))'./diag(full(squeeze(Kb_List.Value(p,temp_designer_stats_dna_p{w}{4},temp_designer_stats_dna_p{w}{3}))))');
            end
            parsave_partial_designer_stats([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsDNA_Info_probebatch' num2str(batch_nums_to_check3(v)) '.mat'],temp_designer_stats_dna_p);
            progress(wb);
        end
        wb.delete();
        fprintf('\n')
        fprintf('\n')
        Nvec_DNAsingle = struct('Nvec_DNAsingle',cell(1,N_ProbeBatches));
        Nvec_DNAmulti = struct('Nvec_DNAmulti',cell(1,N_ProbeBatches));
        Svec_DNA = struct('Svec_DNA',cell(1,N_ProbeBatches));
        Tvec_DNA = struct('Tvec_DNA',cell(1,N_ProbeBatches));
        Tvec_logKOFF_DNA = struct('Tvec_logKOFF_DNA',cell(1,N_ProbeBatches));
        Tvec_logKOFFdivON_DNA =struct('Tvec_logKOFFdivON_DNA',cell(1,N_ProbeBatches));
        Tvec_logKONdivOFF_DNA = struct('Tvec_logKONdivOFF_DNA',cell(1,N_ProbeBatches));
        Tvec_logKOFFdivCOMP_DNA = struct('Tvec_logKOFFdivCOMP_DNA',cell(1,N_ProbeBatches));
        Tvec_logKCOMPdivOFF_DNA = struct('Tvec_logKCOMPdivOFF_DNA',cell(1,N_ProbeBatches));
        fprintf("Aggregating probe batch DNA off-target statistics")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_ProbeBatches+9,'WaitMessage','Aggregating');
        for v = 1:N_ProbeBatches
            if isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsDNA_Info_probebatch' num2str(v) '.mat'])
                partial_designer_stats_dna_p_tmp = load([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsDNA_Info_probebatch' num2str(v) '.mat']).partial_designer_stats_tmp;
                Nvec_DNAsingle(v).Nvec_DNAsingle = cell2mat(cellfun(@(x) x{1},partial_designer_stats_dna_p_tmp,'Un',0));
                Nvec_DNAmulti(v).Nvec_DNAmulti = cell2mat(cellfun(@(x) x{2},partial_designer_stats_dna_p_tmp,'Un',0));
                Svec_DNA(v).Svec_DNA = cellfun(@(x) x{3},partial_designer_stats_dna_p_tmp,'Un',0);
                Tvec_DNA(v).Tvec_DNA = cellfun(@(x) x{4},partial_designer_stats_dna_p_tmp,'Un',0);
                Tvec_logKOFF_DNA(v).Tvec_logKOFF_DNA = cellfun(@(x) x{5},partial_designer_stats_dna_p_tmp,'Un',0);
                Tvec_logKOFFdivON_DNA(v).Tvec_logKOFFdivON_DNA = cellfun(@(x) x{6},partial_designer_stats_dna_p_tmp,'Un',0);
                Tvec_logKONdivOFF_DNA(v).Tvec_logKONdivOFF_DNA = cellfun(@(x) x{7},partial_designer_stats_dna_p_tmp,'Un',0);
                Tvec_logKOFFdivCOMP_DNA(v).Tvec_logKOFFdivCOMP_DNA = cellfun(@(x) x{8},partial_designer_stats_dna_p_tmp,'Un',0);
                Tvec_logKCOMPdivOFF_DNA(v).Tvec_logKCOMPdivOFF_DNA = cellfun(@(x) x{9},partial_designer_stats_dna_p_tmp,'Un',0);
            end
            progress(wb);
        end
        Nvec_DNAsingle = vertcat(Nvec_DNAsingle(:).Nvec_DNAsingle);progress(wb);
        Nvec_DNAmulti = vertcat(Nvec_DNAmulti(:).Nvec_DNAmulti);progress(wb);
        Svec_DNA = vertcat(Svec_DNA(:).Svec_DNA);progress(wb);
        Tvec_DNA = vertcat(Tvec_DNA(:).Tvec_DNA);progress(wb);
        Tvec_logKOFF_DNA = vertcat(Tvec_logKOFF_DNA(:).Tvec_logKOFF_DNA);progress(wb);
        Tvec_logKOFFdivON_DNA = vertcat(Tvec_logKOFFdivON_DNA(:).Tvec_logKOFFdivON_DNA);progress(wb);
        Tvec_logKONdivOFF_DNA = vertcat(Tvec_logKONdivOFF_DNA(:).Tvec_logKONdivOFF_DNA);progress(wb);
        Tvec_logKOFFdivCOMP_DNA = vertcat(Tvec_logKOFFdivCOMP_DNA(:).Tvec_logKOFFdivCOMP_DNA);progress(wb);
        Tvec_logKCOMPdivOFF_DNA = vertcat(Tvec_logKCOMPdivOFF_DNA(:).Tvec_logKCOMPdivOFF_DNA);progress(wb);
        wb.delete();
        fprintf('\n')
        fprintf('\n')
        save([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_BasicStats_DNA_ByProbe_Info' settings.designerName '.mat'],...
            'Nvec_DNAsingle','Nvec_DNAmulti','Svec_DNA','Tvec_DNA','Tvec_logKOFF_DNA','Tvec_logKOFFdivON_DNA','Tvec_logKONdivOFF_DNA','Tvec_logKOFFdivCOMP_DNA','Tvec_logKCOMPdivOFF_DNA','-v7.3')
        fprintf('\n')
        fprintf('\n')
        fprintf("Deleting temporary probe batch RNA off-target statistics files")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_ProbeBatches,'WaitMessage', 'Deleting');
        parfor v = 1:N_ProbeBatches
            if exist([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsDNA_Info_probebatch' num2str(v) '.mat'],'file')        %delete temp mat file if already exists
                delete([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsDNA_Info_probebatch' num2str(v) '.mat'])
            end
            progress(wb);
        end
        wb.delete();
        fprintf('\n')
        fprintf('\n')
    else
        load([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_BasicStats_DNA_ByProbe_Info' settings.designerName '.mat'],...
            'Nvec_DNAsingle','Nvec_DNAmulti','Svec_DNA','Tvec_DNA','Tvec_logKOFF_DNA','Tvec_logKOFFdivON_DNA','Tvec_logKONdivOFF_DNA','Tvec_logKOFFdivCOMP_DNA','Tvec_logKCOMPdivOFF_DNA')
    end
    if (~isfile([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_BasicStats_DNA_ByTarget_Info' settings.designerName '.mat']))%check if temp file exists
        N_DNATargetBatches = ceil(length(DNA_IDs)/targetBatchSize);
        R = mod(length(DNA_IDs),targetBatchSize);
        dnaTargetBatches = cell(1,N_DNATargetBatches);
        if (R==0)
            for k = 1:N_DNATargetBatches
                dnaTargetBatches{k} = targetBatchSize*(k-1)+1:targetBatchSize*k;
            end
        else
            for k = 1:N_DNATargetBatches-1
                dnaTargetBatches{k} = targetBatchSize*(k-1)+1:targetBatchSize*k;
            end
            dnaTargetBatches{N_DNATargetBatches} = targetBatchSize*(N_DNATargetBatches-1)+1:targetBatchSize*(N_DNATargetBatches-1)+R;
        end
        ResultsExist4 = zeros(1,N_DNATargetBatches);
        ResultsDate4 = cell(1,N_DNATargetBatches);
        fprintf("Check if DNA off-target by probe statistic files exist")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_DNATargetBatches,'WaitMessage','Checking');
        parfor i = 1:N_DNATargetBatches
            if (isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsDNA_Info_targetbatch' num2str(i) '.mat']))%check if temp file exists
                d = dir([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsDNA_Info_targetbatch' num2str(i) '.mat']);
                if (d.bytes>0)%check size greater than zero
                    ResultsExist4(i) = 1;
                end
                ResultsDate4{i} = datetime(d.date);
            end
            progress(wb);
        end
        wb.delete();
        fprintf('\n')
        fprintf('\n')
        Results_NotMade4 = find(ResultsExist4==0);
        Results_Made4 = find(ResultsExist4==1);
        %Sort get most recent ResultsMade GeneHitsMade and GeneHitsTable Made and add to probe_check_list
        if (length(Results_Made4)<=most_recent_num)
            results_check4 = Results_Made4;
        else
            Results_RecentMade_Dates4(:,1) = ResultsDate4(Results_Made4);
            Results_RecentMade_Dates4(:,2) = num2cell(Results_Made4);
            Results_RecentMade_Dates4 = table2timetable(cell2table(Results_RecentMade_Dates4));
            Results_RecentMade_Dates4 = sortrows(Results_RecentMade_Dates4,1,'descend');
            Results_RecentMade_Dates4.Properties.VariableNames = {'ID'};
            results_check4 = Results_RecentMade_Dates4.ID(1:most_recent_num).';
            clear Results_RecentMade_Dates
        end
        batch_nums_to_check4 = union(Results_NotMade4,results_check4);
        Kb_List = parallel.pool.Constant(Kb);
        Kon_List = parallel.pool.Constant(Kon);
        dnaTargetBatches_List = parallel.pool.Constant(dnaTargetBatches);
        DNA_IDs_List = parallel.pool.Constant(DNA_IDs);
        DoesProbeBindSite_List = parallel.pool.Constant(DoesProbeBindSite);
        fprintf(strcat("Number of DNA Target Batches: ",string(N_DNATargetBatches)))
        fprintf('\n')
        fprintf('\n')
        fprintf("Batches to Compute DNA Target Statistics: ")
        fprintf(num2str(batch_nums_to_check4))
        fprintf('\n')
        fprintf('\n')
        fprintf("Computing batch DNA off-target by probe statistics")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(length(batch_nums_to_check4),'WaitMessage','Computing');
        parfor v = 1:length(batch_nums_to_check4)
            Tp =  @(x) find(sum(squeeze(DoesProbeBindSite_List.Value(:,x,:)),2)>0);
            Tx2 = @(y,Z) arrayfun(@(x) find(squeeze(DoesProbeBindSite_List.Value(x,y,:))==1)',Z,'Un',0);
            temp_designer_stats_dna_t = cell(1,length(dnaTargetBatches_List.Value{batch_nums_to_check4(v)}));
            for w = 1:length(dnaTargetBatches_List.Value{batch_nums_to_check4(v)})
                temp_DNAV0b = Tp(DNA_IDs_List.Value(dnaTargetBatches_List.Value{batch_nums_to_check4(v)}(w)))';%does not have multiplicity
                temp_DNAV1b = Tx2(DNA_IDs_List.Value(dnaTargetBatches_List.Value{batch_nums_to_check4(v)}(w)),temp_DNAV0b);
                temp_DNAV2b = cellfun(@length,temp_DNAV1b);
                temp_designer_stats_dna_t{w}{1} = length(temp_DNAV0b);
                temp_designer_stats_dna_t{w}{2} = sum(temp_DNAV2b);
                temp_designer_stats_dna_t{w}{3} = cell2mat(arrayfun(@(x) repmat(temp_DNAV0b(x),[1 temp_DNAV2b(x)]),1:length(temp_DNAV0b),'Un',0));
                temp_designer_stats_dna_t{w}{4} = cell2mat(temp_DNAV1b);%site locations
                temp_designer_stats_dna_t{w}{5} = log10(diag(full(squeeze(Kb_List.Value(temp_designer_stats_dna_t{w}{3},DNA_IDs_List.Value(dnaTargetBatches_List.Value{batch_nums_to_check4(v)}(w)),temp_designer_stats_dna_t{w}{4}))))');
                temp_designer_stats_dna_t{w}{6} = log10(diag(full(squeeze(Kb_List.Value(temp_designer_stats_dna_t{w}{3},DNA_IDs_List.Value(dnaTargetBatches_List.Value{batch_nums_to_check4(v)}(w)),temp_designer_stats_dna_t{w}{4}))))'./Kon_List.Value(temp_designer_stats_dna_t{w}{3}));
                temp_designer_stats_dna_t{w}{7} = log10(Kon_List.Value(temp_designer_stats_dna_t{w}{3})./diag(full(squeeze(Kb_List.Value(temp_designer_stats_dna_t{w}{3},DNA_IDs_List.Value(dnaTargetBatches_List.Value{batch_nums_to_check4(v)}(w)),temp_designer_stats_dna_t{w}{4}))))');
                temp_designer_stats_dna_t{w}{8} = log10(diag(full(squeeze(Kb_List.Value(temp_designer_stats_dna_t{w}{3},DNA_IDs_List.Value(dnaTargetBatches_List.Value{batch_nums_to_check4(v)}(w)),temp_designer_stats_dna_t{w}{4}))))'./diag(full(squeeze(Kb_Complement_List.Value(DNA_IDs_List.Value(dnaTargetBatches_List.Value{batch_nums_to_check4(v)}(w)),temp_designer_stats_dna_t{w}{4}))))');
                temp_designer_stats_dna_t{w}{9} = log10(diag(full(squeeze(Kb_Complement_List.Value(DNA_IDs_List.Value(dnaTargetBatches_List.Value{batch_nums_to_check4(v)}(w)),temp_designer_stats_dna_t{w}{4}))))'./diag(full(squeeze(Kb_List.Value(temp_designer_stats_dna_t{w}{3},DNA_IDs_List.Value(dnaTargetBatches_List.Value{batch_nums_to_check4(v)}(w)),temp_designer_stats_dna_t{w}{4}))))');
            end
            parsave_partial_designer_stats([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsDNA_Info_targetbatch' num2str(batch_nums_to_check4(v)) '.mat'],temp_designer_stats_dna_t);
            progress(wb);
        end
        wb.delete();
        fprintf('\n')
        fprintf('\n')
        NTPvec_DNAsingle = struct('NTPvec_DNAsingle',cell(1,N_DNATargetBatches));
        NTPvec_DNAmulti = struct('NTPvec_DNAmulti',cell(1,N_DNATargetBatches));
        TPvec_DNA = struct('TPvec_DNA',cell(1,N_DNATargetBatches));
        TSvec_DNA = struct('TSvec_DNA',cell(1,N_DNATargetBatches));
        TPvec_logKOFF_DNA = struct('TPvec_logKOFF_DNA',cell(1,N_DNATargetBatches));
        TPvec_logKOFFdivON_DNA =struct('TPvec_logKOFFdivON_DNA',cell(1,N_DNATargetBatches));
        TPvec_logKONdivOFF_DNA = struct('TPvec_logKONdivOFF_DNA',cell(1,N_DNATargetBatches));
        TPvec_logKOFFdivCOMP_DNA = struct('TPvec_logKOFFdivCOMP_DNA',cell(1,N_DNATargetBatches));
        TPvec_logKCOMPdivOFF_DNA = struct('TPvec_logKCOMPdivOFF_DNA',cell(1,N_DNATargetBatches));
        fprintf("Aggregating batch RNA off-target by target statistics")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_DNATargetBatches+9,'WaitMessage','Aggregating');
        parfor v = 1:N_DNATargetBatches
            if isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsDNA_Info_targetbatch' num2str(v) '.mat'])
                partial_designer_stats_dna_t_tmp = load([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsDNA_Info_targetbatch' num2str(v) '.mat']).partial_designer_stats_tmp;
                NTPvec_DNAsingle(v).NTPvec_DNAsingle = cell2mat(cellfun(@(x) x{1},partial_designer_stats_dna_t_tmp,'Un',0));
                NTPvec_DNAmulti(v).NTPvec_DNAmulti = cell2mat(cellfun(@(x) x{2},partial_designer_stats_dna_t_tmp,'Un',0));
                TPvec_DNA(v).TPvec_DNA = cellfun(@(x) x{3},partial_designer_stats_dna_t_tmp,'Un',0);
                TSvec_DNA(v).TSvec_DNA = cellfun(@(x) x{4},partial_designer_stats_dna_t_tmp,'Un',0);
                TPvec_logKOFF_DNA(v).TPvec_logKOFF_DNA = cellfun(@(x) x{5},partial_designer_stats_dna_t_tmp,'Un',0);
                TPvec_logKOFFdivON_DNA(v).TPvec_logKOFFdivON_DNA = cellfun(@(x) x{6},partial_designer_stats_dna_t_tmp,'Un',0);
                TPvec_logKONdivOFF_DNA(v).TPvec_logKONdivOFF_DNA = cellfun(@(x) x{7},partial_designer_stats_dna_t_tmp,'Un',0);
                TPvec_logKOFFdivCOMP_DNA(v).TPvec_logKOFFdivCOMP_DNA = cellfun(@(x) x{8},partial_designer_stats_dna_t_tmp,'Un',0);
                TPvec_logKCOMPdivOFF_DNA(v).TPvec_logKCOMPdivOFF_DNA= cellfun(@(x) x{9},partial_designer_stats_dna_t_tmp,'Un',0);
            end
            progress(wb);
        end
        NTPvec_DNAsingle = vertcat(NTPvec_DNAsingle(:).NTPvec_DNAsingle);progress(wb);
        NTPvec_DNAmulti = vertcat(NTPvec_DNAmulti(:).NTPvec_DNAmulti);progress(wb);
        TSvec_DNA = vertcat(TSvec_DNA(:).TSvec_DNA);progress(wb);
        TPvec_DNA = vertcat(TPvec_DNA(:).TPvec_DNA);progress(wb);
        TPvec_logKOFF_DNA = vertcat(TPvec_logKOFF_DNA(:).TPvec_logKOFF_DNA);progress(wb);
        TPvec_logKOFFdivON_DNA = vertcat(TPvec_logKOFFdivON_DNA(:).TPvec_logKOFFdivON_DNA);progress(wb);
        TPvec_logKONdivOFF_DNA = vertcat(TPvec_logKONdivOFF_DNA(:).TPvec_logKONdivOFF_DNA);progress(wb);
        TPvec_logKOFFdivCOMP_DNA = vertcat(TPvec_logKOFFdivCOMP_DNA(:).TPvec_logKOFFdivCOMP_DNA);progress(wb);
        TPvec_logKCOMPdivOFF_DNA = vertcat(TPvec_logKCOMPdivOFF_DNA(:).TPvec_logKCOMPdivOFF_DNA);progress(wb);
        wb.delete();
        fprintf('\n')
        fprintf('\n')
        save([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_BasicStats_DNA_ByTarget_Info' settings.designerName '.mat'],...
            'NTPvec_DNAsingle','NTPvec_DNAmulti','TSvec_DNA','TPvec_DNA','TPvec_logKOFF_DNA','TPvec_logKOFFdivON_DNA','TPvec_logKONdivOFF_DNA','TPvec_logKOFFdivCOMP_DNA','TPvec_logKCOMPdivOFF_DNA','-v7.3')
        fprintf('\n')
        fprintf('\n')
        fprintf("Deleting temporary target batch RNA off-target statistics files")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_DNATargetBatches,'WaitMessage', 'Deleting');
        parfor v = 1:N_DNATargetBatches
            if exist([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsDNA_Info_targetbatch' num2str(v) '.mat'],'file')        %delete temp mat file if already exists
                delete([FolderRootName filesep '(' TranscriptName ')' designerName '_basicStatsDNA_Info_targetbatch' num2str(v) '.mat']);
            end
            progress(wb);
        end
        wb.delete();
        fprintf('\n')
        fprintf('\n')
    else
        load([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_BasicStats_DNA_ByTarget_Info' settings.designerName '.mat'],...
            'NTPvec_DNAsingle','NTPvec_DNAmulti','TSvec_DNA','TPvec_DNA','TPvec_logKOFF_DNA','TPvec_logKOFFdivON_DNA','TPvec_logKONdivOFF_DNA','TPvec_logKOFFdivCOMP_DNA','TPvec_logKCOMPdivOFF_DNA')
    end
    Cout{2}{1} = Tvec_DNA;
    Cout{2}{2} = Svec_DNA;
    Cout{2}{3} = TPvec_DNA;
    Cout{2}{4} = TSvec_DNA;
    Cout{2}{5} = TPvec_logKOFF_DNA;
    Cout{2}{6} = TPvec_logKOFFdivON_DNA;
    Cout{2}{7} = TPvec_logKONdivOFF_DNA;
    Cout{2}{8} = TPvec_logKOFFdivCOMP_DNA;
    Cout{2}{9} = TPvec_logKCOMPdivOFF_DNA;
    Probes_With_DNAOFF = AllowableProbes(Nvec_DNAmulti(AllowableProbes)>0);
    NumDNAOffTargetOptions = unique(Nvec_DNAmulti(AllowableProbes));
    Probes_WithNDNAOFF = arrayfun(@(i) AllowableProbes(Nvec_DNAmulti(AllowableProbes)==NumDNAOffTargetOptions(i)),1:length(NumDNAOffTargetOptions),'Un',0);
    fprintf("Converting DNA Statistics into Probe Specificity and OFF-target Scores")
    fprintf('\n')
    fprintf('\n')
    wb = parwaitbar(length(Probes_With_DNAOFF),'WaitMessage','Converting');
    for v = 1:length(Probes_With_DNAOFF)
        temp_T = Tvec_DNA{Probes_With_DNAOFF(v)};
        temp_KOFF = Tvec_logKOFF_DNA{Probes_With_DNAOFF(v)};
        temp_KOFFdivON = Tvec_logKOFFdivON_DNA{Probes_With_DNAOFF(v)};
        DNASpecificity_Score(Probes_With_DNAOFF(v)) = dot(EKernel(temp_T)',temp_KOFFdivON);
        DNAOFF_Score(Probes_With_DNAOFF(v)) = dot(EKernel(temp_T)',temp_KOFF);
        progress(wb);
    end
    wb.delete();
    fprintf('\n')
    fprintf('\n')
end
Off_Score = RNAOFF_Score+DNAOFF_Score;
Specificity_Score = RNASpecificity_Score+DNASpecificity_Score;
if (~isfile([settings.FolderRootName filesep '(' TranscriptName ')_' settings.rootName '_Tm' num2str(T_hybrid) '_BasicDesignerStats' designerName '.mat']))
    save([settings.FolderRootName filesep '(' TranscriptName ')_' settings.rootName '_Tm' num2str(T_hybrid) '_BasicDesignerStats' designerName '.mat'],...
        'Nvec_RNAsingle','Nvec_RNAmulti','Svec_RNA','Tvec_RNA','Tvec_logKOFF_RNA','Tvec_logKOFFdivON_RNA','Tvec_logKONdivOFF_RNA',...
        'Nvec_DNAsingle','Nvec_DNAmulti','Svec_DNA','Tvec_DNA','Tvec_logKOFF_DNA','Tvec_logKOFFdivON_DNA','Tvec_logKONdivOFF_DNA','Tvec_logKOFFdivCOMP_DNA','Tvec_logKCOMPdivOFF_DNA',...
        'NTPvec_RNAsingle','NTPvec_RNAmulti','TSvec_RNA','TPvec_RNA','TPvec_logKOFF_RNA','TPvec_logKOFFdivON_RNA','TPvec_logKONdivOFF_RNA',...
        'NTPvec_DNAsingle','NTPvec_DNAmulti','TSvec_DNA','TPvec_DNA','TPvec_logKOFF_DNA','TPvec_logKOFFdivON_DNA','TPvec_logKONdivOFF_DNA','TPvec_logKOFFdivCOMP_DNA','TPvec_logKCOMPdivOFF_DNA',...
        'Off_Score','Specificity_Score','NumRNAOffTargetOptions','Probes_WithNRNAOFF','NumDNAOffTargetOptions','Probes_WithNDNAOFF','-v7.3')
end
end