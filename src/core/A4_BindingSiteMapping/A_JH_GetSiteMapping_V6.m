function [Kb_mod,Kb_Complement,DoesProbeBindSite2,Num_of_Molecule_Sites,MolProbesAtEvents,Mol_ProbesAtEventsID,nascentInfo,dHeq_mod,dSeq_mod,dHf_mod,dSf_mod,dHr_mod,dSr_mod,Tm_mod,dCp_mod,dHeq_Complement,dSeq_Complement,dHf_Complement,dSf_Complement,dHr_Complement,dSr_Complement,dCp_Complement] = A_JH_GetSiteMapping_V6(probes,settings,gene_table,Kb_Match,dHeq_Match,dSeq_Match,dHf_Match,dSf_Match,dHr_Match,dSr_Match,dCpeq_Match,Tm_Match)
N_methods = 8;
N_methods2 = 3;
N_methods3 = 9;
nascentInfo = [];
kb = 0.001987204259;%boltzman constant
%Jason Hughes code to parsing gene_table to get sites where probes bind
%For RNA, DNA, and complementary binding for double stranded DNA
%Code also computes nascent transcription sites given expression profile
%get sites for on-target DNA, RNA, and Nascent.
%or or getting mapped.
%Also bug in Koff giving off-target score for on-target ID
Lpmin = min(cell2mat(cellfun(@length,{probes{:,2}},'UniformOutput',false)));
most_recent_num_local = settings.num_parpool_local;
T_hybrid = settings.HybridizationTemperature;
SaltConcentration = settings.SaltConcentration;
PrimerConcentration = settings.PrimerConcentration;

Organism = settings.Organism;
TranscriptName = settings.GeneName;
designerName = settings.designerName;
FolderRootName = settings.FolderRootName;
targetBatchSize = settings.TargetBatchSize;
if (settings.clusterStatus)
    most_recent_num = str2num(getenv('SLURM_JOB_CPUS_PER_NODE'));
else
    most_recent_num = most_recent_num_local;
end
if (isfile(strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_probetarget_flanking_sequences', settings.designerName,'.mat')));%save flanking sequence file if it does not exist already   
load([strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_probetarget_flanking_sequences', settings.designerName,'.mat')],'probetarget_flanking_info');%save flanking sequence file if it does not exist already
gene_table(:,'Alignment') = array2table(arrayfun(@(n) strcat(probetarget_flanking_info.ProbeRevCompSequence_5primeTo3prime{n}',newline,'|',newline,probetarget_flanking_info.TargetSequence_5primeTo3prime{n}')',1:size(probetarget_flanking_info,1),'Un',0)');
%(link 1 is probe, line 3 is target)
end
gene_table = sortrows(gene_table,[7 6],'ascend');
gene_table = gene_table(gene_table.Match>=settings.MinHomologySearchTargetSize,:);
MinusStrandedHits = find(contains(gene_table.Strand,'Minus'));
gene_table_NamesZ = convertCharsToStrings(gene_table.Name);
contains_RNA = find(ismember(gene_table_NamesZ,settings.RNAdbParser));
RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);
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
%% Parse Gene_Hits to know which molecules to actually look exp values for.
calcSiteMap = 0;
try
    load([settings.FolderRootName filesep '(' TranscriptName ')_binding_hits_map' settings.designerName '.mat'],'DoesProbeBindSite2','MolN_ProbesAtEvents','Num_of_Molecule_Sites','Mol_ProbesAtEventsID','MolProbesAtEvents')
    calcSiteMap = calcSiteMap + 0;
catch
    calcSiteMap = calcSiteMap + 1;
end
try
    load([settings.FolderRootName filesep '(' TranscriptName ')_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' settings.designerName '.mat'],'Kb_mod')
    load([settings.FolderRootName filesep '(' TranscriptName ')_BindingMatrices' settings.designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','dCp_mod','Tm_mod')
    calcSiteMap = calcSiteMap + 0;
catch
    calcSiteMap = calcSiteMap + 1;
end
try
    load([settings.FolderRootName filesep '(' TranscriptName ')_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' settings.designerName '.mat'],'Kb_Complement')
    load([settings.FolderRootName filesep '(' TranscriptName ')_BindingMatrices2' settings.designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','Tm__Complement')
    calcEnergyMatrix2 = 0;
catch
    calcEnergyMatrix2 = 1;
end
if (calcSiteMap > 0)
    N_siteMappingBatches = ceil(length(Names)/targetBatchSize);
    R = mod(length(Names),targetBatchSize);
    Batch_siteMapping = cell(1,N_siteMappingBatches);
    if (R==0)
        for k = 1:N_siteMappingBatches
            Batch_siteMapping{k} = [targetBatchSize*(k-1)+1:targetBatchSize*k];
        end
    else
        for k = 1:N_siteMappingBatches-1
            Batch_siteMapping{k} = [targetBatchSize*(k-1)+1:targetBatchSize*k];
        end
        Batch_siteMapping{N_siteMappingBatches} = [targetBatchSize*(N_siteMappingBatches-1)+1:targetBatchSize*(N_siteMappingBatches-1)+R];
    end
    ResultsExist = zeros(1,N_siteMappingBatches);
    ResultsDate = cell(1,N_siteMappingBatches);
    fprintf("Check if probe target batch binding site map files exist")
    fprintf('\n')
    fprintf('\n')
    wb = parwaitbar(N_siteMappingBatches,'WaitMessage','Checking');
    parfor i = 1:N_siteMappingBatches
        if (isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_BindingSiteMapInfo_batch' num2str(i) '.mat']))%check if temp file exists
            d = dir([FolderRootName filesep '(' TranscriptName ')' designerName '_BindingSiteMapInfo_batch' num2str(i) '.mat']);
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
    MolProbesAtEvents = cell(1,length(Names));
    Num_of_Molecule_Sites = zeros(1,length(Names));
    Mol_ProbesAtEventsID = cell(1,length(Names));
    MolN_ProbesAtEvents = cell(1,length(Names));
    %MolN_SitesBoundaryMax = cell(1,length(Names));
    %MolN_SitesBoundaryMin = cell(1,length(Names));
    Event_Rates_In_Site_at_Target= cell(1,length(Names));
    dHeq_In_Site_at_Target= cell(1,length(Names));
    dSeq_In_Site_at_Target= cell(1,length(Names));
    dHf_In_Site_at_Target= cell(1,length(Names));
    dSf_In_Site_at_Target= cell(1,length(Names));
    dHr_In_Site_at_Target= cell(1,length(Names));
    dSr_In_Site_at_Target= cell(1,length(Names));
    Tm_In_Site_at_Target= cell(1,length(Names));
    dCp_In_Site_at_Target= cell(1,length(Names));
    P_ic_In_Site_at_Target = cell(1,length(Names));
    %Test Examples
    %An = [1 2 4 6].'; Bn = [3 5 7 9].';  12 23 34 [Works] Sx_o 12 123 23
    %An = [1 2 3 5 6 9 11 13].'; Bn = [3 4 7 8 12 10 14 15].'; 123 345 56 57 78
    %An = [1 2 3 5 6 9 11 13].'; Bn = [3 4 7 8 10 12 14 15].';  123 345 56 67 78
    %An = [1 2 4 5 7 9 12 15].'; Bn = [3 6 8 11 10 13 14 16].'; 12 234 345 456 67 8
    %An = [1 3 5 6 9].'; Bn = [2 4 6 8 10].'; 1 2 34 5
    %An = [1 3 5 7 9].'; Bn = [2 4 6 8 10].'; 1 2 3 4 5
    %An = [1 2 5 8].'; Bn = [3 9 6 10].'; 12 23 24
    %An = [1 6 12 17 22 34 41 51].'; Bn = [9 14 20 25 30 40 49 60].'; 12 23 34 45 6 7 8
    %An = [1 2 3 4 5 6 8 11 15 17].'; Bn = [7 9 10 12 13 14 16 17 18 19].';%123456 234567 45678 789 8910
    %An = [1 2 3 4 5 6 8 12 15 17].'; Bn = [7 9 10 11 13 14 16 17 18 19].';%123456 234567 5678 789 8910
    gene_names = gene_table.Name;
    BLASTrna_status = settings.BLASTrna;
    BLASTdna_status = settings.BLASTdna;
    TargetNameList = parallel.pool.Constant(Names);
    AnEventList = parallel.pool.Constant(gene_table.Ax);
    BnEventList = parallel.pool.Constant(gene_table.Bx);
    PnEventList = parallel.pool.Constant(gene_table.ProbeNum);
    Kb_Match_List = parallel.pool.Constant(Kb_Match);
    dHeq_Match_List = parallel.pool.Constant(dHeq_Match);
    dSeq_Match_List = parallel.pool.Constant(dSeq_Match);
    dHf_Match_List = parallel.pool.Constant(dHf_Match);
    dSf_Match_List = parallel.pool.Constant(dSf_Match);
    dHr_Match_List = parallel.pool.Constant(dHr_Match);
    dSr_Match_List = parallel.pool.Constant(dSr_Match);
    Tm_Match_List = parallel.pool.Constant(Tm_Match);
    dCpeq_Match_List = parallel.pool.Constant(dCpeq_Match);
    Batch_siteMapping_Constant = parallel.pool.Constant(Batch_siteMapping);
    if (~isempty(batch_nums_to_check))
        fprintf("Generating probe target binding site map in batches:")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(length(batch_nums_to_check),'WaitMessage', 'Mapping');
        parfor w =1:length(batch_nums_to_check)
            pause(0.1);
            P_ic_In_Site_at_Target{w} = cell(1,length(Batch_siteMapping_Constant.Value{batch_nums_to_check(w)}));
            MolProbesAtEvents{w} = cell(1,length(Batch_siteMapping_Constant.Value{batch_nums_to_check(w)}));
            MolN_ProbesAtEvents{w} = cell(1,length(Batch_siteMapping_Constant.Value{batch_nums_to_check(w)}));
            Event_Rates_In_Site_at_Target{w} = cell(1,length(Batch_siteMapping_Constant.Value{batch_nums_to_check(w)}));
            dHeq_In_Site_at_Target{w} = cell(1,length(Batch_siteMapping_Constant.Value{batch_nums_to_check(w)}));
            dSeq_In_Site_at_Target{w} = cell(1,length(Batch_siteMapping_Constant.Value{batch_nums_to_check(w)}));
            dHf_In_Site_at_Target{w} = cell(1,length(Batch_siteMapping_Constant.Value{batch_nums_to_check(w)}));
            dSf_In_Site_at_Target{w} = cell(1,length(Batch_siteMapping_Constant.Value{batch_nums_to_check(w)}));
            dHr_In_Site_at_Target{w} = cell(1,length(Batch_siteMapping_Constant.Value{batch_nums_to_check(w)}));
            dSr_In_Site_at_Target{w} = cell(1,length(Batch_siteMapping_Constant.Value{batch_nums_to_check(w)}));
            Tm_In_Site_at_Target{w} = cell(1,length(Batch_siteMapping_Constant.Value{batch_nums_to_check(w)}));
            dCp_In_Site_at_Target{w} = cell(1,length(Batch_siteMapping_Constant.Value{batch_nums_to_check(w)}));
            partial_binding_site_map_info_tmp{w} = cell(1,length(Batch_siteMapping_Constant.Value{batch_nums_to_check(w)}));
            for w_sub = 1:length(Batch_siteMapping_Constant.Value{batch_nums_to_check(w)})
                u = Batch_siteMapping_Constant.Value{batch_nums_to_check(w)}(w_sub);
                checkRun = 0;
                isDNA = double(ismember(u,DNA_IDs));
                isRNA = double(ismember(u,NonDNA_IDs));
                if (isRNA == BLASTrna_status)
                    checkRun = 1;
                end
                if (isDNA == BLASTdna_status)
                    checkRun = 1;
                end
                if (checkRun)
                    rowz = find(strcmp(gene_names,TargetNameList.Value{u}));
                    Pn = PnEventList.Value(rowz);%Note: Might contain probes with repeats if probe binds site multiple times
                    An = AnEventList.Value(rowz);
                    Bn = BnEventList.Value(rowz);
                    % get unique start and end regions where probes exist
                    Cn = [unique([An Bn]).'];
                    %get all intervals and assign which probes are inside and finds for each unique region
                    % which events on a target are above the bottom end of it and length is below the size of the segment
                    probesInInterval = cell(1,length(Cn));
                    for x = 1:length(Cn)
                        probesInInterval{x} = find(ge(Cn(x)-An,0).*ge(Bn-An,Cn(x)-An)).';
                    end
                    %Gets sites with unique sets of probes
                    charIntervalArray = cellfun(@num2str, probesInInterval,'Un',0);
                    [charIntervalArray_Unique,ia,~] = unique(charIntervalArray,'stable');
                    probesInInterval_Unique = cellfun(@str2num,charIntervalArray_Unique,'Un',0);
                    %Lists number of probes binding in a site
                    InInterval_order = cellfun(@length,probesInInterval_Unique);
                    %Gets boundary of probe binding sites
                    InInterval_Cmax = cell2mat(cellfun(@(x) max([An(x).' Bn(x).']),probesInInterval_Unique,'UniformOutput',false));
                    InInterval_Cmin = cell2mat(cellfun(@(x) min([An(x).' Bn(x).']),probesInInterval_Unique,'UniformOutput',false));
                    IsSubSet = ndSparse.build([length(probesInInterval_Unique),length(probesInInterval_Unique)],0);
                    %get list of proposed sites which are a subset of another sites
                    points_belowMax = find(InInterval_order<max(InInterval_order)); %Finds points below maximum size
                    for v = points_belowMax
                        %subset to a larger interval
                        points_req1 = find(InInterval_order>InInterval_order(v));
                        %and contained within that interval
                        %i.e. bounds of subset I within larger set J
                        points_req = points_req1(find(ge(InInterval_Cmax(points_req1),InInterval_Cmax(v)).*ge(InInterval_Cmin(v),InInterval_Cmin(points_req1))));
                        if (~isempty(points_req))
                            IsSubSet(v,points_req) = cell2mat(cellfun(@(x) all(ismember(probesInInterval_Unique{v},x)),{probesInInterval_Unique{points_req}},'UniformOutput',false));
                        end
                    end
                    IsSubSet1D = sum(IsSubSet,2);
                    nonSubset = find(IsSubSet1D==0);
                    %get unique sites without subsets
                    Sx = {probesInInterval_Unique{nonSubset}};
                    try
                        Mol_ProbesAtEventsID{w}{w_sub} = cellfun(@(x) rowz(x),Sx,'UniformOutput',false);% or rowz(1:length(An));
                    catch
                        Mol_ProbesAtEventsID{w}{w_sub} = [];
                    end
                    Kb_Sub = Kb_Match_List.Value(rowz,:);
                    dHeq_Sub = dHeq_Match_List.Value(rowz,:);
                    dSeq_Sub = dSeq_Match_List.Value(rowz,:);
                    dHf_Sub = dHf_Match_List.Value(rowz,:);
                    dSf_Sub = dSf_Match_List.Value(rowz,:);
                    dHr_Sub = dHr_Match_List.Value(rowz,:);
                    dSr_Sub = dSr_Match_List.Value(rowz,:);
                    Tm_Sub  = Tm_Match_List.Value(rowz,:);
                    dCp_Sub = dCpeq_Match_List.Value(rowz,:);
                    for x = 1:length(Sx)
                        [P_unique,~,P_ic] = unique(Pn(Sx{x}));
                        P_ic_In_Site_at_Target{w}{w_sub}{x} = P_ic;
                        MolProbesAtEvents{w}{w_sub}{x} = P_unique;
                        MolN_ProbesAtEvents{w}{w_sub}(x) = length(P_unique);
                        Event_Rates_In_Site_at_Target{w}{w_sub}{x}  = Kb_Sub(Sx{x},:);
                        dHeq_In_Site_at_Target{w}{w_sub}{x} = dHeq_Sub(Sx{x},:);
                        dSeq_In_Site_at_Target{w}{w_sub}{x}  = dSeq_Sub(Sx{x},:);
                        dHf_In_Site_at_Target{w}{w_sub}{x}  = dHf_Sub(Sx{x},:);
                        dSf_In_Site_at_Target{w}{w_sub}{x}  =  dSf_Sub(Sx{x},:);
                        dHr_In_Site_at_Target{w}{w_sub}{x} = dHr_Sub(Sx{x},:);
                        dSr_In_Site_at_Target{w}{w_sub}{x}  = dSr_Sub(Sx{x},:);
                        Tm_In_Site_at_Target{w}{w_sub}{x}  = Tm_Sub(Sx{x},:);
                        dCp_In_Site_at_Target{w}{w_sub}{x}  = dCp_Sub(Sx{x},:);
                    end
                    partial_binding_site_map_info_tmp{w}{w_sub} = {Sx,Mol_ProbesAtEventsID{w}{w_sub},MolProbesAtEvents{w}{w_sub},...
                        P_ic_In_Site_at_Target{w}{w_sub},Event_Rates_In_Site_at_Target{w}{w_sub},...
                        dHeq_In_Site_at_Target{w}{w_sub},dSeq_In_Site_at_Target{w}{w_sub},...
                        dHf_In_Site_at_Target{w}{w_sub},dSf_In_Site_at_Target{w}{w_sub},...
                        dHr_In_Site_at_Target{w}{w_sub},dSr_In_Site_at_Target{w}{w_sub},...
                        Tm_In_Site_at_Target{w}{w_sub},dCp_In_Site_at_Target{w}{w_sub}};
                end
            end
            parsave_partial_binding_site_map_info([FolderRootName filesep '(' TranscriptName ')' designerName '_BindingSiteMapInfo_batch' num2str(batch_nums_to_check(w)) '.mat'],partial_binding_site_map_info_tmp{w})
            P_ic_In_Site_at_Target{w}=[];
            MolProbesAtEvents{w}=[];
            MolN_ProbesAtEvents{w}=[];
            Event_Rates_In_Site_at_Target{w}=[];
            dHeq_In_Site_at_Target{w}=[];
            dSeq_In_Site_at_Target{w}=[];
            dHf_In_Site_at_Target{w}=[];
            dSf_In_Site_at_Target{w}=[];
            dHr_In_Site_at_Target{w}=[];
            dSr_In_Site_at_Target{w}=[];
            Tm_In_Site_at_Target{w}=[];
            dCp_In_Site_at_Target{w}=[];
            progress(wb);
        end
        wb.delete();
        fprintf('\n')
        fprintf('\n')
    end
    fprintf("Finding Largest Number of Target Binding Sites Across All Batches")
    fprintf('\n')
    fprintf('\n')
    MaxSitesInBatch = zeros(1,N_siteMappingBatches);
    wb = parwaitbar(N_siteMappingBatches,'WaitMessage','Checking');
    parfor w = 1:N_siteMappingBatches
        if isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_BindingSiteMapInfo_batch' num2str(w) '.mat'])
            partial_binding_site_map_info_tmp = load([FolderRootName filesep '(' TranscriptName ')' designerName '_BindingSiteMapInfo_batch' num2str(w) '.mat']).partial_binding_site_map_info_tmp;
            MaxSitesInBatch(w) = max(cell2mat(cellfun(@(x) length(x{1}),partial_binding_site_map_info_tmp,'Un',0)));
        end
        progress(wb);
    end
    wb.delete();
    fprintf('\n')
    fprintf('\n')
    fprintf("Aggregating probe target binding site map batches")
    fprintf('\n')
    fprintf('\n')
    Num_of_Molecule_Sites_batch = cell(1,N_siteMappingBatches);
    Mol_ProbesAtEventsID_batch = cell(1,N_siteMappingBatches);
    MolProbesAtEvents_batch = cell(1,N_siteMappingBatches);
    MolN_ProbesAtEvents_batch = cell(1,N_siteMappingBatches);
    DoesProbeBindSite_batch = cell(1,N_siteMappingBatches);
    Kb_mod_batch = cell(1,N_siteMappingBatches);
    dHeq_mod_batch = cell(1,N_siteMappingBatches);
    dSeq_mod_batch = cell(1,N_siteMappingBatches);
    dCp_mod_batch = cell(1,N_siteMappingBatches);
    dHf_mod_batch = cell(1,N_siteMappingBatches);
    dHr_mod_batch = cell(1,N_siteMappingBatches);
    dSf_mod_batch = cell(1,N_siteMappingBatches);
    dSr_mod_batch = cell(1,N_siteMappingBatches);
    Tm_mod_batch = cell(1,N_siteMappingBatches);
    wb = parwaitbar(N_siteMappingBatches+14,'WaitMessage','Aggregating');
    parfor w = 1:N_siteMappingBatches
        pause(0.1);
        if isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_BindingSiteMapInfo_batch' num2str(w) '.mat'])
            partial_binding_site_map_info_tmp = load([FolderRootName filesep '(' TranscriptName ')' designerName '_BindingSiteMapInfo_batch' num2str(w) '.mat']).partial_binding_site_map_info_tmp;
            PTS_DPS_unique_vector = CATnWrapper(arrayfun(@(w_sub) [cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'Un', 0))   ...
                w_sub*ones(length(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'UniformOutput', 0))),1) ...Prog
                cell2mat(arrayfun(@(x) x*ones(1,length(partial_binding_site_map_info_tmp{w_sub}{3}{x})),1:length(partial_binding_site_map_info_tmp{w_sub}{1}),'Un',0))' ],...
                1:length(partial_binding_site_map_info_tmp),'Un',0),1);
            PTSM_DPS_eq_unique_vector = CATnWrapper(arrayfun(@(w_sub) ...
                [repelem(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'Un', 0)),N_methods,1) ...
                w_sub*repelem(ones(length(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'UniformOutput', 0))),1),N_methods,1) ...
                repelem(cell2mat(arrayfun(@(x) x*ones(1,length(partial_binding_site_map_info_tmp{w_sub}{3}{x})),1:length(partial_binding_site_map_info_tmp{w_sub}{1}),'Un',0))',N_methods,1) ....
                repmat([1:N_methods]',[length(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'Un', 0))) 1])],...
                1:length(partial_binding_site_map_info_tmp),'Un',0),1);
            PTSM_DPS_fr_unique_vector = CATnWrapper(arrayfun(@(w_sub) ...
                [repelem(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'Un', 0)),N_methods2,1) ...
                w_sub*repelem(ones(length(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'UniformOutput', 0))),1),N_methods2,1) ...
                repelem(cell2mat(arrayfun(@(x) x*ones(1,length(partial_binding_site_map_info_tmp{w_sub}{3}{x})),1:length(partial_binding_site_map_info_tmp{w_sub}{1}),'Un',0))',N_methods2,1) ....
                repmat([1:N_methods2]',[length(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'Un', 0))) 1])],...
                1:length(partial_binding_site_map_info_tmp),'Un',0),1);
            PTSM_DPS_Tm_unique_vector = CATnWrapper(arrayfun(@(w_sub) ...
                [repelem(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'Un', 0)),N_methods3,1) ...
                w_sub*repelem(ones(length(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'UniformOutput', 0))),1),N_methods3,1) ...
                repelem(cell2mat(arrayfun(@(x) x*ones(1,length(partial_binding_site_map_info_tmp{w_sub}{3}{x})),1:length(partial_binding_site_map_info_tmp{w_sub}{1}),'Un',0))',N_methods3,1) ....
                repmat([1:N_methods3]',[length(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'Un', 0))) 1])],...
                1:length(partial_binding_site_map_info_tmp),'Un',0),1);
            Kb_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) sum(y{5}{x}(y{4}{x}==p,:),1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
            dHeq_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{6}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
            dSeq_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{7}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
            dHf_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{8}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
            dSf_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{9}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
            dHr_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{10}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
            dSr_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{11}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
            Tm_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{12}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
            dCp_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{13}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
            Num_of_Molecule_Sites_batch{w} = cellfun(@(x) length(x{1}),partial_binding_site_map_info_tmp);
            Mol_ProbesAtEventsID_batch{w} = cellfun(@(x) x{2},partial_binding_site_map_info_tmp,'Un',0);
            MolProbesAtEvents_batch{w} = cellfun(@(x) x{3},partial_binding_site_map_info_tmp,'Un',0);
            MolN_ProbesAtEvents_batch{w} = cellfun(@(x) cellfun(@length,x{3}),partial_binding_site_map_info_tmp,'Un',0);
            DoesProbeBindSite_batch{w} = ndSparse.build(PTS_DPS_unique_vector,ones(size(PTS_DPS_unique_vector,1),1),[size(probes,1) length(partial_binding_site_map_info_tmp) max(MaxSitesInBatch)]);
            Kb_mod_batch{w} = ndSparse.build(PTSM_DPS_eq_unique_vector,Kb_mod_vector,[size(probes,1) length(partial_binding_site_map_info_tmp) max(MaxSitesInBatch) N_methods]);
            dHeq_mod_batch{w} = ndSparse.build(PTSM_DPS_eq_unique_vector,dHeq_mod_vector,[size(probes,1) length(partial_binding_site_map_info_tmp) max(MaxSitesInBatch) N_methods]);
            dSeq_mod_batch{w} = ndSparse.build(PTSM_DPS_eq_unique_vector,dSeq_mod_vector,[size(probes,1) length(partial_binding_site_map_info_tmp) max(MaxSitesInBatch) N_methods]);
            dHf_mod_batch{w} = ndSparse.build(PTSM_DPS_fr_unique_vector,dHf_mod_vector,[size(probes,1) length(partial_binding_site_map_info_tmp) max(MaxSitesInBatch) N_methods2]);
            dSf_mod_batch{w} = ndSparse.build(PTSM_DPS_fr_unique_vector,dSf_mod_vector,[size(probes,1) length(partial_binding_site_map_info_tmp) max(MaxSitesInBatch) N_methods2]);
            dHr_mod_batch{w} = ndSparse.build(PTSM_DPS_fr_unique_vector,dHr_mod_vector,[size(probes,1) length(partial_binding_site_map_info_tmp) max(MaxSitesInBatch) N_methods2]);
            dSr_mod_batch{w} = ndSparse.build(PTSM_DPS_fr_unique_vector,dSr_mod_vector,[size(probes,1) length(partial_binding_site_map_info_tmp) max(MaxSitesInBatch) N_methods2]);
            Tm_mod_batch{w} = ndSparse.build(PTSM_DPS_Tm_unique_vector,Tm_mod_vector,[size(probes,1) length(partial_binding_site_map_info_tmp) max(MaxSitesInBatch) N_methods3]);
            dCp_mod_batch{w} = ndSparse.build(PTSM_DPS_eq_unique_vector,dCp_mod_vector,[size(probes,1) length(partial_binding_site_map_info_tmp) max(MaxSitesInBatch) N_methods]);
        end
        progress(wb);
    end
    Num_of_Molecule_Sites = horzcat(Num_of_Molecule_Sites_batch{:});progress(wb);
    Mol_ProbesAtEventsID = horzcat(Mol_ProbesAtEventsID_batch{:});progress(wb);
    MolProbesAtEvents = horzcat(MolProbesAtEvents_batch{:});progress(wb);
    MolN_ProbesAtEvents = horzcat(MolN_ProbesAtEvents_batch{:});progress(wb);
    DoesProbeBindSite = horzcat(DoesProbeBindSite_batch{:});progress(wb);
    Kb_mod = horzcat(Kb_mod_batch{:});progress(wb);
    dHeq_mod =horzcat(dHeq_mod_batch{:});progress(wb);
    dSeq_mod = horzcat(dSeq_mod_batch{:});progress(wb);
    dHf_mod = horzcat(dHf_mod_batch{:});progress(wb);
    dSf_mod = horzcat(dSf_mod_batch{:});progress(wb);
    dHr_mod = horzcat(dHr_mod_batch{:});progress(wb);
    dSr_mod = horzcat(dSr_mod_batch{:});progress(wb);
    dCp_mod = horzcat(dCp_mod_batch{:});progress(wb);
    Tm_mod = horzcat(Tm_mod_batch{:});progress(wb);
    wb.delete();
    fprintf('\n')
    fprintf('\n')
    %MolN_SitesBoundaryMax{u} = InInterval_Cmax(nonSubset);
    %MolN_SitesBoundaryMin{u} = InInterval_Cmin(nonSubset);
    %Num_of_Molecule_Sites(u) = length(Sx);
    % one directory (Pi,Ti).
    % 1 probe binds proposed site
    % 2+ probes bind proposed site
    %find(x,u,:) == 1. fixed u, changing x.
    % sites Cmin and Cmax.; type filter set some to zero.
    %finding first one then second.
    %overlap of Cmin and Cmax
    %  (P's, Specific Target, S's)
    %needs to be careful and not losen restrictions on sites that cannot
    %both be bound (if probes at distant regions used in the same probe set
    %appear close to each other on off-target within site mapping)/
    %An = [1 2 4 6].'; Bn = [3 5 7 9].';  12 23 34
    %An = [1 2 3 5 6 9 11 13].'; Bn = [3 4 7 8 12 10 14 15].'; 123 345 56 57 78
    %An = [1 2 3 5 6 9 11 13].'; Bn = [3 4 7 8 10 12 14 15].';  123 345 56 67 78
    %An = [1 2 4 5 7 9 12 15].'; Bn = [3 6 8 11 10 13 14 16].'; 12 234 345 456 67 8
    %An = [1 3 5 6 9].'; Bn = [2 4 6 8 10].'; 1 2 34 5
    %An = [1 3 5 7 9].'; Bn = [2 4 6 8 10].'; 1 2 3 4 5
    %An = [1 2 5 8].'; Bn = [3 9 6 10].'; 12 23 24
    %An = [1 6 12 17 22 34 41 51].'; Bn = [9 14 20 25 30 40 49 60].'; 12 23 34 45 6 7 8
    %An = [1 2 3 4 5 6 8 11 15 17].'; Bn = [7 9 10 12 13 14 16 17 18 19].';%123456 234567 45678 789 8910
    %An = [1 2 3 4 5 6 8 12 15 17].'; Bn = [7 9 10 11 13 14 16 17 18 19].';%123456 234567 5678 789 8910
    %Cmin  Cmax
    %for proposed p at target j
    %find all sites bind,  ordered.
    %if only one probe binds site (filter to get unique)
    %but if multiple probes bind site then grouped ordering of sites
    %turning to zero.
    %P_S
    fprintf("Filtering Binding Site Map to be non-overlapping in adjacent binding site regions between probes on individual targets.")
    fprintf('\n')
    fprintf('\n')
    DoesProbeBindSite2 = DoesProbeBindSite;
    wb = parwaitbar(size(probes,1),'WaitMessage', 'Filtering');
    for p=1:size(probes,1)
        for i=find(sum(DoesProbeBindSite2(p,:,:),3)>0)%molecules where probe hits
            Iz = find(diff([0 full(reshape(DoesProbeBindSite2(p,i,:),[1 size(DoesProbeBindSite2,3)])) 0])>0);
            if (~isempty(Iz))
                DoesProbeBindSite2(p,i,setdiff(1:size(DoesProbeBindSite2,3),Iz)) = 0;
            end
        end
        progress(wb);
    end
    wb.delete();
    save([settings.FolderRootName filesep '(' TranscriptName ')_binding_hits_map' settings.designerName '.mat'],'DoesProbeBindSite','DoesProbeBindSite2','MolN_ProbesAtEvents','Num_of_Molecule_Sites','Mol_ProbesAtEventsID','MolProbesAtEvents','-v7.3')
    save([settings.FolderRootName filesep '(' TranscriptName  ')_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' settings.designerName '.mat'],'Kb_mod','-v7.3')
    save([settings.FolderRootName filesep '(' TranscriptName ')_BindingMatrices' settings.designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod','dCp_mod','-v7.3')
    fprintf('\n')
    fprintf('\n')
    fprintf("Deleting temporary probe-target batch binding site map files")
    fprintf('\n')
    fprintf('\n')
    wb = parwaitbar(N_siteMappingBatches,'WaitMessage', 'Deleting');
    parfor i = 1:N_siteMappingBatches
        if exist([FolderRootName filesep '(' TranscriptName ')' designerName '_BindingSiteMapInfo_batch' num2str(i) '.mat'],'file')        %delete temp mat file if already exists
            delete([FolderRootName filesep '(' TranscriptName ')' designerName '_BindingSiteMapInfo_batch' num2str(i) '.mat'])
        end
        progress(wb);
    end
    wb.delete();
    fprintf('\n')
    fprintf('\n')
end
%MolProbesAtEvents tells you probes at each event on target for each site;=
%  for T tells P and S
%An and Bn probably rederive SiteLoc  but has to be for DPS2
%MolN_Probes at event? number of probes at site S on target T
%  targetMatch = arrayfun(@(x) strrep(gene_table.Alignment{x}(3,:),'-','N'),1:size(gene_table,1),'UniformOutput',false);%slow not so slow
%     for i=DNA_IDs
%         for site=1:length(MolN_ProbesAtEvents{i})
%             for l=1:MolN_ProbesAtEvents{i}(site)
%                 PI = MolProbesAtEvents{i}{site}(l);
%                 currentEvent = Mol_ProbesAtEventsID{i}{site}(l);
%                 POGmod_Complement(PI,i,site) = F_DeltaGibson(targetMatch{currentEvent},seqrcomplement(lower(targetMatch{currentEvent})),SaltConcentration,T_hybrid,RemoveMisMatches);
%       rowz = find(strcmp(gene_table.Name,Names{u}));
%         Kb_Sub = Kb_Match(rowz);
%         Pn = gene_table.ProbeNum(rowz);%Might contain probes with repeats if probe binds site multiple times
%         An = gene_table.Ax(rowz);
%         Bn = gene_table.Bx(rowz);
%         Cn = [unique([An Bn]).'];
if (calcEnergyMatrix2)
    POGmod_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)-Lpmin+1,N_methods],0);
    Kb_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)-Lpmin+1,N_methods],0);
    dCp_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)-Lpmin+1,N_methods],0);
    Tm_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)-Lpmin+1,N_methods+1],0);
    dHeq_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)-Lpmin+1,N_methods],0);
    dSeq_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)-Lpmin+1,N_methods],0);
    dHf_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)-Lpmin+1,N_methods2],0);
    dSf_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)-Lpmin+1,N_methods2],0);
    dHr_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)-Lpmin+1,N_methods2],0);
    dSr_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)-Lpmin+1,N_methods2],0);
    targetMatch = arrayfun(@(x) strrep(gene_table.Alignment{x}(3,:),'-','N'),1:size(gene_table,1),'UniformOutput',false);%slow not so slow
    if (settings.BLASTdna)
     sequence_duplexes_thermo_generator_struct_Multi = struct();
    sequence_duplexes_thermo_generator_struct_Multi.Model{1} = F_NearestNeighbors_Parser('Bres86','src/thirdparty/VarGibbs-4.1/P-BS86.par',[]);
    sequence_duplexes_thermo_generator_struct_Multi.Model{2}  = F_NearestNeighbors_Parser('Sant96','src/thirdparty/VarGibbs-4.1/AOP-SL96.par',[]);
    sequence_duplexes_thermo_generator_struct_Multi.Model{3}  = F_NearestNeighbors_Parser('Sant98','src/thirdparty/VarGibbs-4.1/AOP-SL98.par',[]);
    sequence_duplexes_thermo_generator_struct_Multi.Model{4}   = F_NearestNeighbors_Parser('Sugi96','src/thirdparty/VarGibbs-4.1/P-SG96.par',[]);
    sequence_duplexes_thermo_generator_struct_Multi.Model{5}   = F_NearestNeighbors_Parser('Sant04','src/thirdparty/VarGibbs-4.1/P-SL04.par',[]);
    sequence_duplexes_thermo_generator_struct_Multi.Model{6}   = F_NearestNeighbors_Parser('Allawi97','src/thirdparty/VarGibbs-4.1/P-AL97.par',[]);
    sequence_duplexes_thermo_generator_struct_Multi.Model{7}   = F_NearestNeighbors_Parser('Rejali21','src/thirdparty/VarGibbs-4.1/AOP-RJ21KE.par',[]);
    sequence_duplexes_thermo_generator_struct_Multi.Model{8}   = F_NearestNeighbors_Parser('Martins24','src/thirdparty/VarGibbs-4.1/AOP-OW04-69.par','src/thirdparty/VarGibbs-4.1/AOP-MM-60.par');
    sequence_duplexes_thermo_generator_structure = struct2table([sequence_duplexes_thermo_generator_struct_Multi.Model{:}]);

        fprintf("Getting binding affinity of DNA probe targets complementary reactions")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(length(DNA_IDs),'WaitMessage', 'Computing');
        for i=DNA_IDs
            for site=1:length(MolN_ProbesAtEvents{i})
                for l=1:MolN_ProbesAtEvents{i}(site)
                    PI = MolProbesAtEvents{i}{site}(l);
                    currentEvent = Mol_ProbesAtEventsID{i}{site}(l);
                    [dHeq, dSeq, dGeq, dHf, dSf, ~, dHr, dSr, ~,dCpeq, dTm] = F_DeltaGibson_V3(targetMatch{currentEvent},seqrcomplement(lower(targetMatch{currentEvent})),SaltConcentration,T_hybrid,PrimerConcentration,sequence_duplexes_thermo_generator_structure);
                    POGmod_Complement(PI,i,site,:) = dGeq;
                    dHeq_Complement(PI,i,site,:) = dHeq;
                    dSeq_Complement(PI,i,site,:) = dSeq;
                    dCp_Complement(PI,i,site,:) = dCpeq;
                    dHf_Complement(PI,i,site,:) = dHf;
                    dSf_Complement(PI,i,site,:) = dSf;
                    dHr_Complement(PI,i,site,:) = dHr;
                    dSr_Complement(PI,i,site,:) = dSr;
                    Tm_Complement(PI,i,site,:) = dTm;
                    Kb_Complement(PI,i,site,:) = exp(-full(POGmod_Complement(PI,i,site,:))/(kb*(T_hybrid+273.15)));
                end
            end
            progress(wb);
        end
        delete(wb);
        fprintf('\n')
        fprintf('\n')
    end
    Kb_Complement = squeeze(max(Kb_Complement,[],1));
    dHeq_Complement = squeeze(max(dHeq_Complement,[],1));
    dSeq_Complement = squeeze(max(dSeq_Complement,[],1));
    dHf_Complement = squeeze(max(dHf_Complement,[],1));
    dSf_Complement = squeeze(max(dSf_Complement,[],1));
    dHr_Complement = squeeze(max(dHr_Complement,[],1));
    dSr_Complement = squeeze(max(dSr_Complement,[],1));
    dCp_Complement = squeeze(max(dCp_Complement,[],1));
    Tm_Complement = squeeze(max(Tm_Complement,[],1));
    save([settings.FolderRootName filesep '(' TranscriptName  ')_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' settings.designerName '.mat'],'POGmod_Complement','Kb_Complement','-v7.3')
    save([settings.FolderRootName filesep '(' TranscriptName ')_BindingMatrices2' settings.designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','Tm_Complement','dCp_Complement','-v7.3')
end
end

