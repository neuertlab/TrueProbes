function [Kb_mod,Kb_Complement,DoesProbeBindSite2,Num_of_Molecule_Sites,MolProbesAtEvents,Mol_ProbesAtEventsID,nascentInfo,dHeq_mod,dSeq_mod,dHf_mod,dSf_mod,dHr_mod,dSr_mod,Tm_mod,dCp_mod,dHeq_Complement,dSeq_Complement,dHf_Complement,dSf_Complement,dHr_Complement,dSr_Complement,dCp_Complement] = A_JH_GetSiteMapping_NonTile_V6(probes,settings,gene_table,Kb_Match,dHeq_Match,dSeq_Match,dHf_Match,dSf_Match,dHr_Match,dSr_Match,dCpeq_Match,Tm_Match)
N_methods = 8;
N_methods2 = 3;
N_methods3 = 9;
%Jason Hughes code to parsing gene_table to get sites where probes bind
%For RNA, DNA, and complementary binding for double stranded DNA
%Code also computes nascent transcription sites given expression profile
%get sites for on-target DNA, RNA, and Nascent.
%Fix Bug in first and last probes in tiling not entering DoesProbeBindSite
%or or getting mapped.
%Also bug in Koff giving off-target score for on-target ID
Lpmin = min(cell2mat(cellfun(@length,{probes{:,2}},'UniformOutput',false)));
RemoveMisMatches = settings.RemoveMisMatches;

kb = 0.001987204259;%boltzman constant
T_hybrid = settings.HybridizationTemperature;
SaltConcentration = settings.SaltConcentration;
PrimerConcentration = settings.PrimerConcentration;

most_recent_num_local = settings.num_parpool_local;
Organism = settings.Organism;
ChrNum = settings.ChrNum;
GeneChr = settings.GeneChr;
TranscriptName = settings.GeneName;
transcriptID = settings.transcript_IDs;
FolderRootName = settings.FolderRootName;
withNascent = settings.withNascent;
probeBatchSize = settings.BLASTbatchSize;
targetBatchSize = settings.TargetBatchSize;
designerName = settings.designerName;
if (settings.clusterStatus)
    most_recent_num = str2num(getenv('SLURM_JOB_CPUS_PER_NODE'));
else
    most_recent_num = most_recent_num_local;
end
gene_table = sortrows(gene_table,[7 6],'ascend');
gene_table = gene_table(gene_table.Match>=settings.MinHomologySearchTargetSize,:);
MinusStrandedHits = find(contains(gene_table.Strand,'Minus'));
gene_table_NamesZ = convertCharsToStrings(gene_table.Names);
contains_RNA = find(ismember(gene_table_NamesZ,settings.RNAdbParser));clear gene_table_NamesZ
RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);clear contains_RNA
gene_table = gene_table(setdiff(1:size(gene_table,1),RNA_MissedFilteredHits),:);clear RNA_MissedFilteredHits
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
%% Identify location of on-target DNA and RNA molecules in list of DNA/RNA molecules
%generate ID vector to identify molecules that count as on-target hits

%% Parse Gene_Hits to know which molecules to actually look exp values for.
% Ix_DNA = find(strcmp(uniNames,extractBefore(AN_ON_Chr,'.')));
% if (iscell(transcriptID))
%     for v = 1:length(transcriptID)
%         Ix_RNA(v) = find(strcmp(uniNames,extractBefore(transcriptID{v},'.')));
%     end
% else
%     Ix_RNA = find(strcmp(uniNames,extractBefore(transcriptID,'.')));
% end
% Ix = union(Ix_RNA,Ix_DNA);
% IxTypes{1} = Ix_DNA;
% IxTypes{2} = Ix_RNA;
% IxTypes{3} = Ix;
% cTypes = {'DNA','RNA','ALL'};%which molecules are in each type of category
% MaxcTypes = length(cTypes);
% cTypeIndx{1} = DNA_IDs;
% cTypeIndx{2} = NonDNA_IDs;
% cTypeIndx{3} = 1:length(Names);
% cTypeIndx{1} = setdiff(cTypeIndx{1},Ix);
% cTypeIndx{2} = setdiff(cTypeIndx{2},Ix);
% cTypeIndx{3} = setdiff(cTypeIndx{3},Ix);
% IDon = ismember(cTypeIndx{3},Ix);
calcSiteMap = 0;

try
    load([settings.FolderRootName filesep settings.GeneName '_binding_hits_map' settings.designerName '.mat'],'DoesProbeBindSite2','MolN_ProbesAtEvents','Num_of_Molecule_Sites','Mol_ProbesAtEventsID','MolProbesAtEvents')
    calcSiteMap = calcSiteMap + 0;
catch
    calcSiteMap = calcSiteMap + 1;
end

try
    load([settings.FolderRootName filesep settings.GeneName  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' settings.designerName '.mat'],'Kb_mod')
    load([settings.FolderRootName filesep settings.GeneName '_BindingMatrices' settings.designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod')
    calcSiteMap = calcSiteMap + 0;
catch
    calcSiteMap = calcSiteMap + 1;
end

try
    load([settings.FolderRootName filesep settings.GeneName  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' settings.designerName '.mat'],'Kb_Complement')
    load([settings.FolderRootName filesep settings.GeneName '_BindingMatrices2' settings.designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','Tm__Complement')
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
        if (isfile([FolderRootName filesep '(' TranscriptName ')'  designerName '_BindingSiteMapInfo_batch' num2str(i) '.mat']))%check if temp file exists
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
                    %MolN_SitesBoundaryMax{w} = InInterval_Cmax(nonSubset);
                    %MolN_SitesBoundaryMin{w} = InInterval_Cmin(nonSubset);
                    try
                        Mol_ProbesAtEventsID{w} = cellfun(@(x) rowz(x),Sx,'UniformOutput',false);% or rowz(1:length(An));
                    catch
                        Mol_ProbesAtEventsID{w} = [];
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
    Num_of_Molecule_Sites = struct('Num_of_Molecule_Sites', cell(1, N_siteMappingBatches));
    Mol_ProbesAtEventsID = struct('Mol_ProbesAtEventsID',cell(1,N_siteMappingBatches));
    MolProbesAtEvents = struct('MolProbesAtEvents',cell(1,N_siteMappingBatches));
    MolN_ProbesAtEvents = struct('MolN_ProbesAtEvents',cell(1,N_siteMappingBatches));
    PTS_DPS_unique_vector = struct('PTS_DPS_unique_vector',cell(1,N_siteMappingBatches));
    PTSM_DPS_eq_unique_vector = struct('PTSM_DPS_eq_unique_vector',cell(1,N_siteMappingBatches));
    PTSM_DPS_fr_unique_vector = struct('PTSM_DPS_fr_unique_vector',cell(1,N_siteMappingBatches));
    PTSM_DPS_Tm_unique_vector = struct('PTSM_DPS_Tm_unique_vector',cell(1,N_siteMappingBatches));
    Kb_mod_vector =struct('Kb_mod_vector',cell(1,N_siteMappingBatches));
    dHeq_mod_vector = struct('dHeq_mod_vector',cell(1,N_siteMappingBatches));
    dSeq_mod_vector = struct('dSeq_mod_vector',cell(1,N_siteMappingBatches));
    dHf_mod_vector = struct('dHf_mod_vector',cell(1,N_siteMappingBatches));
    dSf_mod_vector = struct('dSf_mod_vector',cell(1,N_siteMappingBatches));
    dHr_mod_vector = struct('dHr_mod_vector',cell(1,N_siteMappingBatches));
    dSr_mod_vector =struct('dSr_mod_vector',cell(1,N_siteMappingBatches));
    Tm_mod_vector = struct('Tm_mod_vector',cell(1,N_siteMappingBatches));
    dCp_mod_vector = struct('dCp_mod_vector',cell(1,N_siteMappingBatches));
    wb = parwaitbar(N_siteMappingBatches+14,'WaitMessage','Aggregating');
    parfor w = 1:N_siteMappingBatches
        pause(0.1);
        if isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_BindingSiteMapInfo_batch' num2str(w) '.mat'])
            if (load_files)
                partial_binding_site_map_info_tmp = load([FolderRootName filesep '(' TranscriptName ')' designerName '_BindingSiteMapInfo_batch' num2str(w) '.mat']).partial_binding_site_map_info_tmp;
                PTS_DPS_unique_vector(w).PTS_DPS_unique_vector = CATnWrapper(arrayfun(@(w_sub) [cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}', 'Un', 0))   ...
                    Batch_siteMapping_Constant.Value{w}(w_sub)*ones(length(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}', 'UniformOutput', 0))),1) ...
                    cell2mat(arrayfun(@(x) x*ones(1,length(partial_binding_site_map_info_tmp{w_sub}{3}{x})),1:length(partial_binding_site_map_info_tmp{w_sub}{1}),'Un',0))' ],...
                    1:length(partial_binding_site_map_info_tmp),'Un',0),1);
                PTSM_DPS_eq_unique_vector(w).PTSM_DPS_eq_unique_vector = CATnWrapper(arrayfun(@(w_sub) ...
                    [repelem(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'Un', 0)),N_methods,1) ...
                    Batch_siteMapping_Constant.Value{w}(w_sub)*repelem(ones(length(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'UniformOutput', 0))),1),N_methods,1) ...
                    repelem(cell2mat(arrayfun(@(x) x*ones(1,length(partial_binding_site_map_info_tmp{w_sub}{3}{x})),1:length(partial_binding_site_map_info_tmp{w_sub}{1}),'Un',0))',N_methods,1) ....
                    repmat([1:N_methods]',[length(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'Un', 0))) 1])],...
                    1:length(partial_binding_site_map_info_tmp),'Un',0),1);
                PTSM_DPS_fr_unique_vector(w).PTSM_DPS_fr_unique_vector = CATnWrapper(arrayfun(@(w_sub) ...
                    [repelem(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'Un', 0)),N_methods2,1) ...
                    Batch_siteMapping_Constant.Value{w}(w_sub)*repelem(ones(length(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'UniformOutput', 0))),1),N_methods2,1) ...
                    repelem(cell2mat(arrayfun(@(x) x*ones(1,length(partial_binding_site_map_info_tmp{w_sub}{3}{x})),1:length(partial_binding_site_map_info_tmp{w_sub}{1}),'Un',0))',N_methods2,1) ....
                    repmat([1:N_methods2]',[length(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'Un', 0))) 1])],...
                    1:length(partial_binding_site_map_info_tmp),'Un',0),1);
                PTSM_DPS_Tm_unique_vector(w).PTSM_DPS_Tm_unique_vector = CATnWrapper(arrayfun(@(w_sub) ...
                    [repelem(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'Un', 0)),N_methods3,1) ...
                    Batch_siteMapping_Constant.Value{w}(w_sub)*repelem(ones(length(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'UniformOutput', 0))),1),N_methods3,1) ...
                    repelem(cell2mat(arrayfun(@(x) x*ones(1,length(partial_binding_site_map_info_tmp{w_sub}{3}{x})),1:length(partial_binding_site_map_info_tmp{w_sub}{1}),'Un',0))',N_methods3,1) ....
                    repmat([1:N_methods3]',[length(cell2mat(cellfun(@(x) x, partial_binding_site_map_info_tmp{w_sub}{3}(:), 'Un', 0))) 1])],...
                    1:length(partial_binding_site_map_info_tmp),'Un',0),1);            
                Num_of_Molecule_Sites(w).Num_of_Molecule_Sites = cellfun(@(x) length(x{1}),partial_binding_site_map_info_tmp);
                Mol_ProbesAtEventsID(w).Mol_ProbesAtEventsID = cellfun(@(x) x{2},partial_binding_site_map_info_tmp,'Un',0);
                MolProbesAtEvents(w).MolProbesAtEvents = cellfun(@(x) x{3},partial_binding_site_map_info_tmp,'Un',0);
                MolN_ProbesAtEvents(w).MolN_ProbesAtEvents = cellfun(@(x) cellfun(@length,x{3}),partial_binding_site_map_info_tmp,'Un',0);
                Kb_mod_vector(w).Kb_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) sum(y{5}{x}(y{4}{x}==p,:),1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
                dHeq_mod_vector(w).dHeq_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{6}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
                dSeq_mod_vector(w).dSeq_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{7}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
                dHf_mod_vector(w).dHf_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{8}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
                dSf_mod_vector(w).dSf_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{9}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
                dHr_mod_vector(w).dHr_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{10}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
                dSr_mod_vector(w).dSr_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{11}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
                Tm_mod_vector(w).Tm_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{12}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
                dCp_mod_vector(w).dCp_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{13}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),partial_binding_site_map_info_tmp,'Un',0))';
            else
                tmp_file = matfile([FolderRootName filesep '(' TranscriptName ')' designerName '_BindingSiteMapInfo_batch' num2str(w) '.mat']);
                PTS_DPS_unique_vector(w).PTS_DPS_unique_vector = CATnWrapper(cellfun(@(individual_target_binding_site_map_info_tmp,w_sub)...
                    [cell2mat(cellfun(@(x) x, individual_target_binding_site_map_info_tmp{3}(:), 'Un', 0))   ...
                    w_sub*ones(length(cell2mat(cellfun(@(x) x, individual_target_binding_site_map_info_tmp{3}(:), 'UniformOutput', 0))),1) ...
                    cell2mat(arrayfun(@(x) x*ones(1,length(individual_target_binding_site_map_info_tmp{3}{x})),1:length(individual_target_binding_site_map_info_tmp{1}),'Un',0))' ],...
                    tmp_file.partial_binding_site_map_info_tmp,num2cell(Batch_siteMapping_Constant.Value{w}),'Un',0),1);
                PTSM_DPS_eq_unique_vector(w).PTSM_DPS_eq_unique_vector = CATnWrapper(cellfun(@(individual_target_binding_site_map_info_tmp,w_sub)...
                    [repelem(cell2mat(cellfun(@(x) x, individual_target_binding_site_map_info_tmp{3}(:), 'Un', 0)),N_methods,1) ...(tmp_file.partial_binding_site_map_info_tmp)
                    w_sub*repelem(ones(length(cell2mat(cellfun(@(x) x, individual_target_binding_site_map_info_tmp{3}(:), 'UniformOutput', 0))),1),N_methods,1) ...
                    repelem(cell2mat(arrayfun(@(x) x*ones(1,length(individual_target_binding_site_map_info_tmp{3}{x})),1:length(individual_target_binding_site_map_info_tmp{1}),'Un',0))',N_methods,1) ....
                    repmat([1:N_methods]',[length(cell2mat(cellfun(@(x) x, individual_target_binding_site_map_info_tmp{3}(:), 'Un', 0))) 1])],...
                    tmp_file.partial_binding_site_map_info_tmp,num2cell(Batch_siteMapping_Constant.Value{w}),'Un',0),1);
                PTSM_DPS_fr_unique_vector(w).PTSM_DPS_fr_unique_vector = CATnWrapper(cellfun(@(individual_target_binding_site_map_info_tmp,w_sub)...
                    [repelem(cell2mat(cellfun(@(x) x,individual_target_binding_site_map_info_tmp{3}(:), 'Un', 0)),N_methods2,1) ...
                    w_sub*repelem(ones(length(cell2mat(cellfun(@(x) x, individual_target_binding_site_map_info_tmp{3}(:), 'UniformOutput', 0))),1),N_methods2,1) ...
                    repelem(cell2mat(arrayfun(@(x) x*ones(1,length(individual_target_binding_site_map_info_tmp{3}{x})),1:length(individual_target_binding_site_map_info_tmp{1}),'Un',0))',N_methods2,1) ....
                    repmat([1:N_methods2]',[length(cell2mat(cellfun(@(x) x, individual_target_binding_site_map_info_tmp{3}(:), 'Un', 0))) 1])],...
                    tmp_file.partial_binding_site_map_info_tmp,num2cell(Batch_siteMapping_Constant.Value{w}),'Un',0),1);
                PTSM_DPS_Tm_unique_vector(w).PTSM_DPS_Tm_unique_vector = CATnWrapper(cellfun(@(individual_target_binding_site_map_info_tmp,w_sub)...
                    [repelem(cell2mat(cellfun(@(x) x, individual_target_binding_site_map_info_tmp{3}(:), 'Un', 0)),N_methods3,1) ...
                    w_sub*repelem(ones(length(cell2mat(cellfun(@(x) x, individual_target_binding_site_map_info_tmp{3}(:), 'UniformOutput', 0))),1),N_methods3,1) ...
                    repelem(cell2mat(arrayfun(@(x) x*ones(1,length(individual_target_binding_site_map_info_tmp{3}{x})),1:length(individual_target_binding_site_map_info_tmp{1}),'Un',0))',N_methods3,1) ....
                    repmat([1:N_methods3]',[length(cell2mat(cellfun(@(x) x, individual_target_binding_site_map_info_tmp{3}(:), 'Un', 0))) 1])],...
                    tmp_file.partial_binding_site_map_info_tmp,num2cell(Batch_siteMapping_Constant.Value{w}),'Un',0),1);
                Num_of_Molecule_Sites(w).Num_of_Molecule_Sites = cellfun(@(x) length(x{1}),tmp_file.partial_binding_site_map_info_tmp);
                Mol_ProbesAtEventsID(w).Mol_ProbesAtEventsID = cellfun(@(x) x{2},tmp_file.partial_binding_site_map_info_tmp,'Un',0);
                MolProbesAtEvents(w).MolProbesAtEvents = cellfun(@(x) x{3},tmp_file.partial_binding_site_map_info_tmp,'Un',0);
                MolN_ProbesAtEvents(w).MolN_ProbesAtEvents = cellfun(@(x) cellfun(@length,x{3}),tmp_file.partial_binding_site_map_info_tmp,'Un',0);
                Kb_mod_vector(w).Kb_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) sum(y{5}{x}(y{4}{x}==p,:),1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),tmp_file.partial_binding_site_map_info_tmp,'Un',0))';
                dHeq_mod_vector(w).dHeq_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{6}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),tmp_file.partial_binding_site_map_info_tmp,'Un',0))';
                dSeq_mod_vector(w).dSeq_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{7}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),tmp_file.partial_binding_site_map_info_tmp,'Un',0))';
                dHf_mod_vector(w).dHf_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{8}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),tmp_file.partial_binding_site_map_info_tmp,'Un',0))';
                dSf_mod_vector(w).dSf_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{9}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),tmp_file.partial_binding_site_map_info_tmp,'Un',0))';
                dHr_mod_vector(w).dHr_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{10}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),tmp_file.partial_binding_site_map_info_tmp,'Un',0))';
                dSr_mod_vector(w).dSr_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{11}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),tmp_file.partial_binding_site_map_info_tmp,'Un',0))';
                Tm_mod_vector(w).Tm_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{12}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),tmp_file.partial_binding_site_map_info_tmp,'Un',0))';
                dCp_mod_vector(w).dCp_mod_vector = cell2mat(cellfun(@(y) cell2mat(arrayfun(@(x) cell2mat(arrayfun(@(p) min(y{13}{x}(y{4}{x}==p,:),[],1),1:length(y{3}{x}),'Un',0)), 1:length(y{1}),'Un',0)),tmp_file.partial_binding_site_map_info_tmp,'Un',0))';
            end
        end
        progress(wb);
    end
    PTS_DPS_unique_vector = vertcat(PTS_DPS_unique_vector(:).PTS_DPS_unique_vector);
    PTSM_DPS_eq_unique_vector = vertcat(PTSM_DPS_eq_unique_vector(:).PTSM_DPS_eq_unique_vector);
    PTSM_DPS_fr_unique_vector = vertcat(PTSM_DPS_fr_unique_vector(:).PTSM_DPS_fr_unique_vector);
    PTSM_DPS_Tm_unique_vector = vertcat(PTSM_DPS_Tm_unique_vector(:).PTSM_DPS_Tm_unique_vector);
    Kb_mod_vector = vertcat(Kb_mod_vector(:).Kb_mod_vector);
    dHeq_mod_vector =vertcat(dHeq_mod_vector(:).dHeq_mod_vector);
    dSeq_mod_vector = vertcat(dSeq_mod_vector(:).dSeq_mod_vector);
    dHf_mod_vector = vertcat(dHf_mod_vector(:).dHf_mod_vector);
    dSf_mod_vector = vertcat(dSf_mod_vector(:).dSf_mod_vector);
    dHr_mod_vector = vertcat(dHr_mod_vector(:).dHr_mod_vector);
    dSr_mod_vector = vertcat(dSr_mod_vector(:).dSr_mod_vector);
    dCp_mod_vector = vertcat(dCp_mod_vector(:).dCp_mod_vector);
    Tm_mod_vector = vertcat(Tm_mod_vector(:).Tm_mod_vector);
    % dGeq_mod_vector = dHeq_mod_vector - (T_hybrid+273.15)*dSeq_mod_vector;
    % dGf_mod_vector = dHf_mod_vector - (T_hybrid+273.15)*dSf_mod_vector;
    % dGr_mod_vector = dHr_mod_vector - (T_hybrid+273.15)*dSr_mod_vector;
    % Keq_mod_vector = exp(-dGeq_mod_vector/(kb*(T_hybrid+273.15)));
    % Kf_mod_vector = exp(-dGf_mod_vector/(kb*(T_hybrid+273.15)));
    % Kr_mod_vector = exp(-dGr_mod_vector/(kb*(T_hybrid+273.15)));
    Num_of_Molecule_Sites = horzcat(Num_of_Molecule_Sites(:).Num_of_Molecule_Sites);progress(wb);
    Mol_ProbesAtEventsID = horzcat(Mol_ProbesAtEventsID(:).Mol_ProbesAtEventsID);progress(wb);
    MolProbesAtEvents = horzcat(MolProbesAtEvents(:).MolProbesAtEvents);progress(wb);
    MolN_ProbesAtEvents = horzcat(MolN_ProbesAtEvents(:).MolN_ProbesAtEvents);progress(wb);
    DoesProbeBindSite = ndSparse.build(PTS_DPS_unique_vector,ones(size(PTS_DPS_unique_vector,1),1),[size(probes,1) numNames max(MaxSitesInBatch)]);progress(wb);clear PTS_DPS_unique_vector
    dHf_mod = ndSparse.build(PTSM_DPS_fr_unique_vector,dHf_mod_vector,[size(probes,1) numNames max(MaxSitesInBatch) N_methods2]);progress(wb);clear dHf_mod_vector
    dSf_mod = ndSparse.build(PTSM_DPS_fr_unique_vector,dSf_mod_vector,[size(probes,1) numNames max(MaxSitesInBatch) N_methods2]);progress(wb);clear dSf_mod_vector
    dHr_mod = ndSparse.build(PTSM_DPS_fr_unique_vector,dHr_mod_vector,[size(probes,1) numNames max(MaxSitesInBatch) N_methods2]);progress(wb);clear dHr_mod_vector
    dSr_mod = ndSparse.build(PTSM_DPS_fr_unique_vector,dSr_mod_vector,[size(probes,1) numNames max(MaxSitesInBatch) N_methods2]);progress(wb);clear dSr_mod_vector PTSM_DPS_fr_unique_vector
    Tm_mod = ndSparse.build(PTSM_DPS_Tm_unique_vector,Tm_mod_vector,[size(probes,1) numNames max(MaxSitesInBatch) N_methods3]);progress(wb);clear Tm_mod_vector PTSM_DPS_Tm_unique_vector
    Kb_mod = ndSparse.build(PTSM_DPS_eq_unique_vector,Kb_mod_vector,[size(probes,1) numNames max(MaxSitesInBatch) N_methods]);progress(wb);clear Kb_mod_vector 
    dHeq_mod = ndSparse.build(PTSM_DPS_eq_unique_vector,dHeq_mod_vector,[size(probes,1) numNames max(MaxSitesInBatch) N_methods]);progress(wb);clear dHeq_mod_vector
    dSeq_mod = ndSparse.build(PTSM_DPS_eq_unique_vector,dSeq_mod_vector,[size(probes,1) numNames max(MaxSitesInBatch) N_methods]);progress(wb);clear dSeq_mod_vector
    dCp_mod = ndSparse.build(PTSM_DPS_eq_unique_vector,dCp_mod_vector,[size(probes,1) numNames max(MaxSitesInBatch) N_methods]);progress(wb);clear dCp_mod_vector PTSM_DPS_eq_unique_vector
    wb.delete();
    fprintf('\n')
    fprintf('\n')
    fprintf("Filtering Binding Site Map to be non-overlapping in adjacent binding site regions between probes on individual targets.")
    fprintf('\n')
    fprintf('\n')
    DoesProbeBindSite2 = DoesProbeBindSite;
       wb = parwaitbar(size(probes,1),'WaitMessage', 'Filtering');
    for p=1:size(probes,1)
        for i=find(sum(DoesProbeBindSite2(p,:,:),3)>0)%molecules where probe hits
            %check for overlap j>i
            %for j=i
            %    (MolN_SitesBoundaryMax{u}(i)-MolN_SitesBoundaryMin{u}(j))*DoesProbeBindSite2(p,u,i)
            %
            %end
            Iz = find(diff([0 full(reshape(DoesProbeBindSite2(p,i,:),[1 size(DoesProbeBindSite2,3)])) 0])>0);
            if (~isempty(Iz))
                DoesProbeBindSite2(p,i,setdiff(1:size(DoesProbeBindSite2,3),Iz)) = 0;
            end
        end
             progress(wb);
    end
    wb.delete();
    save([settings.FolderRootName filesep settings.GeneName '_binding_hits_map' settings.designerName '.mat'],'DoesProbeBindSite','DoesProbeBindSite2','MolN_ProbesAtEvents','Num_of_Molecule_Sites','Mol_ProbesAtEventsID','MolProbesAtEvents','-v7.3')
    save([settings.FolderRootName filesep settings.GeneName  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' settings.designerName '.mat'],'Kb_mod','-v7.3')
    save([settings.FolderRootName filesep settings.GeneName '_BindingMatrices' settings.designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod','dCp_mod','-v7.3')
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
%  targetMatch = arrayfun(@(x) strrep(gene_table3.Alignment{x}(3,:),'-','N'),1:size(gene_table3,1),'UniformOutput',false);%slow not so slow
%     for i=DNA_IDs
%         for site=1:length(MolN_ProbesAtEvents{i})
%             for l=1:MolN_ProbesAtEvents{i}(site)
%                 PI = MolProbesAtEvents{i}{site}(l);
%                 currentEvent = Mol_ProbesAtEventsID{i}{site}(l);
%                 POGmod_Complement(PI,i,site) = F_DeltaGibson(targetMatch{currentEvent},seqrcomplement(lower(targetMatch{currentEvent})),SaltConcentration,T_hybrid,RemoveMisMatches);
%       rowz = find(strcmp(gene_table3.Name,Names{u}));
%         Kb_Sub = Kb_Match(rowz);
%         Pn = gene_table3.ProbeNum(rowz);%Might contain probes with repeats if probe binds site multiple times
%         An = gene_table3.Ax(rowz);
%         Bn = gene_table3.Bx(rowz);
%         Cn = [unique([An Bn]).'];

if (calcEnergyMatrix2)
    POGmod_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods],0);
    Kb_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods],0);
    dCp_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods],0);
    Tm_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods+1],0);
    dHeq_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods],0);
    dSeq_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods],0);
    dHf_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods2],0);
    dSf_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods2],0);
    dHr_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods2],0);
    dSr_Complement = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods2],0);
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
        unique_target_seqs = unique(targetMatch);
        N_DNA_target_seqs = ceil(length(unique_target_seqs)/targetBatchSize);
        R = mod(length(unique_target_seqs),targetBatchSize);
        dnaTargetSeqBatches = cell(1,N_DNA_target_seqs);
        if (R==0)
            for k = 1:N_DNA_target_seqs
                dnaTargetSeqBatches{k} = targetBatchSize*(k-1)+1:targetBatchSize*k;
            end
        else
            for k = 1:N_DNA_target_seqs-1
                dnaTargetSeqBatches{k} = targetBatchSize*(k-1)+1:targetBatchSize*k;
            end
            dnaTargetSeqBatches{N_DNA_target_seqs} = targetBatchSize*(N_DNA_target_seqs-1)+1:targetBatchSize*(N_DNA_target_seqs-1)+R;
        end

        ResultsExist_DNA = zeros(1,N_DNA_target_seqs);
        ResultsDate_DNA = cell(1,N_DNA_target_seqs);
        fprintf("Check if DNA target site mapping probe batch files exist")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_DNA_target_seqs,'WaitMessage','Checking');
        parfor i = 1:N_DNA_target_seqs
            if (isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_TargetComplementBindingSiteMapInfo_batch' num2str(i) '.mat']))%check if temp file exists
                d = dir([FolderRootName filesep '(' TranscriptName ')' designerName '_TargetComplementBindingSiteMapInfo_batch' num2str(i) '.mat']);
                if (d.bytes>0)%check size greater than zero
                    ResultsExist_DNA(i) = 1;
                end
                ResultsDate_DNA{i} = datetime(d.date);
            end
            progress(wb);
        end
        wb.delete();
        fprintf('\n')
        fprintf('\n')
        Results_NotMade_DNA = find(ResultsExist_DNA==0);
        Results_Made_DNA = find(ResultsExist_DNA==1);
        %Sort get most recent ResultsMade GeneHitsMade and GeneHitsTable Made and add to probe_check_list
        if (length(Results_Made_DNA)<=most_recent_num)
            results_check_DNA = Results_Made_DNA;
        else
            Results_RecentMade_Dates_DNA(:,1) = ResultsDate_DNA(Results_Made_DNA);
            Results_RecentMade_Dates_DNA(:,2) = num2cell(Results_Made_DNA);
            Results_RecentMade_Dates_DNA = table2timetable(cell2table(Results_RecentMade_Dates_DNA));
            Results_RecentMade_Dates_DNA = sortrows(Results_RecentMade_Dates_DNA,1,'descend');
            Results_RecentMade_Dates_DNA.Properties.VariableNames = {'ID'};
            results_check_DNA = Results_RecentMade_Dates_DNA.ID(1:most_recent_num).';
            clear Results_RecentMade_Dates_DNA
        end
        batch_nums_to_check_DNA = union(Results_NotMade_DNA,results_check_DNA);
        
CompEQ_I_vector =struct('CrossEQ_I_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
CompEQ_T_vector =struct('CrossEQ_J_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
CompEQ_S_vector =struct('CrossEQ_K_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
CompEQ_M_vector =struct('CrossEQ_M_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
CompFR_P_vector =struct('CrossFR_I_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
CompFR_T_vector =struct('CrossFR_J_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
CompFR_S_vector =struct('CrossFR_K_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
CompFR_M_vector =struct('CrossFR_M_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
Kd_eq_vector =struct('Kd_eq_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dH_Complement_vector = struct('dHd_eq_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dS_Complement_vector = struct('dSd_eq_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dH_Complement_vector = struct('dHd_f_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dSd_f_vector = struct('dSd_f_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dH_Complement_vector =struct('dHd_r_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dSd_Complement_vector = struct('dSd_r_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dCp_Complement_vector = struct('dCpd_eq_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));




unique_ordered_binding_paired_input_sequences_List = parallel.pool.Constant(unique_ordered_binding_paired_input_sequences);
unique_secondary_structure_pair_to_nonunique_entries_List = parallel.pool.Constant(unique_secondary_structure_pair_to_nonunique_entries);
FinalProbeSet_List = parallel.pool.Constant(FinalProbeSet);
targetMatch_List = parallel.pool.Constant(targetMatch);
dnaTargetSeqBatches_List = parallel.pool.Constant(dnaTargetSeqBatches);
targetMatch_List = parallel.pool.Constant(targetMatch);
unique_target_seqs_List  = parallel.pool.Constant(unique_target_seqs);


for unique_calc = 1:length(batch_nums_to_check_DNA)
            % for i=DNA_IDs
            % for site=1:length(MolN_ProbesAtEvents{i})
            %     for l=1:MolN_ProbesAtEvents{i}(site)
            %         PI = MolProbesAtEvents{i}{site}(l);
            %         currentEvent = Mol_ProbesAtEventsID{i}{site}(l);
        for w = 1:length(dnaTargetSeqBatches_List.Value{batch_nums_to_check_DNA(unique_calc)})
               target_seq = unique_target_seqs_List.Value{dnaTargetSeqBatches_List.Value{batch_nums_to_check_DNA(unique_calc)}(w)};
                    [dHeq, dSeq, dGeq, dHf, dSf, ~, dHr, dSr, ~,dCpeq, dTm] = F_DeltaGibson_V3(target_seq,seqrcomplement(lower(target_seq)),SaltConcentration,T_hybrid,PrimerConcentration,sequence_duplexes_thermo_generator_structure);
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
                    cross_locs = find(double(unique_secondary_structure_pair_to_nonunique_entries_List.Value(:,7)==unique_calc).*...
                        double(unique_secondary_structure_pair_to_nonunique_entries_List.Value(:,1)==0));
                    if (~isempty(cross_locs))
                        V_vector = unique_secondary_structure_pair_to_nonunique_entries_List.Value(cross_locs,2);
                        W_vector = unique_secondary_structure_pair_to_nonunique_entries_List.Value(cross_locs,3);
                        K_vector = unique_secondary_structure_pair_to_nonunique_entries_List.Value(cross_locs,4);
                        CompEQ_P_vector(unique_calc).CrossEQ_I_vector = repmat(FinalProbeSet_List.Value(V_vector),[1 N_methods])';
                        CompEQ_T_vector(unique_calc).CrossEQ_J_vector = repmat(FinalProbeSet_List.Value(W_vector),[1 N_methods])';
                        CompEQ_S_vector(unique_calc).CrossEQ_K_vector = repmat(reshape(K_vector,1,[]),[1 N_methods])';
                        CompEQ_M_vector(unique_calc).CrossEQ_M_vector = repelem((1:N_methods)',length(cross_locs),1);
                        CompFR_P_vector(unique_calc).CrossFR_I_vector = repmat(FinalProbeSet_List.Value(V_vector),[1 N_methods2])';
                        CompFR_T_vector(unique_calc).CrossFR_J_vector = repmat(FinalProbeSet_List.Value(W_vector),[1 N_methods2])';
                        CompFR_S_vector(unique_calc).CrossFR_K_vector = repmat(reshape(K_vector,1,[]),[1 N_methods2])';
                        CompFR_M_vector(unique_calc).CrossFR_M_vector = repelem((1:N_methods2)',length(cross_locs),1);
                        Kd_eq_vector(unique_calc).Kd_eq_vector = repelem(exp(-temp_dGeq/(kb*(T_hybrid+273.15))),length(cross_locs),1);
                        dHc_eq_vector(unique_calc).dHd_eq_vector = repelem(temp_dHeq,length(cross_locs),1);
                        dSc_eq_vector(unique_calc).dSd_eq_vector = repelem(temp_dSeq,length(cross_locs),1);
                        dHc_f_vector(unique_calc).dHd_f_vector = repelem(temp_dHf,length(cross_locs),1);
                        dSc_f_vector(unique_calc).dSd_f_vector = repelem(temp_dSf,length(cross_locs),1);
                        dHc_r_vector(unique_calc).dHd_r_vector = repelem(temp_dHr,length(cross_locs),1);
                        dSc_r_vector(unique_calc).dSd_r_vector = repelem(temp_dSr,length(cross_locs),1);
                        dCc_eq_vector(unique_calc).dCpd_eq_vector = repelem(temp_dCpeq,length(cross_locs),1);
                    end

        end




end
ComplementEQ_P_vector =struct('ComplementEQ_P_vector',cell(1,N_DNA_target_seqs));
ComplementEQ_T_vector =struct('ComplementEQ_T_vector',cell(1,N_DNA_target_seqs));
ComplementEQ_S_vector =struct('ComplementEQ_S_vector',cell(1,N_DNA_target_seqs));
ComplementEQ_M_vector =struct('ComplementEQ_M_vector',cell(1,N_DNA_target_seqs));
ComplementFR_P_vector =struct('ComplementFR_P_vector',cell(1,N_DNA_target_seqs));
ComplementFR_T_vector =struct('ComplementFR_T_vector',cell(1,N_DNA_target_seqs));
ComplementFR_S_vector =struct('ComplementFR_S_vector',cell(1,N_DNA_target_seqs));
ComplementFR_M_vector =struct('ComplementFR_M_vector',cell(1,N_DNA_target_seqs));
POGmod_Complement_vector =struct('POGmod_Complement_vector',cell(1,N_DNA_target_seqs));
Kb_Complement_vector =struct('Kb_Complement_vector',cell(1,N_DNA_target_seqs));
dHeq_Complement_vector = struct('dHeq_Complement_vector',cell(1,N_DNA_target_seqs));
dSeq_Complement_vector = struct('dSeq_Complement_vector',cell(1,N_DNA_target_seqs));
dHf_Complement_vector = struct('dHf_Complement_vector',cell(1,N_DNA_target_seqs));
dSf_Complement_vector = struct('dSf_Complement_vector',cell(1,N_DNA_target_seqs));
dHr_Complement_vector =struct('dHr_Complement_vector',cell(1,N_DNA_target_seqs));
dSr_Complement_vector = struct('dSr_Complement_vector',cell(1,N_DNA_target_seqs));
dCp_Complement_vector = struct('dCp_Complement_vector',cell(1,N_DNA_target_seqs));



ComplementEQ_P_vector = vertcat(ComplementEQ_P_vector(:).ComplementEQ_P_vector);
ComplementEQ_T_vector = vertcat(ComplementEQ_T_vector(:).ComplementEQ_T_vector);
ComplementEQ_S_vector= vertcat(ComplementEQ_S_vector(:).ComplementEQ_S_vector);
ComplementEQ_M_vector = vertcat(ComplementEQ_M_vector(:).ComplementEQ_M_vector);
ComplementFR_P_vector = vertcat(ComplementFR_P_vector(:).ComplementEQ_P_vector);
ComplementFR_T_vector = vertcat(ComplementFR_T_vector(:).ComplementEQ_T_vector);
ComplementFR_S_vector= vertcat(ComplementFR_S_vector(:).ComplementEQ_S_vector);
ComplementFR_M_vector = vertcat(ComplementFR_M_vector(:).ComplementEQ_M_vector);
Kb_Complement_vector = vertcat(Kb_Complement_vector(:).Kb_Complement_vector);
POGmod_Complement_vector = vertcat(POGmod_Complement_vector(:).POGmod_Complement_vector);
dHeq_Complement_vector = vertcat(dHeq_Complement_vector(:).dHeq_Complement_vector);
dSeq_Complement_vector = vertcat(dSeq_Complement_vector(:).dSeq_Complement_vector);
dCp_Complement_vector = vertcat(dCp_Complement_vector(:).dCp_Complement_vector);
dHf_Complement_vector = vertcat(dHf_Complement_vector(:).dHf_Complement_vector);
dSf_Complement_vector = vertcat(dSf_Complement_vector(:).dSf_Complement_vector);
dHr_Complement_vector = vertcat(dHr_Complement_vector(:).dHr_Complement_vector);
dSr_Complement_vector = vertcat(dSr_Complement_vector(:).dSr_Complement_vector);
ComplementEQ_PTSM_vector = [ComplementEQ_P_vector ComplementEQ_T_vector ComplementEQ_S_vector ComplementEQ_M_vector];
ComplementFR_PTSM_vector = [ComplementFR_P_vector ComplementFR_T_vector ComplementFR_S_vector ComplementFR_M_vector];
Kb_Complement = ndSparse.build(ComplementEQ_PTSM_vector,Kb_Complement_vector,[size(probes,1),size(probes,1),max(MaxSitesInBatch),N_methods]);
POGmod_Complement = ndSparse.build(ComplementEQ_PTSM_vector,POGmod_Complement_vector,[size(probes,1),size(probes,1),max(MaxSitesInBatch),N_methods]);
dHeq_Complement = ndSparse.build(ComplementEQ_PTSM_vector,dHeq_Complement_vector,[size(probes,1),size(probes,1),max(MaxSitesInBatch),N_methods]);
dSeq_Complement = ndSparse.build(ComplementEQ_PTSM_vector,dSeq_Complement_vector,[size(probes,1),size(probes,1),max(MaxSitesInBatch),N_methods]);
dCp_Complement = ndSparse.build(ComplementEQ_PTSM_vector,dCp_Complement_vector,[size(probes,1),size(probes,1),max(MaxSitesInBatch),N_methods]);
dHf_Complement = ndSparse.build(ComplementFR_PTSM_vector,dHf_Complement_vector,[size(probes,1),size(probes,1),max(MaxSitesInBatch),N_methods2]);
dSf_Complement = ndSparse.build(ComplementFR_PTSM_vector,dSf_Complement_vector,[size(probes,1),size(probes,1),max(MaxSitesInBatch),N_methods2]);
dHr_Complement = ndSparse.build(ComplementFR_PTSM_vector,dHr_Complement_vector,[size(probes,1),size(probes,1),max(MaxSitesInBatch),N_methods2]);
dSr_Complement = ndSparse.build(ComplementFR_PTSM_vector,dSr_Complement_vector,[size(probes,1),size(probes,1),max(MaxSitesInBatch),N_methods2]);




        




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
    save([settings.FolderRootName filesep settings.GeneName  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix2' settings.designerName '.mat'],'POGmod_Complement','Kb_Complement','-v7.3')
    save([settings.FolderRootName filesep settings.GeneName '_BindingMatrices2' settings.designerName '.mat'],'dHeq_Complement','dSeq_Complement','dHf_Complement','dSf_Complement','dHr_Complement','dSr_Complement','Tm_Complement','dCp_Complement','-v7.3')
end
%filter out repeats for DoesProbeBindSite2-3
nascentInfo = [];
% try
%     if (withNascent)
%         if (strcmp(settings.Organism,'Human'))
%             Expression_Null = settings.DoAllGenesHaveSameExpression; %do genes have same expression
%             CellTypeID = settings.CellType_ExprID;
%             HumanOntology = settings.HumanSpecific.Ontology;%Normal/Cancer
%             HumanTissueOnly = settings.HumanSpecific.TissueOrTissueAndCellType; % 0 (TissueOnssue/Cell Type)
%             HumanExpGeneOrTransc = settings.HumanSpecific.HumanExpGeneOrTransc; % 1 (Gene/EMBL GENEID) , 0 (Transcript/EMBL Transcript ID)
%             Human_SelTrack = settings.HumanSpecific.SCOutputTrack;
%             expValType = settings.expressionValType;% 1-4 (expCounts,expValues,mean(CellTypeExpValues),one cell types CellTypeExpValues)
%             NascentRecordID = strcat(HumanOntology,num2str(Expression_Null),num2str(CellTypeID),num2str(HumanTissueOnly),...
%                 num2str(HumanExpGeneOrTransc),num2str(Human_SelTrack),num2str(expValType));
%         else
%             NascentRecordID = '1';
%         end
%         cTypes = {'DNA','Mature RNA','DNA & Mature RNA','Nascent RNA','All'};
%         MaxcTypes = length(cTypes);
%         try
%             load([settings.FolderRootName filesep settings.GeneName '_Tm' num2str(T_hybrid) '_Nascent' NascentRecordID '' settings.designerName '.mat'],'DoesProbeBindSite_TS','Kb_TS','OnTargetTS',...
%                 'NumPolymeraseOnTS','TS_SenseInfo','TS_ChrInfo','AN_ListTS','Num_of_Molecule_SitesTS','maxPolymeraseOnTS',...
%                 'dHeq_TS','dSeq_TS','dHf_TS','dSf_TS','dHr_TS','dSr_TS')
%             calcNascent = 0;
%         catch
%             calcNascent = 1;
%         end
%         if (calcNascent)
%             [DoesProbeBindSite_TS,Kb_TS,OnTargetTS,NumPolymeraseOnTS,maxPolymeraseOnTS,TS_SenseInfo,TS_ChrInfo,Names_TS,Num_of_Molecule_SitesTS,dHeq_TS,dSeq_TS,dHf_TS,dSf_TS,dHr_TS,dSr_TS,dCp_TS] = ...
%                 getNascentBindingInfo_JH3(NumberOfPolymerase,settings,probes,transcript_ID,DoesProbeBindSite2,Kb_mod,cTypeIndx,Ix_DNA,T_hybrid,Names,MolProbesAtEvents,Mol_ProbesAtEventsID,Event_Name,Event_Strands,Event_LocStart,Event_ID,dHeq_mod,dSeq_mod,dHf_mod,dSf_mod,dHr_mod,dSr_mod,dCp_mod);
%             save([settings.FolderRootName filesep settings.GeneName '_Tm' num2str(T_hybrid) '_Nascent' NascentRecordID settings.designerName '.mat'],'DoesProbeBindSite_TS','Kb_TS','OnTargetTS',...
%                 'NumPolymeraseOnTS','TS_SenseInfo','TS_ChrInfo','AN_ListTS','Num_of_Molecule_SitesTS','maxPolymeraseOnTS',...
%                 'dHeq_TS','dSeq_TS','dHf_TS','dSf_TS','dHr_TS','dSr_TS','-v7.3');
%         end
%         Exp_TS = 2*ones(1,length(OnTargetTS));
%         NamesV2 = [Names(:).' Names_TS(:).'];
%         Ix_Nascent = length(Names) + find(OnTargetTS);
%         Ix_All = [Ix Ix_Nascent];
%         IxTypes{1} = Ix_DNA;IxTypes{2} = Ix_RNA;IxTypes{3} = Ix;IxTypes{4} = Ix_Nascent;IxTypes{5} = Ix_All;
%         cTypeIndx{1} = DNA_IDs;cTypeIndx{2} = NonDNA_IDs;cTypeIndx{3} = 1:length(Names);
%         cTypeIndx{4} = length(Names)+1:length(OnTargetTS);cTypeIndx{5} = 1:length(OnTargetTS);
%         cTypeIndx{1} = setdiff(cTypeIndx{1},Ix_All);
%         cTypeIndx{2} = setdiff(cTypeIndx{2},Ix_All);
%         cTypeIndx{3} = setdiff(cTypeIndx{3},Ix_All);
%         cTypeIndx{4} = setdiff(cTypeIndx{4},Ix_All);
%         cTypeIndx{5} = setdiff(cTypeIndx{5},Ix_All);
%         Names = NamesV2;
%         if (size(DoesProbeBindSite2,3)<size(DoesProbeBindSite_TS,3))
%             DoesProbeBindSite2(1,1,size(DoesProbeBindSite_TS,3))=0;
%             Kb_mod(1,1,size(Kb_mod,3),1)=0;
%             dHeq_mod(1,1,size(dHeq_mod,3),1)=0;
%             dSeq_mod(1,1,size(dSeq_mod,3),1)=0;
%             dHf_mod(1,1,size(dHf_mod,3),1)=0;
%             dSf_mod(1,1,size(dSf_mod,3),1)=0;
%             dHr_mod(1,1,size(dHr_mod,3),1)=0;
%             dSr_mod(1,1,size(dSr_mod,3),1)=0;
%             dCp_mod(1,1,size(dCp_mod,3),1)=0;
%         end
%         DoesProbeBindSite2(1:size(probes,1),...
%             size(DoesProbeBindSite2,2)+1:size(DoesProbeBindSite2,2)+size(DoesProbeBindSite_TS,2),...
%             1:size(DoesProbeBindSite_TS,3)) = DoesProbeBindSite_TS;
%         Kb_mod(1:size(probes,1),...
%             size(Kb_mod,2)+1:size(Kb_mod,2)+size(Kb_TS,2),...
%             1:size(Kb_TS,3),:) = Kb_TS;
%         dHeq_mod(1:size(probes,1),...
%             size(dHeq_mod,2)+1:size(dHeq_mod,2)+size(dHeq_TS,2),...
%             1:size(dHeq_TS,3),:) = dHeq_TS;
%         dSeq_mod(1:size(probes,1),...
%             size(dSeq_mod,2)+1:size(dSeq_mod,2)+size(dSeq_TS,2),...
%             1:size(dSeq_TS,3),:) = dSeq_TS;
%         dHf_mod(1:size(probes,1),...
%             size(dHf_mod,2)+1:size(dHf_mod,2)+size(dHf_TS,2),...
%             1:size(dHf_TS,3),:) = dHf_TS;
%         dSf_mod(1:size(probes,1),...
%             size(dSf_mod,2)+1:size(dSf_mod,2)+size(dSf_TS,2),...
%             1:size(dSf_TS,3),:) = dSf_TS;
%         dHr_mod(1:size(probes,1),...
%             size(dHr_mod,2)+1:size(dHr_mod,2)+size(dHr_TS,2),...
%             1:size(dHr_TS,3),:) = dHr_TS;
%         dSr_mod(1:size(probes,1),...
%             size(dSr_mod,2)+1:size(dSr_mod,2)+size(dSr_TS,2),...
%             1:size(dSr_TS,3),:) = dSr_TS;
%         dCp_mod(1:size(probes,1),...
%             size(dCp_mod,2)+1:size(dCp_mod,2)+size(dCp_TS,2),...
%             1:size(dCp_TS,3),:) = dCp_TS;
%         Num_of_Molecule_Sites = [Num_of_Molecule_Sites Num_of_Molecule_SitesTS];
%         IDon(Ix_Nascent) = 1;
%         nascentInfo.Updated.IDon = IDon;
%         nascentInfo.Updated.NumMolSites = Num_of_Molecule_Sites;
%         nascentInfo.Updated.Names = Names;
%         nascentInfo.Expression_TS = Exp_TS;
%         nascentInfo.Updated.cMaxcTypes =  MaxcTypes;
%         nascentInfo.Updated.cTypeIndx =  cTypeIndx;
%         nascentInfo.Updated.ixTypeIndx = IxTypeIndx;
%         nascentInfo.Updated.Kbmod = Kb_mod;
%         nascentInfo.Updated.dHeqmod = dHeq_mod;
%         nascentInfo.Updated.dSeqmod = dSeq_mod;
%         nascentInfo.Updated.dHfmod = dHf_mod;
%         nascentInfo.Updated.dSfmod = dSf_mod;
%         nascentInfo.Updated.dHrmod = dHr_mod;
%         nascentInfo.Updated.dSrmod = dSr_mod;
%         nascentInfo.Updated.DoesProbeBindSite = DoesProbeBindSite2;
%         try
%             nascentInfo.maxPolymerase = maxPolymeraseOnTS;
%         catch ME
%             disp(ME.message)
%         end
%         try
%             nascentInfo.numPolymerase = NumPolymeraseOnTS;
%             nascentInfo.sense_info = TS_SenseInfo;
%             nascentInfo.chr_info = TS_ChrInfo;
%         catch ME
%             disp(ME.message)
%         end
%         save([settings.FolderRootName filesep settings.GeneName '_Tm' num2str(T_hybrid) '_NascentInfo' NascentRecordID settings.designerName '.mat'],'nascentInfo','-v7.3','-append');
% 
%     end
% catch ME
%     disp(ME.message)
% end
end

%math for site conflict and probs

%Overmap(P,I,J) = DPS(P,T,I)*DPS(P,T,J)*double(Cmax{T}(I)>=Cmin{T}(J));
%%J>I
%Tree of overlapping sites.

%decision process for update rule

%store and save branches.

% The first six lines of the code create empty ndSparse matrices and cell arrays to store the output of the function.
%
% The following commented lines are test examples that can be used to verify the function.
%
% The for loop that iterates over each element of the Names cell array. For each Name, it finds the corresponding row in gene_table3, which contains information about the probe's location and the target gene's name.
%
% The Ax and Bx values in the gene_table3 table represent the start and end positions of the probes that target the gene. The code first determines the unique start and end regions where probes exist (Cn). It then finds all intervals and assigns which probes are inside each interval.
%
% The probesInInterval variable is a cell array that contains the indices of the probes that are inside each interval. For example, probesInInterval{1} contains the indices of the probes that are inside the first interval.
%
% The code then finds all sites with unique sets of probes and assigns a number of probes binding to each site (InInterval_order). It also finds the boundary of each probe site (InInterval_Cmax and InInterval_Cmin).
%
% The IsSubSet variable is an ndSparse matrix that contains information about the probe sites that are subsets of another site.
%
% The remaining code inside the for loop checks which sites are a subset of another site and returns the output.

%The second part of the code is modifying the "DoesProbeBindSite2" variable based on the results of the first part.

% The purpose of this code is to ensure that each probe only binds to one site on a molecule. The code accomplishes this by iterating over each probe and each molecule it binds to, and then checking if the probe binds to multiple sites on that molecule. If so, it updates the "DoesProbeBindSite2" variable to only indicate the site where the probe has the highest binding affinity.
%
% Specifically, for each molecule where a probe binds, the code first finds the indices where the "DoesProbeBindSite2" variable is equal to 1 (i.e., where the probe binds to that molecule). It then identifies any discontinuities in these indices, which correspond to locations where the probe binds to multiple sites on the molecule. Finally, it updates the "DoesProbeBindSite2" variable to set all but the highest-affinity binding site to 0.
%
% This is done to ensure that the analysis only considers the highest-affinity binding site for each probe on each molecule. Without this step, the results could be skewed by probes that bind to multiple sites on a single molecule.

%molecule, only one of them is kept. This happens because the code uses the diff function to find the locations of the binding events and keep only the first one. As a result, if the binding events occur sequentially and do not overlap, the code will keep only the first event and discard the rest, leading to an incorrect result.