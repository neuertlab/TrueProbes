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

Organism = settings.Organism;
ChrNum = settings.ChrNum;
GeneChr = settings.GeneChr;
TranscriptName = settings.GeneName;
transcriptID = settings.transcript_IDs;
FolderRootName = settings.FolderRootName;
withNascent = settings.withNascent;
probeBatchSize = settings.BLASTbatchSize;
targetBatchSize = settings.TargetBatchSize;


gene_table = sortrows(gene_table,[7 6],'ascend');
gene_table = gene_table(gene_table.Match>=settings.MinHomologySearchTargetSize,:);
MinusStrandedHits = find(contains(gene_table.Strand,'Minus'));


RNA_IDs_1 = find(contains(gene_table.Name,'NM_'));
RNA_IDs_2 = find(contains(gene_table.Name,'NR_'));
RNA_IDs_3 = find(contains(gene_table.Name,'XM_'));
RNA_IDs_4 = find(contains(gene_table.Name,'XR_'));
contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);


RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);
gene_table = gene_table(setdiff(1:size(gene_table,1),RNA_MissedFilteredHits),:);
gene_table.Ax = min(gene_table.SubjectIndices,[],2);
gene_table.Bx = max(gene_table.SubjectIndices,[],2);
gene_table = sortrows(gene_table,[7 13],'ascend');
Names = unique(gene_table.Name);
Names = convertCharsToStrings(Names);
uniNames = extractBefore(Names,'.');

DNA_IDs_1 = find(contains(uniNames,'NC_'));%IDs
DNA_IDs_2 = find(contains(uniNames,'NT_'));%IDs
DNA_IDs_3 = find(contains(uniNames,'NW_'));%IDs
NonDNA_IDs_1 = find(~contains(uniNames,'NC_'));%IDs
NonDNA_IDs_2 = find(~contains(uniNames,'NT_'));%IDs
NonDNA_IDs_3 = find(~contains(uniNames,'NW_'));%

DNA_IDs =union(union(DNA_IDs_1,DNA_IDs_2),DNA_IDs_3).';
NonDNA_IDs = intersect(intersect(NonDNA_IDs_1,NonDNA_IDs_2),NonDNA_IDs_3).';



%% Identify location of on-target DNA and RNA molecules in list of DNA/RNA molecules
%generate ID vector to identify molecules that count as on-target hits
if (strcmp(Organism,'Human'))
    S8 = readtable(settings.HumanGenomeAssemblyReportFile);
    AN_Chr = S8.RefSeq_Accn(1:639);%%chromosomes 1-22, X (23), Y(24), patchs/scaffolds 25-638, chrM MT (639)
    if (strcmpi(ChrNum,'X'))
        AN_ON_Chr = AN_Chr{23};
    elseif (strcmpi(ChrNum,'Y'))
        AN_ON_Chr = AN_Chr{24};
    elseif (strcmpi(ChrNum,'M'))
        AN_ON_Chr = 'NC_012920.1';
    else
        AN_ON_Chr = AN_Chr{str2double(ChrNum)};
    end
elseif strcmp(Organism,'Mouse')
    S16 = readtable(settings.MouseGenomeAssemblyReportFile);
    AN_Chr = S16.RefSeq_Accn;%chromosomes 1-19, X (20), Y(21), patche scaffolds 22-60, chrM MT (61)
    if (strcmpi(ChrNum,'X'))
        AN_ON_Chr = AN_Chr{20};
    elseif (strcmpi(ChrNum,'Y'))
        AN_ON_Chr = AN_Chr{21};
    elseif (strcmpi(GeneChr,'MT'))
        AN_ON_Chr = 'NC_005089.1';
    else
        AN_ON_Chr = AN_Chr{str2double(ChrNum)};
    end
elseif strcmp(Organism,'Yeast')
    S21 = readtable(settings.YeastGenomeAssemblyReportFile);
    AN_Chr = S21.RefSeq_Accn;%chromosome 1-16, chromosome M (Mitocondria) (17)
    if (strcmpi(GeneChr,'MT'))
        AN_ON_Chr = AN_Chr{17};
    else
        AN_ON_Chr = AN_Chr{roman2num(GeneChr)};
    end
else
    S22 = readtable(settings.CustomGenomeAssemblyReportFile);
    AN_Chr = S22.RefSeq_Accn;%chromosome 1-16, chromosome M (Mitocondria) (17)
    if (strcmpi(ChrNum,'X'))
        AN_ON_Chr = AN_Chr{settings.Custom_X_ChromNumber};
    elseif (strcmpi(ChrNum,'Y'))
        AN_ON_Chr = AN_Chr{settings.Custom_Y_ChromNumber};
    elseif (strcmpi(GeneChr,'MT'))
        AN_ON_Chr = AN_Chr{settings.Custom_MT_ChromNumber};
    else
        AN_ON_Chr = AN_Chr{str2double(ChrNum)};
    end
end
%% Parse Gene_Hits to know which molecules to actually look exp values for.
Ix_DNA = find(strcmp(uniNames,extractBefore(AN_ON_Chr,'.')));
if (iscell(transcriptID))
    for v = 1:length(transcriptID)
        Ix_RNA(v) = find(strcmp(uniNames,extractBefore(transcriptID{v},'.')));
    end
else
    Ix_RNA = find(strcmp(uniNames,extractBefore(transcriptID,'.')));
end
Ix = union(Ix_RNA,Ix_DNA);
IxTypes{1} = Ix_DNA;
IxTypes{2} = Ix_RNA;
IxTypes{3} = Ix;
cTypes = {'DNA','RNA','ALL'};%which molecules are in each type of category
MaxcTypes = length(cTypes);
cTypeIndx{1} = DNA_IDs;
cTypeIndx{2} = NonDNA_IDs;
cTypeIndx{3} = 1:length(Names);
cTypeIndx{1} = setdiff(cTypeIndx{1},Ix);
cTypeIndx{2} = setdiff(cTypeIndx{2},Ix);
cTypeIndx{3} = setdiff(cTypeIndx{3},Ix);
IDon = ismember(cTypeIndx{3},Ix);
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
    ResultsSize = zeros(1,N_siteMappingBatches);
    ResultsDate = cell(1,N_siteMappingBatches);
    for i = 1:N_siteMappingBatches
        if (isfile([settings.FolderRootName filesep TranscriptName  settings.designerName '_BindingSiteMapInfo_batch' num2str(i) '.mat']))%check if temp file exists
            d = dir([settings.FolderRootName filesep TranscriptName settings.designerName '_BindingSiteMapInfo_batch' num2str(i) '.mat']);
            if (d.bytes>0)%check size greater than zero
                ResultsExist(i) = 1;
            end
            ResultsSize(i) = d.bytes;
            ResultsDate{i} = datetime(d.date);
            clear d
        end
    end
    Results_NotMade = find(ResultsExist==0);
    Results_Made = find(ResultsExist==1);
    if (settings.clusterStatus)
        most_recent_num = str2num(getenv('SLURM_JOB_CPUS_PER_NODE'));
    else
        most_recent_num = 8;
    end
    %Sort get most 8 recent ResultsMade GeneHitsMade and GeneHitsTable Made and
    %add to probe_check_list
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

    DoesProbeBindSite = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1],0);
    Kb_mod = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods],0);
    dHeq_mod = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods],0);
    dSeq_mod = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods],0);
    dCp_mod = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods],0);
    dHf_mod = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods2],0);
    dHr_mod = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods2],0);
    dSf_mod = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods2],0);
    dSr_mod = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods2],0);
    Tm_mod = ndSparse.build([size(probes,1),length(Names),size(probes,1)+1,N_methods+1],0);

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
        parsave_partial_binding_site_map_info([FolderRootName filesep TranscriptName designerName '_BindingSiteMapInfo_batch' num2str(batch_nums_to_check(w)) '.mat'],partial_binding_site_map_info_tmp{w})
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
    end

    for w = 1:N_siteMappingBatches
        if isfile([FolderRootName filesep TranscriptName designerName '_BindingSiteMapInfo_batch' num2str(w) '.mat'])
            load([settings.FolderRootName filesep TranscriptName designerName '_BindingSiteMapInfo_batch' num2str(w) '.mat'],'partial_binding_site_map_info_tmp');
            for w_sub = 1:length(Batch_siteMapping{w})
                u = Batch_siteMapping{w}(w_sub);
                Sx = partial_binding_site_map_info_tmp{w_sub}{1};
                Num_of_Molecule_Sites(u) = length(Sx);
                Mol_ProbesAtEventsID{u} = partial_binding_site_map_info_tmp{w_sub}{2};
                MolProbesAtEvents{u} = partial_binding_site_map_info_tmp{w_sub}{3};
                MolN_ProbesAtEvents{u} = cellfun(@length,partial_binding_site_map_info_tmp{w_sub}{3});
                P_DPS_unique_list = CATnWrapper(cellfun(@(x) x',partial_binding_site_map_info_tmp{w_sub}{3},'Un',0),2);
                S_DPS_unique_list = CATnWrapper(arrayfun(@(x) x*ones(size(partial_binding_site_map_info_tmp{w_sub}{3}{x}))',1:length(partial_binding_site_map_info_tmp{w_sub}{3}),'Un',0),2);
                T_DPS_unique_list = u*ones(size(CATnWrapper(arrayfun(@(x) x*ones(size(partial_binding_site_map_info_tmp{w_sub}{3}{x}))',1:length(partial_binding_site_map_info_tmp{w_sub}{3}),'Un',0),2)));
                peq_linear_list = CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(p) partial_binding_site_map_info_tmp{w_sub}{3}{x}(p)*ones(1,N_methods),1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                pfr_linear_list = CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(p) partial_binding_site_map_info_tmp{w_sub}{3}{x}(p)*ones(1,N_methods2),1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                pTm_linear_list = CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(p) partial_binding_site_map_info_tmp{w_sub}{3}{x}(p)*ones(1,N_methods3),1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                xeq_linear_list = CATnWrapper(arrayfun(@(x) x*CATnWrapper(arrayfun(@(p) ones(1,N_methods),1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                xfr_linear_list = CATnWrapper(arrayfun(@(x) x*CATnWrapper(arrayfun(@(p) ones(1,N_methods2),1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                xTm_linear_list = CATnWrapper(arrayfun(@(x) x*CATnWrapper(arrayfun(@(p) ones(1,N_methods3),1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                meq_linear_list = CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(p) 1:N_methods,1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                mfr_linear_list = CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(p) 1:N_methods2,1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                mTm_linear_list = CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(p) 1:N_methods3,1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                if (length(Sx)>size(DoesProbeBindSite,3))
                    DoesProbeBindSite(1,u,length(Sx)) = 0;
                    Kb_mod(1,u,length(Sx),1) = 0;
                    dHeq_mod(1,u,length(Sx),1) = 0;
                    dSeq_mod(1,u,length(Sx),1) = 0;
                    dHf_mod(1,u,length(Sx),1) = 0;
                    dSf_mod(1,u,length(Sx),1) = 0;
                    dHr_mod(1,u,length(Sx),1) = 0;
                    dSr_mod(1,u,length(Sx),1) = 0;
                    Tm_mod(1,u,length(Sx),1) = 0;
                    dCp_mod(1,u,length(Sx),1) = 0;
                end
                DoesProbeBindSite(sub2ind(size(DoesProbeBindSite),P_DPS_unique_list,T_DPS_unique_list,S_DPS_unique_list)) = 1;
                Kb_mod(sub2ind(size(Kb_mod),peq_linear_list,u*ones(size(peq_linear_list)),xeq_linear_list,meq_linear_list)) = ...
                    CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(p) sum(partial_binding_site_map_info_tmp{w_sub}{5}{x}(partial_binding_site_map_info_tmp{w_sub}{4}{x}==p,:),1),1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                dHeq_mod(sub2ind(size(dHeq_mod),peq_linear_list,u*ones(size(peq_linear_list)),xeq_linear_list,meq_linear_list)) = ...
                    CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(p) min(partial_binding_site_map_info_tmp{w_sub}{6}{x}(partial_binding_site_map_info_tmp{w_sub}{4}{x}==p,:),[],1),1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                dSeq_mod(sub2ind(size(dSeq_mod),peq_linear_list,u*ones(size(peq_linear_list)),xeq_linear_list,meq_linear_list)) = ...
                    CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(p) min(partial_binding_site_map_info_tmp{w_sub}{7}{x}(partial_binding_site_map_info_tmp{w_sub}{4}{x}==p,:),[],1),1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                dHf_mod(sub2ind(size(dHf_mod),pfr_linear_list,u*ones(size(pfr_linear_list)),xfr_linear_list,mfr_linear_list)) = ...
                    CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(p) min(partial_binding_site_map_info_tmp{w_sub}{8}{x}(partial_binding_site_map_info_tmp{w_sub}{4}{x}==p,:),[],1),1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                dSf_mod(sub2ind(size(dSf_mod),pfr_linear_list,u*ones(size(pfr_linear_list)),xfr_linear_list,mfr_linear_list)) = ...
                    CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(p) min(partial_binding_site_map_info_tmp{w_sub}{9}{x}(partial_binding_site_map_info_tmp{w_sub}{4}{x}==p,:),[],1),1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                dHr_mod(sub2ind(size(dHr_mod),pfr_linear_list,u*ones(size(pfr_linear_list)),xfr_linear_list,mfr_linear_list)) = ...
                    CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(p) min(partial_binding_site_map_info_tmp{w_sub}{10}{x}(partial_binding_site_map_info_tmp{w_sub}{4}{x}==p,:),[],1),1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                dSr_mod(sub2ind(size(dSr_mod),pfr_linear_list,u*ones(size(pfr_linear_list)),xfr_linear_list,mfr_linear_list)) = ...
                    CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(p) min(partial_binding_site_map_info_tmp{w_sub}{11}{x}(partial_binding_site_map_info_tmp{w_sub}{4}{x}==p,:),[],1),1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                Tm_mod(sub2ind(size(Tm_mod),pTm_linear_list,u*ones(size(pTm_linear_list)),xTm_linear_list,mTm_linear_list)) = ...
                    CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(p) min(partial_binding_site_map_info_tmp{w_sub}{12}{x}(partial_binding_site_map_info_tmp{w_sub}{4}{x}==p,:),[],1),1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
                dCp_mod(sub2ind(size(dCp_mod),peq_linear_list,u*ones(size(peq_linear_list)),xeq_linear_list,meq_linear_list)) = ...
                    CATnWrapper(arrayfun(@(x) CATnWrapper(arrayfun(@(p) min(partial_binding_site_map_info_tmp{w_sub}{13}{x}(partial_binding_site_map_info_tmp{w_sub}{4}{x}==p,:),[],1),1:length(partial_binding_site_map_info_tmp{w_sub}{3}{x}),'Un',0),2),1:length(Sx),'Un',0),2);
            end
            clear partial_binding_site_map_info_tmp
        end
    end

    DoesProbeBindSite2 = DoesProbeBindSite;
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
    end

    save([settings.FolderRootName filesep settings.GeneName '_binding_hits_map' settings.designerName '.mat'],'DoesProbeBindSite','DoesProbeBindSite2','MolN_ProbesAtEvents','Num_of_Molecule_Sites','Mol_ProbesAtEventsID','MolProbesAtEvents','-v7.3')
    save([settings.FolderRootName filesep settings.GeneName  '_Tm' num2str(T_hybrid) '_BindingEnergyMatrix' settings.designerName '.mat'],'Kb_mod','-v7.3')
    save([settings.FolderRootName filesep settings.GeneName '_BindingMatrices' settings.designerName '.mat'],'dHeq_mod','dSeq_mod','dHf_mod','dSf_mod','dHr_mod','dSr_mod','Tm_mod','dCp_mod','-v7.3')
    for i = 1:N_siteMappingBatches
        if exist([settings.FolderRootName filesep TranscriptName settings.designerName '_BindingSiteMapInfo_batch' num2str(i) '.mat'],'file')        %delete temp mat file if already exists
            delete([settings.FolderRootName filesep TranscriptName settings.designerName '_BindingSiteMapInfo_batch' num2str(i) '.mat'])
        end
    end

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
        for i=DNA_IDs
            for site=1:length(MolN_ProbesAtEvents{i})
                for l=1:MolN_ProbesAtEvents{i}(site)
                    PI = MolProbesAtEvents{i}{site}(l);
                    currentEvent = Mol_ProbesAtEventsID{i}{site}(l);
                    [dHeq, dSeq, dGeq, dHf, dSf, ~, dHr, dSr, ~,dCpeq, dTm] = F_DeltaGibson_V3(targetMatch{currentEvent},seqrcomplement(lower(targetMatch{currentEvent})),SaltConcentration,T_hybrid);
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
        end
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
try
    if (withNascent)
        if (strcmp(settings.Organism,'Human'))
            Expression_Null = settings.DoAllGenesHaveSameExpression; %do genes have same expression
            CellTypeID = settings.CellType_ExprID;
            HumanOntology = settings.HumanSpecific.Ontology;%Normal/Cancer
            HumanTissueOnly = settings.HumanSpecific.TissueOrTissueAndCellType; % 0 (TissueOnssue/Cell Type)
            HumanExpGeneOrTransc = settings.HumanSpecific.HumanExpGeneOrTransc; % 1 (Gene/EMBL GENEID) , 0 (Transcript/EMBL Transcript ID)
            Human_SelTrack = settings.HumanSpecific.SCOutputTrack;
            expValType = settings.expressionValType;% 1-4 (expCounts,expValues,mean(CellTypeExpValues),one cell types CellTypeExpValues)
            NascentRecordID = strcat(HumanOntology,num2str(Expression_Null),num2str(CellTypeID),num2str(HumanTissueOnly),...
                num2str(HumanExpGeneOrTransc),num2str(Human_SelTrack),num2str(expValType));
        else
            NascentRecordID = '1';
        end
        cTypes = {'DNA','Mature RNA','DNA & Mature RNA','Nascent RNA','All'};
        MaxcTypes = length(cTypes);
        try
            load([settings.FolderRootName filesep settings.GeneName '_Tm' num2str(T_hybrid) '_Nascent' NascentRecordID '' settings.designerName '.mat'],'DoesProbeBindSite_TS','Kb_TS','OnTargetTS',...
                'NumPolymeraseOnTS','TS_SenseInfo','TS_ChrInfo','AN_ListTS','Num_of_Molecule_SitesTS','maxPolymeraseOnTS',...
                'dHeq_TS','dSeq_TS','dHf_TS','dSf_TS','dHr_TS','dSr_TS')
            calcNascent = 0;
        catch
            calcNascent = 1;
        end
        if (calcNascent)
            [DoesProbeBindSite_TS,Kb_TS,OnTargetTS,NumPolymeraseOnTS,maxPolymeraseOnTS,TS_SenseInfo,TS_ChrInfo,Names_TS,Num_of_Molecule_SitesTS,dHeq_TS,dSeq_TS,dHf_TS,dSf_TS,dHr_TS,dSr_TS,dCp_TS] = ...
                getNascentBindingInfo_JH3(NumberOfPolymerase,settings,probes,transcript_ID,DoesProbeBindSite2,Kb_mod,cTypeIndx,Ix_DNA,T_hybrid,Names,MolProbesAtEvents,Mol_ProbesAtEventsID,Event_Name,Event_Strands,Event_LocStart,Event_ID,dHeq_mod,dSeq_mod,dHf_mod,dSf_mod,dHr_mod,dSr_mod,dCp_mod);
            save([settings.FolderRootName filesep settings.GeneName '_Tm' num2str(T_hybrid) '_Nascent' NascentRecordID settings.designerName '.mat'],'DoesProbeBindSite_TS','Kb_TS','OnTargetTS',...
                'NumPolymeraseOnTS','TS_SenseInfo','TS_ChrInfo','AN_ListTS','Num_of_Molecule_SitesTS','maxPolymeraseOnTS',...
                'dHeq_TS','dSeq_TS','dHf_TS','dSf_TS','dHr_TS','dSr_TS','-v7.3');
        end
        Exp_TS = 2*ones(1,length(OnTargetTS));
        NamesV2 = [Names(:).' Names_TS(:).'];
        Ix_Nascent = length(Names) + find(OnTargetTS);
        Ix_All = [Ix Ix_Nascent];
        IxTypes{1} = Ix_DNA;IxTypes{2} = Ix_RNA;IxTypes{3} = Ix;IxTypes{4} = Ix_Nascent;IxTypes{5} = Ix_All;
        cTypeIndx{1} = DNA_IDs;cTypeIndx{2} = NonDNA_IDs;cTypeIndx{3} = 1:length(Names);
        cTypeIndx{4} = length(Names)+1:length(OnTargetTS);cTypeIndx{5} = 1:length(OnTargetTS);
        cTypeIndx{1} = setdiff(cTypeIndx{1},Ix_All);
        cTypeIndx{2} = setdiff(cTypeIndx{2},Ix_All);
        cTypeIndx{3} = setdiff(cTypeIndx{3},Ix_All);
        cTypeIndx{4} = setdiff(cTypeIndx{4},Ix_All);
        cTypeIndx{5} = setdiff(cTypeIndx{5},Ix_All);
        Names = NamesV2;
        if (size(DoesProbeBindSite2,3)<size(DoesProbeBindSite_TS,3))
            DoesProbeBindSite2(1,1,size(DoesProbeBindSite_TS,3))=0;
            Kb_mod(1,1,size(Kb_mod,3),1)=0;
            dHeq_mod(1,1,size(dHeq_mod,3),1)=0;
            dSeq_mod(1,1,size(dSeq_mod,3),1)=0;
            dHf_mod(1,1,size(dHf_mod,3),1)=0;
            dSf_mod(1,1,size(dSf_mod,3),1)=0;
            dHr_mod(1,1,size(dHr_mod,3),1)=0;
            dSr_mod(1,1,size(dSr_mod,3),1)=0;
            dCp_mod(1,1,size(dCp_mod,3),1)=0;
        end
        DoesProbeBindSite2(1:size(probes,1),...
            size(DoesProbeBindSite2,2)+1:size(DoesProbeBindSite2,2)+size(DoesProbeBindSite_TS,2),...
            1:size(DoesProbeBindSite_TS,3)) = DoesProbeBindSite_TS;
        Kb_mod(1:size(probes,1),...
            size(Kb_mod,2)+1:size(Kb_mod,2)+size(Kb_TS,2),...
            1:size(Kb_TS,3),:) = Kb_TS;
        dHeq_mod(1:size(probes,1),...
            size(dHeq_mod,2)+1:size(dHeq_mod,2)+size(dHeq_TS,2),...
            1:size(dHeq_TS,3),:) = dHeq_TS;
        dSeq_mod(1:size(probes,1),...
            size(dSeq_mod,2)+1:size(dSeq_mod,2)+size(dSeq_TS,2),...
            1:size(dSeq_TS,3),:) = dSeq_TS;
        dHf_mod(1:size(probes,1),...
            size(dHf_mod,2)+1:size(dHf_mod,2)+size(dHf_TS,2),...
            1:size(dHf_TS,3),:) = dHf_TS;
        dSf_mod(1:size(probes,1),...
            size(dSf_mod,2)+1:size(dSf_mod,2)+size(dSf_TS,2),...
            1:size(dSf_TS,3),:) = dSf_TS;
        dHr_mod(1:size(probes,1),...
            size(dHr_mod,2)+1:size(dHr_mod,2)+size(dHr_TS,2),...
            1:size(dHr_TS,3),:) = dHr_TS;
        dSr_mod(1:size(probes,1),...
            size(dSr_mod,2)+1:size(dSr_mod,2)+size(dSr_TS,2),...
            1:size(dSr_TS,3),:) = dSr_TS;
        dCp_mod(1:size(probes,1),...
            size(dCp_mod,2)+1:size(dCp_mod,2)+size(dCp_TS,2),...
            1:size(dCp_TS,3),:) = dCp_TS;
        Num_of_Molecule_Sites = [Num_of_Molecule_Sites Num_of_Molecule_SitesTS];
        IDon(Ix_Nascent) = 1;
        nascentInfo.Updated.IDon = IDon;
        nascentInfo.Updated.NumMolSites = Num_of_Molecule_Sites;
        nascentInfo.Updated.Names = Names;
        nascentInfo.Expression_TS = Exp_TS;
        nascentInfo.Updated.cMaxcTypes =  MaxcTypes;
        nascentInfo.Updated.cTypeIndx =  cTypeIndx;
        nascentInfo.Updated.ixTypeIndx = IxTypeIndx;
        nascentInfo.Updated.Kbmod = Kb_mod;
        nascentInfo.Updated.dHeqmod = dHeq_mod;
        nascentInfo.Updated.dSeqmod = dSeq_mod;
        nascentInfo.Updated.dHfmod = dHf_mod;
        nascentInfo.Updated.dSfmod = dSf_mod;
        nascentInfo.Updated.dHrmod = dHr_mod;
        nascentInfo.Updated.dSrmod = dSr_mod;
        nascentInfo.Updated.DoesProbeBindSite = DoesProbeBindSite2;
        try
            nascentInfo.maxPolymerase = maxPolymeraseOnTS;
        catch ME
            disp(ME.message)
        end
        try
            nascentInfo.numPolymerase = NumPolymeraseOnTS;
            nascentInfo.sense_info = TS_SenseInfo;
            nascentInfo.chr_info = TS_ChrInfo;
        catch ME
            disp(ME.message)
        end
        save([settings.FolderRootName filesep settings.GeneName '_Tm' num2str(T_hybrid) '_NascentInfo' NascentRecordID settings.designerName '.mat'],'nascentInfo','-v7.3','-append');

    end
catch ME
    disp(ME.message)
end
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