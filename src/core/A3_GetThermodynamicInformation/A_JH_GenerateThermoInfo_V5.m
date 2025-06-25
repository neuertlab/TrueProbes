function [Kb_Match,Kon,Koff,dHeq_Match,dSeq_Match,dHf_Match,dSf_Match,dHr_Match,dSr_Match,dCpeq_Match,dHon_eq,dSon_eq,dHon_f,dSon_f,dHon_r,dSon_r,dCpon_eq,Tm_on,Tm_Match] = A_JH_GenerateThermoInfo_V5(probes,gene_table,TranscriptName,settings)
% This Function takes probe sequences for on/off-target duplexes in in blast table
% and computes binding energy and rates associated with off-target and on-targets
N_methods = 8;
N_methods2 = 3;
N_methods3 = 9;
kb = 0.001987204259;%boltzman constant
T_hybrid = settings.HybridizationTemperature;
SaltConcentration = settings.SaltConcentration;
PrimerConcentration = settings.PrimerConcentration;
RemoveMisMatches = settings.RemoveMisMatches;
chromosome_ID = settings.chromosome_IDs;
transcript_ID = settings.transcript_IDs;
most_recent_num_local = settings.num_parpool_local;
designerName = settings.designerName;
FolderRootName = settings.FolderRootName;
probeBatchSize = settings.BLASTbatchSize;
targetBatchSize = settings.TargetBatchSize;
Organism = settings.Organism;
%% Identify location of on-target DNA and RNA molecules in list of DNA/RNA molecules
if (iscell(chromosome_ID))
    chromosome_IDon = cell(1,length(chromosome_ID));
    for v = 1:length(chromosome_ID)
        chromosome_IDon{v} = chromosome_ID{v};
    end
else
    chromosome_IDon{1} = chromosome_ID;
end
chromosome_IDon_temp = extractBefore(chromosome_IDon,'.');
if (sum(ismissing(chromosome_IDon_temp))>0)
    chromosome_IDon_temp(ismissing(chromosome_IDon_temp)) = extractBefore(chromosome_IDon(ismissing(chromosome_IDon_temp)),' ');
end
chromosome_IDon = chromosome_IDon_temp;clear chromosome_IDon_temp
if (iscell(transcript_ID))
    transcript_IDon = cell(1,length(transcript_ID));
    for v = 1:length(transcript_ID)
        transcript_IDon{v} = transcript_ID{v};
    end
else
    transcript_IDon{1} = transcript_ID;
end
transcript_IDon_temp = extractBefore(transcript_IDon,'.');
if (sum(ismissing(transcript_IDon_temp))>0)
    transcript_IDon_temp(ismissing(transcript_IDon_temp)) = extractBefore(transcript_IDon(ismissing(transcript_IDon_temp)),' ');
end
transcript_IDon = transcript_IDon_temp;clear transcript_IDon_temp
%% Getting Probe Sequences from reverse complement of tile sequences
RV = @(x) (seqrcomplement(x));
pi_seq = cell(1,size(probes,1));
for i=1:size(probes,1)
    pi_seq{i} = RV(probes{i,2});
end
N_Probes = size(probes,1);
N_Batches = ceil(N_Probes/probeBatchSize);
R = mod(N_Probes,probeBatchSize);
Batch = cell(1,N_Batches);
if (R==0)
    for k = 1:N_Batches
        Batch{k} = [probeBatchSize*(k-1)+1:probeBatchSize*k];
    end
else
    for k = 1:N_Batches-1
        Batch{k} = [probeBatchSize*(k-1)+1:probeBatchSize*k];
    end
    Batch{N_Batches} = [probeBatchSize*(N_Batches-1)+1:probeBatchSize*(N_Batches-1)+R];
end
%% Extract Gene Target Names
% if (isfile(strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_probetarget_flanking_sequences', settings.designerName,'.mat')));%save flanking sequence file if it does not exist already   
% load([strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_probetarget_flanking_sequences', settings.designerName,'.mat')],'probetarget_flanking_info');%save flanking sequence file if it does not exist already
% gene_table(:,'Alignment') = array2table(arrayfun(@(n) strcat(probetarget_flanking_info.ProbeRevCompSequence_5primeTo3prime{n}',newline,'|',newline,probetarget_flanking_info.TargetSequence_5primeTo3prime{n}')',1:size(probetarget_flanking_info,1),'Un',0)');
% %(link 1 is probe, line 3 is target)
% end
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
%either have seperate outpit fro mismatches, or switch to using flanking
%sequences if those files exist, thne add step to remove mismatches
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

%% Check if output file already exists
try
    if (settings.clusterStatus == 1)
        task_ID = getenv('SLURM_ARRAY_TASK_ID');
        if (isempty(task_ID))
            task_ID = 0;
        else
            task_ID = str2num(task_ID);
        end
        sleepRate = sum(1:task_ID);
        pause(5*sleepRate)
    end
catch
end
if (settings.clusterStatus)
    most_recent_num = str2num(getenv('SLURM_JOB_CPUS_PER_NODE'));
else
    most_recent_num = most_recent_num_local;
end
if (isfile([settings.FolderRootName filesep '(' TranscriptName ')_Tm' num2str(T_hybrid) '_RM' num2str(RemoveMisMatches) '_OnOffThermoInfo' settings.designerName '.mat']))
    calcOnOff = 0;
else
    calcOnOff = 1;
end
if (calcOnOff)
    fprintf("Get List of Unique Probe and Target Alignment Sequences")
    fprintf('\n')
    fprintf('\n')
    %gets sequences of probe and target in blast records
    targetMatch = arrayfun(@(x) strrep(gene_table.Alignment{x}(3,:),'-','N'),1:size(gene_table,1),'UniformOutput',false);%slow not so slow
    probeMatch = arrayfun(@(x) seqrcomplement(strrep(gene_table.Alignment{x}(1,:),'-','N')),1:size(gene_table,1),'UniformOutput',false);%slow
    %probeMatch is reverse complement
    %N is to deal with insertions in alignment matches, that 
    %Find Unique Pairs  (Do Mapping and Reverse Mapping)
    CalcDictionary = unique([targetMatch(:)' probeMatch(:)']);
    CalcDict2 = convertCharsToStrings(CalcDictionary);%1xN string
    indices = 1:length(CalcDict2);
    idMap = dictionary(CalcDict2,indices);
    CalcPairs(:,1) = idMap(targetMatch);
    CalcPairs(:,2) = idMap(probeMatch);
    uniqueCalcPairs = unique(CalcPairs,'rows');
    targetMatch_Unique = {CalcDictionary{uniqueCalcPairs(:,1)}};
    probeMatch_Unique = {CalcDictionary{uniqueCalcPairs(:,2)}};
    NamesZ = gene_table.Name;
    ProbeNumZ = gene_table.ProbeNum;
    MatchZ = gene_table.Match;
    clear gene_table CalcDictionary CalcDict2 idMap indices
    %Re-Map Energies to non-unique sites
    %Unique_to_Original = @(x) find(strcmpi(targetMatch_Unique{x},targetMatch).*strcmpi(probeMatch_Unique{x},probeMatch));  
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
    

% parameter_file_canonical = 'src/thirdparty/VarGibbs-4.1/AOP-SL96.par';
% parameter_file_canonical = 'src/thirdparty/VarGibbs-4.1/AOP-SL98.par';
% parameter_file_canonical = 'src/thirdparty/VarGibbs-4.1/P-SL04.par';
% parameter_file_canonical = 'src/thirdparty/VarGibbs-4.1/P-SL98.par';
% parameter_file_canonical = 'src/thirdparty/VarGibbs-4.1/P-SL96.par';
% parameter_file_canonical = 'src/thirdparty/VarGibbs-4.1/P-CH12.par';
% %DNA-RNA
% parameter_file_canonical = 'src/thirdparty/VarGibbs-4.1/AOP-DRVH.par';
% parameter_file_canonical = 'src/thirdparty/VarGibbs-4.1/AOP-DRFT.par';
% parameter_file_canonical = 'src/thirdparty/VarGibbs-4.1/AOP-DRLS.par';
% parameter_file_canonical = 'src/thirdparty/VarGibbs-4.1/P-BN21.par';
%parameter_file_noncanonical = 'src/thirdparty/VarGibbs-4.1/AOP-MM-60.par';
%parameter_file_canonical = 'src/thirdparty/VarGibbs-4.1/AOP-OW04-69.par';
     if (isfile(strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_probetarget_flanking_sequences', settings.designerName,'.mat')))
     %   Alignment_Target_Sequences = probetarget_flanking_info.TargetSequence_5primeTo3prime;
      %  Alignment_Probe_Sequences = probetarget_flanking_info.ProbeRevCompSequence_5primeTo3prime;
     end

    Kon = zeros(size(probes,1),N_methods);
    Tm_on = zeros(size(probes,1),N_methods3);
    dHeq_Match = zeros(length(targetMatch),N_methods);
    dSeq_Match = zeros(length(targetMatch),N_methods);
    dCpeq_Match = zeros(length(targetMatch),N_methods);
    dHf_Match = zeros(length(targetMatch),N_methods2);
    dSf_Match = zeros(length(targetMatch),N_methods2);
    dHr_Match = zeros(length(targetMatch),N_methods2);
    dSr_Match = zeros(length(targetMatch),N_methods2);
    dHon_eq = zeros(size(probes,1),N_methods);
    dSon_eq = zeros(size(probes,1),N_methods);
    dCpon_eq = zeros(size(probes,1),N_methods);
    dHon_f = zeros(size(probes,1),3);
    dSon_f = zeros(size(probes,1),3);
    dHon_r = zeros(size(probes,1),3);
    dSon_r = zeros(size(probes,1),3);
    Tm_Match = zeros(length(targetMatch),N_methods3);
    Kb_Match = zeros(length(targetMatch),N_methods);
    try
        n_cores = str2num(getenv('SLURM_JOB_CPUS_PER_NODE'));
        p = parpool(n_cores);
        %p.IdleTimeout = Inf;
        spmd
            warning('off','all')
        end
    catch
        %p = parpool;
        %p.IdleTimeout = Inf;
    end
    fprintf("Compute On-Target Binding Affinity for All Probes")
    fprintf('\n')
    fprintf('\n')
    wb = parwaitbar(size(probes,1),'WaitMessage', 'Computing');
    parfor i=1:size(probes,1)
        pause(0.1);
        [temp_dHeq{i}, temp_dSeq{i}, temp_dGeq{i}, ...
            temp_dHf{i}, temp_dSf{i}, ~, ...
            temp_dHr{i}, temp_dSr{i}, ~,temp_dCpeq{i},temp_Tm{i}] = F_DeltaGibson_V3(lower(pi_seq{i}),seqrcomplement(lower(pi_seq{i})),SaltConcentration,T_hybrid,PrimerConcentration,sequence_duplexes_thermo_generator_structure);
        dHon_eq(i,:) = temp_dHeq{i};
        dSon_eq(i,:) = temp_dSeq{i};
        dHon_f(i,:) = temp_dHf{i};
        dSon_f(i,:) = temp_dSf{i};
        dHon_r(i,:) = temp_dHr{i};
        dSon_r(i,:) = temp_dSr{i};
        dCpon_eq(i,:) = temp_dCpeq{i};
        Tm_on(i,:) = temp_Tm{i};
        Kon(i,:) = exp(-temp_dGeq{i}/(kb*(T_hybrid+273.15))); %Keq rates on-Target binding
        temp_dHeq{i} = [];
        temp_dSeq{i} = [];
        temp_dGeq{i} = [];
        temp_dHf{i} = [];
        temp_dSf{i} = [];
        temp_dHr{i} = [];
        temp_dSr{i} = [];
        temp_dCpeq{i} = [];
        temp_Tm{i} = [];
        progress(wb);
    end
    wb.delete();
    fprintf('\n')
    fprintf('\n')
    clear temp_dHeq temp_dSeq temp_dGeq temp_dHf temp_dSf temp_dHr temp_dSr temp_dCpeq temp_Tm
    %compute binding thermodynamic/kinetic information for on-target hits parallelized
    if (~isfile([settings.FolderRootName filesep '(' TranscriptName ')_' settings.rootName '_dCpInfo' settings.designerName '.mat']))%check if temp file exists
        N_uniqueCalcPairsBatches = ceil(size(uniqueCalcPairs,1)/targetBatchSize);
        R = mod(size(uniqueCalcPairs,1),targetBatchSize);
        Batch_uniqueCalcPairs = cell(1,N_uniqueCalcPairsBatches);
        if (R==0)
            for k = 1:N_uniqueCalcPairsBatches
                Batch_uniqueCalcPairs{k} = [targetBatchSize*(k-1)+1:targetBatchSize*k];
            end
        else
            for k = 1:N_uniqueCalcPairsBatches-1
                Batch_uniqueCalcPairs{k} = [targetBatchSize*(k-1)+1:targetBatchSize*k];
            end
            Batch_uniqueCalcPairs{N_uniqueCalcPairsBatches} = [targetBatchSize*(N_uniqueCalcPairsBatches-1)+1:targetBatchSize*(N_uniqueCalcPairsBatches-1)+R];
        end
        ResultsExist = zeros(1,N_uniqueCalcPairsBatches);
        ResultsDate = cell(1,N_uniqueCalcPairsBatches);
        fprintf("Check if probe batch gibbs free energy files exist")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_uniqueCalcPairsBatches,'WaitMessage', 'Checking');
        parfor i = 1:N_uniqueCalcPairsBatches
            if (isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_deltaGibsonInfo_batch' num2str(i) '.mat']))%check if temp file exists
                d = dir([FolderRootName filesep '(' TranscriptName ')' designerName '_deltaGibsonInfo_batch' num2str(i) '.mat']);
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
        %PARTIALLY SAVED compute binding thermodynamic/kinetic information for unique off-target hit sequences parallelized
        if (~isempty(batch_nums_to_check))
            targetMatch_Unique_Constant = parallel.pool.Constant(targetMatch_Unique);
            probeMatch_Unique_Constant = parallel.pool.Constant(probeMatch_Unique);
            Batch_uniqueCalcPairs_Constant = parallel.pool.Constant(Batch_uniqueCalcPairs);
            fprintf("Computing Binding Affinity for Batches of Unique Probe Target Sequence Pairs")
            fprintf('\n')
            fprintf('\n')
            wb = parwaitbar(length(batch_nums_to_check),'WaitMessage','Computing');
            parfor v = 1:length(batch_nums_to_check)
                pause(0.1)
                uniquePair_temp_dHeq{v} = cell(1,length(Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)}));
                uniquePair_temp_dSeq{v} = cell(1,length(Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)}));
                uniquePair_temp_dGeq{v} = cell(1,length(Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)}));
                uniquePair_temp_dHf{v} = cell(1,length(Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)}));
                uniquePair_temp_dSf{v} = cell(1,length(Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)}));
                uniquePair_temp_dHr{v} = cell(1,length(Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)}));
                uniquePair_temp_dSr{v} = cell(1,length(Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)}));
                uniquePair_temp_dCpeq{v} = cell(1,length(Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)}));
                uniquePair_temp_Tm{v} = cell(1,length(Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)}));
                temp_indx{v} = cell(1,length(Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)}));
                temp_GibsonInfo{v} = cell(1,length(Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)}));
                for w = 1:length(Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)})
                    %i = Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)}(w);
                    [uniquePair_temp_dHeq{v}{w}, uniquePair_temp_dSeq{v}{w}, uniquePair_temp_dGeq{v}{w}, ...
                        uniquePair_temp_dHf{v}{w}, uniquePair_temp_dSf{v}{w}, ~, ...
                        uniquePair_temp_dHr{v}{w}, uniquePair_temp_dSr{v}{w}, ~,uniquePair_temp_dCpeq{v}{w},uniquePair_temp_Tm{v}{w}] = ...
                        F_DeltaGibson_V3(targetMatch_Unique_Constant.Value{Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)}(w)},probeMatch_Unique_Constant.Value{Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)}(w)},SaltConcentration,T_hybrid,PrimerConcentration);
                    temp_indx{v}{w}= find(strcmpi(targetMatch_Unique_Constant.Value{Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)}(w)},targetMatch).*strcmpi(probeMatch_Unique_Constant.Value{Batch_uniqueCalcPairs_Constant.Value{batch_nums_to_check(v)}(w)},probeMatch));
                    temp_GibsonInfo{v}{w} = {temp_indx{v}{w};uniquePair_temp_dGeq{v}{w};uniquePair_temp_dHeq{v}{w}';uniquePair_temp_dSeq{v}{w}';...
                        uniquePair_temp_dHf{v}{w}';uniquePair_temp_dSf{v}{w}';...
                        uniquePair_temp_dHr{v}{w}';uniquePair_temp_dSr{v}{w}';...
                        uniquePair_temp_Tm{v}{w}';uniquePair_temp_dCpeq{v}{w}'};
                end
                parsave_partial_delta_gibson([FolderRootName filesep '(' TranscriptName ')' designerName '_deltaGibsonInfo_batch' num2str(batch_nums_to_check(v)) '.mat'],temp_GibsonInfo{v})
                temp_indx{v} = [];
                uniquePair_temp_dHeq{v} = [];
                uniquePair_temp_dSeq{v} = [];
                uniquePair_temp_dGeq{v} = [];
                uniquePair_temp_dHf{v} = [];
                uniquePair_temp_dSf{v} = [];
                uniquePair_temp_dHr{v} = [];
                uniquePair_temp_dSr{v} = [];
                uniquePair_temp_dCpeq{v} = [];
                uniquePair_temp_Tm{v} = [];
                progress(wb);
            end
            wb.delete();
            fprintf('\n')
            fprintf('\n')
        end
        fprintf("Aggregating probe target batch thermodynamic information")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_uniqueCalcPairsBatches,'WaitMessage','Aggregating');
        for v = 1:N_uniqueCalcPairsBatches
            if isfile([settings.FolderRootName filesep '(' TranscriptName ')' designerName '_deltaGibsonInfo_batch' num2str(v) '.mat'])
                load([settings.FolderRootName filesep '(' TranscriptName ')' designerName '_deltaGibsonInfo_batch' num2str(v) '.mat'],'partial_delta_gibson_tmp');
                for w = 1:length(Batch_uniqueCalcPairs{v})
                    Kb_Match(partial_delta_gibson_tmp{w}{1},:) = repmat(partial_delta_gibson_tmp{w}{2}'/(kb*(T_hybrid+273.15)),[length(partial_delta_gibson_tmp{w}{1}) 1]);
                    dHeq_Match(partial_delta_gibson_tmp{w}{1},:) = repmat(partial_delta_gibson_tmp{w}{3},[length(partial_delta_gibson_tmp{w}{1}) 1]);
                    dSeq_Match(partial_delta_gibson_tmp{w}{1},:) = repmat(partial_delta_gibson_tmp{w}{4},[length(partial_delta_gibson_tmp{w}{1}) 1]);
                    dHf_Match(partial_delta_gibson_tmp{w}{1},:) = repmat(partial_delta_gibson_tmp{w}{5},[length(partial_delta_gibson_tmp{w}{1}) 1]);
                    dSf_Match(partial_delta_gibson_tmp{w}{1},:) = repmat(partial_delta_gibson_tmp{w}{6},[length(partial_delta_gibson_tmp{w}{1}) 1]);
                    dHr_Match(partial_delta_gibson_tmp{w}{1},:) = repmat(partial_delta_gibson_tmp{w}{7},[length(partial_delta_gibson_tmp{w}{1}) 1]);
                    dSr_Match(partial_delta_gibson_tmp{w}{1},:) = repmat(partial_delta_gibson_tmp{w}{8},[length(partial_delta_gibson_tmp{w}{1}) 1]);
                    Tm_Match(partial_delta_gibson_tmp{w}{1},:) = repmat(partial_delta_gibson_tmp{w}{9}',[length(partial_delta_gibson_tmp{w}{1}) 1]);
                    dCpeq_Match(partial_delta_gibson_tmp{w}{1},:) = repmat(partial_delta_gibson_tmp{w}{10},[length(partial_delta_gibson_tmp{w}{1}) 1]);
                end
                clear partial_delta_gibson_tmp
            end
            progress(wb);
        end
        wb.delete();



        
        save([settings.FolderRootName filesep '(' TranscriptName ')_' settings.rootName '_dHInfo' settings.designerName '.mat'],'dHon_f','dHon_r','dHon_eq','dHeq_Match','dHf_Match','dHr_Match','-v7.3');
        save([settings.FolderRootName filesep '(' TranscriptName ')_' settings.rootName '_dSInfo' settings.designerName '.mat'],'dSon_f','dSon_r','dSon_eq','dSeq_Match','dSf_Match','dSr_Match','-v7.3');
        save([settings.FolderRootName filesep '(' TranscriptName ')_' settings.rootName '_TmInfo' settings.designerName '.mat'],'Tm_on','Tm_Match','-v7.3');
        save([settings.FolderRootName filesep '(' TranscriptName ')_' settings.rootName '_dCpInfo' settings.designerName '.mat'],'dCpon_eq','dCpeq_Match','-v7.3');
        save([settings.FolderRootName filesep '(' TranscriptName ')_' settings.rootName '_dKbInfo' settings.designerName '.mat'],'Kb_Match','-v7.3');
        fprintf('\n')
        fprintf('\n')
        fprintf("Deleting temporary gibbs energy files")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_uniqueCalcPairsBatches,'WaitMessage', 'Deleting');
        parfor i = 1:N_uniqueCalcPairsBatches
            if exist([FolderRootName filesep '(' TranscriptName ')' designerName '_deltaGibsonInfo_batch' num2str(i) '.mat'],'file')        %delete temp mat file if already exists
                delete([FolderRootName filesep '(' TranscriptName ')' designerName '_deltaGibsonInfo_batch' num2str(i) '.mat'])
            end
            progress(wb);
        end
        wb.delete();
    else
    end
    load([settings.FolderRootName filesep '(' TranscriptName ')_' settings.rootName '_dKbInfo' settings.designerName '.mat'],'Kb_Match');
    %store results for on/off-target binding information matches to unique sequence pairs
    indices = 1:length(uniNames);
    idMap2 = dictionary(uniNames',indices);
    if (and(strcmp(settings.referenceType,"ENSEMBL"),max(double(contains(extractBefore(Names,' '),'ENS')))==0))
        uniNamesZ = extractBefore(NamesZ,' ');
    else
        uniNamesZ = extractBefore(NamesZ,'.');
        if (sum(ismissing(uniNamesZ))>0)
            uniNamesZ(ismissing(uniNamesZ)) = extractBefore(NamesZ(ismissing(uniNamesZ)),' ');
        end
    end
    if (~isfile([settings.FolderRootName filesep '(' TranscriptName ')_' settings.rootName '_OnOffInfo' settings.designerName '.mat']))%check if temp file exists
        ResultsExist_OnOff_1D = zeros(1,N_Batches);
        ResultsDate_OnOff_1D = cell(1,N_Batches);
        fprintf('\n')
        fprintf('\n')
        fprintf("Check if on/off-target thermo batch files exist")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(N_Batches,'WaitMessage', 'Checking');
        parfor i = 1:N_Batches
            if (isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_OnOffInfo_batch' num2str(i) '.mat']))%check if temp file exists
                d = dir([FolderRootName filesep '(' TranscriptName ')' designerName '_OnOffInfo_batch' num2str(i) '.mat']);
                if (d.bytes>0)%check size greater than zero
                    ResultsExist_OnOff_1D(i) = 1;
                end
                ResultsDate_OnOff_1D{i} = datetime(d.date);
            end
            progress(wb);
        end
        wb.delete();
        fprintf('\n')
        fprintf('\n')
        Results_NotMade_OnOff_1D = find(ResultsExist_OnOff_1D ==0);
        Results_Made_OnOff_1D = find(ResultsExist_OnOff_1D ==1);
        %Sort get most recent ResultsMade GeneHitsMade and GeneHitsTable Made and add to probe_check_list
        if (length(Results_Made_OnOff_1D)<=most_recent_num)
            results_check1_OnOff_1D = Results_Made_OnOff_1D;
        else
            Results_RecentMade_Dates_OnOff_1D(:,1) = ResultsDate_OnOff_1D(Results_Made_OnOff_1D);
            Results_RecentMade_Dates_OnOff_1D(:,2) = num2cell(Results_Made_OnOff_1D);
            Results_RecentMade_Dates_OnOff_1D  = table2timetable(cell2table(Results_RecentMade_Dates_OnOff_1D));
            Results_RecentMade_Dates_OnOff_1D  = sortrows(Results_RecentMade_Dates_OnOff_1D,1,'descend');
            Results_RecentMade_Dates_OnOff_1D.Properties.VariableNames = {'ID'};
            results_check1_OnOff_1D = Results_RecentMade_Dates_OnOff_1D .ID(1:most_recent_num).';
            clear Results_RecentMade_Dates_OnOff_1D
        end
        batch_nums_to_check_OnOff_1D = union(Results_NotMade_OnOff_1D,results_check1_OnOff_1D);
        probes_List = parallel.pool.Constant(probes);
        idMap2_List = parallel.pool.Constant(idMap2);
        uniNames_List = parallel.pool.Constant(uniNames);
        uniNamesZ_List = parallel.pool.Constant(uniNamesZ);
        NamesZ_List = parallel.pool.Constant(NamesZ);
        Batch_List = parallel.pool.Constant(Batch);
        %build structure with on/off-target binding by each target in blast results.
        fprintf("Computing 1D Probe Batch On/Off-Target affinity index values")
        fprintf('\n')
        fprintf('\n')
        wb = parwaitbar(length(batch_nums_to_check_OnOff_1D),'WaitMessage', 'Computing');
        parfor v=1:length(batch_nums_to_check_OnOff_1D)
            pause(0.1);
            temp_on_off{v} = cell(1,length(Batch_List.Value{batch_nums_to_check_OnOff_1D(v)}));
            for w = 1:length(Batch_List.Value{batch_nums_to_check_OnOff_1D(v)})
                %p = Batch_List.Value{batch_nums_to_check_OnOff_1D(v)}(w);
                IDy = find(ProbeNumZ==Batch_List.Value{batch_nums_to_check_OnOff_1D(v)}(w));%add step to remove Off-target from chromosome
                IDs = find(MatchZ<length(probes_List.Value{Batch_List.Value{batch_nums_to_check_OnOff_1D(v)}(w),2}));
                IDf = find(MatchZ==length(probes_List.Value{Batch_List.Value{batch_nums_to_check_OnOff_1D(v)}(w),2}));
                Indx_OFF = intersect(IDy,IDs);
                Indx_ON = intersect(IDy,IDf);%contains off hits
                NamesZ_OFF = NamesZ_List.Value(Indx_OFF);
                trimmed_NamesZ_OFF = extractBefore(NamesZ_OFF,'.');
                if (sum(ismissing(trimmed_NamesZ_OFF))>0)
                    trimmed_NamesZ_OFF(ismissing(trimmed_NamesZ_OFF)) = extractBefore(NamesZ_OFF(ismissing(NamesZ_OFF)),' ');
                end
                Names_OFF = unique(trimmed_NamesZ_OFF);
                NamesZ_ON = NamesZ_List.Value(Indx_ON);
                trimmed_NamesZ_ON = extractBefore(NamesZ_ON,'.');
                if (sum(ismissing(trimmed_NamesZ_ON))>0)
                    trimmed_NamesZ_ON(ismissing(trimmed_NamesZ_ON)) = extractBefore(NamesZ_ON(ismissing(trimmed_NamesZ_ON)),' ');
                end
                Names_ON = unique(trimmed_NamesZ_ON);
                Names_p = union(Names_OFF,Names_ON);
                OFF_uniNamesZ =  uniNamesZ_List.Value(Indx_OFF);
                ON_uniNamesZ =  uniNamesZ_List.Value(Indx_ON);
                IDon_temp = idMap2_List.Value(Names_ON);
                IDoff_temp = idMap2_List.Value(Names_OFF);
                actual_IDon_DNA = cell2mat(arrayfun(@(x) strcmp(Names_ON{x},chromosome_IDon),1:size(Names_ON,1),'UniformOutput',false));
                actual_IDon_RNA = cell2mat(arrayfun(@(x) strcmp(Names_ON{x},transcript_IDon),1:size(Names_ON,1),'UniformOutput',false));
                actual_IDon = find(actual_IDon_DNA+actual_IDon_RNA>0);
                false_IDon = find(actual_IDon_DNA+actual_IDon_RNA==0);
                IDon = IDon_temp(actual_IDon);
                IDoff = union(IDoff_temp,IDon_temp(false_IDon));
                IDshared = intersect(IDon,IDoff);
                SharedName = {uniNames_List.Value{IDshared}};
                falseNames = Names_ON(false_IDon);
                Indx_ON_False = Indx_ON(find(cell2mat(arrayfun(@(x) sum(strcmp(falseNames,ON_uniNamesZ(x))),1:length(ON_uniNamesZ),'Un',0))>0));
                t1_vector = arrayfun(@(t) find(strcmp(uniNames,Names_p{t}))',1:length(Names_p),'Un',0);
                Index_OFF1_vector = arrayfun(@(t) Indx_OFF(find(strcmp(OFF_uniNamesZ,Names_p{t})))',1:length(Names_p),'Un',0);
                Index_OFF2_vector = arrayfun(@(t) Indx_ON_False(find(strcmp(falseNames,Names_p{t})))',1:length(Names_p),'Un',0);
                temp_on_off{v}{w} = {Names_p, IDshared,ON_uniNamesZ,SharedName,Indx_ON,t1_vector,Index_OFF1_vector,Index_OFF2_vector};
            end
            parsave_partial_on_off([FolderRootName filesep '(' TranscriptName ')' designerName '_OnOffInfo_batch' num2str(batch_nums_to_check_OnOff_1D(v)) '.mat'],temp_on_off{v});
            temp_on_off{v} = [];
            progress(wb);
        end
        wb.delete();
        fprintf('\n')
        fprintf('\n')
    end
    fprintf("Aggregating 1D On/Off-Target Simple affinity index values")
    fprintf('\n')
    fprintf('\n')
    Koff = cell(1,N_Batches);
    Kb_Match_List = parallel.pool.Constant(Kb_Match);
    wb = parwaitbar(N_Batches,'WaitMessage', 'Aggregating');
    parfor v=1:N_Batches
        pause(0.1);
        if isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_OnOffInfo_batch' num2str(v) '.mat'])
            partial_on_off_tmp = load([FolderRootName filesep '(' TranscriptName ')' designerName '_OnOffInfo_batch' num2str(v) '.mat']).partial_on_off_tmp;
            Koff{v} = ndSparse.build([length(Batch{v}) length(Names) N_methods],0);
            for w = 1:length(Batch{v})
                Names_p = partial_on_off_tmp{w}{1};
                IDshared = partial_on_off_tmp{w}{2};
                ON_uniNamesZ = partial_on_off_tmp{w}{3};
                SharedName = partial_on_off_tmp{w}{4};
                Indx_ON = partial_on_off_tmp{w}{5};
                t1_vector = partial_on_off_tmp{w}{6};
                Index_OFF1_vector = partial_on_off_tmp{w}{7};
                Index_OFF2_vector = partial_on_off_tmp{w}{8};
                if (isempty(IDshared))
                    Koff{v}(w,cell2mat(t1_vector(1:length(Names_p))),:) = CATnWrapper(arrayfun(@(t) sum(Kb_Match_List.Value(Index_OFF1_vector{t},:),1),1:length(Names_p),'Un',0),1);
                else
                    is_shared = cell2mat(arrayfun(@(t) ~ismember(Names_p{t},uniNames(IDshared)),1:length(Names_p),'Un',0));
                    t_shared = find(is_shared);
                    t_not_shared = find(is_shared==0);
                    if (~isempty(t_shared))
                        Koff{v}(w,cell2mat(t1_vector(t_shared)),:) = CATnWrapper(arrayfun(@(t) sum(Kb_Match_List.Value(Index_OFF1_vector{t},:),1)+sum(Kb_Match_List.Value(Index_OFF2_vector{t},:),1),t_shared,'Un',0),1);
                    end
                    if (~isempty(t_not_shared))
                        Indx_OFF_3 = Indx_ON(find(cell2mat(arrayfun(@(x) sum(strcmp(ON_uniNamesZ,SharedName(x))),1:length(SharedName),'Un',0))>0));
                        % remove first ON event for not shared
                        if (length(Indx_OFF_3)>1)
                            Indx_OFF_3 = Indx_OFF_3(2:end);
                            Koff{v}(w,cell2mat(t1_vector(t_not_shared)),:) = CATnWrapper(arrayfun(@(t) sum(Kb_Match_List.Value(Index_OFF1_vector{t},:),1)+sum(Kb_Match_List.Value(Index_OFF2_vector{t},:),1),t_not_shared,'Un',0),1)+sum(Kb_Match_List.Value(Indx_OFF_3,:),1);
                        else
                            Koff{v}(w,cell2mat(t1_vector(t_not_shared)),:) = CATnWrapper(arrayfun(@(t) sum(Kb_Match_List.Value(Index_OFF1_vector{t},:),1)+sum(Kb_Match_List.Value(Index_OFF2_vector{t},:),1),t_not_shared,'Un',0),1);
                        end
                    end
                end
            end
        end
        progress(wb);
    end
    wb.delete();
    Koff = CATnWrapper(Koff,1);
    save([settings.FolderRootName filesep '(' TranscriptName ')_' settings.rootName '_Tm' num2str(settings.HybridizationTemperature) '_OnOffThermoInfo' settings.designerName '.mat'],'Kon','Koff','Kb_Match','-v7.3');
    fprintf('\n')
    fprintf('\n')
    fprintf("Deleting temporary 1D On/Off-Target Simple affinity index value files")
    fprintf('\n')
    fprintf('\n')
    wb = parwaitbar(N_Batches,'WaitMessage','Deleting');
    parfor i = 1:N_Batches
        if exist([FolderRootName filesep '(' TranscriptName ')' designerName '_OnOffInfo_batch' num2str(i) '.mat'],'file')        %delete temp mat file if already exists
            delete([FolderRootName filesep '(' TranscriptName ')' designerName '_OnOffInfo_batch' num2str(i) '.mat'])
        end
        progress(wb);
    end
    wb.delete();
    fprintf('\n')
    fprintf('\n')
end
end
