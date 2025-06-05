function [probes,gene_table] = Probe_checker_general_JH10(probes,Organism,TranscriptName,TranscriptName2,Chromosome,settings)
%% This searches a target genome for off-target matches to a set of probes
%Use "Import data" in the Home tab, select excel file and the column with the probes you
%want to check, click on "cell array", and then import data. change what
%probes = to match the title (line 26), and then run.

%if running on the cluster, you will need to run lines 27 and 28 to save the probes that will be used. Then upload
%the correct file (TranscriptName_probes) to the cluster where this code is.

%These names will exclude the top fit for the chromosome and every fit
%within the same transcript when determining off-target binding
%tic
% [probes,tent_probes,seqs,gene_names,organisms] = BK_Probe_Generator({},{'AREG.txt'},{});
% %probes = {};
% Organism = 'Human'              %'Mouse' 'Cerevisiae' 'Human'
% TranscriptName = '(AREG)'                %BLAST you transcript name and then match exactly before running
% TranscriptName2 = '(AREG)'                %BLAST you transcript name and then match exactly before running
% Chromosome = '4'                      %Specify which chromosome (Just say the number or letter). Put a comma after if single digit

% probes = PossibleRepAprobes;            %Set this equal to the imported file name
% save([TranscriptName '_probes'],'probes')
db1 = [];
db2 = [];
batchSize = settings.BLASTbatchSize;
simultaneous = settings.BLASTsimultaneousParsingOverSequentialParsing;
runDNA = settings.BLASTdna;
runRNA = settings.BLASTrna;
reward = settings.BlastParameters.reward;
penalty = settings.BlastParameters.penalty;
word_size = settings.BlastParameters.wordsize;
gapopen = settings.BlastParameters.gapopen;
gapextend = settings.BlastParameters.gapextend;
evalue = settings.BlastParameters.evalue;
num_alignments = settings.BlastParameters.num_alignments;
dust = settings.BlastParameters.dust;
runDNA
runRNA
if strcmp(Organism,'Mouse')
    if (runDNA)
        db1 = settings.mLoc;
    end
    if (runRNA)
        db2 = settings.mLoc2;
    end
elseif strcmp(Organism, 'Human')
    if (runDNA)
        db1 = settings.hLoc;
    end
    if (runRNA)
        db2 = settings.hLoc2;
    end
elseif strcmp(Organism, 'Yeast')
    if (runDNA)
        db1 = settings.yLoc;
    end
    if (runRNA)
        db2 = settings.yLoc2;
    end
else
    if (runDNA)
        db1 = settings.specificBlastDatabase;
    end
    if (runRNA)
        db2 = settings.specificBlastDatabase2;
    end
end


%probe_nums = 1:size(probes,1);
% each row is information for a probe, then within the cell array in the second column, rows are genes, columns are number of hits
first_loaded = 0;
fail_nums = zeros(size(probes,1));
FolderRootName = settings.FolderRootName;
designerName = settings.designerName;
MaxProbeSize = settings.MaxProbeSize;

% generate batches of probes to blast
N_Probes = size(probes,1);
N_Batches = ceil(N_Probes/batchSize);
batch_nums = 1:N_Batches;
R = mod(N_Probes,batchSize);
Batch = cell(1,N_Batches);

if (R==0)
    for k = 1:N_Batches
        Batch{k} = [batchSize*(k-1)+1:batchSize*k];
    end
else
    for k = 1:N_Batches-1
        Batch{k} = [batchSize*(k-1)+1:batchSize*k];
    end
    Batch{N_Batches} = [batchSize*(N_Batches-1)+1:batchSize*(N_Batches-1)+R];
end
try
    pc = parcluster('local');
    %pc.JobStorageLocation = strcat(getenv('SCRATCH'),filesep,getenv('SLURM_JOB_ID'));
    parpool(pc,str2num(getenv('SLURM_JOB_CPUS_PER_NODE')));
    spmd
        warning('off','all')
    end
catch
end
ResultsExist = zeros(1,N_Batches);
ResultsSize = zeros(1,N_Batches);
ResultsDate = cell(1,N_Batches);
ResultsExist2 = zeros(1,N_Batches);
ResultsSize2 = zeros(1,N_Batches);
ResultsDate2 = cell(1,N_Batches);
GeneHitsTableExist = zeros(1,N_Batches);
GeneHitsTableSize = zeros(1,N_Batches);
GeneHitsTableDate = cell(1,N_Batches);
%sort on size and date time (remove last 8 in time) if any exist

fprintf('\n')
fprintf('\n')
fprintf("Check if probe batch BLAST report xml files exist")
fprintf('\n')
fprintf('\n')
progBar = ProgressBar(N_Batches);
for i = 1:N_Batches
    if (isfile([FolderRootName filesep TranscriptName settings.designerName 'probebatch' num2str(i) 'resultsDNA.xml']))%check if XML file exists
        d = dir([FolderRootName filesep TranscriptName settings.designerName 'probebatch' num2str(i) 'resultsDNA.xml']);
        if (d.bytes>0)%check size greater than zero
            ResultsExist(i) = 1;
        end
        ResultsSize(i) = d.bytes;
        ResultsDate{i} = datetime(d.date);
        clear d
    end
    if (isfile([FolderRootName filesep TranscriptName settings.designerName 'probebatch' num2str(i) 'resultsRNA.xml']))%check if XML file exists
        d = dir([FolderRootName filesep TranscriptName settings.designerName 'probebatch' num2str(i) 'resultsRNA.xml']);
        if (d.bytes>0)%check size greater than zero
            ResultsExist2(i) = 1;
        end
        ResultsSize2(i) = d.bytes;
        ResultsDate2{i} = datetime(d.date);
        clear d
    end
    if (isfile([FolderRootName filesep TranscriptName settings.designerName 'gene_hits_table_batch' num2str(i) '.mat']))%check if temp file exists
        d = dir([FolderRootName filesep TranscriptName settings.designerName 'gene_hits_table_batch' num2str(i) '.mat']);
        if (d.bytes>0)%check size greater than zero
            GeneHitsTableExist(i) = 1;
        end
        GeneHitsTableSize(i) = d.bytes;
        GeneHitsTableDate{i} = datetime(d.date);
        clear d
    end
    progBar([],[],[]);
end
progBar.release();

Results_NotMade = find(ResultsExist==0);
Results_Made = find(ResultsExist==1);
Results_NotMade2 = find(ResultsExist2==0);
Results_Made2 = find(ResultsExist2==1);
GeneHitsTable_NotMade = find(GeneHitsTableExist==0);
GeneHitsTable_Made = find(GeneHitsTableExist==1);
probes_to_check1 = union(Results_NotMade,Results_NotMade2);
probes_to_check2 = GeneHitsTable_NotMade;
%Sort get most 8 recent ResultsMade GeneHitsMade and GeneHitsTable Made and
%add to probe_check_list
if (settings.clusterStatus)
    most_recent_num = str2num(getenv('SLURM_JOB_CPUS_PER_NODE'));
else
    most_recent_num = 8;
end

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
if (length(GeneHitsTable_Made)<=most_recent_num)
    results_check4 = GeneHitsTable_Made;
else
    GeneHitsTable_RecentMade_Dates(:,1) = GeneHitsTableDate(GeneHitsTable_Made);
    GeneHitsTable_RecentMade_Dates(:,2) = num2cell(GeneHitsTable_Made);
    GeneHitsTable_RecentMade_Dates = table2timetable(cell2table(GeneHitsTable_RecentMade_Dates));
    GeneHitsTable_RecentMade_Dates = sortrows(GeneHitsTable_RecentMade_Dates,1,'descend');
    GeneHitsTable_RecentMade_Dates.Properties.VariableNames = {'ID'};
    results_check4 = GeneHitsTable_RecentMade_Dates.ID(1:most_recent_num).';
    clear GeneHitsTable_RecentMade_Dates
end

probes_to_check3 = union(results_check1,results_check2);
probes_to_check4 = results_check4;
batch_nums_to_check1 = union(probes_to_check1,probes_to_check3);
batch_nums_to_check2 = union(probes_to_check2,probes_to_check4);

%BatchMaxLen = zeros(1,N_Batches);
fprintf('\n')
fprintf('\n')
fprintf("Generating probe batch fasta files")
fprintf('\n')
fprintf('\n')
progBar = ProgressBar(N_Batches);
for i = 1:N_Batches
    data = [];
    %% For local blast on server
    if exist([FolderRootName filesep TranscriptName settings.designerName 'probebatch' num2str(i) '.fa'],'file')        %delete fasta if already exists
        delete([FolderRootName filesep TranscriptName settings.designerName 'probebatch' num2str(i) '.fa'])
    end
    for v = 1:length(Batch{i})
        data(v).Sequence = probes{Batch{i}(v),2};  %gets the sequence
        data(v).Header = ['p' num2str(Batch{i}(v))];
    end
    fastawrite([FolderRootName filesep TranscriptName settings.designerName 'probebatch' num2str(i) '.fa'],data);
    % BatchMaxLen(i) = max(cellfun(@length,probes(Batch{i},2)));
    progBar([],[],[]);
end
progBar.release();

N_Batches
batch_nums_to_check1
MinHomologySize_Constant = parallel.pool.Constant(settings.MinHomologySearchTargetSize);
%if no database can make database
%how to check if
if (~isempty(batch_nums_to_check1))
fprintf('\n')
fprintf('\n')
fprintf("blasting probe batch fasta files")
fprintf('\n')
fprintf('\n')
progBar = ProgressBar(length(batch_nums_to_check1),'IsParallel', true,'WorkerDirectory', pwd(),'Title', 'Blasting');
progBar.setup([],[],[]);
parfor v = 1:length(batch_nums_to_check1)
    pause(0.1);
    i = batch_nums_to_check1(v);
    if (runDNA)
        sequence_data = [FolderRootName filesep TranscriptName designerName 'probebatch' num2str(i) '.fa'];
        outputfile_DNA = [FolderRootName  filesep TranscriptName designerName 'probebatch' num2str(i) 'resultsDNA.xml'];
        database_DNA = db1;
        per_id = 100*MinHomologySize_Constant.Value/MaxProbeSize;
        strand = 'both';
        bnopts_DNA = blastplusoptions("blastn",strcat(" -num_alignments ", num2str(num_alignments)," -evalue ",num2str(evalue)," -word_size ",num2str(word_size)," -gapopen ",num2str(gapopen), ...
            " -gapextend ",num2str(gapextend)," -strand ",strand," -penalty ",num2str(penalty)," -reward ",num2str(reward)," -dust ",dust," -perc_identity ",num2str(per_id)))
        bnopts_DNA.Task = "blastn";
        bnopts_DNA.ReportFormat = "BLASTXML";
        blastplus("blastn",sequence_data,database_DNA ,outputfile_DNA,bnopts_DNA)
    end
    if (runRNA)
        sequence_data = [FolderRootName filesep TranscriptName designerName 'probebatch' num2str(i) '.fa'];
        outputfile_RNA = [FolderRootName  filesep TranscriptName designerName 'probebatch' num2str(i) 'resultsRNA.xml'];
        database_RNA = db2;
        per_id = 100*MinHomologySize_Constant.Value/MaxProbeSize;
        strand = 'plus';
        bnopts_RNA = blastplusoptions("blastn",strcat(" -num_alignments ", num2str(num_alignments)," -evalue ",num2str(evalue)," -word_size ",num2str(word_size)," -gapopen ",num2str(gapopen), ...
            " -gapextend ",num2str(gapextend)," -strand ",strand," -penalty ",num2str(penalty)," -reward ",num2str(reward)," -dust ",dust," -perc_identity ",num2str(per_id)))
        bnopts_RNA.Task = "blastn";
        bnopts_RNA.ReportFormat = "BLASTXML";
        blastplus("blastn",sequence_data,database_RNA,outputfile_RNA,bnopts_RNA)
    end
    updateParallel([],pwd);
end
progBar.release();
end

% process each batch blast report into structure needed for downstream usage in probe design and evaluation
% concantenate each batches results together into a single output file
if (~isempty(batch_nums_to_check2))
Batch_List = parallel.pool.Constant(Batch);
fprintf('\n')
fprintf('\n')
fprintf("Generating MATLAB tables from probe BLAST batch results")
fprintf('\n')
fprintf('\n')
progBar = ProgressBar(length(batch_nums_to_check2),'IsParallel',true,'WorkerDirectory', pwd(),'Title', 'Converting');
progBar.setup([],[],[]);
parfor y = 1:length(batch_nums_to_check2)
    pause(0.1);
    i = batch_nums_to_check2(y);
    total_temp_hits_subnode = [];    %Will store the gene hits information for storage
    %% For local blast on server
    pSeq = fastaread([FolderRootName filesep TranscriptName designerName 'probebatch' num2str(i) '.fa']);
    if (runDNA)
        if (simultaneous)%If All at Once
            temp_hits_subnode = readstruct([FolderRootName filesep TranscriptName designerName 'probebatch' num2str(i) 'resultsDNA.xml'],...
                "StructSelector",strcat('(//BlastOutput//BlastOutput_iterations)'));
            temp_hits_subnode = struct2cell(temp_hits_subnode.Iteration(:));
            probe_szs = cellfun(@(x) length(x.Hit),temp_hits_subnode(5,:));
            temp_hits_subnode = cellfun(@(x) x.Hit,temp_hits_subnode(5,:),'Un',0);
            temp_hits_all = [temp_hits_subnode{:}];
            probe_ids = [];probe_ord = [];
            for v = 1:length(probe_szs)
                probe_ord = [probe_ord v*ones(1,probe_szs(v))];
                probe_ids = [probe_ids Batch_List.Value{i}(v)*ones(1,probe_szs(v))];
            end
            temp_hits_subnode = squeeze(struct2cell([temp_hits_all.Hit_hsps]));
            for v = 1:length(temp_hits_subnode)
                [temp_hits_subnode{v}.align_num] = deal(temp_hits_all(v).Hit_num);
                [temp_hits_subnode{v}.Name] = deal(temp_hits_all(v).Hit_def);
                [temp_hits_subnode{v}.ProbeNum] = deal(probe_ids(v));
                [temp_hits_subnode{v}.ProbeSequence] = deal(pSeq(probe_ord(v)).Sequence);
                S2 = [temp_hits_subnode{v}.Hsp_align_len];
                S2b = [temp_hits_subnode{v}.Hsp_gaps];
                S2c = arrayfun(@(x)S2(x)-S2b(x),1:length(S2),'Un',0);
                [temp_hits_subnode{v}.Match] = S2c{:};
                if (sum(([temp_hits_subnode{v}.Match]>=MinHomologySize_Constant.Value))>0)
                    temp_hits_subnode{v} = temp_hits_subnode{v}([temp_hits_subnode{v}.Match]>=MinHomologySize_Constant.Value);
                end
            end
            temp_hits_subnode  = [temp_hits_subnode{:}];%cat cell array of structs
            S2 = [temp_hits_subnode.Hsp_align_len];
            S1 = [temp_hits_subnode.Hsp_query_frame;temp_hits_subnode.Hsp_hit_frame].';
            S1a = arrayfun(@(x)strplaceStrand(S1(x,:)),1:size(S1,1),'Un',0);
            [temp_hits_subnode.Strand] = S1a{:};
            S_SubjectIndices = [temp_hits_subnode.Hsp_hit_from;temp_hits_subnode.Hsp_hit_to].';
            S_SubjectIndices = arrayfun(@(x) S_SubjectIndices(x,:),1:length(S2),'Un',0).';
            [temp_hits_subnode.SubjectIndices] = S_SubjectIndices{:};
            S_QueryIndices = [temp_hits_subnode.Hsp_query_from;temp_hits_subnode.Hsp_query_to].';
            S_QueryIndices = arrayfun(@(x) S_QueryIndices(x,:),1:length(S2),'Un',0).';
            [temp_hits_subnode.QueryIndices] = S_QueryIndices{:};
            S3a = cellfun(@(x) x(2),S_QueryIndices,'Un',0);
            [temp_hits_subnode.Possible] = S3a{:};
            S4 = arrayfun(@(x) round(100*S2(x)/S3a{x},0),1:length(S2),'Un',0);
            [temp_hits_subnode.Percent] = S4{:};
            S_Score = num2cell([temp_hits_subnode.Hsp_score]);
            [temp_hits_subnode.Score] = S_Score{:};
            S_Expect = num2cell([temp_hits_subnode.Hsp_evalue]);
            [temp_hits_subnode.Expect] = S_Expect{:};
            S_Alignment = arrayfun(@(x) [char(temp_hits_subnode(x).Hsp_qseq);...
                char(temp_hits_subnode(x).Hsp_midline);char(temp_hits_subnode(x).Hsp_hseq)],...
                1:length(S2),'Un',0);
            [temp_hits_subnode.Alignment] = S_Alignment{:};
            temp_hits_subnode = rmfield(temp_hits_subnode,...
                {'Hsp_query_frame','Hsp_hit_frame','Hsp_hit_from','Hsp_hit_to',...
                'Hsp_query_from','Hsp_query_to','Hsp_qseq','Hsp_hseq','Hsp_midline',...
                'Hsp_num','Hsp_bit_score','Hsp_gaps','Hsp_align_len',...
                'Hsp_positive','Hsp_identity','Hsp_score','Hsp_evalue'});
            total_temp_hits_subnode = [total_temp_hits_subnode temp_hits_subnode];
        else%If Individual Probes Sequentially
            for w = 1:length(Batch_List.Value{i})
                temp_hits_subnode = readstruct([FolderRootName filesep TranscriptName designerName 'probebatch' num2str(i) 'resultsDNA.xml'],...
                    "StructSelector",strcat('(//BlastOutput//BlastOutput_iterations//Iteration[',num2str(w),']//Iteration_hits)'));
                sub_neststruc = squeeze(struct2cell([temp_hits_subnode.Hit(:).Hit_hsps]));
                for v = 1:length(sub_neststruc)
                    [sub_neststruc{v}.align_num] = deal(temp_hits_subnode.Hit(v).Hit_num);
                    [sub_neststruc{v}.Name] = deal(temp_hits_subnode.Hit(v).Hit_def);
                    [sub_neststruc{v}.ProbeNum] = deal(Batch_List.Value{i}(w));
                    [sub_neststruc{v}.ProbeSequence] = deal(pSeq(w).Sequence);
                end
                temp_hits_subnode  = [sub_neststruc{:}];%cat cell array of structs
                S2 = [temp_hits_subnode.Hsp_align_len];
                S2b = [temp_hits_subnode.Hsp_gaps];
                S2c = arrayfun(@(x)S2(x)-S2b(x),1:length(S2),'Un',0);
                [temp_hits_subnode.Match] = S2c{:};
                if (sum(([temp_hits_subnode.Match]>=MinHomologySize_Constant.Value))>0)%What if no hits over MinHomologySize_Constant.Value?
                    temp_hits_subnode = temp_hits_subnode([temp_hits_subnode.Match]>=MinHomologySize_Constant.Value);
                    S2 = [temp_hits_subnode.Hsp_align_len];
                    S1 = [temp_hits_subnode.Hsp_query_frame;temp_hits_subnode.Hsp_hit_frame].';
                    S1a = arrayfun(@(x)strplaceStrand(S1(x,:)),1:size(S1,1),'Un',0);
                    [temp_hits_subnode.Strand] = S1a{:};
                    S_SubjectIndices = [temp_hits_subnode.Hsp_hit_from;temp_hits_subnode.Hsp_hit_to].';
                    S_SubjectIndices = arrayfun(@(x) S_SubjectIndices(x,:),1:length(S2),'Un',0).';
                    [temp_hits_subnode.SubjectIndices] = S_SubjectIndices{:};
                    S_QueryIndices = [temp_hits_subnode.Hsp_query_from;temp_hits_subnode.Hsp_query_to].';
                    S_QueryIndices = arrayfun(@(x) S_QueryIndices(x,:),1:length(S2),'Un',0).';
                    [temp_hits_subnode.QueryIndices] = S_QueryIndices{:};
                    S3a = cellfun(@(x) x(2),S_QueryIndices,'Un',0);
                    [temp_hits_subnode.Possible] = S3a{:};
                    S4 = arrayfun(@(x) round(100*S2(x)/S3a{x},0),1:length(S2),'Un',0);
                    [temp_hits_subnode.Percent] = S4{:};
                    S_Score = num2cell([temp_hits_subnode.Hsp_score]);
                    [temp_hits_subnode.Score] = S_Score{:};
                    S_Expect = num2cell([temp_hits_subnode.Hsp_evalue]);
                    [temp_hits_subnode.Expect] = S_Expect{:};
                    S_Alignment = arrayfun(@(x) [char(temp_hits_subnode(x).Hsp_qseq);...
                        char(temp_hits_subnode(x).Hsp_midline);char(temp_hits_subnode(x).Hsp_hseq)],...
                        1:length(S2),'Un',0);
                    [temp_hits_subnode.Alignment] = S_Alignment{:};
                    temp_hits_subnode = rmfield(temp_hits_subnode,...
                        {'Hsp_query_frame','Hsp_hit_frame','Hsp_hit_from','Hsp_hit_to',...
                        'Hsp_query_from','Hsp_query_to','Hsp_qseq','Hsp_hseq','Hsp_midline',...
                        'Hsp_num','Hsp_bit_score','Hsp_gaps','Hsp_align_len',...
                        'Hsp_positive','Hsp_identity','Hsp_score','Hsp_evalue'});
                    total_temp_hits_subnode = [total_temp_hits_subnode temp_hits_subnode];
                end
            end
        end
    end
    if (runRNA)
        if (simultaneous)%If All at Once
            temp_hits_subnode = readstruct([FolderRootName filesep TranscriptName designerName 'probebatch' num2str(i) 'resultsRNA.xml'],...
                "StructSelector",strcat('(//BlastOutput//BlastOutput_iterations)'));
            temp_hits_subnode = struct2cell(temp_hits_subnode.Iteration(:));
            probe_szs = cellfun(@(x) length(x.Hit),temp_hits_subnode(5,:));
            temp_hits_subnode = cellfun(@(x) x.Hit,temp_hits_subnode(5,:),'Un',0);
            temp_hits_all = [temp_hits_subnode{:}];
            probe_ids = [];probe_ord = [];
            for v = 1:length(probe_szs)
                probe_ord = [probe_ord v*ones(1,probe_szs(v))];
                probe_ids = [probe_ids Batch_List.Value{i}(v)*ones(1,probe_szs(v))];
            end
            temp_hits_subnode = squeeze(struct2cell([temp_hits_all.Hit_hsps]));
            for v = 1:length(temp_hits_subnode)
                [temp_hits_subnode{v}.align_num] = deal(temp_hits_all(v).Hit_num);
                [temp_hits_subnode{v}.Name] = deal(temp_hits_all(v).Hit_def);
                [temp_hits_subnode{v}.ProbeNum] = deal(probe_ids(v));
                [temp_hits_subnode{v}.ProbeSequence] = deal(pSeq(probe_ord(v)).Sequence);
                S2 = [temp_hits_subnode{v}.Hsp_align_len];
                S2b = [temp_hits_subnode{v}.Hsp_gaps];
                S2c = arrayfun(@(x)S2(x)-S2b(x),1:length(S2),'Un',0);
                [temp_hits_subnode{v}.Match] = S2c{:};
                if (sum(([temp_hits_subnode{v}.Match]>=MinHomologySize_Constant.Value))>0)
                    temp_hits_subnode{v} = temp_hits_subnode{v}([temp_hits_subnode{v}.Match]>=MinHomologySize_Constant.Value);
                end
            end
            temp_hits_subnode  = [temp_hits_subnode{:}];%cat cell array of structs
            S2 = [temp_hits_subnode.Hsp_align_len];
            S1 = [temp_hits_subnode.Hsp_query_frame;temp_hits_subnode.Hsp_hit_frame].';
            S1a = arrayfun(@(x)strplaceStrand(S1(x,:)),1:size(S1,1),'Un',0);
            [temp_hits_subnode.Strand] = S1a{:};
            S_SubjectIndices = [temp_hits_subnode.Hsp_hit_from;temp_hits_subnode.Hsp_hit_to].';
            S_SubjectIndices = arrayfun(@(x) S_SubjectIndices(x,:),1:length(S2),'Un',0).';
            [temp_hits_subnode.SubjectIndices] = S_SubjectIndices{:};
            S_QueryIndices = [temp_hits_subnode.Hsp_query_from;temp_hits_subnode.Hsp_query_to].';
            S_QueryIndices = arrayfun(@(x) S_QueryIndices(x,:),1:length(S2),'Un',0).';
            [temp_hits_subnode.QueryIndices] = S_QueryIndices{:};
            S3a = cellfun(@(x) x(2),S_QueryIndices,'Un',0);
            [temp_hits_subnode.Possible] = S3a{:};
            S4 = arrayfun(@(x) round(100*S2(x)/S3a{x},0),1:length(S2),'Un',0);
            [temp_hits_subnode.Percent] = S4{:};
            S_Score = num2cell([temp_hits_subnode.Hsp_score]);
            [temp_hits_subnode.Score] = S_Score{:};
            S_Expect = num2cell([temp_hits_subnode.Hsp_evalue]);
            [temp_hits_subnode.Expect] = S_Expect{:};
            S_Alignment = arrayfun(@(x) [char(temp_hits_subnode(x).Hsp_qseq);...
                char(temp_hits_subnode(x).Hsp_midline);char(temp_hits_subnode(x).Hsp_hseq)],...
                1:length(S2),'Un',0);
            [temp_hits_subnode.Alignment] = S_Alignment{:};
            temp_hits_subnode = rmfield(temp_hits_subnode,...
                {'Hsp_query_frame','Hsp_hit_frame','Hsp_hit_from','Hsp_hit_to',...
                'Hsp_query_from','Hsp_query_to','Hsp_qseq','Hsp_hseq','Hsp_midline',...
                'Hsp_num','Hsp_bit_score','Hsp_gaps','Hsp_align_len',...
                'Hsp_positive','Hsp_identity','Hsp_score','Hsp_evalue'});
            total_temp_hits_subnode = [total_temp_hits_subnode temp_hits_subnode];
        else%If Individual Probes Sequentially
            for w = 1:length(Batch_List.Value{i})
                temp_hits_subnode = readstruct([FolderRootName filesep TranscriptName designerName 'probebatch' num2str(i) 'resultsRNA.xml'],...
                    "StructSelector",strcat('(//BlastOutput//BlastOutput_iterations//Iteration[',num2str(w),']//Iteration_hits)'));
                sub_neststruc = squeeze(struct2cell([temp_hits_subnode.Hit(:).Hit_hsps]));
                for v = 1:length(sub_neststruc)
                    [sub_neststruc{v}.align_num] = deal(temp_hits_subnode.Hit(v).Hit_num);
                    [sub_neststruc{v}.Name] = deal(temp_hits_subnode.Hit(v).Hit_def);
                    [sub_neststruc{v}.ProbeNum] = deal(Batch_List.Value{i}(w));
                    [sub_neststruc{v}.ProbeSequence] = deal(pSeq(w).Sequence);
                end
                temp_hits_subnode  = [sub_neststruc{:}];%cat cell array of structs
                S2 = [temp_hits_subnode.Hsp_align_len];
                S2b = [temp_hits_subnode.Hsp_gaps];
                S2c = arrayfun(@(x)S2(x)-S2b(x),1:length(S2),'Un',0);
                [temp_hits_subnode.Match] = S2c{:};
                if (sum(([temp_hits_subnode.Match]>=MinHomologySize_Constant.Value))>0)%What if no hits over MinHomologySize_Constant.Value?
                    temp_hits_subnode = temp_hits_subnode([temp_hits_subnode.Match]>=MinHomologySize_Constant.Value);
                    S2 = [temp_hits_subnode.Hsp_align_len];
                    S1 = [temp_hits_subnode.Hsp_query_frame;temp_hits_subnode.Hsp_hit_frame].';
                    S1a = arrayfun(@(x)strplaceStrand(S1(x,:)),1:size(S1,1),'Un',0);
                    [temp_hits_subnode.Strand] = S1a{:};
                    S_SubjectIndices = [temp_hits_subnode.Hsp_hit_from;temp_hits_subnode.Hsp_hit_to].';
                    S_SubjectIndices = arrayfun(@(x) S_SubjectIndices(x,:),1:length(S2),'Un',0).';
                    [temp_hits_subnode.SubjectIndices] = S_SubjectIndices{:};
                    S_QueryIndices = [temp_hits_subnode.Hsp_query_from;temp_hits_subnode.Hsp_query_to].';
                    S_QueryIndices = arrayfun(@(x) S_QueryIndices(x,:),1:length(S2),'Un',0).';
                    [temp_hits_subnode.QueryIndices] = S_QueryIndices{:};
                    S3a = cellfun(@(x) x(2),S_QueryIndices,'Un',0);
                    [temp_hits_subnode.Possible] = S3a{:};
                    S4 = arrayfun(@(x) round(100*S2(x)/S3a{x},0),1:length(S2),'Un',0);
                    [temp_hits_subnode.Percent] = S4{:};
                    S_Score = num2cell([temp_hits_subnode.Hsp_score]);
                    [temp_hits_subnode.Score] = S_Score{:};
                    S_Expect = num2cell([temp_hits_subnode.Hsp_evalue]);
                    [temp_hits_subnode.Expect] = S_Expect{:};
                    S_Alignment = arrayfun(@(x) [char(temp_hits_subnode(x).Hsp_qseq);...
                        char(temp_hits_subnode(x).Hsp_midline);char(temp_hits_subnode(x).Hsp_hseq)],...
                        1:length(S2),'Un',0);
                    [temp_hits_subnode.Alignment] = S_Alignment{:};
                    temp_hits_subnode = rmfield(temp_hits_subnode,...
                        {'Hsp_query_frame','Hsp_hit_frame','Hsp_hit_from','Hsp_hit_to',...
                        'Hsp_query_from','Hsp_query_to','Hsp_qseq','Hsp_hseq','Hsp_midline',...
                        'Hsp_num','Hsp_bit_score','Hsp_gaps','Hsp_align_len',...
                        'Hsp_positive','Hsp_identity','Hsp_score','Hsp_evalue'});
                    total_temp_hits_subnode = [total_temp_hits_subnode temp_hits_subnode];
                end
            end
        end
    end
    total_temp_hits_subnode = struct2table(total_temp_hits_subnode);
    total_temp_hits_subnode = total_temp_hits_subnode(:,{'Score' 'Expect' 'Strand' 'Alignment' 'QueryIndices' 'SubjectIndices' 'Name' 'ProbeSequence' 'ProbeNum' 'Match' 'Possible' 'Percent'});
    total_temp_hits_subnode = sortrows(total_temp_hits_subnode,10,'descend');%Check is Match
    %filter out hits to Bacterial Artifical Chromosomes, and first on-target match.
    loc_TranscriptName = find(contains(total_temp_hits_subnode.Name,TranscriptName));
    loc_TranscriptName2 = find(contains(total_temp_hits_subnode.Name,TranscriptName2));
    loc_BAC = find(contains(total_temp_hits_subnode.Name,'BAC clone'));
    loc_Chromosome = find(contains(total_temp_hits_subnode.Name,['chromosome ' Chromosome]));%only k<=2 kept, so remove first entry
    loc_ORF = find(contains(total_temp_hits_subnode.Name,['open reading frame' Chromosome]));%only k<=2 kept, so remove first entry
    loc_DNA = union(loc_Chromosome,loc_ORF);
    columns_All = 1:size(total_temp_hits_subnode,1);
    columns_Remaining1 = setdiff(columns_All,loc_BAC);
    if (length(loc_DNA)>1)
        columns_Remaining2 = setdiff(columns_All,union(loc_DNA(2:end),union(loc_BAC,union(loc_TranscriptName,loc_TranscriptName2))));
    else
        columns_Remaining2 = setdiff(columns_All,union(loc_BAC,union(loc_TranscriptName,loc_TranscriptName2)));
    end
    temp2 = total_temp_hits_subnode(columns_Remaining1,:);
    parsave_gene_hits_table([FolderRootName filesep TranscriptName designerName '_gene_hits_table_batch' num2str(i) '.mat'],temp2)
    if (runDNA)
        if exist([FolderRootName filesep TranscriptName designerName 'probebatch' num2str(i) 'resultsDNA.xml'],'file')        %delete xml if already exists
            delete([FolderRootName filesep TranscriptName designerName 'probebatch' num2str(i) 'resultsDNA.xml'])
        end
    end
    if (runRNA)
        if exist([FolderRootName filesep TranscriptName designerName 'probebatch' num2str(i) 'resultsRNA.xml'],'file')        %delete xml if already exists
            delete([FolderRootName filesep TranscriptName designerName 'probebatch' num2str(i) 'resultsRNA.xml'])
        end
    end
    if exist([FolderRootName filesep TranscriptName designerName 'probebatch' num2str(i) '.fa'],'file')        %delete fasta if already exists
        delete([FolderRootName filesep TranscriptName designerName 'probebatch' num2str(i) '.fa'])
    end
    updateParallel([],pwd);
end
progBar.release();
end
%% Load all files and clear
fprintf('\n')
fprintf('\n')
fprintf("Aggregating probe BLAST batch MATLAB tables into a single MATLAB BLAST gene hits table")
fprintf('\n')
fprintf('\n')
% delete temporary files generated for each batches blast report
tic
progBar = ProgressBar(N_Batches);
for i = 1:N_Batches
    if isfile([FolderRootName filesep TranscriptName designerName '_gene_hits_table_batch' num2str(i) '.mat'])
        load([FolderRootName filesep TranscriptName designerName '_gene_hits_table_batch' num2str(i) '.mat'],'gene_hits_table_tmp');
    end
    try
        if (first_loaded == 1)
            gene_table = vertcat(gene_table,gene_hits_table_tmp);
        else
            gene_table = gene_hits_table_tmp;
            first_loaded = 1;
        end
        clear gene_hits_table_tmp
    catch
    end
        progBar([],[],[]);
end
progBar.release();
save([settings.FolderRootName filesep TranscriptName '_' settings.rootName '_hits_table' settings.designerName '.mat'],'gene_table','-v7.3')
save([settings.FolderRootName filesep TranscriptName '_' settings.rootName '_fail_nums' settings.designerName '.mat'],'fail_nums','-v7.3')
tEnd = toc;fprintf('\n')
fprintf("Time elapsed to aggregate BLAST result batch tables %g seconds",round(tEnd,3,"significant"))
fprintf('\n')
fprintf('\n')
fprintf("Deleting temporary probe BLAST batch gene hits tables")
fprintf('\n')
fprintf('\n')
progBar = ProgressBar(N_Batches);
for i = 1:N_Batches
    if exist([settings.FolderRootName filesep TranscriptName settings.designerName '_gene_hits_table_batch' num2str(i) '.mat'],'file')        %delete temp mat file if already exists
        delete([settings.FolderRootName filesep TranscriptName settings.designerName '_gene_hits_table_batch' num2str(i) '.mat'])
    end
    progBar([],[],[]);
end
progBar.release();
end
function o = strplaceStrand(vec)
if (prod(vec)>0)
    o = 'Plus/Plus';
else
    o = 'Plus/Minus';
end
end

