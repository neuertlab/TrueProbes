function [RNAdbParser, DNAdbParser] = A3_BlastDBCMD_JH(settings,gene_table)
%% This takes input data and gets flanking sequence of blast hits in blastdb
doFlankCalc = settings.UseFlankingInfo;
runDNA = settings.BLASTdna;
runRNA = settings.BLASTrna;
most_recent_num_local = settings.num_parpool_local;
designerName = settings.designerName;
FolderRootName = settings.FolderRootName;
TranscriptName = settings.GeneName;
Organism = settings.Organism;
targetBatchSize = settings.TargetBatchSize;
SEQdbRoot = settings.SEQdbRoot;
if (settings.clusterStatus)
    most_recent_num = double(string(getenv('SLURM_JOB_CPUS_PER_NODE')));
else
    most_recent_num = most_recent_num_local;
end
gene_table = sortrows(gene_table,[7 6],'ascend');
gene_table = gene_table(gene_table.Match>=settings.MinHomologySearchTargetSize,:);
gene_table.Ax = min(gene_table.SubjectIndices,[],2);
gene_table.Bx = max(gene_table.SubjectIndices,[],2);
gene_table = sortrows(gene_table,[7 13],'ascend');
gene_table_ANs = extractBefore(gene_table.Name,' ');
NamesZ = gene_table.Name;
A1_MakeBlastDB_JH(settings)
if (isKey(settings.LocRoot_FASTA,Organism))
    if (runDNA)
        db1 = char(settings.LocBLASTDB_DNA(Organism));
    end
    if (runRNA)
        db2 = char(settings.LocBLASTDB_RNA(Organism));
    end
else
    if (runDNA)
        db1 = char(settings.specificBlastDatabase);
    end
    if (runRNA)
        db2 = char(settings.specificBlastDatabase2);
    end
end
if (ismac)
    system('ulimit -n 65536');
end
if (runDNA)
    BLASTdb_DNAexists = isfile(strcat(db1,'.nsq'));
else
    BLASTdb_DNAexists = 1;
end
if (runRNA)
    BLASTdb_RNAexists = isfile(strcat(db2,'.nsq'));
else
    BLASTdb_RNAexists = 1;
end
DNAdbParser = [];
RNAdbParser = [];
if (or(BLASTdb_DNAexists,BLASTdb_RNAexists))
    FNAFiles_FNA = dir([SEQdbRoot '**/*.fna']);
    FNAFiles_FA = dir([SEQdbRoot '**/*.fa']);
    FNAFiles_FASTA = dir([SEQdbRoot '**/*.fasta']);
    FNAFiles_FRN = dir([SEQdbRoot '**/*.frn']);
    if (~isempty(FNAFiles_FNA)&&~isempty(FNAFiles_FA))
        FNAFiles_1 = cat(1,FNAFiles_FNA,FNAFiles_FA);
    elseif (~isempty(FNAFiles_FNA))
        FNAFiles_1 = FNAFiles_FNA;
    elseif (~isempty(FNAFiles_FA))
        FNAFiles_1 = FNAFiles_FA;
    else
        FNAFiles_1 = [];
    end
    if (~isempty(FNAFiles_FASTA)&&~isempty(FNAFiles_FRN))
        FNAFiles_2 = cat(1,FNAFiles_FASTA,FNAFiles_FRN);
    elseif (~isempty(FNAFiles_FASTA))
        FNAFiles_2 = FNAFiles_FNA;
    elseif (~isempty(FNAFiles_FRN))
        FNAFiles_2 = FNAFiles_FRN;
    else
        FNAFiles_2 = [];
    end
    if (~isempty(FNAFiles_1)&&~isempty(FNAFiles_2))
        FNAFiles = cat(1,FNAFiles_1,FNAFiles_2);
    elseif (~isempty(FNAFiles_1))
        FNAFiles = FNAFiles_1;
    elseif (~isempty(FNAFiles_2))
        FNAFiles = FNAFiles_2;
    else
        FNAFiles = [];
    end
    if (~isempty(FNAFiles))
        FNAFileLoc = cell(1,length(FNAFiles));
        BioIFobj = cell(1,length(FNAFiles));
        BioKeys = cell(1,length(FNAFiles));
        for k=1:length(FNAFiles)
            FNAFileLoc{k} = strcat(FNAFiles(k).folder,filesep,FNAFiles(k).name);
            BioIFobj{k} = BioIndexedFile('FASTA',FNAFileLoc{k});
            BioIFobj{k}.Interpreter = @(x)fastaread(x,'TRIMHEADERS',false);
            BioKeys{k}{:} = getKeys(BioIFobj{k});
        end
    else
        msg = strcat('Error. There must be FASTA DNA type files in the Root FASTA Folder location: ',SEQdbRoot);
        error(msg)
    end
    strand_dict = dictionary(["plus" "minus"],1:2);
    if (runDNA)
        if (BLASTdb_DNAexists)
            if (~runRNA)
                DNAdbParser =  unique(convertCharsToStrings(table2array(gene_table(:,'Name'))));
            else
                outputDNAdb_sseqids = strcat(settings.FolderRootName,filesep, 'BLASTdb_DNA_sseqids.txt');
                dbtype ="nucl";
                if (strcmp(filesep,'/'))
                    inputDatabaseDNA = db1;
                else
                    inputDatabaseDNA = strcat(strrep(db1,'/',filesep));
                end
                if (~isfile(outputDNAdb_sseqids))
                    outfmt = '%i | %l | %a | %t';
                    fprintf("Getting DNA BLAST database sseqids and header names")
                    fprintf('\n')
                    fprintf('\n')
                    if (~ismac)
                        args_DNAparser = sprintf('%s -db %s -entry all -dbtype %s -out %s -outfmt %s', "blastdbcmd",inputDatabaseDNA, dbtype,outputDNAdb_sseqids,strcat('"',outfmt,'"'));
                    else
                        args_DNAparser = sprintf('%s -db %s -entry all -dbtype %s -out %s -outfmt %s', "blastdbcmd",inputDatabaseDNA, dbtype,strcat(char(39),outputDNAdb_sseqids,char(39)),strcat('"',outfmt,'"'));
                    end
                    [status,msg] = system(args_DNAparser);
                    if status
                        fprintf(msg)
                        fprintf('\n')
                        fprintf('\n')
                    end
                end
                fprintf('Checked DNA BLAST database sseqids and header names')
                fprintf('\n')
                fprintf('\n')
                DNAdb_sseqids_optsFile = detectImportOptions(outputDNAdb_sseqids,'FileType','delimitedtext','Delimiter',' | ');
                DNAdb_sseqids_optsFile.VariableNames = {'sseqid','seqlength','AccessionNumber','Name'};
                DNAdb_sseqids_Table = readtable(outputDNAdb_sseqids,DNAdb_sseqids_optsFile);
                if (sum(contains(DNAdb_sseqids_Table.sseqid,'BL_ORD_ID'))>0)
                    sseqids_preserved = 0;
                else
                    sseqids_preserved = 1;
                end
                if sseqids_preserved
                    DNAdbParser = join(convertCharsToStrings(table2array(DNAdb_sseqids_Table(:,3:4)))," ");
                else
                    DNAdbParser = convertCharsToStrings(table2array(DNAdb_sseqids_Table(:,'Name')));
                end
            end
            if (doFlankCalc)
                if (~isfile([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_flanking_DNA_sequences' settings.designerName '.mat']))
                    if (~isfile([strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_probetarget_flanking_sequences', settings.designerName,'.mat')]))
                        DNA_entries = find(ismember(NamesZ,DNAdbParser));
                        DNA_Accession_Numbers = extractBefore(gene_table.Name(DNA_entries),' ');
                        if (~isempty(DNA_entries))
                            unique_DNA_ANs = unique(DNA_Accession_Numbers,'stable');
                            Num_Unique_DNA_ANs = length(unique(DNA_Accession_Numbers));
                            Num_Unique_DNA_ANs_Batches = ceil(Num_Unique_DNA_ANs/targetBatchSize);
                            R = mod(Num_Unique_DNA_ANs,targetBatchSize);
                            Unique_DNA_ANs_Batches = cell(1,Num_Unique_DNA_ANs_Batches);
                            if (R==0)
                                for k = 1:Num_Unique_DNA_ANs_Batches
                                    Unique_DNA_ANs_Batches{k} = targetBatchSize*(k-1)+1:targetBatchSize*k;
                                end
                            else
                                for k = 1:Num_Unique_DNA_ANs_Batches-1
                                    Unique_DNA_ANs_Batches{k} = targetBatchSize*(k-1)+1:targetBatchSize*k;
                                end
                                Unique_DNA_ANs_Batches{Num_Unique_DNA_ANs_Batches} = targetBatchSize*(Num_Unique_DNA_ANs_Batches-1)+1:targetBatchSize*(Num_Unique_DNA_ANs_Batches-1)+R;
                            end
                            ResultsExist_DNA_ANs_Batches = zeros(1,Num_Unique_DNA_ANs_Batches);
                            ResultsDate_DNA_ANs_Batches = cell(1,Num_Unique_DNA_ANs_Batches);
                            fprintf("Check if DNA Hit BLAST alignment flanking sequence batch files exist")
                            fprintf('\n')
                            fprintf('\n')
                            wb = parwaitbar(Num_Unique_DNA_ANs_Batches,'WaitMessage','Checking');
                            parfor i = 1:Num_Unique_DNA_ANs_Batches
                                if (isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_flanking_DNA_sequences_targetbatch' num2str(i) '.mat']))%check if temp file exists
                                    d = dir([FolderRootName filesep '(' TranscriptName ')' designerName '_flanking_DNA_sequences_targetbatch' num2str(i) '.mat']);
                                    if (d.bytes>0)%check size greater than zero
                                        ResultsExist_DNA_ANs_Batches(i) = 1;
                                    end
                                    ResultsDate_DNA_ANs_Batches{i} = datetime(d.date);
                                end
                                progress(wb);
                            end
                            wb.delete();
                            fprintf('\n')
                            fprintf('\n')
                            Results_DNA_ANs_Batches_NotMade = find(ResultsExist_DNA_ANs_Batches==0);
                            ResultsExist_DNA_ANs_Batches_Made = find(ResultsExist_DNA_ANs_Batches==1);
                            %Sort get most recent ResultsMade GeneHitsMade and GeneHitsTable Made and add to probe_check_list
                            if (length(ResultsExist_DNA_ANs_Batches_Made)<=most_recent_num)
                                DNA_results_check = ResultsExist_DNA_ANs_Batches_Made;
                            else
                                Results_RecentMade_DNA(:,1) = ResultsDate_DNA_ANs_Batches(ResultsExist_DNA_ANs_Batches_Made);
                                Results_RecentMade_DNA(:,2) = num2cell(ResultsExist_DNA_ANs_Batches_Made);
                                Results_RecentMade_DNA = table2timetable(cell2table(Results_RecentMade_DNA));
                                Results_RecentMade_DNA = sortrows(Results_RecentMade_DNA,1,'descend');
                                Results_RecentMade_DNA.Properties.VariableNames = {'ID'};
                                DNA_results_check = Results_RecentMade_DNA.ID(1:most_recent_num).';
                                clear Results_RecentMade_Dates
                            end
                            batch_nums_to_check_DNA = union(Results_DNA_ANs_Batches_NotMade,DNA_results_check);
                            BioIFobj_List = parallel.pool.Constant(BioIFobj);
                            gene_table_List = parallel.pool.Constant(gene_table);
                            strand_dict_List = parallel.pool.Constant(strand_dict);
                            Unique_DNA_ANs_Batches_List = parallel.pool.Constant(Unique_DNA_ANs_Batches);
                            unique_DNA_ANs_List = parallel.pool.Constant(unique_DNA_ANs);
                            BioKeys_List = parallel.pool.Constant(BioKeys);
                            fprintf("Getting batch DNA Hit BLAST alignment flanking sequences")
                            fprintf('\n')
                            fprintf('\n')
                            wb = parwaitbar(length(batch_nums_to_check_DNA),'WaitMessage','Computing');
                            for nn = 1:length(batch_nums_to_check_DNA)
                                unique_DNA_ANs_in_batch = unique_DNA_ANs_List.Value(Unique_DNA_ANs_Batches_List.Value{batch_nums_to_check_DNA(nn)});
                                AN_DNA_hit_locs_in_batch = find(ismember(gene_table_ANs,unique_DNA_ANs_in_batch));
                                probe_DNA_Target_flanking_info_tmp = table('Size',[length(AN_DNA_hit_locs_in_batch) 4],'VariableNames',{'Name','gene_table_location','TargetSequence_5primeTo3prime','ProbeRevCompSequence_5primeTo3prime'},'VariableTypes',{'string','double','string','string'});
                                for kk = 1:length(unique_DNA_ANs_in_batch)
                                    AN_DNA_hit_locs = find(strcmp(gene_table_ANs,unique_DNA_ANs_in_batch(kk)));
                                    DNA_info = read(BioIFobj_List.Value{find(arrayfun(@(k) sum(contains(BioKeys_List.Value{k}{:},unique_DNA_ANs_in_batch(kk))),1:length(BioKeys)),1)},unique_DNA_ANs_in_batch(kk));
                                    sub_DNA_Target_plus_Sequence = DNA_info.Sequence;
                                    sub_DNA_Target_minus_Sequence = seqrcomplement(DNA_info.Sequence);
                                    sub_DNA_Target_Sequence_Case = {sub_DNA_Target_plus_Sequence sub_DNA_Target_minus_Sequence};
                                    sub_DNA_length_low = gene_table_List.Value(AN_DNA_hit_locs,'SubjectIndices').SubjectIndices(:,1) - gene_table_List.Value(AN_DNA_hit_locs,'QueryIndices').QueryIndices(:,1) + 1;
                                    sub_DNA_length_high = gene_table_List.Value(AN_DNA_hit_locs,'SubjectIndices').SubjectIndices(:,2) - gene_table_List.Value(AN_DNA_hit_locs,'QueryIndices').QueryIndices(:,2) + cellfun(@length,gene_table_List.Value(AN_DNA_hit_locs,'ProbeSequence').ProbeSequence);
                                    sub_DNA_length_max = length(DNA_info.Sequence)*ones(size(sub_DNA_length_low));
                                    sub_DNA_length_lower = sub_DNA_length_low;
                                    sub_DNA_length_lower(sub_DNA_length_lower<1) = 1;
                                    sub_DNA_length_upper = min([sub_DNA_length_high sub_DNA_length_max],[],2);
                                    sub_DNA_Hit_strands = convertCharsToStrings(lower(extractAfter(gene_table_List.Value.Strand(AN_DNA_hit_locs),'/')));
                                    sub_AN_DNA_hit_TargetSequences = arrayfun(@(n_sub) sub_DNA_Target_Sequence_Case{strand_dict_List.Value(sub_DNA_Hit_strands(n_sub))}(sub_DNA_length_lower(n_sub):sub_DNA_length_upper(n_sub)),1:length(AN_DNA_hit_locs),'Un',0);
                                    sub_DNA_target_Block_1 = convertCharsToStrings(arrayfun(@(ni) sub_AN_DNA_hit_TargetSequences{ni}(1:gene_table_List.Value(AN_DNA_hit_locs(ni),:).SubjectIndices(1)-sub_DNA_length_lower(ni)),1:length(AN_DNA_hit_locs),'Un',0));
                                    sub_DNA_target_Block_2 = convertCharsToStrings(arrayfun(@(ni) strrep(gene_table_List.Value.Alignment{ni}(3,:),'-','N'),AN_DNA_hit_locs,'Un',0));
                                    sub_DNA_target_Block_3 = convertCharsToStrings(arrayfun(@(ni) sub_AN_DNA_hit_TargetSequences{ni}(end-sub_DNA_length_upper(ni)+gene_table_List.Value(AN_DNA_hit_locs(ni),:).SubjectIndices(2)+1:end),1:length(AN_DNA_hit_locs),'Un',0));
                                    sub_DNA_target_matching_flank_Trimmed_blocks = [reshape(sub_DNA_target_Block_1,[],1) reshape(sub_DNA_target_Block_2,[],1) reshape(sub_DNA_target_Block_3,[],1)]';
                                    sub_DNA_target_matching_flank_Trimmed = join(sub_DNA_target_matching_flank_Trimmed_blocks,"",1)';
                                    sub_DNA_probe_Block_1 = convertCharsToStrings(arrayfun(@(ni) gene_table_List.Value(AN_DNA_hit_locs(ni),:).ProbeSequence{:}(gene_table_List.Value(AN_DNA_hit_locs(ni),:).QueryIndices(1)-gene_table_List.Value(AN_DNA_hit_locs(ni),:).SubjectIndices(1)+sub_DNA_length_lower(ni):gene_table_List.Value(AN_DNA_hit_locs(ni),:).QueryIndices(1)-1),1:length(AN_DNA_hit_locs),'Un',0));
                                    sub_DNA_probe_Block_2 =  convertCharsToStrings(arrayfun(@(ni) strrep(gene_table_List.Value.Alignment{ni}(1,:),'-','N'),AN_DNA_hit_locs,'Un',0));
                                    sub_DNA_probe_Block_3 = convertCharsToStrings(arrayfun(@(ni) gene_table_List.Value(AN_DNA_hit_locs(ni),:).ProbeSequence{:}(gene_table_List.Value(AN_DNA_hit_locs(ni),:).QueryIndices(2)+1:gene_table_List.Value(AN_DNA_hit_locs(ni),:).QueryIndices(2)+sub_DNA_length_upper(ni)-gene_table_List.Value(AN_DNA_hit_locs(ni),:).SubjectIndices(2)),1:length(AN_DNA_hit_locs),'Un',0));
                                    sub_DNA_probe_matching_flank_Trimmed_blocks = [reshape(sub_DNA_probe_Block_1,[],1) reshape(sub_DNA_probe_Block_2,[],1) reshape(sub_DNA_probe_Block_3,[],1)]';
                                    sub_DNA_probe_matching_flank_Trimmed = join(sub_DNA_probe_matching_flank_Trimmed_blocks,"",1)';
                                    probe_DNA_Target_flanking_info_tmp(ismember(AN_DNA_hit_locs_in_batch,AN_DNA_hit_locs),'Name') = gene_table_List.Value(AN_DNA_hit_locs,'Name');
                                    probe_DNA_Target_flanking_info_tmp(ismember(AN_DNA_hit_locs_in_batch,AN_DNA_hit_locs),'gene_table_location') = array2table(AN_DNA_hit_locs);
                                    probe_DNA_Target_flanking_info_tmp(ismember(AN_DNA_hit_locs_in_batch,AN_DNA_hit_locs),'TargetSequence_5primeTo3prime') = array2table(sub_DNA_target_matching_flank_Trimmed);
                                    probe_DNA_Target_flanking_info_tmp(ismember(AN_DNA_hit_locs_in_batch,AN_DNA_hit_locs),'ProbeRevCompSequence_5primeTo3prime') = array2table(sub_DNA_probe_matching_flank_Trimmed);
                                end
                                parsave_probe_target_flanking_info_tmp([FolderRootName filesep '(' TranscriptName ')' designerName '_flanking_DNA_sequences_targetbatch' num2str(batch_nums_to_check_DNA(nn)) '.mat'],probe_DNA_Target_flanking_info_tmp);
                                progress(wb);
                            end
                            fprintf('Acquired blast DNA hit target flanking sequences')
                            fprintf('\n')
                            fprintf('\n')
                        end
                        fprintf("Aggregating hit blast DNA target flanking sequence batches")
                        fprintf('\n')
                        fprintf('\n')
                        probe_DNA_target_flanking_info = cell(1,Num_Unique_DNA_ANs_Batches);
                        wb = parwaitbar(Num_Unique_DNA_ANs_Batches,'WaitMessage','Loading Batches');
                        parfor i = 1:Num_Unique_DNA_ANs_Batches
                            if isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_flanking_DNA_sequences_targetbatch' num2str(i) '.mat'])
                                probe_DNA_target_flanking_info{i} = load([FolderRootName filesep '(' TranscriptName ')' designerName '_flanking_DNA_sequences_targetbatch' num2str(i) '.mat']).probe_target_flanking_info_tmp;
                            end
                            progress(wb);
                        end
                        wb.delete();
                        probe_DNA_target_flanking_info = vertcat(probe_DNA_target_flanking_info{:});
                        save([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_flanking_DNA_sequences' settings.designerName '.mat'],'probe_DNA_target_flanking_info','-v7.3')
                        fprintf("Deleting temporary batch DNA Hit BLAST alignment flanking sequence files")
                        fprintf('\n')
                        fprintf('\n')
                        wb = parwaitbar(Num_Unique_DNA_ANs_Batches,'WaitMessage','Deleting Batches');
                        parfor i = 1:Num_Unique_DNA_ANs_Batches
                            if exist([FolderRootName filesep '(' TranscriptName ')' designerName '_flanking_DNA_sequences_targetbatch' num2str(i) '.mat'],'file')        %delete temp mat file if already exists
                                delete([FolderRootName filesep '(' TranscriptName ')' designerName '_flanking_DNA_sequences_targetbatch' num2str(i) '.mat'])
                            end
                            progress(wb);
                        end
                        wb.delete();
                        fprintf('\n')
                        fprintf('\n')
                    end
                end
            end
        end
    end
    if (runRNA)
        if (BLASTdb_RNAexists)
            if (~runDNA)
                RNAdbParser =  unique(convertCharsToStrings(table2array(gene_table(:,'Name'))));
            else
                outputRNAdb_sseqids = strcat(settings.FolderRootName,filesep, 'BLASTdb_RNA_sseqids.txt');
                dbtype ="nucl";
                if (strcmp(filesep,'/'))
                    inputDatabaseRNA = db2;
                else
                    inputDatabaseRNA = strcat(strrep(db2,'/',filesep));
                end
                if (~isfile(outputRNAdb_sseqids))
                    outfmt = '%i | %l | %a | %t';
                    fprintf("Getting RNA BLAST database sseqids and header names")
                    fprintf('\n')
                    fprintf('\n')
                    if (~ismac)
                        args_RNAparser = sprintf('%s -db %s -entry all -dbtype %s -out %s -outfmt %s', "blastdbcmd",inputDatabaseRNA, dbtype,outputRNAdb_sseqids,strcat('"',outfmt,'"'));
                    else
                        args_RNAparser = sprintf('%s -db %s -entry all -dbtype %s -out %s -outfmt %s', "blastdbcmd",inputDatabaseRNA, dbtype,strcat(char(39),outputRNAdb_sseqids,char(39)),strcat('"',outfmt,'"'));
                    end
                    [status,msg] = system(args_RNAparser);
                    if status
                        fprintf(msg)
                        fprintf('\n')
                        fprintf('\n')
                    end
                end
                fprintf('Checked RNA BLAST database sseqids and header names')
                fprintf('\n')
                fprintf('\n')
                %message = strcat('blastdbcmd -db', " ",inputDatabaseRNA,' -entry'," ",'''','gnl|BL_ORD_ID|23','''')
                RNAdb_sseqids_optsFile = detectImportOptions(outputRNAdb_sseqids,'FileType','delimitedtext','Delimiter',' | ');
                RNAdb_sseqids_optsFile.VariableNames = {'sseqid','seqlength','AccessionNumber','Name'};
                RNAdb_sseqids_Table = readtable(outputRNAdb_sseqids,RNAdb_sseqids_optsFile);
                if (sum(contains(RNAdb_sseqids_Table.sseqid,'BL_ORD_ID'))>0)
                    sseqids_preserved = 0;
                else
                    sseqids_preserved = 1;
                end
                if sseqids_preserved
                    RNAdbParser = join(convertCharsToStrings(table2array(RNAdb_sseqids_Table(:,3:4)))," ");
                else
                    RNAdbParser = convertCharsToStrings(table2array(RNAdb_sseqids_Table(:,'Name')));
                end
            end
            if (doFlankCalc)
                if (~isfile([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_flanking_RNA_sequences' settings.designerName '.mat']))
                    if (~isfile([strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_probetarget_flanking_sequences', settings.designerName,'.mat')]))
                        RNA_entries = find(ismember(NamesZ,RNAdbParser));
                        RNA_Accession_Numbers = extractBefore(gene_table.Name(RNA_entries),' ');
                        if (~isempty(RNA_entries))
                            unique_RNA_ANs = unique(RNA_Accession_Numbers,'stable');
                            Num_Unique_RNA_ANs = length(unique(RNA_Accession_Numbers));
                            Num_Unique_RNA_ANs_Batches = ceil(Num_Unique_RNA_ANs/targetBatchSize);
                            R = mod(Num_Unique_RNA_ANs,targetBatchSize);
                            Unique_RNA_ANs_Batches = cell(1,Num_Unique_RNA_ANs_Batches);
                            if (R==0)
                                for k = 1:Num_Unique_RNA_ANs_Batches
                                    Unique_RNA_ANs_Batches{k} = targetBatchSize*(k-1)+1:targetBatchSize*k;
                                end
                            else
                                for k = 1:Num_Unique_RNA_ANs_Batches-1
                                    Unique_RNA_ANs_Batches{k} = targetBatchSize*(k-1)+1:targetBatchSize*k;
                                end
                                Unique_RNA_ANs_Batches{Num_Unique_RNA_ANs_Batches} = targetBatchSize*(Num_Unique_RNA_ANs_Batches-1)+1:targetBatchSize*(Num_Unique_RNA_ANs_Batches-1)+R;
                            end
                            ResultsExist_RNA_ANs_Batches = zeros(1,Num_Unique_RNA_ANs_Batches);
                            ResultsDate_RNA_ANs_Batches = cell(1,Num_Unique_RNA_ANs_Batches);
                            fprintf("Check if RNA Hit BLAST alignment flanking sequence batch files exist")
                            fprintf('\n')
                            fprintf('\n')
                            wb = parwaitbar(Num_Unique_RNA_ANs_Batches,'WaitMessage','Checking');
                            parfor i = 1:Num_Unique_RNA_ANs_Batches
                                if (isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_flanking_RNA_sequences_targetbatch' num2str(i) '.mat']))%check if temp file exists
                                    d = dir([FolderRootName filesep '(' TranscriptName ')' designerName '_flanking_RNA_sequences_targetbatch' num2str(i) '.mat']);
                                    if (d.bytes>0)%check size greater than zero
                                        ResultsExist_RNA_ANs_Batches(i) = 1;
                                    end
                                    ResultsDate_RNA_ANs_Batches{i} = datetime(d.date);
                                end
                                progress(wb);
                            end
                            wb.delete();
                            fprintf('\n')
                            fprintf('\n')
                            Results_RNA_ANs_Batches_NotMade = find(ResultsExist_RNA_ANs_Batches==0);
                            ResultsExist_RNA_ANs_Batches_Made = find(ResultsExist_RNA_ANs_Batches==1);
                            %Sort get most recent ResultsMade GeneHitsMade and GeneHitsTable Made and add to probe_check_list
                            if (length(ResultsExist_RNA_ANs_Batches_Made)<=most_recent_num)
                                RNA_results_check = ResultsExist_RNA_ANs_Batches_Made;
                            else
                                Results_RecentMade_RNA(:,1) = ResultsDate_RNA_ANs_Batches(ResultsExist_RNA_ANs_Batches_Made);
                                Results_RecentMade_RNA(:,2) = num2cell(ResultsExist_RNA_ANs_Batches_Made);
                                Results_RecentMade_RNA = table2timetable(cell2table(Results_RecentMade_RNA));
                                Results_RecentMade_RNA = sortrows(Results_RecentMade_RNA,1,'descend');
                                Results_RecentMade_RNA.Properties.VariableNames = {'ID'};
                                RNA_results_check = Results_RecentMade_RNA.ID(1:most_recent_num).';
                                clear Results_RecentMade_Dates
                            end
                            batch_nums_to_check_RNA = union(Results_RNA_ANs_Batches_NotMade,RNA_results_check);
                            BioIFobj_List = parallel.pool.Constant(BioIFobj);
                            gene_table_List = parallel.pool.Constant(gene_table);
                            strand_dict_List = parallel.pool.Constant(strand_dict);
                            Unique_RNA_ANs_Batches_List = parallel.pool.Constant(Unique_RNA_ANs_Batches);
                            unique_RNA_ANs_List = parallel.pool.Constant(unique_RNA_ANs);
                            BioKeys_List = parallel.pool.Constant(BioKeys);
                            fprintf("Getting batch RNA Hit BLAST alignment flanking sequences")
                            fprintf('\n')
                            fprintf('\n')
                            wb = parwaitbar(length(batch_nums_to_check_RNA),'WaitMessage','Computing');
                            for nn = 1:length(batch_nums_to_check_RNA)
                                unique_RNA_ANs_in_batch = unique_RNA_ANs_List.Value(Unique_RNA_ANs_Batches_List.Value{batch_nums_to_check_RNA(nn)});
                                AN_RNA_hit_locs_in_batch = find(ismember(gene_table_ANs,unique_RNA_ANs_in_batch));
                                probe_RNA_Target_flanking_info_tmp = table('Size',[length(AN_RNA_hit_locs_in_batch) 4],'VariableNames',{'Name','gene_table_location','TargetSequence_5primeTo3prime','ProbeRevCompSequence_5primeTo3prime'},'VariableTypes',{'string','double','string','string'});
                                for kk = 1:length(unique_RNA_ANs_in_batch)
                                    AN_RNA_hit_locs = find(strcmp(gene_table_ANs,unique_RNA_ANs_in_batch(kk)));
                                    RNA_info = read(BioIFobj_List.Value{find(arrayfun(@(k) sum(contains(BioKeys_List.Value{k}{:},unique_RNA_ANs_in_batch(kk))),1:length(BioKeys)),1)},unique_RNA_ANs_in_batch(kk));
                                    sub_RNA_Target_plus_Sequence = RNA_info.Sequence;
                                    sub_RNA_Target_minus_Sequence = seqrcomplement(RNA_info.Sequence);
                                    sub_RNA_Target_Sequence_Case = {sub_RNA_Target_plus_Sequence sub_RNA_Target_minus_Sequence};
                                    sub_RNA_length_low = gene_table_List.Value(AN_RNA_hit_locs,'SubjectIndices').SubjectIndices(:,1) - gene_table_List.Value(AN_RNA_hit_locs,'QueryIndices').QueryIndices(:,1) + 1;
                                    sub_RNA_length_high = gene_table_List.Value(AN_RNA_hit_locs,'SubjectIndices').SubjectIndices(:,2) - gene_table_List.Value(AN_RNA_hit_locs,'QueryIndices').QueryIndices(:,2) + cellfun(@length,gene_table_List.Value(AN_RNA_hit_locs,'ProbeSequence').ProbeSequence);
                                    sub_RNA_length_max = length(RNA_info.Sequence)*ones(size(sub_RNA_length_low));
                                    sub_RNA_length_lower = sub_RNA_length_low;
                                    sub_RNA_length_lower(sub_RNA_length_lower<1) = 1;
                                    sub_RNA_length_upper = min([sub_RNA_length_high sub_RNA_length_max],[],2);
                                    sub_RNA_Hit_strands = convertCharsToStrings(lower(extractAfter(gene_table_List.Value.Strand(AN_RNA_hit_locs),'/')));
                                    sub_AN_RNA_hit_TargetSequences = arrayfun(@(n_sub) sub_RNA_Target_Sequence_Case{strand_dict_List.Value(sub_RNA_Hit_strands(n_sub))}(sub_RNA_length_lower(n_sub):sub_RNA_length_upper(n_sub)),1:length(AN_RNA_hit_locs),'Un',0);
                                    sub_RNA_target_Block_1 = convertCharsToStrings(arrayfun(@(ni) sub_AN_RNA_hit_TargetSequences{ni}(1:gene_table_List.Value(AN_RNA_hit_locs(ni),:).SubjectIndices(1)-sub_RNA_length_lower(ni)),1:length(AN_RNA_hit_locs),'Un',0));
                                    sub_RNA_target_Block_2 = convertCharsToStrings(arrayfun(@(ni) strrep(gene_table_List.Value.Alignment{ni}(3,:),'-','N'),AN_RNA_hit_locs,'Un',0));
                                    sub_RNA_target_Block_3 = convertCharsToStrings(arrayfun(@(ni) sub_AN_RNA_hit_TargetSequences{ni}(end-sub_RNA_length_upper(ni)+gene_table_List.Value(AN_RNA_hit_locs(ni),:).SubjectIndices(2)+1:end),1:length(AN_RNA_hit_locs),'Un',0));
                                    sub_RNA_target_matching_flank_Trimmed_blocks = [reshape(sub_RNA_target_Block_1,[],1) reshape(sub_RNA_target_Block_2,[],1) reshape(sub_RNA_target_Block_3,[],1)]';
                                    sub_RNA_target_matching_flank_Trimmed = join(sub_RNA_target_matching_flank_Trimmed_blocks,"",1)';
                                    sub_RNA_probe_Block_1 = convertCharsToStrings(arrayfun(@(ni) gene_table_List.Value(AN_RNA_hit_locs(ni),:).ProbeSequence{:}(gene_table_List.Value(AN_RNA_hit_locs(ni),:).QueryIndices(1)-gene_table_List.Value(AN_RNA_hit_locs(ni),:).SubjectIndices(1)+sub_RNA_length_lower(ni):gene_table_List.Value(AN_RNA_hit_locs(ni),:).QueryIndices(1)-1),1:length(AN_RNA_hit_locs),'Un',0));
                                    sub_RNA_probe_Block_2 =  convertCharsToStrings(arrayfun(@(ni) strrep(gene_table_List.Value.Alignment{ni}(1,:),'-','N'),AN_RNA_hit_locs,'Un',0));
                                    sub_RNA_probe_Block_3 = convertCharsToStrings(arrayfun(@(ni) gene_table_List.Value(AN_RNA_hit_locs(ni),:).ProbeSequence{:}(gene_table_List.Value(AN_RNA_hit_locs(ni),:).QueryIndices(2)+1:gene_table_List.Value(AN_RNA_hit_locs(ni),:).QueryIndices(2)+sub_RNA_length_upper(ni)-gene_table_List.Value(AN_RNA_hit_locs(ni),:).SubjectIndices(2)),1:length(AN_RNA_hit_locs),'Un',0));
                                    sub_RNA_probe_matching_flank_Trimmed_blocks = [reshape(sub_RNA_probe_Block_1,[],1) reshape(sub_RNA_probe_Block_2,[],1) reshape(sub_RNA_probe_Block_3,[],1)]';
                                    sub_RNA_probe_matching_flank_Trimmed = join(sub_RNA_probe_matching_flank_Trimmed_blocks,"",1)';
                                    probe_RNA_Target_flanking_info_tmp(ismember(AN_RNA_hit_locs_in_batch,AN_RNA_hit_locs),'Name') = gene_table_List.Value(AN_RNA_hit_locs,'Name');
                                    probe_RNA_Target_flanking_info_tmp(ismember(AN_RNA_hit_locs_in_batch,AN_RNA_hit_locs),'gene_table_location') = array2table(AN_RNA_hit_locs);
                                    probe_RNA_Target_flanking_info_tmp(ismember(AN_RNA_hit_locs_in_batch,AN_RNA_hit_locs),'TargetSequence_5primeTo3prime') = array2table(sub_RNA_target_matching_flank_Trimmed);
                                    probe_RNA_Target_flanking_info_tmp(ismember(AN_RNA_hit_locs_in_batch,AN_RNA_hit_locs),'ProbeRevCompSequence_5primeTo3prime') = array2table(sub_RNA_probe_matching_flank_Trimmed);
                                end
                                parsave_probe_target_flanking_info_tmp([FolderRootName filesep '(' TranscriptName ')' designerName '_flanking_RNA_sequences_targetbatch' num2str(batch_nums_to_check_RNA(nn)) '.mat'],probe_RNA_Target_flanking_info_tmp);
                                progress(wb);
                            end
                            fprintf('Acquired blast RNA hit target flanking sequences')
                            fprintf('\n')
                            fprintf('\n')
                        end
                        fprintf("Aggregating hit blast RNA target flanking sequence batches")
                        fprintf('\n')
                        fprintf('\n')
                        probe_RNA_target_flanking_info = cell(1,Num_Unique_RNA_ANs_Batches);
                        wb = parwaitbar(Num_Unique_RNA_ANs_Batches,'WaitMessage','Loading Batches');
                        parfor i = 1:Num_Unique_RNA_ANs_Batches
                            if isfile([FolderRootName filesep '(' TranscriptName ')' designerName '_flanking_RNA_sequences_targetbatch' num2str(i) '.mat'])
                                probe_RNA_target_flanking_info{i} = load([FolderRootName filesep '(' TranscriptName ')' designerName '_flanking_RNA_sequences_targetbatch' num2str(i) '.mat']).probe_target_flanking_info_tmp;
                            end
                            progress(wb);
                        end
                        wb.delete();
                        probe_RNA_target_flanking_info = vertcat(probe_RNA_target_flanking_info{:});
                        save([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_flanking_RNA_sequences' settings.designerName '.mat'],'probe_RNA_target_flanking_info','-v7.3')
                        fprintf("Deleting temporary batch RNA Hit BLAST alignment flanking sequence files")
                        fprintf('\n')
                        fprintf('\n')
                        wb = parwaitbar(Num_Unique_RNA_ANs_Batches,'WaitMessage','Deleting Batches');
                        parfor i = 1:Num_Unique_RNA_ANs_Batches
                            if exist([FolderRootName filesep '(' TranscriptName ')' designerName '_flanking_RNA_sequences_targetbatch' num2str(i) '.mat'],'file')        %delete temp mat file if already exists
                                delete([FolderRootName filesep '(' TranscriptName ')' designerName '_flanking_RNA_sequences_targetbatch' num2str(i) '.mat'])
                            end
                            progress(wb);
                        end
                        wb.delete();
                        fprintf('\n')
                        fprintf('\n')
                    end
                end
            end
        end
    end
end
if (doFlankCalc)
    if (~isfile([strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_probetarget_flanking_sequences', settings.designerName,'.mat')]))
        if (and(runDNA,runRNA))
            probetarget_flanking_info = cell(1,2);
            probetarget_flanking_info{1} = load([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_flanking_RNA_sequences' settings.designerName '.mat']).probe_RNA_target_flanking_info;
            probetarget_flanking_info{2} = load([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_flanking_DNA_sequences' settings.designerName '.mat']).probe_DNA_target_flanking_info;
            probetarget_flanking_info = vertcat(probetarget_flanking_info{:});
            probetarget_flanking_info = sortrows(probetarget_flanking_info,'gene_table_location','descend');
            save([strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_probetarget_flanking_sequences', settings.designerName,'.mat')],'probetarget_flanking_info','-v7.3');%save flanking sequence file if it does not exist already
            if exist([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_flanking_RNA_sequences' settings.designerName '.mat'],'file')        %delete temp mat file if already exists
                delete([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_flanking_RNA_sequences' settings.designerName '.mat'])
            end
            if exist([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_flanking_DNA_sequences' settings.designerName '.mat'],'file')        %delete temp mat file if already exists
                delete([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_flanking_DNA_sequences' settings.designerName '.mat'])
            end
        elseif (runRNA)
            probetarget_flanking_info = load([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_flanking_RNA_sequences' settings.designerName '.mat']).probe_RNA_target_flanking_info;
            save([strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_probetarget_flanking_sequences', settings.designerName,'.mat')],'probetarget_flanking_info','-v7.3');%save flanking sequence file if it does not exist already
            if exist([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_flanking_RNA_sequences' settings.designerName '.mat'],'file')        %delete temp mat file if already exists
                delete([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_flanking_RNA_sequences' settings.designerName '.mat'])
            end
        elseif (runDNA)
            probetarget_flanking_info = load([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_flanking_DNA_sequences' settings.designerName '.mat']).probe_RNA_target_flanking_info;
            save([strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_probetarget_flanking_sequences', settings.designerName,'.mat')],'probetarget_flanking_info','-v7.3');%save flanking sequence file if it does not exist already
            if exist([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_flanking_DNA_sequences' settings.designerName '.mat'],'file')        %delete temp mat file if already exists
                delete([settings.FolderRootName filesep '(' TranscriptName ')' '_' settings.rootName '_flanking_DNA_sequences' settings.designerName '.mat'])
            end
        end
    end
end
end


