function [RNAdbParser, DNAdbParser,probetarget_flanking_info] = A3_BlastDBCMD_JH(settings,gene_table)
%% This takes input data and gets flanking sequence of blast hits in blastdb
Organism = settings.Organism;
runDNA = settings.BLASTdna;
runRNA = settings.BLASTrna;
NamesZ = gene_table.Name;
A1_MakeBlastDB_JH(settings)
% if (settings.BLASTdna)
% DNA_IDs = find(~ismember(Names,settings.dbDNAparser));%IDs
% else
% DNA_IDs = [];
% end
% if (settings.BLASTrna)
% NonDNA_IDs = find(ismember(Names,settings.dbRNAParser));%IDs
% else
% NonDNA_IDs =[];
% end
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
sseqids_preserved = 1;
% gene_table = sortrows(gene_table,[7 6],'ascend');
% gene_table = gene_table(gene_table.Match>=settings.MinHomologySearchTargetSize,:);
% MinusStrandedHits = find(contains(gene_table.Strand,'Minus'));
% gene_table_NamesZ = convertCharsToStrings(gene_table.Name);
% gene_table_uniNamesZ = extractBefore(gene_table_NamesZ,'.');
% if (sum(ismissing(gene_table_uniNamesZ))>0)
%     gene_table_uniNamesZ(ismissing(gene_table_uniNamesZ)) = extractBefore(gene_table_NamesZ(ismissing(gene_table_uniNamesZ)),' ');
% end
% contains_RNA = find(ismember(gene_table_uniNamesZ,settings.RNAparser));
% RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);
% gene_table = gene_table(setdiff(1:size(gene_table,1),RNA_MissedFilteredHits),:);
% gene_table.Ax = min(gene_table.SubjectIndices,[],2);
% gene_table.Bx = max(gene_table.SubjectIndices,[],2);
% gene_table = sortrows(gene_table,[7 13],'ascend');
probetarget_flanking_info = table('Size',[size(gene_table,1) 3],'VariableNames',{'Name','TargetSequence_5primeTo3prime','ProbeRevCompSequence_5primeTo3prime'},'VariableTypes',{'string','string','string'});
output_PairwiseTrimmedFlankingFile = strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_ProbeTarget_SequencesWithFlanking', settings.designerName,'.mat');
if (or(BLASTdb_DNAexists,BLASTdb_RNAexists))
    if (runDNA)
        if (BLASTdb_DNAexists)
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
                    error(msg)
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
            if (~isfile([strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_ProbeTarget_SequencesWithFlanking', settings.designerName,'.mat')]))
                dna_entries = find(ismember(NamesZ,DNAdbParser));
                if (~isempty(dna_entries))
                    if sseqids_preserved
                    DNA_sseq_id_dictionary = dictionary(DNAdbParser',1:size(DNAdb_sseqids_Table,1)); 
                    else
                    DNA_sseq_id_dictionary = dictionary(convertCharsToStrings(DNAdb_sseqids_Table.Name)',1:size(DNAdb_sseqids_Table,1));
                    end
                    corresponding_rows_DNA = DNA_sseq_id_dictionary(table2array(gene_table(dna_entries,{'Name'})));
                    corresponding_rows_DNA_i = dna_entries(1:length(corresponding_rows_DNA));
                    if sseqids_preserved
                        DNA_sseq_ids = convertCharsToStrings(DNAdb_sseqids_Table(corresponding_rows_DNA,'AccessionNumber').AccessionNumber);
                    else
                        DNA_sseq_ids = convertCharsToStrings(DNAdb_sseqids_Table(corresponding_rows_DNA,'sseqid').sseqid);
                    end
                    DNA_length_low = gene_table(corresponding_rows_DNA_i,'SubjectIndices').SubjectIndices(:,1) - gene_table(corresponding_rows_DNA_i,'QueryIndices').QueryIndices(:,1) + 1;
                    DNA_length_high = gene_table(corresponding_rows_DNA_i,'SubjectIndices').SubjectIndices(:,2) - gene_table(corresponding_rows_DNA_i,'QueryIndices').QueryIndices(:,2) + cellfun(@length,gene_table(corresponding_rows_DNA_i,'ProbeSequence').ProbeSequence);
                    DNA_length_max = DNAdb_sseqids_Table(corresponding_rows_DNA,'seqlength').seqlength;
                    DNA_length_lower = DNA_length_low;
                    DNA_length_lower(DNA_length_lower<1) = 1;
                    DNA_length_upper = min([DNA_length_high DNA_length_max],[],2);
                    DNA_Hit_length_ranges = join(string([DNA_length_low DNA_length_upper]),'-');
                    DNA_Hit_strands = convertCharsToStrings(lower(extractAfter(gene_table.Strand(corresponding_rows_DNA_i),'/')));
                    DNA_getflank_command_blocks = [DNA_sseq_ids'; DNA_Hit_length_ranges'; DNA_Hit_strands'];
                    DNA_getflank_command = join(DNA_getflank_command_blocks," ",1)';
                    inputDNAFile = strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName, '_getflanksequencesDNA',settings.designerName,'.csv');
                    outputDNAFile = strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_dnaFlanks', settings.designerName,'.fasta');
                    if (~isfile(inputDNAFile))
                        writematrix(DNA_getflank_command,inputDNAFile);
                    end
                    if (~isfile(outputDNAFile))
                        outfmt_flank = '%f';
                        if (~ismac)
                        args_DNA_flanking = sprintf('%s -db %s -dbtype %s -entry_batch %s -out %s -outfmt %s', "blastdbcmd",inputDatabaseDNA, dbtype, inputDNAFile, outputDNAFile,strcat('"',outfmt_flank,'"'));
                        else
                        args_DNA_flanking = sprintf('%s -db %s -dbtype %s -entry_batch %s -out %s -outfmt %s', "blastdbcmd",inputDatabaseDNA, dbtype,strcat(char(39),inputDNAFile,char(39)), strcat(char(39),outputDNAFile,char(39)),strcat('"',outfmt_flank,'"'));
                        end
                        [status,msg] = system(args_DNA_flanking);
                        if status
                            error(msg)
                            fprintf('\n')
                            fprintf('\n')
                        end
                    end
                    flanking_info_DNA = struct2table(fastaread(outputDNAFile));
                    flanking_info_DNA.Header = extractAfter(flanking_info_DNA.Header," ");
                    flanking_info_DNA = renamevars(flanking_info_DNA,'Header','Name');
                    flanking_info_DNA = renamevars(flanking_info_DNA,'Sequence','TargetSequence');
                    %DNA_length_lower  gene_table(ni,:).SubjectIndices(1)  gene_table(ni,:).SubjectIndices(2) DNA_length_upper
                    DNA_target_Block_1 = convertCharsToStrings(arrayfun(@(ni) flanking_info_DNA(ni,:).TargetSequence{:}(1:gene_table(dna_entries(ni),:).SubjectIndices(1)-DNA_length_lower(ni)),1:length(dna_entries),'Un',0));
                    DNA_target_Block_2 = convertCharsToStrings(arrayfun(@(ni) strrep(gene_table.Alignment{ni}(3,:),'-','N'),dna_entries,'Un',0));
                    DNA_target_Block_3 = convertCharsToStrings(arrayfun(@(ni) flanking_info_DNA(ni,:).TargetSequence{:}(end-DNA_length_upper(ni)+gene_table(dna_entries(ni),:).SubjectIndices(2)+1:end),1:length(dna_entries),'Un',0));
                    DNA_target_matching_flank_Trimmed_blocks = [reshape(DNA_target_Block_1,[],1) reshape(DNA_target_Block_2,[],1) reshape(DNA_target_Block_3,[],1)]';
                    DNA_target_matching_flank_Trimmed = join(DNA_target_matching_flank_Trimmed_blocks,"",1)';
                    DNA_probe_Block_1 = convertCharsToStrings(arrayfun(@(ni) gene_table(dna_entries(ni),:).ProbeSequence{:}(gene_table(dna_entries(ni),:).QueryIndices(1)-gene_table(dna_entries(ni),:).SubjectIndices(1)+DNA_length_lower(ni):gene_table(dna_entries(ni),:).QueryIndices(1)-1),1:length(dna_entries),'Un',0));
                    DNA_probe_Block_2 =  convertCharsToStrings(arrayfun(@(ni) strrep(gene_table.Alignment{ni}(1,:),'-','N'),dna_entries,'Un',0));
                    DNA_probe_Block_3 = convertCharsToStrings(arrayfun(@(ni) gene_table(dna_entries(ni),:).ProbeSequence{:}(gene_table(ni,:).QueryIndices(2)+1:gene_table(ni,:).QueryIndices(2)+DNA_length_upper(ni)-gene_table(dna_entries(ni),:).SubjectIndices(2)),1:length(dna_entries),'Un',0));
                    DNA_probe_matching_flank_Trimmed_blocks = [reshape(DNA_probe_Block_1,[],1) reshape(DNA_probe_Block_2,[],1) reshape(DNA_probe_Block_3,[],1)]';
                    DNA_probe_matching_flank_Trimmed = join(DNA_probe_matching_flank_Trimmed_blocks,"",1)';
                    if sseqids_preserved
                    probetarget_flanking_info(dna_entries,"Name") = array2table(DNAdbParser(corresponding_rows_DNA));
                    else
                    probetarget_flanking_info(dna_entries,"Name") = flanking_info_DNA(:,"Name");
                    end
                    probetarget_flanking_info(dna_entries,'TargetSequence_5primeTo3prime') = array2table(DNA_target_matching_flank_Trimmed);
                    probetarget_flanking_info(dna_entries,'ProbeRevCompSequence_5primeTo3prime') = array2table(DNA_probe_matching_flank_Trimmed);
                    fprintf('Acquired blast DNA hit target flanking sequences')
                    fprintf('\n')
                    fprintf('\n')
                end
            end
        end
    end
    if (runRNA)
        if (BLASTdb_RNAexists)
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
                    error(msg)
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
            %RNAdb_sseqids_Table(:,'AccessionNumber') = array2table(extractBefore(table2array(RNAdb_sseqids_Table(:,'Name')),' '));
            end
            if (~isfile([strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_ProbeTarget_SequencesWithFlanking', settings.designerName,'.mat')]))
                rna_entries = find(ismember(NamesZ,RNAdbParser));
                if (~isempty(rna_entries))
                    if sseqids_preserved
                    RNA_sseq_id_dictionary = dictionary(RNAdbParser',1:size(RNAdb_sseqids_Table,1)); 
                    else
                    RNA_sseq_id_dictionary = dictionary(convertCharsToStrings(RNAdb_sseqids_Table.Name)',1:size(RNAdb_sseqids_Table,1));
                    end
                    corresponding_rows_RNA = RNA_sseq_id_dictionary(table2array(gene_table(rna_entries,{'Name'})));
                    corresponding_rows_RNA_i = rna_entries(1:length(corresponding_rows_RNA));
                    if sseqids_preserved
                        RNA_sseq_ids = convertCharsToStrings(RNAdb_sseqids_Table(corresponding_rows_RNA,'AccessionNumber').AccessionNumber);
                    else
                        RNA_sseq_ids = convertCharsToStrings(RNAdb_sseqids_Table(corresponding_rows_RNA,'sseqid').sseqid);
                    end
                    RNA_length_low = gene_table(corresponding_rows_RNA_i,'SubjectIndices').SubjectIndices(:,1) - gene_table(corresponding_rows_RNA_i,'QueryIndices').QueryIndices(:,1) + 1;
                    RNA_length_high = gene_table(corresponding_rows_RNA_i,'SubjectIndices').SubjectIndices(:,2) - gene_table(corresponding_rows_RNA_i,'QueryIndices').QueryIndices(:,2) + cellfun(@length,gene_table(corresponding_rows_RNA_i,'ProbeSequence').ProbeSequence);
                    RNA_length_max = RNAdb_sseqids_Table(corresponding_rows_RNA,'seqlength').seqlength;
                    RNA_length_lower = RNA_length_low;
                    RNA_length_lower(RNA_length_lower<1) = 1;
                    RNA_length_upper = min([RNA_length_high RNA_length_max],[],2);
                    RNA_Hit_length_ranges = join(string([RNA_length_lower RNA_length_upper]),'-');
                    RNA_Hit_strands = convertCharsToStrings(lower(extractAfter(gene_table.Strand(corresponding_rows_RNA_i),'/')));
                    RNA_getflank_command_blocks = [RNA_sseq_ids'; RNA_Hit_length_ranges'; RNA_Hit_strands'];
                    RNA_getflank_command = join(RNA_getflank_command_blocks," ",1)';
                    inputRNAFile = strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName, '_getflanksequencesRNA',settings.designerName,'.csv');
                    outputRNAFile = strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_rnaFlanks', settings.designerName,'.fasta');
                    if (~isfile(inputRNAFile))
                        writematrix(RNA_getflank_command,inputRNAFile);
                    end
                    if (~isfile(outputRNAFile))
                        outfmt_flank = '%f';
                        if (~ismac)
                        args_RNA_flanking = sprintf('%s -db %s -dbtype %s -entry_batch %s -out %s', "blastdbcmd",inputDatabaseRNA, dbtype, inputRNAFile, outputRNAFile);
                        else
                        args_RNA_flanking = sprintf('%s -db %s -dbtype %s -entry_batch %s -out %s -outfmt %s', "blastdbcmd",inputDatabaseRNA, dbtype, strcat(char(39),inputRNAFile,char(39)), strcat(char(39),outputRNAFile,char(39)),strcat('"',outfmt_flank,'"'));
                        end
                        [status,msg] = system(args_RNA_flanking);
                        if status
                            error(msg)
                            fprintf('\n')
                            fprintf('\n')
                        end
                    end
                    flanking_info_RNA = struct2table(fastaread(outputRNAFile));
                    flanking_info_RNA.Header = extractAfter(flanking_info_RNA.Header," ");
                    flanking_info_RNA = renamevars(flanking_info_RNA,'Header','Name');
                    flanking_info_RNA = renamevars(flanking_info_RNA,'Sequence','TargetSequence');
                    %RNA_length_lower  gene_table(ni,:).SubjectIndices(1)  gene_table(ni,:).SubjectIndices(2) RNA_length_upper
                    RNA_target_Block_1 = convertCharsToStrings(arrayfun(@(ni) flanking_info_RNA(ni,:).TargetSequence{:}(1:gene_table(rna_entries(ni),:).SubjectIndices(1)-RNA_length_lower(ni)),1:length(rna_entries),'Un',0));
                    RNA_target_Block_2 = convertCharsToStrings(arrayfun(@(ni) strrep(gene_table.Alignment{ni}(3,:),'-','N'),rna_entries,'Un',0));
                    RNA_target_Block_3 = convertCharsToStrings(arrayfun(@(ni) flanking_info_RNA(ni,:).TargetSequence{:}(end-RNA_length_upper(ni)+gene_table(rna_entries(ni),:).SubjectIndices(2)+1:end),1:length(rna_entries),'Un',0));
                    RNA_target_matching_flank_Trimmed_blocks = [reshape(RNA_target_Block_1,[],1) reshape(RNA_target_Block_2,[],1) reshape(RNA_target_Block_3,[],1)]';
                    RNA_target_matching_flank_Trimmed = join(RNA_target_matching_flank_Trimmed_blocks,"",1)';
                    RNA_probe_Block_1 = convertCharsToStrings(arrayfun(@(ni) gene_table(rna_entries(ni),:).ProbeSequence{:}(gene_table(rna_entries(ni),:).QueryIndices(1)-gene_table(rna_entries(ni),:).SubjectIndices(1)+RNA_length_lower(ni):gene_table(rna_entries(ni),:).QueryIndices(1)-1),1:length(rna_entries),'Un',0));
                    RNA_probe_Block_2 =  convertCharsToStrings(arrayfun(@(ni) strrep(gene_table.Alignment{ni}(1,:),'-','N'),rna_entries,'Un',0));
                    RNA_probe_Block_3 = convertCharsToStrings(arrayfun(@(ni) gene_table(rna_entries(ni),:).ProbeSequence{:}(gene_table(rna_entries(ni),:).QueryIndices(2)+1:gene_table(rna_entries(ni),:).QueryIndices(2)+RNA_length_upper(ni)-gene_table(rna_entries(ni),:).SubjectIndices(2)),1:length(rna_entries),'Un',0));
                    RNA_probe_matching_flank_Trimmed_blocks = [reshape(RNA_probe_Block_1,[],1) reshape(RNA_probe_Block_2,[],1) reshape(RNA_probe_Block_3,[],1)]';
                    RNA_probe_matching_flank_Trimmed = join(RNA_probe_matching_flank_Trimmed_blocks,"",1)';
                    if sseqids_preserved
                    probetarget_flanking_info(rna_entries,"Name") = array2table(RNAdbParser(corresponding_rows_RNA));
                    else
                    probetarget_flanking_info(rna_entries,"Name") = flanking_info_RNA(:,"Name");
                    end
                    probetarget_flanking_info(rna_entries,'TargetSequence_5primeTo3prime') = array2table(RNA_target_matching_flank_Trimmed);
                    probetarget_flanking_info(rna_entries,'ProbeRevCompSequence_5primeTo3prime') = array2table(RNA_probe_matching_flank_Trimmed);
                    fprintf('Acquired blast RNA hit target flanking sequences')
                    fprintf('\n')
                    fprintf('\n')
                end
            end
        end
    end
end
if (~isfile([strcat(settings.FolderRootName,filesep,'(',  settings.GeneName ,')_', settings.rootName,'_ProbeTarget_SequencesWithFlanking', settings.designerName,'.mat')]))
    save(output_PairwiseTrimmedFlankingFile,'probetarget_flanking_info','-v7.3');%save flanking sequence file if it does not exist already
end
end


