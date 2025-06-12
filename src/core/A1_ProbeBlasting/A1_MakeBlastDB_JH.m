function A1_MakeBlastDB_JH(settings)
%% This takes input data and makes blastdb if they do not exist
Organism = settings.Organism;
runDNA = settings.BLASTdna;
runRNA = settings.BLASTrna;
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
fprintf('\n')
fprintf("Checking if BLAST Database exists")
fprintf('\n')
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
if (or(~BLASTdb_DNAexists,~BLASTdb_RNAexists))
    FNAFiles_FNA = dir([settings.SEQdbRoot '**/*.fna']);
    FNAFiles_FA = dir([settings.SEQdbRoot '**/*.fa']);
    if (~isempty(FNAFiles_FNA)&&~isempty(FNAFiles_FA))
        FNAFiles = cat(1,FNAFiles_FNA,FNAFiles_FA);
    elseif (~isempty(FNAFiles_FNA))
        FNAFiles = FNAFiles_FNA;
    elseif (~isempty(FNAFiles_FA))
        FNAFiles = FNAFiles_FA;
    end
    if (runDNA)
        if (~BLASTdb_DNAexists)
            fprintf('\n')
            fprintf("Making DNA BLAST database")
            fprintf('\n')
            dbtype = extractBefore("nucleotide", 5); % "nucl" or "prot" shorthand
            BLASTdb_DNAoptions = bioinfo.blastplus.MakeDatabaseOptions;
            if (strcmp(settings.referenceType,'RefSeq'))
                BLASTdb_DNAoptions.Title = strcat(settings.Organism,'_NCBI_genomic');
                dna_files = find(~contains({FNAFiles.name},'rna').*contains({FNAFiles.name},'genomic'),1);
            elseif (strcmp(settings.referenceType,'ENSEMBL'))
                BLASTdb_DNAoptions.Title = strcat(settings.Organism,'_ENSEMBL_genomic');
                dna_files = find(contains({FNAFiles.name},'.dna'));
            end
            if (~isempty(dna_files))
                DNAFile = strcat('"',strjoin(string(cellfun(@(x) strcat(settings.SEQdbRoot,x),{FNAFiles(dna_files).name},'Un',0))," "),'"');
            elseif (isempty(dna_files))
                msg = 'Error. No genome file is in the target blast database directory for making the DNA BLASTDB.';
                error(msg);
            end
            filesizeDNA = sum([FNAFiles(dna_files).bytes]);
            inputFileDNA = DNAFile;
            if (strcmp(filesep,'/'))
                outputDatabaseDNA = strcat(pwd,filesep,db1);
                outputRoot = strcat(pwd,filesep,settings.SEQdbRoot);
            else
                outputDatabaseDNA = strcat(pwd,filesep,strrep(db1,'/',filesep));
                outputRoot = strcat(pwd,filesep,strrep(settings.SEQdbRoot,'/',filesep));
            end
            outputRoot_temp = strrep(outputRoot,' ','_');
            outputDatabaseDNA_temp = strrep(outputDatabaseDNA,' ','_');
            opts = bioinfo.internal.CLIOptions.create('bioinfo.blastplus.MakeDatabaseOptions',BLASTdb_DNAoptions);
            args = sprintf('%s -dbtype %s -in %s -out %s %s -max_file_sz %d', "makeblastdb",dbtype, inputFileDNA, outputDatabaseDNA_temp, opts.getCommand(),filesizeDNA);
            [status, result] = system(args,'-echo');
            if status
                error(message('bioinfo:blastplus:blastplus:NativeErrorOrWarning','blastplusdatabase', result));
            end
        end
    end
    if (runRNA)
        if (~BLASTdb_RNAexists)
            fprintf('\n')
            fprintf("Making RNA BLAST database")
            fprintf('\n')
            dbtype = extractBefore("nucleotide", 5); % "nucl" or "prot" shorthand
            BLASTdb_RNAoptions = bioinfo.blastplus.MakeDatabaseOptions;
            if (strcmp(settings.referenceType,'RefSeq'))
                BLASTdb_RNAoptions.Title = strcat(settings.Organism,'_NCBI_transcript');
                rna_files = find(contains({FNAFiles.name},'rna').*~contains({FNAFiles.name},'genomic'),1);
            elseif (strcmp(settings.referenceType,'ENSEMBL'))
                BLASTdb_RNAoptions.Title = strcat(settings.Organism,'_ENSEMBL_transcript');
                rna_files = find(~contains({FNAFiles.name},'.dna'));
            end
            if (~isempty(rna_files))
                RNAFile = strcat('"',strjoin(string(cellfun(@(x) strcat(settings.SEQdbRoot,x),{FNAFiles(rna_files).name},'Un',0))," "),'"');
            elseif (isempty(rna_files))
                msg = 'Error. No transcriptome file is in the target blast database directory for making the RNA BLASTDB.';
                error(msg);
            end
            filesizeRNA = sum([FNAFiles(rna_files).bytes]);
            inputFileRNA = RNAFile;
            if (strcmp(filesep,'/'))
                outputDatabaseRNA = strcat(pwd,filesep,db2);
                outputRoot = strcat(pwd,filesep,settings.SEQdbRoot);
            else
                outputDatabaseRNA = strcat(pwd,filesep,strrep(db2,'/',filesep));
                outputRoot = strcat(pwd,filesep,strrep(settings.SEQdbRoot,'/',filesep));
            end
            outputRoot_temp = strrep(outputRoot,' ','_');
            outputDatabaseRNA_temp = strrep(outputDatabaseRNA,' ','_');
            opts = bioinfo.internal.CLIOptions.create('bioinfo.blastplus.MakeDatabaseOptions',BLASTdb_RNAoptions);
            args = sprintf('%s -dbtype %s -in %s -out %s %s -max_file_sz %d', "makeblastdb",dbtype, inputFileRNA, outputDatabaseRNA_temp, opts.getCommand(),filesizeRNA);
            [status, result] = system(args,'-echo');
            if status
                error(message('bioinfo:blastplus:blastplus:NativeErrorOrWarning','blastplusdatabase', result));
            end
        end
        sep_locations = strfind(outputRoot,filesep);
        space_locations = find(isspace(outputRoot));
        outputRoot_temp_char =  char(outputRoot_temp);
        if (~isempty(union(space_locations,sep_locations)))
            files_in_temp_folder = dir(outputRoot_temp);
            for k = find([files_in_temp_folder.isdir]==0)
                [status,msg] = movefile(strcat(outputRoot_temp,files_in_temp_folder(k).name), strcat(outputRoot,files_in_temp_folder(k).name));
            end
            rmdir(outputRoot_temp_char(1:sep_locations(find(sep_locations>space_locations(1),1))-1),'s')
        end
    end
end
end


