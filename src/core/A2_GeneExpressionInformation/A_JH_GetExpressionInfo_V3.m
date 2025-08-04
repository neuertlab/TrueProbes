function [ExpressionMatrix,get_expression_time] = A_JH_GetExpressionInfo_V3(gene_table,settings,input_gene_expression_file_locations)
% This function gets the gene expression information for the on/off-targets that probes bind genome-wide.
Organism = settings.Organism;
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
%% Check if Expression IDs are in terms of NCBI or ENSEMBL
%% Get gene expression levels in all tissues and cell types
Null_DNAcopynum = settings.DNAPloidy;
Null_RNAcopynum = settings.UniformRNAExpression;
Expression_Null = settings.DoAllGenesHaveSameExpression; %do genes have same expressipon
get_exp_start = tic;
if (Expression_Null==1)
    DNAcopynum = Null_DNAcopynum;
    ExpressionMatrix_Gene = zeros(length(uniNames),1);
    ExpressionMatrix_Gene(DNA_IDs) = DNAcopynum;
    ExpressionMatrix_Gene(NonDNA_IDs) = Null_RNAcopynum;
    ExpressionMatrix_Transcript = zeros(length(uniNames),1);
    ExpressionMatrix_Transcript(DNA_IDs) = DNAcopynum;
    ExpressionMatrix_Transcript(NonDNA_IDs) = Null_RNAcopynum;
else
    if (isKey(settings.LocRoot_FASTA,Organism))
        % JC = fileread('homo_sapiens.json');
        % JCi = jsondecode(JC);
        EMBLtoNCBI_opts = detectImportOptions(settings.EMBLtoNCBI(Organism),'FileType','delimitedtext','Delimiter','\t');
        EMBLtoNCBI_stableIDs = readtable(settings.EMBLtoNCBI(Organism),EMBLtoNCBI_opts);
        %filter out source identity and xref identity NaN?
        EMBLtoNCBI_mappings = EMBLtoNCBI_stableIDs(sum(cell2mat(cellfun(@(x) contains(EMBLtoNCBI_stableIDs.xref,x),{'NM','NR','XM','XR'},'Un',0)),2)==1,{'gene_stable_id','transcript_stable_id','xref','source_identity','xref_identity','info_type','db_name'});
        EMBLtoNCBI_mappings_Matched_P0 = EMBLtoNCBI_mappings(~or(isnan(EMBLtoNCBI_mappings.source_identity),isnan(EMBLtoNCBI_mappings.xref_identity)),:);
        EMBLtoNCBI_mappings_Matched_P0 = sortrows(EMBLtoNCBI_mappings_Matched_P0,[2 3 1 4 5 6 7],{'ascend','ascend','ascend','descend','descend','ascend','ascend'});
        unique_stable_ids = unique(EMBLtoNCBI_mappings_Matched_P0.transcript_stable_id,'stable');
        temp_loc_stable_dict = dictionary(convertCharsToStrings(unique_stable_ids'),1:length(unique_stable_ids));
        temp_inverse_stable_dict = dictionary(1:length(unique_stable_ids),convertCharsToStrings(unique_stable_ids'));
        num_matches_loc_stable = temp_loc_stable_dict(EMBLtoNCBI_mappings_Matched_P0.transcript_stable_id);
        repeat_locs_stable = find(groupcounts(num_matches_loc_stable)>1)';
        EMBLtoNCBI_mappings_Matched_Single_P0 = EMBLtoNCBI_mappings_Matched_P0(~ismember(EMBLtoNCBI_mappings_Matched_P0.transcript_stable_id,temp_inverse_stable_dict(repeat_locs_stable)),:);
        EMBLtoNCBI_mappings_Matched_Multi_P0 = EMBLtoNCBI_mappings_Matched_P0(ismember(EMBLtoNCBI_mappings_Matched_P0.transcript_stable_id,temp_inverse_stable_dict(repeat_locs_stable)),:);
        multi_mapped_EMBL_transcripts = unique(EMBLtoNCBI_mappings_Matched_Multi_P0.transcript_stable_id);
        for kk = 1:length(multi_mapped_EMBL_transcripts)
            rows = find(ismember(EMBLtoNCBI_mappings_Matched_Multi_P0.transcript_stable_id,multi_mapped_EMBL_transcripts{kk}));
            identity = mean([EMBLtoNCBI_mappings_Matched_Multi_P0(rows,:).source_identity EMBLtoNCBI_mappings_Matched_Multi_P0(rows,:).xref_identity],2);
            EMBLtoNCBI_mappings_Matched_Multi_P0(rows(identity<max(identity)),:) = [];
        end
        EMBLtoNCBI_mappings_Matched_Multi_P0 = sortrows(EMBLtoNCBI_mappings_Matched_Multi_P0,"xref","ascend");
        EMBLtoNCBI_mappings_Matched_P1 = [EMBLtoNCBI_mappings_Matched_Multi_P0;EMBLtoNCBI_mappings_Matched_Single_P0];
        unique_xref = unique(EMBLtoNCBI_mappings_Matched_P1.xref,'stable');
        temp_loc_dict_xref = dictionary(convertCharsToStrings(unique_xref'),1:length(unique_xref));
        temp_inverse_dict_xref = dictionary(1:length(unique_xref),convertCharsToStrings(unique_xref'));
        num_matches_loc_xref = temp_loc_dict_xref(EMBLtoNCBI_mappings_Matched_P1.xref);
        repeat_locs_xref = find(groupcounts(num_matches_loc_xref)>1)';
        EMBLtoNCBI_mappings_Matched_Single_P1 = EMBLtoNCBI_mappings_Matched_P1(~ismember(EMBLtoNCBI_mappings_Matched_P1.xref,temp_inverse_dict_xref(repeat_locs_xref)),:);
        EMBLtoNCBI_mappings_Matched_Multi_P1 = EMBLtoNCBI_mappings_Matched_P1(ismember(EMBLtoNCBI_mappings_Matched_P1.xref,temp_inverse_dict_xref(repeat_locs_xref)),:);
        multi_mapped_xref = unique(EMBLtoNCBI_mappings_Matched_Multi_P1.xref);
        for kk = 1:length(multi_mapped_xref)
            rows = find(ismember(EMBLtoNCBI_mappings_Matched_Multi_P1.xref,multi_mapped_xref{kk}));
            identity = mean([EMBLtoNCBI_mappings_Matched_Multi_P1(rows,:).source_identity EMBLtoNCBI_mappings_Matched_Multi_P1(rows,:).xref_identity],2);
            EMBLtoNCBI_mappings_Matched_Multi_P1(rows(identity<max(identity)),:) = [];
        end
        EMBLtoNCBI_mappings_Matched_Multi_P1 = sortrows(EMBLtoNCBI_mappings_Matched_Multi_P1,"xref","ascend");
        EMBLtoNCBI_mappings_Matched = [EMBLtoNCBI_mappings_Matched_Multi_P1;EMBLtoNCBI_mappings_Matched_Single_P1];
        EMBLgeneIds = unique(EMBLtoNCBI_mappings_Matched.gene_stable_id);
        EMBLtranscriptIds = unique(EMBLtoNCBI_mappings_Matched.transcript_stable_id);
        NCBItranscriptIds = unique(EMBLtoNCBI_mappings_Matched.xref);
        EMBLgeneIds_dict = dictionary(string(EMBLgeneIds)',1:length(EMBLgeneIds));
        EMBLtranscriptIds_dict = dictionary(string(EMBLtranscriptIds)',1:length(EMBLtranscriptIds));
        NCBItranscriptIds_dict = dictionary(string(NCBItranscriptIds)',1:length(NCBItranscriptIds));
        inverseEMBLgeneIds_dict = dictionary(1:length(EMBLgeneIds),string(EMBLgeneIds)');
        inverseEMBLtranscriptIds_dict = dictionary(1:length(EMBLtranscriptIds),string(EMBLtranscriptIds)');
        inverseNCBItranscriptIds_dict = dictionary(1:length(NCBItranscriptIds),string(NCBItranscriptIds)');
        maps = [EMBLgeneIds_dict(EMBLtoNCBI_mappings_Matched.gene_stable_id) ...
            EMBLtranscriptIds_dict(EMBLtoNCBI_mappings_Matched.transcript_stable_id) ...
            NCBItranscriptIds_dict(EMBLtoNCBI_mappings_Matched.xref)];
        geneExpressionDataLocations = struct2table(readstruct(input_gene_expression_file_locations).row);
        geneExpressionDataLocationsVars = struct2table(readstruct(input_gene_expression_file_locations).tracks);
        values = geneExpressionDataLocations(strcmp(geneExpressionDataLocations.Organism,Organism),:);
        values = values(:,~ismissing(values));
        GeneExpressionFiles = [table2cell(values(:,values.Properties.VariableNames(~strcmp(values.Properties.VariableNames,'Organism'))))];
        GeneExpressionFiles = [GeneExpressionFiles{:}];
        GeneExpressionFiles(ismissing(GeneExpressionFiles)) = [];
        if (~isempty(GeneExpressionFiles))
            optsFile = cell(1,length(GeneExpressionFiles));
            outFormat = cell(1,length(GeneExpressionFiles));
            expFileVals = cell(1,length(GeneExpressionFiles));
            expFileVals_table = table('Size', [length(GeneExpressionFiles), 4], 'VariableTypes', ["string", "string", "string", "cell"],'VariableNames',{'source','resolvedLevel','referenceType','id_column_values'});
            for ii = 1:length(GeneExpressionFiles)
                optsFile{ii} = detectImportOptions(GeneExpressionFiles{ii},'FileType','delimitedtext','Delimiter','\t');
                if (contains(GeneExpressionFiles{ii},'bed'))
                    optsFile{ii}.VariableNames = strsplit(geneExpressionDataLocationsVars.bed.columns,',');
                    outFormat{ii} = 'tab';
                else
                end
                expFileVals{ii} = readtable(GeneExpressionFiles{ii},optsFile{ii});
                tmpPosition1 = table2cell(expFileVals{ii}(:,4));
                tmpPosition1(contains(tmpPosition1,'.')) = extractBefore(tmpPosition1(contains(tmpPosition1,'.')),'.');
                tmpPosition2 = table2cell(expFileVals{ii}(:,7));
                tmpPosition2(contains(tmpPosition2,'.')) = extractBefore(tmpPosition2(contains(tmpPosition2,'.')),'.');
                if (sum(ismember(tmpPosition1,NCBItranscriptIds))>0)
                    expFileVals_table(ii, :) = {GeneExpressionFiles{ii},'T','RefSeq', string(tmpPosition1)};
                end
                if (sum(ismember(tmpPosition2,NCBItranscriptIds))>0)
                    expFileVals_table(ii, :) = {GeneExpressionFiles{ii},'T','RefSeq', string(tmpPosition2)};
                end
                if (sum(ismember(tmpPosition1,EMBLgeneIds))>0)
                    expFileVals_table(ii, :) = {GeneExpressionFiles{ii},'G','ENSEMBL', string(tmpPosition1)};
                end
                if (sum(ismember(tmpPosition2,EMBLgeneIds))>0)
                    expFileVals_table(ii, :) = {GeneExpressionFiles{ii},'G','ENSEMBL', string(tmpPosition2)};
                end
                if (sum(ismember(tmpPosition1,EMBLtranscriptIds))>0)
                    expFileVals_table(ii, :) = {GeneExpressionFiles{ii},'T','ENSEMBL', string(tmpPosition1)};
                end
                if (sum(ismember(tmpPosition2,EMBLtranscriptIds))>0)
                    expFileVals_table(ii, :) = {GeneExpressionFiles{ii},'T','ENSEMBL', string(tmpPosition2)};
                end
            end
            positions_ref = cell(1,length(GeneExpressionFiles));
            for ii = 1:length(GeneExpressionFiles)
                %(multiple xref to transcirpt or vis versa need to share values)
                order_out = uniNames(NonDNA_IDs);
                positions_ref{ii} = zeros(length(NonDNA_IDs),1);
                if (strcmp(settings.referenceType,"ENSEMBL")*isequal(table2cell(expFileVals_table(ii, 2:3)),{'G','ENSEMBL'})==1)
                    %G-RES EXPR EMBL    & REF EMBL
                    exp_geneIds = expFileVals_table.id_column_values{ii};
                    mapping_geneIds_to_geneIds = dictionary(exp_geneIds',1:length(exp_geneIds));
                    gtf_transcript_to_gene_names = dictionary(settings.pairedGTF_TranscriptIDs,settings.pairedGTF_GeneIDs);
                    hit_in_gtf = isKey(gtf_transcript_to_gene_names,order_out);
                    positions_ref{ii}(~hit_in_gtf) = NaN;
                    hit_in_gtf_pos = find(hit_in_gtf);
                    order_out_gene_names_in_gtf = gtf_transcript_to_gene_names(order_out(hit_in_gtf));
                    isInExpressionReference = isKey(mapping_geneIds_to_geneIds,order_out_gene_names_in_gtf);
                    positions_ref{ii}(hit_in_gtf_pos(isInExpressionReference)) = mapping_geneIds_to_geneIds(order_out_gene_names_in_gtf(isInExpressionReference));
                    positions_ref{ii}(hit_in_gtf_pos(~isInExpressionReference)) = NaN;
                end
                if (strcmp(settings.referenceType,"RefSeq")*isequal(table2cell(expFileVals_table(ii, 2:3)),{'G','RefSeq'})==1)
                    %G-RES EXPR NCBI    & REF NCBI
                    exp_geneIds = expFileVals_table.id_column_values{ii};
                    mapping_geneIds_to_geneIds = dictionary(exp_geneIds',1:length(exp_geneIds));
                    gtf_transcript_to_gene_names = dictionary(settings.pairedGTF_TranscriptIDs,settings.pairedGTF_GeneIDs);
                    hit_in_gtf = isKey(gtf_transcript_to_gene_names,order_out);
                    positions_ref{ii}(~hit_in_gtf) = NaN;
                    hit_in_gtf_pos = find(hit_in_gtf);
                    order_out_gene_names_in_gtf = gtf_transcript_to_gene_names(order_out(hit_in_gtf));
                    isInExpressionReference = isKey(mapping_geneIds_to_geneIds,order_out_gene_names_in_gtf);
                    positions_ref{ii}(hit_in_gtf_pos(isInExpressionReference)) = mapping_geneIds_to_geneIds(order_out_gene_names_in_gtf(isInExpressionReference));
                    positions_ref{ii}(hit_in_gtf_pos(~isInExpressionReference)) = NaN;
                end
                if (strcmp(settings.referenceType,"RefSeq")*isequal(table2cell(expFileVals_table(ii, 2:3)),{'G','ENSEMBL'})==1)
                    %G-RES EXPR EMBL  & REF  NCBI
                    maps = sortrows(maps,[3 1 2]);
                    maps_unique_NCBItranscriptId_to_unique_EMBLgeneId = unique(maps(:,[3 1]),'rows');
                    exp_geneIds = expFileVals_table.id_column_values{ii};
                    mapping_geneIds_to_geneIds = dictionary(exp_geneIds',1:length(exp_geneIds));
                    isInEMBL_Mapping = isKey(NCBItranscriptIds_dict,order_out);
                    positions_ref{ii}(~isInEMBL_Mapping) = NaN;
                    pos_in_mapping = find(isInEMBL_Mapping);
                    gene_ids_in_EMBL_mapping = inverseEMBLgeneIds_dict(maps_unique_NCBItranscriptId_to_unique_EMBLgeneId(arrayfun(@(zz) find(maps_unique_NCBItranscriptId_to_unique_EMBLgeneId(:,1)==zz,1),NCBItranscriptIds_dict(order_out(isInEMBL_Mapping))),2));
                    isInExpressionReference = isKey(mapping_geneIds_to_geneIds,gene_ids_in_EMBL_mapping);
                    positions_ref{ii}(pos_in_mapping(isInExpressionReference)) = mapping_geneIds_to_geneIds(gene_ids_in_EMBL_mapping(isInExpressionReference));
                    positions_ref{ii}(pos_in_mapping(~isInExpressionReference)) = NaN;
                end
                if (strcmp(settings.referenceType,"ENSEMBL")*isequal(table2cell(expFileVals_table(ii, 2:3)),{'G','RefSeq'})==1)
                    %G-RES EXPR NCBI    & REF EMBL
                    exp_geneIds = expFileVals_table.id_column_values{ii};
                    mapping_geneIds_to_geneIds = dictionary(exp_geneIds',1:length(exp_geneIds));
                    gtf_transcript_to_gene_names = dictionary(settings.pairedGTF_TranscriptIDs,settings.pairedGTF_GeneNames_EMBLonly);
                    hit_in_gtf = isKey(gtf_transcript_to_gene_names,order_out);
                    positions_ref{ii}(~hit_in_gtf) = NaN;
                    hit_in_gtf_pos = find(hit_in_gtf);
                    order_out_gene_names_in_gtf = gtf_transcript_to_gene_names(order_out(hit_in_gtf));
                    isInExpressionReference = isKey(mapping_geneIds_to_geneIds,order_out_gene_names_in_gtf);
                    positions_ref{ii}(hit_in_gtf_pos(isInExpressionReference)) = mapping_geneIds_to_geneIds(order_out_gene_names_in_gtf(isInExpressionReference));
                    positions_ref{ii}(hit_in_gtf_pos(~isInExpressionReference)) = NaN;
                end
                if (strcmp(settings.referenceType,"ENSEMBL")*isequal(table2cell(expFileVals_table(ii, 2:3)),{'T','ENSEMBL'})==1)
                    %T-RES EXPR EMBL    & REF EMBL
                    exp_transcriptIds = expFileVals_table.id_column_values{ii};
                    mapping_transcriptIds_to_transcriptIds = dictionary(exp_transcriptIds',1:length(exp_transcriptIds));
                    isInExpressionReference = isKey(mapping_transcriptIds_to_transcriptIds,order_out);
                    positions_ref{ii}(isInExpressionReference) = mapping_transcriptIds_to_transcriptIds(order_out(isInExpressionReference));
                    positions_ref{ii}(~isInExpressionReference) = NaN;
                end
                if (strcmp(settings.referenceType,"RefSeq")*isequal(table2cell(expFileVals_table(ii, 2:3)),{'T','RefSeq'})==1)
                    %T-RES EXPR NCBI    & REF NCBI
                    exp_transcriptIds = expFileVals_table.id_column_values{ii};
                    mapping_transcriptIds_to_transcriptIds = dictionary(exp_transcriptIds',1:length(exp_transcriptIds));
                    isInExpressionReference = isKey(mapping_transcriptIds_to_transcriptIds,order_out);
                    positions_ref{ii}(isInExpressionReference) = mapping_transcriptIds_to_transcriptIds(order_out(isInExpressionReference));
                    positions_ref{ii}(~isInExpressionReference) = NaN;
                end
                if (strcmp(settings.referenceType,"RefSeq")*isequal(table2cell(expFileVals_table(ii, 2:3)),{'T','ENSEMBL'})==1)
                    %T-RES EXPR EMBL  & REF  NCBI
                    maps = sortrows(maps,[3 2 1]);
                    maps_unique_NCBItranscriptId_to_unique_EMBLtranscriptId = unique(maps(:,[3 2]),'rows');
                    exp_transcriptIds = expFileVals_table.id_column_values{ii};
                    mapping_transcriptIds_to_transcriptIds = dictionary(exp_transcriptIds',1:length(exp_transcriptIds));
                    isInEMBL_Mapping = isKey(NCBItranscriptIds_dict,order_out);
                    positions_ref{ii}(~isInEMBL_Mapping) = NaN;
                    pos_in_mapping = find(isInEMBL_Mapping);
                    transcript_ids_in_EMBL_mapping = inverseEMBLtranscriptIds_dict(maps_unique_NCBItranscriptId_to_unique_EMBLtranscriptId(arrayfun(@(zz) find(maps_unique_NCBItranscriptId_to_unique_EMBLtranscriptId(:,1)==zz,1),NCBItranscriptIds_dict(order_out(isInEMBL_Mapping))),2));
                    isInExpressionReference = isKey(mapping_transcriptIds_to_transcriptIds,transcript_ids_in_EMBL_mapping);
                    positions_ref{ii}(pos_in_mapping(isInExpressionReference)) = mapping_transcriptIds_to_transcriptIds(transcript_ids_in_EMBL_mapping(isInExpressionReference));
                    positions_ref{ii}(pos_in_mapping(~isInExpressionReference)) = NaN;
                end
                if (strcmp(settings.referenceType,"ENSEMBL")*isequal(table2cell(expFileVals_table(ii, 2:3)),{'T','RefSeq'})==1)
                    %T-RES EXPR NCBI    & REF EMBL
                    maps = sortrows(maps,[2 3 1]);
                    maps_unique_EMBLtranscriptId_to_unique_NCBItranscriptId = unique(maps(:,[2 3]),'rows');
                    exp_transcriptIds = expFileVals_table.id_column_values{ii};
                    mapping_transcriptIds_to_transcriptIds = dictionary(exp_transcriptIds',1:length(exp_transcriptIds));
                    isInEMBL_Mapping = isKey(EMBLtranscriptIds_dict,order_out);
                    positions_ref{ii}(~isInEMBL_Mapping) = NaN;
                    pos_in_mapping = find(isInEMBL_Mapping);
                    transcript_ids_in_NCBI_mapping = inverseNCBItranscriptIds_dict(maps_unique_EMBLtranscriptId_to_unique_NCBItranscriptId(arrayfun(@(zz) find(maps_unique_EMBLtranscriptId_to_unique_NCBItranscriptId(:,1)==zz,1),EMBLtranscriptIds_dict(order_out(isInEMBL_Mapping))),2));
                    isInExpressionReference = isKey(mapping_transcriptIds_to_transcriptIds,transcript_ids_in_NCBI_mapping);
                    positions_ref{ii}(pos_in_mapping(isInExpressionReference)) = mapping_transcriptIds_to_transcriptIds(transcript_ids_in_NCBI_mapping(isInExpressionReference));
                    positions_ref{ii}(pos_in_mapping(~isInExpressionReference)) = NaN;
                end
            end
            outs_G = find(strcmp(expFileVals_table.resolvedLevel,'G'));
            outs_T = find(strcmp(expFileVals_table.resolvedLevel,'T'));
            outputExpressionG = cell(1,sum(strcmp(expFileVals_table.resolvedLevel,'G')));
            outputExpressionT = cell(1,sum(strcmp(expFileVals_table.resolvedLevel,'T')));
            for ii = 1:length(GeneExpressionFiles)
                if (strcmp(outFormat{ii},'tab'))
                    if (strcmp(expFileVals_table.resolvedLevel(ii),'G'))
                        outputExpressionG{outs_G==ii} = ndSparse.build([length(uniNames) max(expFileVals{ii}.expCount)],0);
                        outputExpressionG{outs_G==ii}(NonDNA_IDs(~isnan(positions_ref{ii})),:) = cellfun(@str2double,split(expFileVals{ii}.expScores(positions_ref{ii}(~isnan(positions_ref{ii}))),','));
                    end
                    if (strcmp(expFileVals_table.resolvedLevel(ii),'T'))
                        outputExpressionT{outs_T==ii} = ndSparse.build([length(uniNames) max(expFileVals{ii}.expCount)],0);
                        outputExpressionT{outs_T==ii}(NonDNA_IDs(~isnan(positions_ref{ii})),:) = cellfun(@str2double,split(expFileVals{ii}.expScores(positions_ref{ii}(~isnan(positions_ref{ii}))),','));
                    end
                end
            end
            ExpressionMatrix_Transcript = CATnWrapper(outputExpressionT(~cellfun(@isempty,outputExpressionT)),2);
            ExpressionMatrix_Gene = CATnWrapper(outputExpressionG(~cellfun(@isempty,outputExpressionG)),2);
        else
            DNAcopynum = Null_DNAcopynum;
            ExpressionMatrix_Gene = zeros(length(uniNames),1);
            ExpressionMatrix_Gene(DNA_IDs) = DNAcopynum;
            ExpressionMatrix_Gene(NonDNA_IDs) = Null_RNAcopynum;
            ExpressionMatrix_Transcript = zeros(length(uniNames),1);
            ExpressionMatrix_Transcript(DNA_IDs) = DNAcopynum;
            ExpressionMatrix_Transcript(NonDNA_IDs) = Null_RNAcopynum;
        end
    else
        DNAcopynum = Null_DNAcopynum;
        ExpressionMatrix_Gene = zeros(length(uniNames),1);
        ExpressionMatrix_Gene(DNA_IDs) = DNAcopynum;
        ExpressionMatrix_Gene(NonDNA_IDs) = Null_RNAcopynum;
        ExpressionMatrix_Transcript = zeros(length(uniNames),1);
        ExpressionMatrix_Transcript(DNA_IDs) = DNAcopynum;
        ExpressionMatrix_Transcript(NonDNA_IDs) = Null_RNAcopynum;
    end
end
%% check for custom
%% Organism not included default to equal expression profile
unique_gene_names = unique(settings.pairedGTF_GeneIDs,'stable');
unique_isoform_count = cellfun(@(V) sum(ismember(settings.pairedGTF_GeneIDs,V)),unique_gene_names);
iso_count_dictionary = dictionary(unique_gene_names,unique_isoform_count);
gtf_transcript_to_gene_names = dictionary(settings.pairedGTF_TranscriptIDs,settings.pairedGTF_GeneIDs);
has_name_mapping = isKey(gtf_transcript_to_gene_names,uniNames(NonDNA_IDs));
has_name_mapping_pos = find(has_name_mapping);
genes_in_gtf = isKey(iso_count_dictionary,gtf_transcript_to_gene_names(uniNames(NonDNA_IDs(has_name_mapping))));
NonDNA_IDs_in_mapping = NonDNA_IDs(has_name_mapping_pos(genes_in_gtf));
Num_Isoforms_Vector = iso_count_dictionary(gtf_transcript_to_gene_names(uniNames(NonDNA_IDs_in_mapping)));
if (and(~isempty(ExpressionMatrix_Gene),~isempty(ExpressionMatrix_Transcript)))
    ExpressionMatrix_Gene(NonDNA_IDs_in_mapping,:) = ExpressionMatrix_Gene(NonDNA_IDs_in_mapping,:)./repmat(Num_Isoforms_Vector,[1 size(ExpressionMatrix_Gene,2)]);
elseif (and(isempty(ExpressionMatrix_Gene),~isempty(ExpressionMatrix_Transcript)))
    ExpressionMatrix_Gene = ExpressionMatrix_Transcript;
    uniNames_GeneID_in_mapping = gtf_transcript_to_gene_names(uniNames(NonDNA_IDs_in_mapping));
    unique_names_in_mapping = unique(uniNames_GeneID_in_mapping);
    iso_total_count_expression = cellfun(@(V) sum(ExpressionMatrix_Transcript(NonDNA_IDs_in_mapping(ismember(uniNames_GeneID_in_mapping,V)),:),1),unique_names_in_mapping,'Un',0);
    for kk = 1:length(unique_names_in_mapping)
        ExpressionMatrix_Gene(NonDNA_IDs_in_mapping(ismember(uniNames_GeneID_in_mapping,unique_names_in_mapping{kk})),:) = repmat(iso_total_count_expression{kk},[sum(ismember(uniNames_GeneID_in_mapping,unique_names_in_mapping{kk})) 1]);
    end
    ExpressionMatrix_Gene(NonDNA_IDs_in_mapping,:) = ExpressionMatrix_Gene(NonDNA_IDs_in_mapping,:)./repmat(Num_Isoforms_Vector,[1 size(ExpressionMatrix_Gene,2)]);
elseif (and(~isempty(ExpressionMatrix_Gene),isempty(ExpressionMatrix_Transcript)))
    ExpressionMatrix_Gene(NonDNA_IDs_in_mapping,:) = ExpressionMatrix_Gene(NonDNA_IDs_in_mapping,:)./repmat(Num_Isoforms_Vector,[1 size(ExpressionMatrix_Gene,2)]);
    ExpressionMatrix_Transcript = ExpressionMatrix_Gene;
elseif (and(isempty(ExpressionMatrix_Gene),isempty(ExpressionMatrix_Transcript)))
    DNAcopynum = Null_DNAcopynum;
    ExpressionMatrix_Gene = zeros(length(uniNames),1);
    ExpressionMatrix_Gene(DNA_IDs) = DNAcopynum;
    ExpressionMatrix_Gene(NonDNA_IDs) = Null_RNAcopynum;
    ExpressionMatrix_Gene(NonDNA_IDs_in_mapping,:) = ExpressionMatrix_Gene(NonDNA_IDs_in_mapping,:)./repmat(Num_Isoforms_Vector,[1 size(ExpressionMatrix_Gene,2)]);
    ExpressionMatrix_Transcript = zeros(length(uniNames),1);
    ExpressionMatrix_Transcript(DNA_IDs) = DNAcopynum;
    ExpressionMatrix_Transcript(NonDNA_IDs) = Null_RNAcopynum;
end
ExpressionMatrix = [ExpressionMatrix_Gene ExpressionMatrix_Transcript];
get_expression_time = toc(get_exp_start);
end

