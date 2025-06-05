function [ExpressionMatrix,get_expression_time] = A_JH_GetExpressionInfo_V2(gene_table,settings)  
% This function gets the gene expression information for the on/off-targets that probes bind genome-wide.
% This code can get the expression of targets that probes bind in any cell line or cells.
% This information is initially obtained for human normal tissue, cancer tissue, and cell-lines. 

% This code also includes expression data for mice but this might not be stable 
% (i.e. error's might occur due to differences in file structures not included in code).

Null_DNAcopynum = settings.DNAPloidy;
Null_RNAcopynum = settings.UniformRNAExpression;
Expression_Null = settings.DoAllGenesHaveSameExpression; %do genes have same expressipon
Organism = settings.Organism;
HumanOntology = settings.HumanSpecific.Ontology;%Normal/Cancer'
HumanTissueOnly = settings.HumanSpecific.TissueOrTissueAndCellType; % 0 (TissueOnly) 1 (Tissue/Cell Type)
HumanExpGeneOrTransc = settings.HumanSpecific.HumanExpGeneOrTransc; % 1 (Gene/EMBL GENEID) , 0 (Transcript/EMBL Transcript ID)
Human_SelTrack = settings.HumanSpecific.SCOutputTrack;
HumanSCTracks = settings.HumanSpecific.SCTracks;
%CellTypeID = settings.CellType_ExprID;
%expValType = settings.expressionValType;% 1-4 (expCounts,expValues,mean(CellTypeExpValues),one cell types CellTypeExpValues)
gene_table = sortrows(gene_table,[7 6],'ascend');
gene_table = gene_table(gene_table.Match>=settings.MinHomologySearchTargetSize,:);
MinusStrandedHits = find(contains(gene_table.Strand,'Minus'));
if (strcmp(settings.referenceType,'RefSeq'))
RNA_IDs_1 = find(contains(gene_table.Name,'NM_'));
RNA_IDs_2 = find(contains(gene_table.Name,'NR_'));
RNA_IDs_3 = find(contains(gene_table.Name,'XM_'));
RNA_IDs_4 = find(contains(gene_table.Name,'XR_'));
contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
elseif (strcmp(settings.referenceType,'ENSEMBL'))
RNA_IDs_1 = find(contains(gene_table.Name,'EN'));
RNA_IDs_2 = find(contains(gene_table.Name,'EN'));
RNA_IDs_3 = find(contains(gene_table.Name,'EN'));
RNA_IDs_4 = find(contains(gene_table.Name,'EN'));
contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
else
    ss = 1;
% RNA_IDs_1 = find(contains(gene_table.Name,'EN'));
% RNA_IDs_2 = find(contains(gene_table.Name,'EN'));
% RNA_IDs_3 = find(contains(gene_table.Name,'EN'));
% RNA_IDs_4 = find(contains(gene_table.Name,'EN'));
% contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
end
RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);
gene_table = gene_table(setdiff(1:size(gene_table,1),RNA_MissedFilteredHits),:);
gene_table.Ax = min(gene_table.SubjectIndices,[],2);
gene_table.Bx = max(gene_table.SubjectIndices,[],2);
gene_table = sortrows(gene_table,[7 13],'ascend');
Names = unique(gene_table.Name);
Names = convertCharsToStrings(Names);
uniNames = extractBefore(Names,'.');
if (sum(ismissing(uniNames))>0)
uniNames(ismissing(uniNames)) = extractBefore(Names(ismissing(uniNames)),' ');
end
if (strcmp(settings.referenceType,'RefSeq'))
DNA_IDs_1 = find(contains(uniNames,'NC_'));%IDs
DNA_IDs_2 = find(contains(uniNames,'NT_'));%IDs
DNA_IDs_3 = find(contains(uniNames,'NW_'));%IDs
NonDNA_IDs_1 = find(~contains(uniNames,'NC_'));%IDs
NonDNA_IDs_2 = find(~contains(uniNames,'NT_'));%IDs
NonDNA_IDs_3 = find(~contains(uniNames,'NW_'));%
else
ss = 1;
DNA_IDs_1 = find(contains(uniNames,'NC_'));%IDs
DNA_IDs_2 = find(contains(uniNames,'NT_'));%IDs
DNA_IDs_3 = find(contains(uniNames,'NW_'));%IDs
NonDNA_IDs_1 = find(~contains(uniNames,'NC_'));%IDs
NonDNA_IDs_2 = find(~contains(uniNames,'NT_'));%IDs
NonDNA_IDs_3 = find(~contains(uniNames,'NW_'));%
end
DNA_IDs =union(union(DNA_IDs_1,DNA_IDs_2),DNA_IDs_3).';
NonDNA_IDs = intersect(intersect(NonDNA_IDs_1,NonDNA_IDs_2),NonDNA_IDs_3).';
%% Load version of databases     
    if (strcmp(Organism,'Yeast'))
        optsTabulaYeast = detectImportOptions(settings.YeastExpressionFile,'FileType','delimitedtext');
        optsTabulaYeast.VariableNames = settings.ExpressionVariableNames;
        TabulaYeast = tdfread(settings.YeastExpressionFile,'\t',optsTabulaYeast);
        TabulaYeast = cell2struct(struct2cell(TabulaYeast),optsTabulaYeast.VariableNames);
        TabulaYeast = struct2table(TabulaYeast);
        TabulaYeast_AN.tabulamurisBarChart = table2struct(TabulaYeast);
        defined = 1;
    end
    if (strcmp(Organism,'Mouse'))
        optsTabulaMuris = detectImportOptions(settings.MouseExpressionFile,'FileType','delimitedtext');
        optsTabulaMuris.VariableNames = settings.ExpressionVariableNames;
        TabulaMuris = tdfread(settings.MouseExpressionFile,'\t',optsTabulaMuris);
        TabulaMuris = cell2struct(struct2cell(TabulaMuris),optsTabulaMuris.VariableNames);
        TabulaMuris = struct2table(TabulaMuris);
        TabulaMuris_AN.tabulamurisBarChart = table2struct(TabulaMuris);
        defined = 1;
    end
    if (strcmp(Organism,'Human'))
        if (strcmp(HumanOntology,'Cancer'))
            optsTcgaTransc= detectImportOptions(settings.HumanTCGA_TranscExpressionFile,'FileType','delimitedtext');
            optsTcgaGene = detectImportOptions(settings.HumanTCGA_GeneExpressionFile,'FileType','delimitedtext');
            optsTcgaTransc.VariableNames = settings.ExpressionVariableNames;
            optsTcgaGene.VariableNames = settings.ExpressionVariableNames;
            TcgaTransc = tdfread(settings.HumanTCGA_TranscExpressionFile,'\t',optsTcgaTransc);
            TcgaGene = tdfread(settings.HumanTCGA_GeneExpressionFile,'\t',optsTcgaGene);
            TcgaGene = cell2struct(struct2cell(TcgaGene),optsTcgaGene.VariableNames);
            TcgaGene = struct2table(TcgaGene);
            TCGA_Gene_AN.tcgaGeneExpr = table2struct(TcgaGene);
            TcgaTransc = cell2struct(struct2cell(TcgaTransc),optsTcgaTransc.VariableNames);
            TcgaTransc = struct2table(TcgaTransc);
            TCGA_Transc_AN.tcgaTranscExpr = table2struct(TcgaTransc);
        end
        if (strcmp(HumanOntology,'Normal'))
            if (HumanTissueOnly==0)
                optsGtexTransc = detectImportOptions(settings.HumanGTEX_TranscExpressionFile,'FileType','delimitedtext');
                optsGtexGeneV8 = detectImportOptions(settings.HumanGTEX_GeneExpressionFile,'FileType','delimitedtext');
                optsGtexTransc.VariableNames = settings.ExpressionVariableNames;
                optsGtexGeneV8.VariableNames = settings.ExpressionVariableNames2;
                GtexTransc = tdfread(settings.HumanGTEX_TranscExpressionFile,'\t',optsGtexTransc);
                GtexGeneV8 = tdfread(settings.HumanGTEX_GeneExpressionFile,'\t',optsGtexGeneV8);
                GtexGeneV8 = cell2struct(struct2cell(GtexGeneV8),optsGtexGeneV8.VariableNames);
                GtexGeneV8 = struct2table(GtexGeneV8);
                GTEx_GeneV8_AN.gtexGeneV8 = table2struct(GtexGeneV8);
                GtexTransc = cell2struct(struct2cell(GtexTransc),optsGtexTransc.VariableNames);
                GtexTransc = struct2table(GtexTransc);
                GTEx_Transc_AN.gtexTranscExpr = table2struct(GtexTransc);
            else
                scTracksBedFiles = settings.scTracksBedFiles;
                scTracks = settings.scTracks;
                optsHumanSS.VariableNames = settings.ExpressionVariableNames2;
                optsGeneSS = cell(1,length(HumanSCTracks));
                ssGene = cell(1,length(HumanSCTracks));
                scTrack_Gene_AN = cell(1,length(HumanSCTracks));
                for j=HumanSCTracks
                    optsGeneSS{j} = detectImportOptions(strcat('DatabaseData/',scTracksBedFiles{j},'.bed'),'FileType','delimitedtext');
                    optsGeneSS{j}.VariableNames = optsHumanSS.VariableNames;
                    ssGene{j} = tdfread(strcat('DatabaseData/',scTracksBedFiles{j},'.bed'),'\t',optsGeneSS{j});
                    ssGene{j} = cell2struct(struct2cell(ssGene{j}),optsGeneSS{j}.VariableNames);
                    ssGene{j} = struct2table(ssGene{j});
                    scTrack_Gene_AN{j}.(scTracks{j}) = table2struct(ssGene{j});        
                end
            end 
            defined = 1;
        end
        if (defined ==0)
            try
                optsCustomTransc = detectImportOptions(settings.Custom_TranscExpressionFile,'FileType','delimitedtext');
                optsCustomTransc.VariableNames = settings.CustomTranscExpressionVariableNames;
                CustomTransc = tdfread(settings.Custom_TranscExpressionFile,'\t',optsCustomTransc);
                CustomTransc = cell2struct(struct2cell(CustomTransc),optsCustomTransc.VariableNames);
                CustomTransc = struct2table(CustomTransc);
                CustomTransc_AN.TranscExpr = table2struct(CustomTransc);
            catch
            end
            try
                optsCustomGene = detectImportOptions(settings.Custom_GeneExpressionFile,'FileType','delimitedtext');
                optsCustomGene.VariableNames = settings.CustomGeneExpressionVariableNames;
                CustomGene = tdfread(settings.Custom_GeneExpressionFile,'\t',optsCustomGene);
                CustomGene = cell2struct(struct2cell(CustomGene),optsCustomGene.VariableNames);
                CustomGene = struct2table(CustomGene);
                CustomGene_AN.GeneExpr = table2struct(CustomGene);
            catch        
            end    
        end
    end
    if (strcmp(Organism,'Human'))
        optsUCSC = detectImportOptions(settings.Human_wgEncodeGencodeRefSeqFile);
        optsUCSC.VariableNames = {'Var1','Var2','Var3'};
        %EMBLRefSeqAlign_db = readtable(settings.Human_wgEncodeGencodeRefSeqFile,optsUCSC); 
        Gencode_db = readtable(settings.Human_GencodeRefSeqMetadataFile,optsUCSC); 
        optsUCSC2 = detectImportOptions(settings.Human_wgEncodeGencodeAttributesFile);
        optsUCSC2.VariableNames = settings.wgEncodeGencodeAttributesVariableNames;
        EMBLAttrAlign_db = readtable(settings.Human_wgEncodeGencodeAttributesFile,optsUCSC2); 
    elseif (strcmp(Organism,'Mouse'))
        optsUCSC = detectImportOptions(settings.Mouse_wgEncodeGencodeRefSeqFile);
        optsUCSC.VariableNames = {'Var1','Var2','Var3'};
        %EMBLRefSeqAlign_db = readtable(settings.Mouse_wgEncodeGencodeRefSeqFile,optsUCSC); 
        Gencode_db = readtable(settings.Mouse_GencodeRefSeqMetadataFile,optsUCSC); 
        optsUCSC2 = detectImportOptions(settings.Mouse_wgEncodeGencodeAttributesFile);
        optsUCSC2.VariableNames =settings.wgEncodeGencodeAttributesVariableNames;
        EMBLAttrAlign_db = readtable(settings.Mouse_wgEncodeGencodeCompFile,optsUCSC2); 
    elseif (strcmp(Organism,'Yeast'))
        optsUCSC = detectImportOptions(settings.Yeast_wgEncodeGencodeRefSeqFile);
        optsUCSC.VariableNames = {'Var1','Var2','Var3'};
        %EMBLRefSeqAlign_db = readtable(settings.Yeast_wgEncodeGencodeRefSeqFile,optsUCSC); 
        Gencode_db = readtable(settings.Yeast_GencodeRefSeqMetadataFile,optsUCSC); 
        optsUCSC2 = detectImportOptions(settings.Yeast_wgEncodeGencodeAttributesFile);
        optsUCSC2.VariableNames = settings.wgEncodeGencodeAttributesVariableNames;
        EMBLAttrAlign_db = readtable(settings.Yeast_wgEncodeGencodeAttributesFile,optsUCSC2);
    else  %custom organism 
        optsUCSC = detectImportOptions(settings.Custom_wgEncodeGencodeRefSeqFile);
        optsUCSC.VariableNames = {'Var1','Var2','Var3'};
        %EMBLRefSeqAlign_db = readtable(settings.Custom_wgEncodeGencodeRefSeqFile,optsUCSC); 
        Gencode_db = readtable(settings.Custom_GencodeRefSeqMetadataFile,optsUCSC); 
        optsUCSC2 = detectImportOptions(settings.Custom_wgEncodeGencodeAttributesFile);
        optsUCSC2.VariableNames = settings.wgEncodeGencodeAttributesVariableNames;
        EMBLAttrAlign_db = readtable(settings.Custom_wgEncodeGencodeAttributesFile,optsUCSC2); 
    end
%% Get gene expression levels in all tissues and cell types
get_exp_start = tic;
if (Expression_Null==1)  
   DNAcopynum = Null_DNAcopynum;
   ExpressionMatrix = zeros(length(uniNames),1);
   ExpressionMatrix(DNA_IDs) = DNAcopynum;
   ExpressionMatrix(NonDNA_IDs) = Null_RNAcopynum;
else                  
   if (strcmp(Organism,'Yeast'))
       k = 7;
   elseif (strcmp(Organism,'Mouse'))
       k = 6;
   elseif (strcmp(Organism,'Human'))
       if (strcmp(HumanOntology,'Normal'))
          if (HumanTissueOnly==0)
              if (HumanExpGeneOrTransc)%GeneID
                  k = 2;
              else    %TranscID
                  k = 1;
              end
          else
              k = 5;
          end
       end
       if (strcmp(HumanOntology,'Cancer'))
          if (HumanExpGeneOrTransc)%GeneID
              k = 3;
          else    %TranscID
              k = 4;
          end 
       end
   else
       if (settings.CustomExpGeneOrTransc)
           k = 8;
       else
           k = 9;
       end
   end  
   switch k 
        case 1
            ExprStruct = GTEx_Transc_AN.gtexTranscExpr;
            Expr_Names = {ExprStruct.name};
        case 2
            ExprStruct = GTEx_GeneV8_AN.gtexGeneV8;
            Expr_Names = {ExprStruct.geneId};
        case 3
            ExprStruct = TCGA_Gene_AN.tcgaGeneExpr;
            Expr_Names = {ExprStruct.name};
        case 4
            ExprStruct = TCGA_Transc_AN.tcgaTranscExpr;
            Expr_Names = {ExprStruct.name};
        case 5
            ExprStruct =  scTrack_Gene_AN(Human_SelTrack).tcgaTranscExpr;
            Expr_Names = {ExprStruct.name};    %is it name?
        case 6
            ExprStruct = TabulaMuris_AN.tabulamurisBarChart;
            Expr_Names = {ExprStruct.name}; %is it name?   
        case 7
            ExprStruct = TabulaYeast_AN;
            Expr_Names = {ExprStruct.name}; %is it name?     
        case 8
            ExprStruct = CustomGene_AN.GeneExpr ;
            Expr_Names = {ExprStruct.name}; %is it name?    
        case 9
            ExprStruct = CustomTransc_AN.TranscExpr;
            Expr_Names = {ExprStruct.name}; %is it name?    
   end
%    RiboGenes1 = find(strcmp(EMBLAttrAlign_db.transcriptType,'rRNA'));
%    RiboGenes2 = find(strcmp(EMBLAttrAlign_db.transcriptType,'rRNA_pseudogene'));
%    rRNA_Transcript_ID_EMBL = EMBLAttrAlign_db.transcriptId(RiboGenes1);
%    rRNA_pseudogene_Transcript_ID_EMBL = EMBLAttrAlign_db.transcriptId(RiboGenes2);
%    RiboGenes_Transcript_ID_EMBL = {rRNA_Transcript_ID_EMBL{:} rRNA_pseudogene_Transcript_ID_EMBL{:}};
%    RiboGenes_Transcript_ID_EMBL = extractBefore(RiboGenes_Transcript_ID_EMBL,'.');
%    Expr_Struct_RiboGeneIDs = cell(1,length(RiboGenes_Transcript_ID_EMBL));
%    for v = 1:length(RiboGenes_Transcript_ID_EMBL)
%       Expr_Struct_RiboGeneIDs{v} = find(strcmp(Expr_Names,RiboGenes_Transcript_ID_EMBL{v}));
%    end
%    Expr_Struct_RiboGeneIDs2 = horzcat(Expr_Struct_RiboGeneIDs{:});
%    OtherGenes = setdiff(1:length(Expr_Names),Expr_Struct_RiboGeneIDs2);
%    Original_rRNA_Expr = cell2mat(arrayfun(@(x) ExprStruct(x).expScores.',Expr_Struct_RiboGeneIDs2,'UniformOutput',false)).';
%    OtherExpr = cell2mat(arrayfun(@(x) ExprStruct(x).expScores.',OtherGenes,'UniformOutput',false)).';
%    rRNA_Tissue = sum(Original_rRNA_Expr,1);
%    Other_Tissue = sum(OtherExpr,1);
   
   Expr_Names = convertCharsToStrings(Expr_Names);
   GencodeList_Transcript = extractBefore(Gencode_db.Var1,'.');%Transcript_ID
   GencodeList_Transcript = convertCharsToStrings(GencodeList_Transcript);%Transcript_ID          
   GencodeList_RefSeq = extractBefore(Gencode_db.Var2,'.');
   GencodeList_RefSeq = convertCharsToStrings(GencodeList_RefSeq);
   paired_EMBLgeneId = cellfun(@(x) extractBefore(x,'.'),EMBLAttrAlign_db.geneId,'UniformOutput',false);
   paired_EMBLtranscId = cellfun(@(x) extractBefore(x,'.'),EMBLAttrAlign_db.transcriptId,'UniformOutput',false);
   paired_EMBLgeneId = convertCharsToStrings(paired_EMBLgeneId);
   paired_EMBLtranscId = convertCharsToStrings(paired_EMBLtranscId);
   [C_GList_T,GList_T_ia,GList_T_ic] = unique(GencodeList_Transcript);
   [C_GList_R,GList_R_ia,GList_R_ic] = unique(GencodeList_RefSeq);
   [C_EList_T,EList_T_ia,EList_T_ic] = unique(paired_EMBLtranscId);
   [C_EList_G,EList_G_ia,EList_G_ic] = unique(paired_EMBLgeneId);
   indx_GList_T = @(x) (find(find(strcmp(x,C_GList_R))==GList_R_ic));
   indx_EList_G = @(x) (find(find(strcmp(x,C_EList_T))==EList_T_ic));
   DNA_IDs_1 = find(contains(uniNames,'NC_'));%IDs
   DNA_IDs_2 = find(contains(uniNames,'NT_'));%IDs
   DNA_IDs =union(DNA_IDs_1,DNA_IDs_2);
   RNA_IDs_1 = find(~contains(uniNames,'NC_'));%IDs
   RNA_IDs_2 = find(~contains(uniNames,'NT_'));%IDs
   RNA_IDs = intersect(RNA_IDs_1,RNA_IDs_2);
   y = @(u) (GencodeList_Transcript{indx_GList_T(uniNames{u})});
   z = @(u) (paired_EMBLgeneId{indx_EList_G(u)});
   try
        n_cores = str2num(getenv('SLURM_JOB_CPUS_PER_NODE')); 
        p = parpool(n_cores);
        spmd
           warning('off','all')
        end
   catch
   end
   parfor u = 1:length(RNA_IDs)
      v = RNA_IDs(u);
      try
         Corresponding_transcID{u} = GencodeList_Transcript{indx_GList_T(uniNames{v})};
      catch
      end 
      try
         Corresponding_geneID{u} = paired_EMBLgeneId{indx_EList_G(GencodeList_Transcript{indx_GList_T(uniNames{v})})};
      catch
      end
   end  
   nTissue = length(ExprStruct(1).expScores);
   ExpressionMatrix = zeros(length(uniNames),nTissue+43);
   ExpressionMatrix(DNA_IDs,1:nTissue) = 2;
   if (ismember(k,[2 3]))
       v = find(~cellfun(@isempty,Corresponding_geneID));%RNA_IDs
       parfor u=1:length(v)
          if (sum(strcmp(Expr_Names,Corresponding_geneID{v(u)}))>0)
             corrRow(u) = find(strcmp(Expr_Names,Corresponding_geneID{v(u)}));
          end
       end
   else
       v = find(~cellfun(@isempty,Corresponding_transcID));
       parfor u=1:length(v)
          if (sum(strcmp(Expr_Names,Corresponding_transcID{v(u)}))>0)
             corrRow(u) = find(strcmp(Expr_Names,Corresponding_transcID{v(u)}));
          end
       end
   end
   RNA_IDs2 = RNA_IDs(v);
   corrIdx = find(corrRow>0);
   ExpressionMatrix(RNA_IDs2(corrIdx),1:nTissue) = cell2mat(arrayfun(@(x) ExprStruct(x).expScores.',corrRow(corrIdx),'UniformOutput',false)).';           
end   


if (strcmp(Organism,'Human'))
%NonDNA_IDs   , RNA_IDs 
sRefSeq_List = uniNames(NonDNA_IDs);



HPA_Expression = readtable(settings.HumanHPAcellLine_TranscExpressionFile,'FileType','text','Delimiter','\t');
HPA_Expression_CellLineInfo = HPA_Expression.CellLine;

% Get cell-lines to keep gene expression information for.
KeepInfo = cell(1,43);
KeepInfo{1} = {'AF22','F','Brain','NC','adherent','stomatic stem cell','neuralepithelial','isenetbiobanking.com/af22','doubling time NA'}; 
KeepInfo{2} = {'ASC52telo','F','Adipose','NC','adherent','mesenchymal stem cell','Adipose Tissue','ATCC:SCRC-4000','doubling time 36-48hr'};
KeepInfo{3} = {'BJ [Human fibroblast]','M','Skin','NC','adherent','fibroblast','skin','ATCC:CRL-2522','doubling time NA'};
KeepInfo{4} = {'BJ1-hTERT','M','Skin','NC','adherent','fibroblast','skin','ATCC:CRL-4001','doubling time 36hrs'};
KeepInfo{5} = {'fHDF/TERT166','M','Skin','NC','adherent','fibroblast','skin','Evercyte:fHDF/TERT166','doubling time 48hrs'};
KeepInfo{6} = {'HaCaT','M','Skin','NC','adherent','keratinocyte','skin','CLS:300493','doubling time 28hrs'};
KeepInfo{7} = {'HBEC3-KT','F','Lung','NC','adherent','epithelial cell','lung;bronchus','ATCC:CRL-4051','doubling time 36-48hrs'};
KeepInfo{8} = {'HEK293','F','Kidney','NC','adherent','epithelial cell','kidney','ATCC:CRL-1573','doubling time 20-30hrs'};
KeepInfo{9} = {'hTCEpi','M','Eye','NC','adherent','epithelial cell','corneal','Evercyte:','doubling time 24-32hrs'};
KeepInfo{10} = {'hTEC/SVTERT24-B','F','Lymphoid','NC','adherent','epithelial cell','thymus','Evercyte:hTEC/SVTERT24-B','doubling time 30hrs'};
KeepInfo{11} = {'hTERT-HME1','F','Breast','NC','adherent','epithelial cell','breast;mammary gland','ATCC:CRL-4010','doubling time 25 or 48-64 hrs'};
KeepInfo{12} = {'hTERT-RPE1','F','Eye','NC','adherent','epithelial cell','eye;retina;retina pigment epithelium','ATCC:CRL-4000','doubling time NA'};
KeepInfo{13} = {'HUVEC/TERT2','F','Endothelial','NC','adherent','endothelial cell','umbilical cord;vascular endothelium','ATCC:CRL-4053','doubling time 28hrs or 2.5 days'};
KeepInfo{14} = {'LHCN-M2','M','Muscle','NC','adherent','satellite cell','pectoralis major muscle','Evercyte:LHCN-M2','doubling time 35-40 hrs'};
KeepInfo{15} = {'MCF-10A','F','Breast','NC','adherent','epithelial cell','breast;mammary gland','ATCC:CRL-10317','doubling time NA'};
KeepInfo{16} = {'PODO/SVTERT152','M','Kidney','NC','adherent','podocytes;visceral epithelial cells','male kidney','Evercyte:podo-svtert152','doubling time NA'};
KeepInfo{17} = {'PODO/TERT256','F','Kidney','NC','adherent','podocytes;visceral epithelial cells','female kidney','Evercyte:podo-tert256','doubling time NA'};
KeepInfo{18} = {'RPTEC/TERT1','M','Kidney','NC','adherent','renal proximal tubular epithelial cell','kidney cortex','Evercyte:','doubling time 72-96 hrs'};
KeepInfo{19} = {'TIME','M','Skin','NC','adherent','microvascular endothelial cell','skin','ATCC:CRL-4025','doubling time NA'};
KeepInfo{20} = {'A-431','F','Skin','C','adherent','epithelial cell','skin;epidermis','ATCC:CRL-1555','doubling time 80-100hrs'};
KeepInfo{21} = {'A-549','M','Lung','C','adherent','epithelial cell','lung','ATCC:CCL-185','doubling time 18-40hrs'};
KeepInfo{22} = {'CACO-2','M','GI Tract','C','adherent','epithelial cell','Colon;Large Intestines','ATCC:HTB-37','doubling time 32-80hrs'};
KeepInfo{23} = {'EFO-21','F','Ovarian','C','adherent','cystadenocarcinoma','','DSMZ:ACC-235','doubling time 40-60hrs'};
KeepInfo{24} = {'GAMG','F','Brain','C','adherent','glioma','','DSMZ:ACC-242','doubling time 30-50hrs'};
KeepInfo{25} = {'HAP1','M','Myeloid','C','adherent','neoplasm','','HorizonDiscovery:C631','doubling time NA'};
KeepInfo{26} = {'HDLM-2','M','Lymphoid','C','suspension','B Lymphocyte','lung;pleural effusion','DSMZ:ACC-17','doubling time 70-75hrs'};
KeepInfo{27} = {'HEL','M','Myeloid','C','suspension','erythroleukemia cell','peripheral blood','DSMZ:ACC-11','doubling time 24hrs'};
KeepInfo{28} = {'HeLa','F','Cervix','C','adherent','epithelial cell','uterus;cervix','ATCC:CCL-2','doubling time 48hrs'};
KeepInfo{29} = {'JURKAT','M','Lymphoid','C','adherent','T cell','peripheral blood','many','doubling time 25-35hrs'};
KeepInfo{30} = {'K-562','F','Myeloid','C','suspension','lymphoblast','Bone;Marrow','ATCC:CCL-243','doubling time 18-47hrs'};
KeepInfo{31} = {'MCF-7','F','Breast','C','adherent','epithelial cell','breast;mammary gland','ATCC:HTB-22','doubling time 25-80hrs'};
KeepInfo{32} = {'NB4','F','Bone Marrow','C','suspension','B Lymphocyte','bone marrow','four,CLS:300299','doubling time 30-45hrs'};
KeepInfo{33} = {'OE19','M','Proximal Digestive Track','C','adherent','epithelial cell','oesophagus','two,DSMZ:ACC-700','doubling time 50-60hrs'};
KeepInfo{34} = {'PC-3','M','Bone Marrow','C','adherent','epithelial cell','prostate','many,ATCC:CRL-1435','doubling time 25-50hrs'};
KeepInfo{35} = {'REH','F','Blood','C','suspension','lymphoblast','blood','many,ATCC:CRL-8286','doubling time 30-70hrs'};
KeepInfo{36} = {'Rh30','M','Bone Marrow','C','adherent','fibroblast','skeletal muscle','many,ATCC:CRL-2061','doubling time 35-37hrs'};
KeepInfo{37} = {'RT-4','M','Urinary Bladder','C','adherent','epithelial cell','urinary bladder','many,ATCC:HTB-2','doubling time 37-80hrs'};
KeepInfo{38} = {'SiHa','F','Cervix','C','adherent','epithelial cell','uterus;cervix','many,ATCC:HTB-35','doubling time 17-21hrs'};
KeepInfo{39} = {'SK-MEL-30','M','Skin','C','adherent','epithelial cell','skin','DSMZ:ACC-151','doubling time 30hrs'};
KeepInfo{40} = {'SuSa','M','Testis','C','adherent','germ cell','testis','DSMZ:ACC-747','doubling time 27-72hrs'};
KeepInfo{41} = {'THP-1','M','Myeloid','C','suspension','monocyte','periheral blood','many,ATCC:TIB-202','doubling time 26-70hrs'};
KeepInfo{42} = {'U-251MG','M','Brain','C','adherent','epithelial cell','brain','many,CLS:000375','doubling time 23-24hrs'};
KeepInfo{43} = {'U2OS','F','Bone Marrow','C','adherent','epithelial cell','bone','many,ATCC:HTB-96','doubling time 25-36hrs'};
KeepInfo_Region = cellfun(@(x) x{3},KeepInfo,'Un',0);
[~,~,RegionNum] = unique(KeepInfo_Region);
[~,newOrder] = sort(RegionNum,'ascend');
KeepOrder_Name = arrayfun(@(x) KeepInfo{newOrder(x)}{1},1:length(newOrder),'Un',0);
N_CellLines = length(KeepOrder_Name);
newHPA_Table = cell(1,N_CellLines);
    for i = 1:N_CellLines 
        newHPA_Table{i} = HPA_Expression(strcmp(HPA_Expression_CellLineInfo,KeepOrder_Name{i}),[1 2 4 5 6]); 
    end
HPAMatrix_TPM = cell2mat(arrayfun(@(v) table2array(newHPA_Table{v}(:,3)),1:N_CellLines,'Un',0));
HPAList_GeneID = newHPA_Table{1}.Gene;
geneID_ENST = extractBefore(EMBLAttrAlign_db.geneId,'.');
transID_ENST = extractBefore(EMBLAttrAlign_db.transcriptId,'.');
transID_Gencode = extractBefore(Gencode_db.Var1,'.');
RefSeqID_Gencode = extractBefore(Gencode_db.Var2,'.');
HPA_ENSTcorr = arrayfun(@(x) find(strcmp(geneID_ENST,HPAList_GeneID{x})),1:size(HPAMatrix_TPM,1),'Un',0);
HPA_TRANS = arrayfun(@(x) unique(transID_ENST(HPA_ENSTcorr{x})),1:size(HPAMatrix_TPM,1),'Un',0);
HPA_TRANScorr = arrayfun(@(x) cell2mat(arrayfun(@(y) find(strcmp(transID_Gencode,HPA_TRANS{x}{y}))',1:length(HPA_TRANS{x}),'Un',0)),1:size(HPAMatrix_TPM,1),'Un',0);
HPA_RefSeq = arrayfun(@(x) unique(RefSeqID_Gencode(HPA_TRANScorr{x})),1:size(HPAMatrix_TPM,1),'Un',0);
HPA_sTabelCorr = arrayfun(@(x) cell2mat(arrayfun(@(y) find(strcmp(sRefSeq_List,HPA_RefSeq{x}{y}))',1:length(HPA_RefSeq{x}),'Un',0)),1:size(HPAMatrix_TPM,1),'Un',0);
CellLine_TPM = zeros(length(sRefSeq_List),N_CellLines);
    for v = 1:length(HPA_sTabelCorr)
        if (~isempty(HPA_sTabelCorr{v}))
            CellLine_TPM(HPA_sTabelCorr{v},:) =  repmat(HPAMatrix_TPM(v,:),[length(HPA_sTabelCorr{v}) 1]);
        end
    end  
ExpressionMatrix(NonDNA_IDs,nTissue+1:nTissue+N_CellLines) = CellLine_TPM;
end


get_expression_time = toc(get_exp_start);
end