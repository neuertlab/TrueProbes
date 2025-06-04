function [AN, AN_Chr, AN_Start, AN_End] = GetAccensionNumber(Organism,targetSymbol,nucType)
%nType (DNA 1 or RNA 2), consider additionally returning start and stop
%location for evaluating nascent better, also not vectorized for getting
%multiple sequence accensionNumbers
%Nuc => Nucleotide, Pro=>Protein, %Accession Number Gene(DNA), RNA, as well%protein
 %S1 (NCBI GeneID, ENSEMBL GeneID, HGNC ID, Symbol, RefSeq Nuc/Pro, ENSEMBL RNA/Protein ID, Chr, Start, End)
 %S2 (NCBI GeneID, Symbol, RSG Refseq Nuc (NG), LRG ID, RNA RefSeq Nuc/Pro
 %S3 (NCBI GeneID, Symbol, RSG RefSeq Nuc (NG),
 %S4 (NCBI GeneID, Symbol,
 %S5 (ENSEMBL RNA ID, NCBI RefSeq Nuc/Pro)
 %S6 (ENSEMBL RNA ID, Pubmed ID)
 %S7 (ENSEMBL RNA ID, Entrez ID)
 %comment throughout
nType = 0;AN_Start = []; AN_End = []; AN_Chr = [];
if (strcmp(nucType,'DNA'))
    nType = 1;
elseif (strcmp(nucType,'RNA'))
    nType = 2;
end
if (strcmp(Organism,'Human'))
    optsHuman = detectImportOptions('Metadata/Human/gencode.v39.metadata.RefSeq.txt');
    optsHuman.DataLines=[1,Inf];
    optsHuman.VariableNamesLine = 0;
    optsHuman.VariableNames = {'Var1','Var2','Var3'};
    S1 = tdfread('Metadata/Human/MANE.GRCh38.v0.95.summary.txt');
    S2 = tdfread('Metadata/Human/Human_LRG_RefSeqGene.txt');
    S3 = tdfread('Metadata/Human/Human_gene_RefSeqGene.txt');
    S5 = readtable('Metadata/Human/gencode.v39.metadata.RefSeq.txt',optsHuman);
    EMBLID = S5.Var1; %ENST
    AN_RNA = S5.Var2; %85199
    ChrNum = cellstr(S1.GRCh38_chr);
    ChrStart = S1.chr_start;
    ChrEnd = S1.chr_end;
    MANE_Symbol = cellstr(S1.symbol);
    MANE_NCBIGeneID = cellstr(S1.x0x23NCBI_GeneID(:,8:end));
    MANE_ENSEMBnucLID = cellstr(S1.Ensembl_nuc);
    MANE_RefSeqAN_RNA = cellstr(S1.RefSeq_nuc); %18640
    LocusRefGenome_Symbol = cellstr(S2.Symbol);
    LocusRefGenome_GeneID = S2.GeneID;
    LocusRefGenome_LRGAN_DNA = cellstr(S2.RSG); %unique 6845
    LocusRefGenome_LRGAN_RNA = cellstr(S2.RNA); %31184
    RefSeqGene_Symbol = cellstr(S3.Symbol);
    RefSeqGene_GeneID = S3.GeneID;
    RefSeqGene_RSGAN_DNA = cellstr(S3.RSG); %6876, %6842 unique
    FindMANESymbol = cell2mat(cellfun(@(x) strcmp(x,targetSymbol),MANE_Symbol,'UniformOutput',false));
    FindLRGSymbol = cell2mat(cellfun(@(x) strcmp(x,targetSymbol),LocusRefGenome_Symbol,'UniformOutput',false));
    FindRSGSymbol = cell2mat(cellfun(@(x) strcmp(x,targetSymbol),RefSeqGene_Symbol,'UniformOutput',false));
    ID_MANE = find(FindMANESymbol);
    ID_LRG = find(FindLRGSymbol);
    ID_RSG = find(FindRSGSymbol);
    S5_AN_RNA = [];MANE_AN_RNA=[];LRG_AN_RNA=[];FindEMBL=[];
    LRG_AN_DNA=[];RSG_AN_DNA=[];AN=[];
    if (~isempty(ID_MANE))
        EMBL_MANE = MANE_ENSEMBnucLID{ID_MANE};
        if (iscell(EMBL_MANE))
            for i=1:length(EMBL_MANE)
            FinEMBL{i}  = cell2mat(cellfun(@(x) strcmp(x,EMBL_MANE{i}),EMBLID,'UniformOutput',false));
            FindEMBL = FindEMBL + FinEMBL{i}.';
            end        
        else
            FindEMBL = cell2mat(cellfun(@(x) strcmp(x,EMBL_MANE),EMBLID,'UniformOutput',false));
        end
        S5_ID = find(FindEMBL);
        S5_AN_RNA = {AN_RNA{S5_ID}};
        MANE_AN_RNA = {MANE_RefSeqAN_RNA{ID_MANE}};
        AN_Chr = ChrNum{ID_MANE};
        AN_Start = ChrStart(ID_MANE);
        AN_End = ChrEnd(ID_MANE);
    end 
    if (~isempty(ID_LRG))
        LRG_AN_DNA = {LocusRefGenome_LRGAN_DNA{ID_LRG}};
        LRG_AN_RNA = {LocusRefGenome_LRGAN_RNA{ID_LRG}};
    end 
    if (~isempty(ID_RSG))
        RSG_AN_DNA = {RefSeqGene_RSGAN_DNA{ID_RSG}};
    end 
    try
    if (nType==1)  
        AN = unique({LRG_AN_DNA{:} RSG_AN_DNA{:}});
    elseif (nType==2)
        AN = unique({S5_AN_RNA{:} MANE_AN_RNA{:} LRG_AN_RNA{:}});
    end
    catch
    end
  
elseif strcmp(Organism,'Mouse') 
    optsMouse = detectImportOptions('Metadata/Mouse/gencode.vM28.metadata.RefSeq.txt');
    optsMouse.DataLines=[1,Inf];
    optsMouse.VariableNamesLine = 0;
    optsMouse.VariableNames = {'Var1','Var2','Var3'};
    S12 = tdfread('Metadata/Mouse/Mus_musculus.gene_info.txt');
    S13 = readtable('Metadata/Mouse/gencode.vM28.metadata.RefSeq.txt',optsMouse);
    S15 = readtable('Metadata/Mouse/gencode.vM28.metadata.EntrezGene.txt');
    SM = tdfread('Metadata/Mouse/Mouse_ncbi_data.tsv');
    EMBLID1 = S13.Var1;
    AN_RNA = S13.Var2;
    EMBLID3 = S15.Var1;
    EntrezID = S15.Var2;
    GeneEntrezID = S12.GeneID;%Entrez
    Symbol = cellstr(S12.Symbol);%upper case
    M_Symbol = cellstr(SM.Symbol);%not upper case
    M_GeneID = SM.NCBI_GeneID;
    AN_RNA = cellstr(SM.Transcript_Accession);
    AN_DNA = cellstr(SM.Transcript_Genomic_Accession);
    M_Chr = cellstr(SM.Chromosomes);
    M_Start = SM.Transcript_Genomic_Start;
    M_End = SM.Transcript_Genomic_Stop;
    FindSymbol = cell2mat(cellfun(@(x) strcmp(upper(x),upper(targetSymbol)),Symbol,'UniformOutput',false));
    IDs = find(FindSymbol);
    FindMSymbol = cell2mat(cellfun(@(x) strcmp(upper(x),upper(targetSymbol)),M_Symbol,'UniformOutput',false));
    ID_M = find(FindMSymbol);
    EntrezIDs = GeneEntrezID(IDs);
    EntrezLoc = find(EntrezID==EntrezIDs);
    EMBLID3s = EMBLID3(EntrezLoc);
    FoundEMBL = zeros(1,length(EMBLID1));
    S13_AN_RNA = [];M_AN_RNA=[];AN = [];
    if(iscell(EMBLID3s))
        for i=1:length(EMBLID3s)
        FinEMBL{i} = cell2mat(cellfun(@(x) strcmp(x,EMBLID3s{i}),EMBLID1,'UniformOutput',false));
        FoundEMBL = FoundEMBL + FinEMBL{i}.';
        end
    else
        FoundEMBL = cell2mat(cellfun(@(x) strcmp(x,EMBLID3s),EMBLID1,'UniformOutput',false)); 
    end
    S13_ID = find(FoundEMBL);
    S13_AN_RNA = {AN_RNA{S13_ID}};
    if (~isempty(ID_M))
        M_AN_RNA = {AN_RNA{ID_M}};
        M_AN_DNA = {AN_DNA{ID_M}};
    end 
    try
    if (nType==1)
        StartLocs = M_Start(ID_M);
        EndLocs = M_End(ID_M);
        ChrNum = {M_Chr{ID_M}};
        [AN,ia,ic] = unique({M_AN_DNA{:}});
        AN_Chr = ChrNum{ia};
        AN_Start = StartLocs(ia);
        AN_End = EndLocs(ia);
    elseif (nType==2)
        AN = unique({S13_AN_RNA{:} M_AN_RNA{:}});
    end
    catch
    end
elseif strcmp(Organism,'Yeast')
    S22 = readtable('Metadata/Yeast/S288C/GCF_000146045.2_R64_feature_table.txt'); 
    %Paired_Feature = S22.x_Feature;
    Paired_Class = S22.class;
    Paired_DNA_ANs = S22.genomic_accession;
    Paired_ProductAccession = S22.product_accession;%includes rna & protein
    Paired_Symbol = S22.symbol;
    Paired_Start = S22.start;
    Paired_End = S22.xEnd;
    Paired_Chr = S22.chromosome;
    FindSymbol = cell2mat(cellfun(@(x) strcmp(x,targetSymbol),Paired_Symbol,'UniformOutput',false));
    AN = [];
    if (max(FindSymbol)==1)
        IDs = find(FindSymbol);
        try
        if (nType==1)
            StartLocs = Paired_Start(IDs);
            EndLocs = Paired_End(IDs);
            ChrNum = {Paired_Chr{IDs}};
            [AN,ia,ic] = unique({Paired_DNA_ANs{IDs}});
            AN_Chr = ChrNum{ia};
            AN_Start = StartLocs(ia);
            AN_End = EndLocs(ia);
        elseif (nType==2)
            FindRNA = cell2mat(cellfun(@(x) strcmp(x,''),Paired_Class,'UniformOutput',false));
            RNAID = find(FindRNA);
            ID2s = intersect(IDs,RNAID);
            AN = unique({Paired_ProductAccession{ID2s}});
        end
        catch
        end
    end 
end

