function globalExpressionMatrix = A_ProbeDesignerCellLineExpression(settings)       
try
    load('DatabaseData/globalExpressionMatrix.mat','globalExpressionMatrix') 
    calc = 0;
catch
    calc = 1;
end
if (calc)    
    KeepInfo = cell(1,44);
    optsUCSC = detectImportOptions(settings.Human_wgEncodeGencodeRefSeqFile);
    optsUCSC.VariableNames = {'Var1','Var2','Var3'};
    Gencode_db = readtable(settings.Human_GencodeRefSeqMetadataFile,optsUCSC); 
    optsUCSC2 = detectImportOptions(settings.Human_wgEncodeGencodeAttributesFile);
    optsUCSC2.VariableNames = settings.wgEncodeGencodeAttributesVariableNames;
    EMBLAttrAlign_db = readtable(settings.Human_wgEncodeGencodeAttributesFile,optsUCSC2);
    geneID_ENST = extractBefore(EMBLAttrAlign_db.geneId,'.');
    transID_ENST = extractBefore(EMBLAttrAlign_db.transcriptId,'.');
    transID_Gencode = extractBefore(Gencode_db.Var1,'.');
    RefSeqID_Gencode = extractBefore(Gencode_db.Var2,'.');
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
    KeepInfo{26} = {'HBF/TERT88','F','Brain','C','adherent','','','Evercyte:','doubling time NA'};
    KeepInfo{27} = {'HDLM-2','M','Lymphoid','C','suspension','B Lymphocyte','lung;pleural effusion','DSMZ:ACC-17','doubling time 70-75hrs'};
    KeepInfo{28} = {'HEL','M','Myeloid','C','suspension','erythroleukemia cell','peripheral blood','DSMZ:ACC-11','doubling time 24hrs'};
    KeepInfo{29} = {'HeLa','F','Cervix','C','adherent','epithelial cell','uterus;cervix','ATCC:CCL-2','doubling time 48hrs'};
    KeepInfo{30} = {'JURKAT','M','Lymphoid','C','adherent','T cell','peripheral blood','many','doubling time 25-35hrs'};
    KeepInfo{31} = {'K-562','F','Myeloid','C','suspension','lymphoblast','Bone;Marrow','ATCC:CCL-243','doubling time 18-47hrs'};
    KeepInfo{32} = {'MCF-7','F','Breast','C','adherent','epithelial cell','breast;mammary gland','ATCC:HTB-22','doubling time 25-80hrs'};
    KeepInfo{33} = {'NB4','F','Bone Marrow','C','suspension','B Lymphocyte','bone marrow','four,CLS:300299','doubling time 30-45hrs'};
    KeepInfo{34} = {'OE19','M','Proximal Digestive Track','C','adherent','epithelial cell','oesophagus','two,DSMZ:ACC-700','doubling time 50-60hrs'};
    KeepInfo{35} = {'PC-3','M','Bone Marrow','C','adherent','epithelial cell','prostate','many,ATCC:CRL-1435','doubling time 25-50hrs'};
    KeepInfo{36} = {'REH','F','Blood','C','suspension','lymphoblast','blood','many,ATCC:CRL-8286','doubling time 30-70hrs'};
    KeepInfo{37} = {'Rh30','M','Bone Marrow','C','adherent','fibroblast','skeletal muscle','many,ATCC:CRL-2061','doubling time 35-37hrs'};
    KeepInfo{38} = {'RT-4','M','Urinary Bladder','C','adherent','epithelial cell','urinary bladder','many,ATCC:HTB-2','doubling time 37-80hrs'};
    KeepInfo{39} = {'SiHa','F','Cervix','C','adherent','epithelial cell','uterus;cervix','many,ATCC:HTB-35','doubling time 17-21hrs'};
    KeepInfo{40} = {'SK-MEL-30','M','Skin','C','adherent','epithelial cell','skin','DSMZ:ACC-151','doubling time 30hrs'};
    KeepInfo{41} = {'SuSa','M','Testis','C','adherent','germ cell','testis','DSMZ:ACC-747','doubling time 27-72hrs'};
    KeepInfo{42} = {'THP-1','M','Myeloid','C','suspension','monocyte','periheral blood','many,ATCC:TIB-202','doubling time 26-70hrs'};
    KeepInfo{43} = {'U-251MG','M','Brain','C','adherent','epithelial cell','brain','many,CLS:000375','doubling time 23-24hrs'};
    KeepInfo{44} = {'U2OS','F','Bone Marrow','C','adherent','epithelial cell','bone','many,ATCC:HTB-96','doubling time 25-36hrs'};
    KeepInfo_Name = cellfun(@(x) x{1},KeepInfo,'Un',0);
    KeepInfo_Gender = cellfun(@(x) x{2},KeepInfo,'Un',0);
    KeepInfo_Region = cellfun(@(x) x{3},KeepInfo,'Un',0);
    KeepInfo_CNC = cellfun(@(x) x{4},KeepInfo,'Un',0);
    KeepInfo_Growth = cellfun(@(x) x{5},KeepInfo,'Un',0);
    KeepInfo_CT = cellfun(@(x) x{6},KeepInfo,'Un',0);
    [~,~,RegionNum] = unique(KeepInfo_Region);
    [~,newOrder] = sort(RegionNum,'ascend');
    KeepOrder_Name = arrayfun(@(x) KeepInfo{newOrder(x)}{1},1:length(newOrder),'Un',0);
    reOrderInfo = cell2mat(cellfun(@(x) find(strcmp(KeepInfo_Name,x)),KeepOrder_Name,'Un',0));
    KeepInfo_Gender1 = KeepInfo_Gender(reOrderInfo);
    KeepInfo_Region1 = KeepInfo_Region(reOrderInfo);
    KeepInfo_CNC1 = KeepInfo_CNC(reOrderInfo);
    KeepInfo_Growth1 = KeepInfo_Growth(reOrderInfo);
    KeepInfo_CT1 = KeepInfo_CT(reOrderInfo);
    Info_Table = table(KeepOrder_Name',KeepInfo_Gender1',KeepInfo_Region1',KeepInfo_CNC1',KeepInfo_Growth1',KeepInfo_CT1');
    HPA_Schema = readtable("DatabaseData/rna_celline_description.tsv",'FileType','text','Delimiter','\t');
    CellLineNames = HPA_Schema.CellLine;
    CellLineType = HPA_Schema.Disease;
    CellLinesToInclude = find(ismember(CellLineType,{'Non-cancerous','Uncategorized'}));
    CellLinesToExclude1 = find(ismember(CellLineNames,{'BEWO','DM-3','HBF/TERT88','HBL-100','HD-MY-Z','HLF-a','NTERA-2','RS-5','T1-73','TE 125.T','TE 159.T','TO 175.T','HEK TE'...
        ,'TIG-3 TD','SALE','PrEC LH','OELE','NHAHTDD','HMEL','HHSteC','HSkMC',...
       'ASC2telo differentiated','BJ hTERT+ SV40 Large T+','BJ hTERT+ SV40 Large T+ RasG12V'}));
    CellLinesToExclude2 = find(contains(CellLineNames,'Hs'));
    CellLinesToExclude = union(CellLinesToExclude1,CellLinesToExclude2);
    CellLinesToKeep = setdiff(CellLinesToInclude,CellLinesToExclude);
    AddBackIn = find(ismember(CellLineNames,{...
        'U-251MG','GAMG',...
        'OE19','CACO-2',...
        'RT-4','SuSa','PC-3',...
        'SK-MEL-30','A-431',...
        'A-549','HeLa','MCF-7',...
        'EFO-21','SiHa','Rh30',...
        'U2OS','HBF/TERT88',...
        'REH','HDLM-2','JURKAT',...
        'HAP1','HEL','NB4',...
        'THP-1','K-562'}));   
    HPA_Expression = readtable("DatabaseData/rna_celline.tsv",'FileType','text','Delimiter','\t');    
    %PotentialSchema = HPA_Schema([CellLinesToKeep;AddBackIn],:);
    HPA_Expression_CellLineInfo = HPA_Expression.CellLine;
    N_CellLines = length(KeepOrder_Name);
    newHPA_Table = cell(1,N_CellLines);
    for i = 1:N_CellLines 
       newHPA_Table{i} = HPA_Expression(strcmp(HPA_Expression_CellLineInfo,KeepOrder_Name{i}),[1 2 4 5 6]); 
    end
    HPAMatrix_TPM = cell2mat(arrayfun(@(v) table2array(newHPA_Table{v}(:,3)),1:N_CellLines,'Un',0));
    HPAList_GeneID = newHPA_Table{1}.Gene;
    load('DatabaseData/HumanGeneProbeDesignSelectionCriteria.mat','sTabel') 
    RefSeq_List = sTabel.RefSeq_ID;
    sRefSeq_List = extractBefore(RefSeq_List,'.');
    geneName_List = sTabel.geneName;
    transcriptType_List = sTabel.transcriptType;
    Num_Isoforms_List = sTabel.Num_Isoforms;
    OtherIso_RefSeqIDs_InSeqDB = sTabel.OtherIso_RefSeqIDs_InSeqDB;
    OtherIso_RefSeqIDs_InSeqDB_Len = cellfun(@length,OtherIso_RefSeqIDs_InSeqDB);
    Num_Isoform_InSeqDB = Num_Isoforms_List - OtherIso_RefSeqIDs_InSeqDB_Len-1;
    HPA_ENSTcorr = arrayfun(@(x) find(strcmp(geneID_ENST,HPAList_GeneID{x})),1:size(HPAMatrix_TPM,1),'Un',0);
    HPA_TRANS = arrayfun(@(x) unique(transID_ENST(HPA_ENSTcorr{x})),1:size(HPAMatrix_TPM,1),'Un',0);
    HPA_TRANScorr = arrayfun(@(x) cell2mat(arrayfun(@(y) find(strcmp(transID_Gencode,HPA_TRANS{x}{y}))',1:length(HPA_TRANS{x}),'Un',0)),1:size(HPAMatrix_TPM,1),'Un',0);
    HPA_RefSeq = arrayfun(@(x) unique(RefSeqID_Gencode(HPA_TRANScorr{x})),1:size(HPAMatrix_TPM,1),'Un',0);
    HPA_sTabelCorr = arrayfun(@(x) cell2mat(arrayfun(@(y) find(strcmp(sRefSeq_List,HPA_RefSeq{x}{y}))',1:length(HPA_RefSeq{x}),'Un',0)),1:size(HPAMatrix_TPM,1),'Un',0);
    sTabel_CellLine_TPM = zeros(length(RefSeq_List),N_CellLines);
    for v = 1:length(HPA_sTabelCorr)
       if (~isempty(HPA_sTabelCorr{v}))
            sTabel_CellLine_TPM(HPA_sTabelCorr{v},:) =  repmat(HPAMatrix_TPM(v,:),[length(HPA_sTabelCorr{v}) 1]);
       end
    end
    CellLines_Keep = KeepOrder_Name;
    CellLines_FinFilt = find(~ismember(CellLines_Keep,{'HBF/TERT88'}));
    CellLines_Keep_FinFilt = CellLines_Keep(CellLines_FinFilt);
    CellL_Table = Info_Table(cell2mat(cellfun(@(x) find(ismember(KeepOrder_Name,x)),CellLines_Keep_FinFilt,'Un',0)),:);
    CellTableColumns = {'Cell-Line','Gender','Tissue','CancerStatus','Growth','Cell-Type'};
    CellLcolOrder = [3 6 2 4 5];
    kCellLines_Keep = CellLines_Keep_FinFilt;
    CellLineSortOrder = [3 6];
    CellL_Table = Info_Table(cell2mat(cellfun(@(x) find(ismember(KeepOrder_Name,x)),kCellLines_Keep,'Un',0)),:);
    CellL_TableSort = sortrows(CellL_Table,CellLineSortOrder);
    CellL_TableSort.Var4(strcmp(CellL_TableSort.Var4,'C')) = {'Cancer'};
    CellL_TableSort.Var4(strcmp(CellL_TableSort.Var4,'NC')) = {'Non-Cancer'};
    CellL_TableSort.Var2(strcmp(CellL_TableSort.Var2,'F')) = {'Female'};
    CellL_TableSort.Var2(strcmp(CellL_TableSort.Var2,'M')) = {'Male'};
    reOrderCellLines = cell2mat(cellfun(@(x) find(ismember(kCellLines_Keep,x)),CellL_TableSort.Var1,'Un',0));
    Tissues = {'Adipose','Adrenal Gland','Artery','Bladder','Brain','Breast','Lymphoid','Myeloid',...
        'Cervix','Colon','Esophagus','Fallopian Tube','Heart','Kidney','Liver','Lung','Minor Salivary Gland','Muscle',...
        'Tibial Nerve','Ovarian','Pancreas', 'Pituitary Gland','Prostate','Skin','Small Intestines','Spleen','Stomach','Testis',...
        'Thyroid','Uterus','Vagina','Blood'};
    CopyN = [2 1 3 1 13 1 1 1 2 2 3 1 2 2 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1];
    TissueGTEX1 = arrayfun(@(x) repmat(string(Tissues{x}),1,CopyN(x)),1:length(Tissues),'Un',0);
    TissueGTEX = [TissueGTEX1{:}];
    corrTissues = arrayfun(@(x) find(strcmp(TissueGTEX,CellL_TableSort(x,3).Var3{:})),1:size(CellL_TableSort,1),'Un',0);
    corrTissues{3} = find(strcmp(TissueGTEX,'Lymphoid'));%Bone Marrow, Lymphocyte
    corrTissues{4} = find(strcmp(TissueGTEX,'Lymphoid'));%Bone Marrow, Epithelial
    corrTissues{5} = find(strcmp(TissueGTEX,'Lymphoid'));%Bone Marrow, Epithelial
    corrTissues{6} = find(strcmp(TissueGTEX,'Myeloid'));%Bone Marrow, Fibroblast
    corrTissues{15} = find(strcmp(TissueGTEX,'Artery'));%Endothelial, Endothelial cell, vascular endiotheliam 
    corrTissues{16} = [find(strcmp(TissueGTEX,'Lung')) find(strcmp(TissueGTEX,'Skin')) find(strcmp(TissueGTEX,'Small Intestines'))];%Eye, Epithelial cell (most common epithelial tissue is in lung, skin, and intestines
    corrTissues{17} = [find(strcmp(TissueGTEX,'Lung')) find(strcmp(TissueGTEX,'Skin')) find(strcmp(TissueGTEX,'Small Intestines'))];%Eye, Epithelial cell
    corrTissues{18} = [find(strcmp(TissueGTEX,'Esophagus')) find(strcmp(TissueGTEX,'Stomach')) find(strcmp(TissueGTEX,'Small Intestines'))];%GI Tract
    corrTissues{34} = find(strcmp(TissueGTEX,'Esophagus'));%Proximal Digestive Track
    corrTissues{43} = find(strcmp(TissueGTEX,'Bladder'));%Urinary Bladder, epithelial cell
    %cancer cell override
    corrTissues2 = cell(1,32);
    %REH is B acute lymphoblastic leukemia/lymphoma is 14
    corrTissues2{2} = 14;
    %NB4 is Acute promyelocytic leukemia  is 14
    corrTissues2{3} = 14;
    %U2OS Osteosarcoma   is 31 scarcoma
    corrTissues2{4} = 31;
    %PC-3  Prostate carcinoma is 19
    corrTissues2{5} = 19;
    %rh30 Alveolar rhabdomyosarcoma
    %U-251MG Astrocytoma
    %GAMG  Glioblastoma is 3 & 7
    corrTissues2{8} = [3 7];
    %MCF-7 breast carcinoma is 4
    corrTissues2{12} = 4;
    %Hela  Cervix/Uterus Adenocarcinoma is 5
    corrTissues2{13} = 5;
    %Siha  cervical squamous cell carcinoma is 28
    corrTissues2{14} = 28;
    %hTERT-RPE1 retinal pigement epithelium is best by 32 uveal melanoma
    corrTissues{16} = 32;
    %hTCEpi corneal eiphelial immortalized is best by 32 uveal melanoma
    corrTissues{17} = 32;
    %caco-2 colorectal adenocarcinoma.  6 IN TCGA
    corrTissues2{18} = 6;
    %A-549 Lung adenocarcinoma  is 12
    corrTissues2{24} = 12;
    %HDLM-2  Hodgkin lymphoma
    %Jurkat T-cell acute lymphoblastic leukemia  is 14
    corrTissues2{26} = 14;
    %HEL Acute erythroid leukemia Erythroleukemia  is 14
    corrTissues2{29} = 14;
    %K-562 Chronic myeloid leukemia  is 14
    corrTissues2{30} = 14;
    %THP-1 Acute monoblastic/monocytic leukemia  is 14
    corrTissues2{31} = 14;
    %HAP1 Chronic myeloid leukemia  is 14
    corrTissues2{32} = 14;
    %EFO-21 Ovarian serous cystadenocarcinoma. 16 in TCGA
    corrTissues2{33} = 16;
    %OE19 Esophageal adenocarcinoma. 30 is TCGA
    corrTissues2{34} = 30;
    %SK-MEL-30 Cutaneous melanoma. is 21 in TCGA
    corrTissues2{36} = 21;
    %A-431 Skin squamous cell carcinoma is  in TCGA use head and neck squamous%cell carcinoma
    corrTissues2{35} = 30;
    %SuSa Testicular teratoma . is 23 in TCGA 
    corrTissues2{42} = 23;
    %RT-4 Bladder carcinoma . is 2 in TCGA
    corrTissues2{43} = 2;
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
    GeneNormalTissue = GTEx_GeneV8_AN.gtexGeneV8;
    TranNormalTissue = GTEx_Transc_AN.gtexTranscExpr;
    GeneCancerTissue = TCGA_Gene_AN.tcgaGeneExpr;
    TranCancerTissue = TCGA_Transc_AN.tcgaTranscExpr;
    
    
    strcmp(sTabel.geneType
    unique(sTabel.transcriptType)
    
    optsInfo = detectImportOptions('DatabaseData/Homo_sapiens.gene_info','FileType','delimitedtext');
    geneInfo_db = tdfread('DatabaseData/Homo_sapiens.gene_info','\t',optsInfo); 
    geneInfoName = arrayfun(@(x) strtrim(geneInfo_db.Symbol(x,:)),1:size(geneInfo_db.GeneID),'Un',0);
    geneInfoSyno = arrayfun(@(x) strsplit(strtrim(geneInfo_db.Synonyms(x,:)),'|'),1:size(geneInfo_db.GeneID),'Un',0);
    geneInfoType = arrayfun(@(x) strtrim(geneInfo_db.type_of_gene(x,:)),1:size(geneInfo_db.GeneID),'Un',0);
    
    geneInfoNamesList = [geneInfoSyno{:}];
    allNames = unique(geneInfoNamesList);
    NamesMatchMatrix = CATnWrapper(cellfun(@(x) ismember(allNames,['*' x])',geneInfoSyno,'Un',0),2);

    GNT_Type = strtrim({GeneNormalTissue.geneType});
    GNT_Name = strtrim({GeneNormalTissue.name});
    GNT_ID = strtrim({GeneNormalTissue.geneId});
    Is_3prime_overlapping_ncRNA = find(strcmp(GNT_Type,'3prime_overlapping_ncRNA'));
    Is_antisense = find(strcmp(GNT_Type,'antisense'));
    Is_bidirectional_promoter_lncRNA = find(strcmp(GNT_Type,'bidirectional_promoter_lncRNA'));
    Is_lincRNA = find(strcmp(GNT_Type,'lincRNA'));
    Is_macro_lncRNA = find(strcmp(GNT_Type,'macro_lncRNA'));
    Is_sense_intronic = find(strcmp(GNT_Type,'sense_intronic'));
    is_sense_overlapping = find(strcmp(GNT_Type,'sense_overlapping'));
    Is_non_coding = find(strcmp(GNT_Type,'non_coding'));
    Check_ForPairs = [Is_3prime_overlapping_ncRNA Is_antisense Is_bidirectional_promoter_lncRNA Is_lincRNA Is_macro_lncRNA Is_sense_intronic is_sense_overlapping Is_non_coding];
    
    %take all lncRNA in RefSeq set and look for in ensembl set? (have diff
    %policy for no matches
    
    
    find(strcmp(allNames,GNT_Name{Check_ForPairs(1)}))
    
    NamesMatchMatrix(:,1);
    
    
    
    strcmp(geneInfoNamesList
    %list all gene names/transcripts for both reference.
%matches and paired sense/anti-sense, any non-paired anti-sense
    EC = find(~cellfun(@isempty,{GeneNormalTissue.geneId})); 
    GCT_Name = strtrim({GeneCancerTissue.name2});
    GCT_ID = strtrim({GeneCancerTissue.name});
    TNT_Name = strtrim({TranNormalTissue.name2});
    TNT_ID = strtrim({TranNormalTissue.name});
    TCT_Name = strtrim({TranCancerTissue.name2});
    TCT_ID = strtrim({TranCancerTissue.name});

    transcriptOverlap = @(x,S) (ge([GeneNormalTissue.chromStart],GeneNormalTissue(x).chromStart-S).*le([GeneNormalTissue.chromStart],GeneNormalTissue(x).chromEnd+S)+...
                             ge([GeneNormalTissue.chromEnd],GeneNormalTissue(x).chromStart-S).*le([GeneNormalTissue.chromEnd],GeneNormalTissue(x).chromEnd+S) + ...
                             le([GeneNormalTissue.chromStart],GeneNormalTissue(x).chromStart-S).*ge([GeneNormalTissue.chromEnd],GeneNormalTissue(x).chromEnd+S) + ...
                             ge([GeneNormalTissue.chromStart],GeneNormalTissue(x).chromStart-S).*le([GeneNormalTissue.chromEnd],GeneNormalTissue(x).chromEnd+S)).*...
                             strcmp({GeneNormalTissue.chrom},GeneNormalTissue(x).chrom).*~strcmp({GeneNormalTissue.name},GeneNormalTissue(x).name); %does not show ones where not overlapping
    transcriptOverlap2 = @(x,S) (ge([TranNormalTissue.chromStart],GeneNormalTissue(x).chromStart-S).*le([TranNormalTissue.chromStart],GeneNormalTissue(x).chromEnd+S)+...
                             ge([TranNormalTissue.chromEnd],GeneNormalTissue(x).chromStart-S).*le([TranNormalTissue.chromEnd],GeneNormalTissue(x).chromEnd+S) + ...
                             le([TranNormalTissue.chromStart],GeneNormalTissue(x).chromStart-S).*ge([TranNormalTissue.chromEnd],GeneNormalTissue(x).chromEnd+S) + ...
                             ge([TranNormalTissue.chromStart],GeneNormalTissue(x).chromStart-S).*le([TranNormalTissue.chromEnd],GeneNormalTissue(x).chromEnd+S)).*...
                             strcmp({TranNormalTissue.chrom},GeneNormalTissue(x).chrom); %does not show ones where not overlapping
                         
    fCheck_ForPairs_0 = arrayfun(@(x) find(transcriptOverlap(x,0)>0),Check_ForPairs,'Un',0);%13749
    Dist0_Pairs_Status = ~cellfun(@isempty,fCheck_ForPairs_0);
    Found_Pairs_Dist0 = find(Dist0_Pairs_Status); 
    NotFound_Pairs_Dist0 = find(~Dist0_Pairs_Status);
    fCheck_ForPairs_250 = arrayfun(@(x) find(transcriptOverlap(x,250)>0),NotFound_Pairs_Dist0,'Un',0);%6170
    Dist250_Pairs_Status = ~cellfun(@isempty,fCheck_ForPairs_250);
    Found_Pairs_Dist250 = find(Dist250_Pairs_Status); 
    NotFound_Pairs_Dist250 = find(~Dist250_Pairs_Status);
    fCheck_ForPairs_500 = arrayfun(@(x) find(transcriptOverlap(x,500)>0),NotFound_Pairs_Dist0(NotFound_Pairs_Dist250),'Un',0);%2581
    Dist500_Pairs_Status = ~cellfun(@isempty,fCheck_ForPairs_500);
    Found_Pairs_Dist500 = find(Dist500_Pairs_Status); 
    NotFound_Pairs_Dist500 = find(~Dist500_Pairs_Status);
    fCheck_ForPairs_1000 = arrayfun(@(x) find(transcriptOverlap(x,1000)>0),NotFound_Pairs_Dist0(NotFound_Pairs_Dist250(NotFound_Pairs_Dist500)),'Un',0);%2431
    Dist1000_Pairs_Status = ~cellfun(@isempty,fCheck_ForPairs_1000);
    Found_Pairs_Dist1000 = find(Dist1000_Pairs_Status); 
    NotFound_Pairs_Dist1000 = find(~Dist1000_Pairs_Status);
    fCheck_ForPairs_2000 = arrayfun(@(x) find(transcriptOverlap(x,2000)>0),NotFound_Pairs_Dist0(NotFound_Pairs_Dist250(NotFound_Pairs_Dist500(NotFound_Pairs_Dist1000))),'Un',0);%2257
    Dist2000_Pairs_Status = ~cellfun(@isempty,fCheck_ForPairs_2000);
    Found_Pairs_Dist2000 = find(Dist2000_Pairs_Status); %239
    NotFound_Pairs_Dist2000 = find(~Dist2000_Pairs_Status);
    Pairs = cell(1,length(Check_ForPairs));
    Pairs{Found_Pairs_Dist0} = fCheck_ForPairs_0{Found_Pairs_Dist0}; %closest + anti-sense to mRNA diff strands
    Pairs{NotFound_Pairs_Dist0(Found_Pairs_Dist250)} = fCheck_ForPairs_250{Found_Pairs_Dist250}; %closest + anti-sense to mRNA diff strands
    Pairs{NotFound_Pairs_Dist0(NotFound_Pairs_Dist250(Found_Pairs_Dist500))} = fCheck_ForPairs_500{Found_Pairs_Dist500}; %closest + anti-sense to mRNA diff strands
    Pairs{NotFound_Pairs_Dist0(NotFound_Pairs_Dist250(NotFound_Pairs_Dist500(Fount_Pairs_Dist1000)))} = fCheck_ForPairs_1000{Found_Pairs_Dist1000}; %closest + anti-sense to mRNA diff strands
    Pairs{NotFound_Pairs_Dist0(NotFound_Pairs_Dist250(NotFound_Pairs_Dist500(NotFount_Pairs_Dist1000(Fount_Pairs_Dist2000))))} = fCheck_ForPairs_2000{Found_Pairs_Dist2000}; %closest + anti-sense to mRNA diff strands
    
    %get alias
    %zinc finger
    Check_ForPairs(117)
    fCheck_ForPairs_0{117}
    %need is alias data, help make quick matches
    
        


    %All Pairs Biotypes be mRNA, opposite strand, if multiple closest, must
    %not be lncRNA? 
    
    %Last find go from RefSeq Gene ID <-> ENSEMBL Gene ID
    %Alias names (different).
    
    %can you do check pairs for refseq reference.
    %how to get pairs  (or is this already in the sTabel
    
    %find pairs (mRNA-lncRNA, ENSEMBL/RefSeq)
    %translation IDs mRNA , ENSEMBL-RefSeq
    %translation IDs lncRNA , ENSEMBL-RefSeq
    %Add in Alias info.
    
    
%Not in: what about ncRNA not with anti-sense just use average tissue level
%Not in: with pair use ratio
%keep the same
%mention antisense


TargetIsLncRNA = find(strcmp(transcriptType_List,'lnc_RNA'));
TargetIsAntiSenseRNA = find(strcmp(transcriptType_List,'antisense_RNA'));
GI = GFFAnnotation('DatabaseData/Homo_sapiens.GRCh38.111.gff3');   
lncRNA_Feature = find(strcmp(GI.Feature,'lnc_RNA'));
GI_gene = arrayfun(@(x) extractBetween(GI.Attributes{x},'gene:',';'),1:length(GI.Attributes),'Un',0);
GI_transcript = arrayfun(@(x)  extractBetween(GI.Attributes{x},'transcript:',';'),1:length(GI.Attributes),'Un',0);
GI_Name = arrayfun(@(x) extractBetween(GI.Attributes{x},'Name=',';'),1:length(GI.Attributes),'Un',0);
GI_biotype = arrayfun(@(x) extractBetween(GI.Attributes{x},'biotype=',';'),1:length(GI.Attributes),'Un',0);
GI_transcriptid =  arrayfun(@(x) extractBetween(GI.Attributes{x},'transcript_id=',';'),1:length(GI.Attributes),'Un',0);
GI_geneid =  arrayfun(@(x) extractBetween(GI.Attributes{x},'gene_id=',';'),1:length(GI.Attributes),'Un',0);
%two things one names have -AS  / sense/pairs dont, aliasses
EC = find(~cellfun(@isempty,GI_geneid));
GI_transcript{lncRNA_Feature(1)}
    %some lncRNA do not have sense pair?  
fid = fopen('DatabaseData/Homo_sapiens.GRCh38.111.gtf','rt');
GJ = textscan(fid,'%s%s%s%d%d%s%s%s%s%*[^\n]',...
    'delimiter','\t','CommentStyle','#');
fclose(fid);
p5 = find(strcmp(GJ{1,3},'five_prime_utr'));
p3 = find(strcmp(GJ{1,3},'three_prime_utr'));
%4 5 start end
GJ_geneid =  arrayfun(@(x) extractBetween(x,'gene_id ',';'),GJ{1,9},'Un',0);
GJ_gene =  arrayfun(@(x) extractBetween(x,'gene_name ',';'),GJ{1,9},'Un',0);
find(contains([GJ_gene{:}],'HELLPAR'))
test=cellfun(@(x) contains(x,'HELLPAR'),[GJ_gene{:}])
  


%mRNA and lncRNA have similar gene_id


    TargetIsLncRNA = find(arrayfun(@(id) contains(inputs1{id,5},'-AS'),1:size(inputs1,1)));
    TissueGeneLevel_Name_mRNA_Loc1 = arrayfun(@(X) find(strcmp(TissueGeneLevel_Name,extractBefore(inputs1{X,5}(2:end-1),'-AS'))),TargetIsLncRNA,'Un',0);
    TissueGeneLevel_Name_mRNA_Loc2 = arrayfun(@(X) find(strcmp(TissueGeneLevel_Name,extractBefore(inputs1{X,6}(2:end-1),'-AS'))),TargetIsLncRNA,'Un',0);
    TissueGeneLevel_Name_mRNA_Loc = arrayfun(@(x) unique([TissueGeneLevel_Name_mRNA_Loc1{x} TissueGeneLevel_Name_mRNA_Loc2{x}]),1:length(TargetIsLncRNA),'Un',0);
    TissueGeneLevel_Name_lRNA_Loc1 = arrayfun(@(X) find(strcmp(TissueGeneLevel_Name,inputs1{X,6}(2:end-1))),TargetIsLncRNA,'Un',0);
    TissueGeneLevel_Name_lRNA_Loc2 = arrayfun(@(X) find(strcmp(TissueGeneLevel_geneId,inputs1{X,6}(2:end-1))),TargetIsLncRNA,'Un',0);
    TissueGeneLevel_Name_lRNA_Loc = arrayfun(@(x) [TissueGeneLevel_Name_lRNA_Loc1{x} TissueGeneLevel_Name_lRNA_Loc2{x}],1:length(TargetIsLncRNA),'Un',0);
    for vi = 1:length(TargetIsLncRNA)
        lncLOC = find(strcmp(RefSeq_List,inputs1{TargetIsLncRNA(vi),1}));
        mrnaLOC = find(strcmp(RefSeq_List,inputs1{find(strcmp({inputs1{:,5}},strcat('(',extractBefore(inputs1{TargetIsLncRNA(vi),5}(2:end-1),'-AS'),')')),1),1}));
        for y = 1:43
        sTabel_CellLine_TPM(lncLOC,CellLines_FinFilt(y)) = sTabel_CellLine_TPM(mrnaLOC,CellLines_FinFilt(y))*mean(GeneNormalTissue(TissueGeneLevel_Name_lRNA_Loc{vi},:).expScores(corrTissues{y})./GeneNormalTissue(TissueGeneLevel_Name_mRNA_Loc{vi},:).expScores(corrTissues{y}));
        end
    end
    TissueGeneLevel2_Name = strtrim({GeneCancerTissue.name2});
    TissueGeneLevel2_geneId = strtrim({GeneCancerTissue.name});
    TissueGeneLevel2_Name_mRNA_Loc1 = arrayfun(@(X) find(strcmp(TissueGeneLevel2_Name,extractBefore(inputs2{X,5}(2:end-1),'-AS'))),TargetIsLncRNA,'Un',0);
    TissueGeneLevel2_Name_mRNA_Loc2 = arrayfun(@(X) find(strcmp(TissueGeneLevel2_Name,extractBefore(inputs2{X,5}(2:end-1),'-AS'))),TargetIsLncRNA,'Un',0);
    TissueGeneLevel2_Name_mRNA_Loc = arrayfun(@(x) unique([TissueGeneLevel2_Name_mRNA_Loc1{x} TissueGeneLevel2_Name_mRNA_Loc2{x}]),1:length(TargetIsLncRNA),'Un',0);
    TissueGeneLevel2_Name_lRNA_Loc1 = arrayfun(@(X) find(strcmp(TissueGeneLevel2_Name,inputs2{X,5}(2:end-1))),TargetIsLncRNA,'Un',0);
    TissueGeneLevel2_Name_lRNA_Loc2 = arrayfun(@(X) find(strcmp(TissueGeneLevel2_geneId,inputs2{X,6}(2:end-1))),TargetIsLncRNA,'Un',0);
    TissueGeneLevel2_Name_lRNA_Loc = arrayfun(@(x) unique([TissueGeneLevel2_Name_lRNA_Loc1{x} TissueGeneLevel2_Name_lRNA_Loc2{x}]),1:length(TargetIsLncRNA),'Un',0);
    for vi = 1:length(TargetIsLncRNA)
        lncLOC = find(strcmp(RefSeq_List,inputs1{TargetIsLncRNA(vi),1}));
        mrnaLOC = find(strcmp(RefSeq_List,inputs1{find(strcmp({inputs1{:,5}},strcat('(',extractBefore(inputs1{TargetIsLncRNA(vi),5}(2:end-1),'-AS'),')')),1),1}));
        for y = find(~cellfun(@isempty,corrTissues2))
            mrnaCancer = str2double(split(GeneCancerTissue(TissueGeneLevel2_Name_mRNA_Loc{vi},:).expScores,'.'));
            lncCancer = str2double(split(GeneCancerTissue(TissueGeneLevel2_Name_lRNA_Loc{vi},:).expScores,','));
        sTabel_CellLine_TPM(lncLOC,CellLines_FinFilt(y)) = sTabel_CellLine_TPM(mrnaLOC,CellLines_FinFilt(y))*mean(lncCancer(corrTissues2{y})./mrnaCancer(corrTissues2{y}));
        end
    end
    %table matches ENST_IDs to ENSG_IDs
    globalExpressionMatrix.CellLineValues_Gene = sTabel_CellLine_TPM;
    globalExpressionMatrix.CellLineValues_inferTran = sTabel_CellLine_TPM;
    globalExpressionMatrix.tcgaValues_Gene = 1;
    globalExpressionMatrix.tcgaValues_Tran = 1;
    globalExpressionMatrix.gtexValues_Gene = 1;
    globalExpressionMatrix.gtexValues_Tran = 1;
    globalExpressionMatrix.RefSeq_TranIDs = 1;
    globalExpressionMatrix.ENSEMBL_GeneIDs = 1;
    globalExpressionMatrix.ENSEMBL_TranIDs = 1;
    globalExpressionMatrix.GeneName = 1;
    globalExpressionMatrix.TranscriptName = 1;
    globalExpressionMatrix.CellLines = 1;
    globalExpressionMatrix.tgcaCancers = 1;
    globalExpressionMatrix.gtexTissues = 1;


    
    %all cancer tissue values
    %all normal tissue values
    %RefSeq_List
    %Gene Name List
    %Cell Line Order List;
    
end

for v = 1:length(outGroup)
settings.Human_wgEncodeGencodeRefSeqFile = 'DatabaseData/wgEncodeGencodeRefSeqV44.txt';
settings.Human_wgEncodeGencodeAttributesFile = 'DatabaseData/wgEncodeGencodeAttrsV44.txt';
settings.Human_wgEncodeGencodeCompFile = 'DatabaseData/wgEncodeGencodeCompV44.txt';
settings.wgEncodeGencodeAttributesVariableNames = {'geneId','geneName','geneType','geneStatus',...
            'transcriptId','transcriptName','transcriptType','transcriptStatus',...
            'havanaGeneId','havanaTranscriptId','ccdsId','level','transcriptClass','proteinId'};
settings.wgEncodeGencodeCompVariableNames = {'nx','transcriptId','chrom','strand','txStart',...
        'txEnd','cdsStart','cdsEnd','exonCount',...
        'exonStarts','exonEnds','score','name2','cdsStartStat','cdsEndStat','exonFrames'}; 
settings.Human_GencodeRefSeqMetadataFile = 'DatabaseData/gencode.v44.metadata.RefSeq';
settings.HumanGenomeAssemblyReportFile = 'Metadata/Human/GCF_000001405.39_GRCh38.p13_assembly_report.txt';    
    ProbeDesignResults3.SoftwareResults(1).ExpressionMatrix
    ProbeDesignResults3.SoftwareResults(1).uniNames
    ProbeDesignResults1.SoftwareResults(1).Names
    




end
end