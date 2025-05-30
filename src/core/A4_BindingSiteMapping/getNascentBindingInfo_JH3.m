function [DoesProbeBindSite_TS,Kb_TS,OnTargetTS,NumPolymeraseOnTS,maxPolymeraseOnTS,TS_SenseInfo,TS_ChrInfo,AN_ListTS,Num_of_Molecule_SitesTS,dHeq_TS,dSeq_TS,dHf_TS,dSf_TS,dHr_TS,dSr_TS,dCp_TS] = getNascentBindingInfo_JH3(NumberOfPolymerase,settings,probes,transcript_ID,DoesProbeBindSite,Kb,cTypeIndx,Ix_DNA,T_hybrid,AN_List,MolProbesAtEvents,Mol_ProbesAtEventsID,Event_Name,Event_Strands,Event_LocStart,Event_ID,dHeq_mod,dSeq_mod,dHf_mod,dSf_mod,dHr_mod,dSr_mod,dCp_mod)
seqType = 'RNA';
if (strcmp(settings.Organism,'Human'))
    PolymeraseSpacing = 90;
    PolymeraseOccupancy = 8;
elseif (strcmp(settings.Organism,'Mouse'))
    PolymeraseSpacing = 90;
    PolymeraseOccupancy = 8;
elseif (strcmp(settings.Organism,'Yeast'))
    PolymeraseSpacing = 90;
    PolymeraseOccupancy = 8;
end

N_methods = size(Kb,4);
N_methods2 = 3;
%All RNA Hits are locations on RNA 
%Exon's database align databases Off-targets values align up.
%want to use LocStart which should align.
%for On-Target use alignment of probes.
%TargetSite_RelativeLocationMap
%DoesProbeBindSite_TS %does probe bind nascent site
%POGmod_TS            %nascent site binding energy
%Kb_TS                %nascent site binding rate
%OnTargetTS           %on-target or off-target  
%NumPolymeraseOnTS    %number of polymerase
%TS_SenseInfo         %sense or antisense
%TS_ChrInfo           %what chromosome it is on

% RPKM,  TPM,,  FPKM
% UCSC reports RPKM,    visual is TPM
%expScores per-tissue median expression levels (RPKM)
%score log transformed and scaled
%legend1 = [p1.mainLine p2.mainLine p3.mainLine p4.mainLine];
%    legend(legend1,[xlabel1 ' with ' ylabel1],[ylabel1 ' with ' xlabel1],['Randomized spots with ' ylabel1],['Randomized spots with ' xlabel1],'Location','northwest')
Lpmin = min(cell2mat(cellfun(@length,{probes{:,2}},'UniformOutput',false)));
Lpmax = max(cell2mat(cellfun(@length,{probes{:,2}},'UniformOutput',false)));
DoesProbeBindSite_TS = ndSparse.build([size(probes,1),1,size(probes,1)-Lpmin+1],0);
Kb_TS = ndSparse.build([size(probes,1),1,size(probes,1)-Lpmin+1,N_methods],0); 
dHeq_TS = ndSparse.build([size(probes,1),1,size(probes,1)-Lpmin+1,N_methods],0); 
dSeq_TS = ndSparse.build([size(probes,1),1,size(probes,1)-Lpmin+1,N_methods],0); 
dCp_TS = ndSparse.build([size(probes,1),1,size(probes,1)-Lpmin+1,N_methods],0); 
dHf_TS = ndSparse.build([size(probes,1),1,size(probes,1)-Lpmin+1,N_methods2],0); 
dSf_TS = ndSparse.build([size(probes,1),1,size(probes,1)-Lpmin+1,N_methods2],0); 
dHr_TS = ndSparse.build([size(probes,1),1,size(probes,1)-Lpmin+1,N_methods2],0); 
dSr_TS = ndSparse.build([size(probes,1),1,size(probes,1)-Lpmin+1,N_methods2],0); 
Expression_Null = settings.DoAllGenesHaveSameExpression; %do genes have same expression
CellTypeID = settings.CellType_ExprID;
Organism = settings.Organism;
%Human Specific Expression Sequence/Level Settings
HumanOntology = settings.HumanSpecific.Ontology;%Normal/Cancer
HumanTissueOnly = settings.HumanSpecific.TissueOrTissueAndCellType; % 0 (TissueOnssue/Cell Type)
HumanExpGeneOrTransc = settings.HumanSpecific.HumanExpGeneOrTransc; % 1 (Gene/EMBL GENEID) , 0 (Transcript/EMBL Transcript ID)
Human_SelTrack = settings.HumanSpecific.SCOutputTrack;
HumanSCTracks = settings.HumanSpecific.SCTracks;
isOffline = settings.IsOffline;
expValType = settings.expressionValType;% 1-4 (expCounts,expValues,mean(CellTypeExpValues),one cell types CellTypeExpValues)
if (isOffline)         
    scTracksBedFiles = {'colonWang_cell_type','ileumWang_cell_type','rectumWang_cell_type',...
                    'bloodHaoL1_cell_type','bloodHaoL2_cell_type','bloodHaoL3_cell_type','cortexVelmeshev_cell_type',...
                    'fetalGeneAtlas_cell_type','fetalGeneAtlas_Organ','fetalGeneAtlas_Organ_cell_lineage',...
                    'heartAtlas_cell_type','heartAtlas_cell_states','kidneyStewart_broad_cell_type',...
                    'kidneyStewart_cell_type','kidneyStewart_compartment','kidneyStewart_detailed_cell_type',...
                    'liverMacParland_broad_cell_type','liverMacParland_cell_type',...
                    'lungTravaglini10x_cell_type','lungTravaglini10x_detailed_cell_type','muscleDeMicheli_cell_type',...
                    'pancreasBaron_detailed_cell_type','pancreasBaron_cell_type',...
                    'placentaVentoTormo10x_cell_type1','placentaVentoTormo10x_cell_type','skinSoleBoldo_cell_type'};
    scTracks = {'colonWangCellType','ileumWangCellType','rectumWangCellType',...
        'bloodHaoCellType','bloodHaoL2','bloodHaoL3','cortexVelmeshevCellType',...
        'fetalGeneAtlasCellType','fetalGeneAtlasOrgan','fetalGeneAtlasOrganCellLineage',...
        'heartAtlasCellType','heartAtlasCellStates','kidneyStewartBroadCellType',...
        'kidneyStewartCellType','kidneyStewartCompartment','kidneyStewartDetailedCellType',...
        'liverMacParlandBroadCellType','liverMacParlandCellType',...
        'lungTravaglini2020CellType10x','lungTravaglini2020DetailedCellType10x','muscleDeMicheliCellType',...
        'pancreasBaronDetailedCellType','pancreasBaronCellType',...
        'placentaVentoTormoCellType10x','placentaVentoTormoCellDetailed10x','skinSoleBoldoCellType'};
    if (strcmp(Organism,'Mouse'))
        optsTabulaMuris = detectImportOptions('DatabaseData/tabulamuris_barChart.bed','FileType','delimitedtext');
        optsTabulaMuris.VariableNames = {'chrom','chromStart','chromEnd','name',...
        'score','strand','name2','expCount','expScores','dataOffset','dataLen'};
        TabulaMuris = tdfread('DatabaseData/tabulamuris_barChart.bed','\t',optsTabulaMuris);
        TabulaMuris = cell2struct(struct2cell(TabulaMuris),optsTabulaMuris.VariableNames);
        TabulaMuris = struct2table(TabulaMuris);
        TabulaMuris_AN.tabulamurisBarChart = table2struct(TabulaMuris);
    end
    if (strcmp(Organism,'Human'))
        if (strcmp(HumanOntology,'Cancer'))
            optsTcgaTransc= detectImportOptions('DatabaseData/tcgaTranscExpr.bed','FileType','delimitedtext');
            optsTcgaGene = detectImportOptions('DatabaseData/tcgaGeneExpr.bed','FileType','delimitedtext');
            optsTcgaTransc.VariableNames = {'chrom','chromStart','chromEnd','name',...
            'score','strand','name2','expCount','expScores','dataOffset','dataLen'};
            optsTcgaGene.VariableNames = {'chrom','chromStart','chromEnd','name',...
            'score','strand','name2','expCount','expScores','dataOffset','dataLen'};
            TcgaTransc = tdfread('DatabaseData/tcgaTranscExpr.bed','\t',optsTcgaTransc);
            TcgaGene = tdfread('DatabaseData/tcgaGeneExpr.bed','\t',optsTcgaGene);
            TcgaGene = cell2struct(struct2cell(TcgaGene),optsTcgaGene.VariableNames);
            TcgaGene = struct2table(TcgaGene);
            TCGA_Gene_AN.tcgaGeneExpr = table2struct(TcgaGene);
            TcgaTransc = cell2struct(struct2cell(TcgaTransc),optsTcgaTransc.VariableNames);
            TcgaTransc = struct2table(TcgaTransc);
            TCGA_Transc_AN.tcgaTranscExpr = table2struct(TcgaTransc);
        end
        if (strcmp(HumanOntology,'Normal'))
            if (HumanTissueOnly==0)
                optsGtexTransc = detectImportOptions('DatabaseData/gtexTranscExpr.bed','FileType','delimitedtext');
                optsGtexGeneV8 = detectImportOptions('DatabaseData/gtexGeneV8.txt','FileType','delimitedtext');
                optsGtexTransc.VariableNames = {'chrom','chromStart','chromEnd','name',...
                'score','strand','name2','expCount','expScores','dataOffset','dataLen'};
                optsGtexGeneV8.VariableNames = {'chrom','chromStart','chromEnd','name',...
                'score','strand','geneId','geneType','expCount','expScores'};
                GtexTransc = tdfread('DatabaseData/gtexTranscExpr.bed','\t',optsGtexTransc);
                GtexGeneV8 = tdfread('DatabaseData/gtexGeneV8.txt','\t',optsGtexGeneV8);
                GtexGeneV8 = cell2struct(struct2cell(GtexGeneV8),optsGtexGeneV8.VariableNames);
                GtexGeneV8 = struct2table(GtexGeneV8);
                GTEx_GeneV8_AN.gtexGeneV8 = table2struct(GtexGeneV8);
                GtexTransc = cell2struct(struct2cell(GtexTransc),optsGtexTransc.VariableNames);
                GtexTransc = struct2table(GtexTransc);
                GTEx_Transc_AN.gtexTranscExpr = table2struct(GtexTransc);
            else
                optsHumanSS.VariableNames = {'chrom','chromStart','chromEnd','name',...
                'score','strand','name2','expCount','expScores'};
                optsGeneSS{1:max(HumanSCTracks)} = [];
                ssGene{1:max(HumanSCTracks)} = [];
                scTrack_Gene_AN{1:max(HumanSCTracks)} = [];
                for j=HumanSCTracks
                    optsGeneSS{j} = detectImportOptions(strcat('DatabaseData/',scTracksBedFiles{j},'.bed'),'FileType','delimitedtext');
                    optsGeneSS{j}.VariableNames = optsHumanSS.VariableNames;
                    ssGene{j} = tdfread(strcat('DatabaseData/',scTracksBedFiles{j},'.bed'),'\t',optsGeneSS);
                    ssGene{j} = cell2struct(struct2cell(ssGene{j}),optsGeneSS{j}.VariableNames);
                    ssGene{j} = struct2table(ssGene{j});
                    scTrack_Gene_AN{j}.(scTracks{j}) = table2struct(ssGene{j});        
                end
            end 
        end
    end
    if (strcmp(Organism,'Human'))
        optsUCSC = detectImportOptions('DatabaseData/wgEncodeGencodeRefSeqV39.txt');
        optsUCSC.VariableNames = {'Var1','Var2','Var3'};
        EMBLRefSeqAlign_db = readtable('DatabaseData/wgEncodeGencodeRefSeqV39.txt',optsUCSC); 
        optsUCSC2 = detectImportOptions('DatabaseData/wgEncodeGencodeAttrsV39.txt');
        optsUCSC2.VariableNames = {'geneId','geneName','geneType','geneStatus',...
            'transcriptId','transcriptName','transcriptType','transcriptStatus',...
            'havanaGeneId','havanaTranscriptId','ccdsId','level','transcriptClass','proteinId'};
        EMBLAttrAlign_db = readtable('DatabaseData/wgEncodeGencodeAttrsV39.txt',optsUCSC2); 
        optsUCSC3 = detectImportOptions('DatabaseData/wgEncodeGencodeCompV39.txt');
        optsUCSC3.VariableNames = {'nx','transcriptId','chrom','strand','txStart',...
            'txEnd','cdsStart','cdsEnd','exonCount',...
            'exonStarts','exonEnds','score','name2','cdsStartStat','cdsEndStat','exonFrames'};
        EMBLComp_db = readtable('DatabaseData/wgEncodeGencodeCompV39.txt',optsUCSC3);
    elseif (strcmp(Organism,'Mouse'))
        optsUCSC = detectImportOptions('DatabaseData/wgEncodeGencodeRefSeqVM25.txt');
        optsUCSC.VariableNames = {'Var1','Var2','Var3'};
        EMBLRefSeqAlign_db = readtable('DatabaseData/wgEncodeGencodeRefSeqVM25.txt',optsUCSC); 
        optsUCSC2 = detectImportOptions('DatabaseData/wgEncodeGencodeAttrsVM25.txt');
        optsUCSC2.VariableNames = {'geneId','geneName','geneType','geneStatus',...
            'transcriptId','transcriptName','transcriptType','transcriptStatus',...
            'havanaGeneId','havanaTranscriptId','ccdsId','level','transcriptClass','proteinId'};
        EMBLAttrAlign_db = readtable('DatabaseData/wgEncodeGencodeCompVM25.txt',optsUCSC2); 
        optsUCSC3 = detectImportOptions('DatabaseData/wgEncodeGencodeCompVM25.txt');
        optsUCSC3.VariableNames = {'nx','transcriptId','chrom','strand','txStart',...
            'txEnd','cdsStart','cdsEnd','exonCount',...
            'exonStarts','exonEnds','score','name2','cdsStartStat','cdsEndStat','exonFrames'};
        EMBLComp_db = readtable('DatabaseData/wgEncodeGencodeCompVM25.txt',optsUCSC3);
    elseif (strcmp(Organism,'Yeast'))
    end
end
           
if (Expression_Null==0)
    if (strcmp(Organism,'Human')) 
        if (strcmp(HumanOntology,'Cancer'))
            if (HumanExpGeneOrTransc)
                ExprGene_AN = TCGA_Gene_AN.tcgaGeneExpr;
                expCase = 0;
            else
                ExprTransc_AN = TCGA_Transc_AN.tcgaTranscExpr;
                expCase = 0;
            end
        end
        if (strcmp(HumanOntology,'Normal'))
            if (HumanTissueOnly==0)
                if (HumanExpGeneOrTransc)
                ExprGene_AN = GTEx_GeneV8_AN.gtexGeneV8;
                expCase = 0;
                else
                ExprTransc_AN = GTEx_Transc_AN.gtexTranscExpr;
                expCase = 1;
                end  
            else
                scExprGene_AN = scTrack_Gene_AN{Human_SelTrack}.(scTracks{Human_SelTrack});
                expCase = 2;
            end   
        end
    elseif strcmp(Organism,'Mouse')
        ExprGene_AN = TabulaMuris_AN.tabulamurisBarChart;
        HumanExpGeneOrTransc = 1;expCase = 0;
    elseif strcmp(Organism,'Yeast')
    end 
else
end
            
%% Nascent RNA Targets
for i=1:length(transcript_ID)
    onEMBL_transcriptID{i} = EMBLRefSeqAlign_db.Var1{find(strcmp(transcript_ID{i},EMBLRefSeqAlign_db.Var2))};
    onEMBL_geneID{i} = EMBLAttrAlign_db.geneId{find(strcmp(onEMBL_transcriptID{i},EMBLAttrAlign_db.transcriptId))};
end
for i = 1:length(transcript_ID)
       if (expCase==2) 
            for j=HumanSCTracks
                tempGene = str2double(split(scExprGene_AN(find(strcmp(extractBefore(onEMBL_geneID{i},'.'),{scExprGene_AN.name}))).expScores,','));
                ExpValues = tempGene(1:end-1);clear tempGene
            end  
       elseif (expCase==1)
           if (HumanExpGeneOrTransc)
                ExpValues =  ExprGene_AN(find(strcmp(extractBefore(onEMBL_geneID{i},'.'),{ExprGene_AN.geneId}))).expScores; 
           else
                ExpValues = ExprTransc_AN(find(strcmp(extractBefore(onEMBL_transcriptID{i},'.'),{ExprTransc_AN.name}))).expScores; 
           end         
       elseif (expCase==0)
           if (HumanExpGeneOrTransc)
                ExpValues = ExprGene_AN(find(strcmp(extractBefore(onEMBL_geneID{i},'.'),{ExprGene_AN.name}))).expScores;
           else     
                ExpValues = ExprTransc_AN(find(strcmp(extractBefore(onEMBL_transcriptID{i},'.'),{ExprTransc_AN.name}))).expScores;
           end
       end
       switch expValType
            case 3
                ExpOn(i) = mean(ExpValues);
            case 4
                ExpOn(i) = ExpValues(CellTypeID);
            otherwise
                ExpOn(i) = mean(ExpValues);
       end 
end

for i = 1:length(onEMBL_transcriptID)
    temp_txStart = EMBLComp_db.txStart(find(strcmp(onEMBL_transcriptID{i},EMBLComp_db.transcriptId)));
    temp_txEnd = EMBLComp_db.txEnd(find(strcmp(onEMBL_transcriptID{i},EMBLComp_db.transcriptId)));
    onEMBL_exonCount = EMBLComp_db.exonCount(find(strcmp(onEMBL_transcriptID{i},EMBLComp_db.transcriptId)));
    temp_exonStarts = EMBLComp_db.exonStarts{find(strcmp(onEMBL_transcriptID{i},EMBLComp_db.transcriptId))};
    temp_exonEnds = EMBLComp_db.exonEnds{find(strcmp(onEMBL_transcriptID{i},EMBLComp_db.transcriptId))};
    temp_exonStarts2 = strsplit(temp_exonStarts,',');
    temp_exonEnds2 = strsplit(temp_exonEnds,',');
    temp_exonStarts3 = cell2mat(cellfun(@(x) str2num(x),temp_exonStarts2,'UniformOutput',false));
    temp_exonEnds3 = cell2mat(cellfun(@(x) str2num(x),temp_exonEnds2,'UniformOutput',false));
    temp_Span = temp_txStart:temp_txEnd;
    temp_counts = [];
    for j = 1:onEMBL_exonCount
    temp_counts = [temp_counts temp_exonStarts3(j):temp_exonEnds3(j)];
    end
    Mapping_Start_to_End{i} =  temp_txStart:temp_txEnd;%
    Mapping_Exonic{i} = temp_counts;
    Mapping_Intronic{i} = setdiff(temp_Span,temp_counts);
    clear temp_txStart temp_txEnd temp_cdsStart temp_cdsEnd temp_exonStarts temp_exonEnds
    clear temp_exonStarts2 temp_exonEnds2 temp_exonStarts3 temp_exonEnds3 temp_Span temp_counts
end
%For current Transcript get location of UTR
%If dealing with RNA they already removed introns

locTranscript_Exonic = unique(cat(1,Mapping_Exonic{:}));
locTranscript_Intronic = unique(cat(1,Mapping_Intronic{:})); 
locTranscript_StartToEnd = unique(cat(1,Mapping_Start_to_End{:})); 
if (strcmp(seqType,'RNA'))
    locTranscript = locTranscript_Exonic;
elseif (strcmp(seqType,'DNA'))
    locTranscript = locTranscript_StartToEnd;
end
for i = 1:size(probes,1)
    probes_location(i) = min(locTranscript(probes{i,3}:probes{i,3}+length(probes{i,2})-1)); 
end
% Sense 
%   0    ______  ___ ___ ___ end
%  end   ______  ___ ___ ___  0
% Anti-sense 
%Off-Target TS's    
if (strcmp(Organism,'Human'))
    S8 = readtable(settings.HumanGenomeAssemblyReportFile);
    AN_Chr = S8.RefSeq_Accn(1:639);%%chromosomes 1-22, X (23), Y(24), patchs/scaffolds 25-638, chrM MT (639)
elseif strcmp(Organism,'Mouse')  
    S16 = readtable(settings.MouseGenomeAssemblyReportFile); 
    AN_Chr = S16.RefSeq_Accn;%chromosomes 1-19, X (20), Y(21), patche scaffolds 22-60, chrM MT (61)
elseif strcmp(Organism,'Yeast')
    S21 = readtable(settings.YeastGenomeAssemblyReportFile); 
    AN_Chr = S21.RefSeq_Accn;%chromosome 1-16, chromosome M (Mitocondria) (17)
else
    S22 = readtable(settings.CustomGenomeAssemblyReportFile); 
    AN_Chr = S22.RefSeq_Accn;%chromosome 1-16, chromosome M (Mitocondria) (17)
end
runningSeperateTS = 1; 
for DNA_OffTargetIDs = cTypeIndx{1}
    Name_Target = AN_List{DNA_OffTargetIDs};
    ChrIdx = find(strcmp(Name_Target,AN_Chr));
    if (strcmp(Organism,'Human'))
        switch ChrIdx 
            case 23
               chr_num = 'chrX';
            case 24
               chr_num = 'chrY';
            case 639
               chr_num = 'chrM';
            otherwise
                chr_num = strcat('chr',num2str(ChrIdx));
        end
    elseif strcmp(Organism,'Mouse') 
        switch ChrIdx 
            case 20
               chr_num = 'chrX';
            case 21
               chr_num = 'chrY';
            case 61
               chr_num = 'chrM';
            otherwise
                chr_num = strcat('chr',num2str(ChrIdx));
        end
    elseif strcmp(Organism,'Yeast')
       switch ChrIdx 
            case 17
               chr_num = 'chrM';
            otherwise
                chr_num = strcat('chr',num2str(ChrIdx));
        end
    end
    gene_strand = EMBLComp_db.strand(find(strcmp(chr_num,EMBLComp_db.chrom)));
    TSS = EMBLComp_db.txStart(find(strcmp(chr_num,EMBLComp_db.chrom)));
    TTS = EMBLComp_db.txEnd(find(strcmp(chr_num,EMBLComp_db.chrom)));
    AssociatedTranscriptIDs = EMBLComp_db.transcriptId(find(strcmp(chr_num,EMBLComp_db.chrom)));
    AssociatedGeneIDs = EMBLAttrAlign_db.geneId(find(strcmp(chr_num,EMBLComp_db.chrom)));
    for i=1:length(onEMBL_transcriptID)
        Associated_OffTx{i} = find(cell2mat(arrayfun(@(x) find(strcmp(onEMBL_transcriptID{i},x)==0),AssociatedTranscriptIDs,'UniformOutput',false)));
    end
    Associated_Off = unique(cat(1,Associated_OffTx{:}));
    AssociatedTranscriptIDs_off = AssociatedTranscriptIDs(Associated_Off);
    AssociatedGeneIDs_off = AssociatedGeneIDs(Associated_Off);
    TSS = TSS(Associated_Off);
    TTS = TTS(Associated_Off);
    gene_strand_off = gene_strand(Associated_Off);
    Idx = find(strcmp(Event_Name,Name_Target));%'NC_000016.10'
    Chr_Event_Strand = Event_Strands(Idx);
    Chr_Event_LocStart = Event_LocStart(Idx);
    Chr_Event_ID = Event_ID(Idx); 
    for u=1:length(AssociatedTranscriptIDs_off)
       Length_Nascent = abs(TTS(u) - TSS(u)) + 1;%Get directionality
       if (strcmp(gene_strand_off{u},'+'))%sense or anti-sense binding
           Nascent_Direction = 1;
       elseif (strcmp(gene_strand_off{u},'-'))
           Nascent_Direction = -1;
       end
       try
       if (expCase==2) 
            for j=HumanSCTracks
                ExpValues = scExprGene_AN(find(strcmp(extractBefore(AssociatedGeneIDs_off{u},'.'),{scExprGene_AN.name}))).expScores;  
            end  
       elseif (expCase==1)
           if (HumanExpGeneOrTransc)
               ExpValues = ExprGene_AN(find(strcmp(extractBefore(AssociatedGeneIDs_off{u},'.'),{ExprGene_AN.geneId}))).expScores;  
           else
               ExpValues = ExprTransc_AN(find(strcmp(extractBefore(AssociatedTranscriptIDs_off{u},'.'),{ExprTransc_AN.name}))).expScores;
           end         
       elseif (expCase==0)
           if (HumanExpGeneOrTransc)
               ExpValues = ExprGene_AN(find(strcmp(extractBefore(AssociatedGeneIDs_off{u},'.'),{ExprGene_AN.name}))).expScores;  
           else     
               ExpValues = ExprTransc_AN(find(strcmp(extractBefore(AssociatedTranscriptIDs_off{u},'.'),{ExprTransc_AN.name}))).expScores;
           end
       end
       switch expValType
            case 3
                ExpOff = mean(ExpValues);
            case 4
                ExpOff = ExpValues(CellTypeID);
            otherwise
                ExpOff = mean(ExpValues);
       end   
       catch
           ExpOff = 0;
       end
       EquivalentNascentPolymerase = floor(NumberOfPolymerase*ExpOff/max(ExpOn));
       maxNascentPolymerase = floor(Length_Nascent/(PolymeraseOccupancy + PolymeraseSpacing));
       if (EquivalentNascentPolymerase>maxNascentPolymerase)
           EquivalentNascentPolymerase = maxNascentPolymerase;
       end
       if (EquivalentNascentPolymerase>0)
            S_running = 1;
            for n = 1:EquivalentNascentPolymerase
                BasePairRegionLength = floor(n*Length_Nascent/EquivalentNascentPolymerase);
                if (Nascent_Direction>0)%Sense
                    StartLoc = min([TSS(u) TTS(u)]);
                    Events_Sx1 = Chr_Event_ID(Chr_Event_LocStart>=StartLoc);
                    Events_Sx2 = Chr_Event_ID(Chr_Event_LocStart<=StartLoc + BasePairRegionLength); 
                    Events_Sx3 = Chr_Event_ID(strcmp(Chr_Event_Strand,'Plus/Plus')); 
                else%Antisense
                    StartLoc = max([TSS(u) TTS(u)]);
                    Events_Sx1 = Chr_Event_ID(Chr_Event_LocStart<=StartLoc);  
                    Events_Sx2 = Chr_Event_ID(Chr_Event_LocStart>=StartLoc - BasePairRegionLength);
                    Events_Sx3 = Chr_Event_ID(strcmp(Chr_Event_Strand,'Plus/Minus')); 
                end   
            Events_Sx = intersect(intersect(Events_Sx1,Events_Sx2),Events_Sx3);
            %Go from Probe and Event_ID to Probe and Sites
            if (~isempty(Events_Sx))
            Sites_Over_Events = find(cell2mat(cellfun(@(x) sum(ismember(x,Events_Sx)),Mol_ProbesAtEventsID{DNA_OffTargetIDs},'UniformOutput',false)));
            Sx = []; Px = [];
            for Site = Sites_Over_Events %get a list of Probes, and corresponding sites
                locSindx = find(ismember(Mol_ProbesAtEventsID{DNA_OffTargetIDs}{Site},Events_Sx));
                if (~isempty(Sx)&&~isempty(locSindx))
                    Sx = [Sx Site*ones(1,length(locSindx))];
                    Px = [Px MolProbesAtEvents{DNA_OffTargetIDs}{Site}(locSindx).'];
                elseif (isempty(Sx)&&~isempty(locSindx))
                    Sx(1:length(locSindx)) = Site*ones(1,length(locSindx));
                    Px(1:length(locSindx)) = MolProbesAtEvents{DNA_OffTargetIDs}{Site}(locSindx).';
                end 
                clear locSindx
            end  
            clear Events_Sx1 Events_Sx2 Events_Sx3 Events_Sx Sites_Over_Events
            Su = unique(Sx);
            Si = length(unique(Sx));
            if (Si>0)%If Probe Binds Sense or AntiSense since site map does not
                tempDoesProbeBindSite = squeeze(DoesProbeBindSite(:,DNA_OffTargetIDs,:));
                idx1 = sub2ind(size(tempDoesProbeBindSite),Px,Sx);
                newAssignment_DoesProbeBindSite = tempDoesProbeBindSite(idx1);
                if (size(DoesProbeBindSite_TS,2)<runningSeperateTS)
                    DoesProbeBindSite_TS(1,runningSeperateTS,1) = 0;
                end
                if (size(DoesProbeBindSite_TS,3)<max(new_Sx))
                    DoesProbeBindSite_TS(1,runningSeperateTS,max(new_Sx)) = 0;
                end
                DoesProbeBindSite_TS(sub2ind(size(DoesProbeBindSite_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx)) = newAssignment_DoesProbeBindSite;
                for m = 1:N_methods
                    tempKb = squeeze(Kb(:,DNA_OffTargetIDs,:,m));
                    tempHeq = squeeze(dHeq_mod(:,DNA_OffTargetIDs,:,m));
                    tempSeq = squeeze(dSeq_mod(:,DNA_OffTargetIDs,:,m));
                    tempCp = squeeze(dCp_mod(:,DNA_OffTargetIDs,:,m));
                    idx2 = sub2ind(size(tempKb),Px,Sx);
                    newAssignment_Kb = tempKb(idx2);
                    newAssignment_dHeq = tempHeq(idx2);
                    newAssignment_dSeq = tempSeq(idx2);
                    newAssignment_dCp = tempCp(idx2);
                    new_Px = Px;
                    new_Sx = cell2mat(arrayfun(@(x) find(x==Su)+S_running-1,Sx,'UniformOutput',false));
                    if (size(Kb_TS,2)<runningSeperateTS)
                        Kb_TS(1,runningSeperateTS,1,m) = 0;
                    end
                    if (size(dHeq_TS,2)<runningSeperateTS)
                        dHeq_TS(1,runningSeperateTS,1,m) = 0;
                    end
                    if (size(dSeq_TS,2)<runningSeperateTS)
                        dSeq_TS(1,runningSeperateTS,1,m) = 0;
                    end
                    if (size(dCp_TS,2)<runningSeperateTS)
                        dCp_TS(1,runningSeperateTS,1,m) = 0;
                    end
                    if (size(Kb_TS,3)<max(new_Sx))
                        Kb_TS(1,runningSeperateTS,max(new_Sx),m) = 0;
                    end
                    if (size(dHeq_TS,3)<max(new_Sx))
                        dHeq_TS(1,runningSeperateTS,max(new_Sx),m) = 0;
                    end
                    if (size(dSeq_TS,3)<max(new_Sx))
                        dSeq_TS(1,runningSeperateTS,max(new_Sx),m) = 0;
                    end
                    if (size(dCp_TS,3)<max(new_Sx))
                        dCp_TS(1,runningSeperateTS,max(new_Sx),m) = 0;
                    end
                    Kb_TS(sub2ind(size(Kb_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx,m*ones(1,length(new_Px)))) = newAssignment_Kb;     
                    dHeq_TS(sub2ind(size(dHeq_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx,m*ones(1,length(new_Px)))) = newAssignment_dHeq;  
                    dSeq_TS(sub2ind(size(dSeq_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx,m*ones(1,length(new_Px)))) = newAssignment_dSeq;    
                    dCp_TS(sub2ind(size(dCp_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx,m*ones(1,length(new_Px)))) = newAssignment_dCp;    
                end 
                for m = 1:N_methods2
                    tempHf = squeeze(dHf_mod(:,DNA_OffTargetIDs,:,m));
                    tempSf = squeeze(dSf_mod(:,DNA_OffTargetIDs,:,m));
                    tempHr = squeeze(dHr_mod(:,DNA_OffTargetIDs,:,m));
                    tempSr = squeeze(dSr_mod(:,DNA_OffTargetIDs,:,m));
                    idx2 = sub2ind(size(tempHf),Px,Sx);
                    newAssignment_Hf = tempHf(idx2);
                    newAssignment_Sf = tempSf(idx2);
                    newAssignment_Hr = tempHr(idx2);
                    newAssignment_Sr = tempSr(idx2); 
                    new_Px = Px;
                    new_Sx = cell2mat(arrayfun(@(x) find(x==Su)+S_running-1,Sx,'UniformOutput',false));
                    if (size(dHf_TS,2)<runningSeperateTS)
                        dHf_TS(1,runningSeperateTS,1,m) = 0;
                    end
                    if (size(dSf_TS,2)<runningSeperateTS)
                        dSf_TS(1,runningSeperateTS,1,m) = 0;
                    end
                    if (size(dHr_TS,2)<runningSeperateTS)
                        dHr_TS(1,runningSeperateTS,1,m) = 0;
                    end
                    if (size(dSr_TS,2)<runningSeperateTS)
                        dSr_TS(1,runningSeperateTS,1,m) = 0;
                    end
                    if (size(dHf_TS,3)<max(new_Sx))
                        dHf_TS(1,runningSeperateTS,max(new_Sx),m) = 0;
                    end
                    if (size(dSf_TS,3)<max(new_Sx))
                        dSf_TS(1,runningSeperateTS,max(new_Sx),m) = 0;
                    end
                    if (size(dHr_TS,3)<max(new_Sx))
                        dHr_TS(1,runningSeperateTS,max(new_Sx),m) = 0;
                    end
                    if (size(dSr_TS,3)<max(new_Sx))
                        dSr_TS(1,runningSeperateTS,max(new_Sx),m) = 0;
                    end
                    dHf_TS(sub2ind(size(dHf_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx,m*ones(1,length(new_Px)))) = newAssignment_Hf;     
                    dSf_TS(sub2ind(size(dSf_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx,m*ones(1,length(new_Px)))) = newAssignment_Sf;   
                    dHr_TS(sub2ind(size(dHr_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx,m*ones(1,length(new_Px)))) = newAssignment_Hr;   
                    dSr_TS(sub2ind(size(dSr_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx,m*ones(1,length(new_Px)))) = newAssignment_Sr;   
                end         
%                 for ck = 1:length(Px)    
%                     DoesProbeBindSite_TS(Px(ck),runningSeperateTS,find(Sx(ck)==Su)+S_running-1) = DoesProbeBindSite(Px(ck),DNA_OffTargetIDs,Sx(ck));
%                     Kb_TS(Px(ck),runningSeperateTS,find(Sx(ck)==Su)+S_running-1) = Kb(Px(ck),DNA_OffTargetIDs,Sx(ck));
%                 end
                S_running = S_running + Si;
                clear BasePairRegionLength Px Sx Si Su idx1 idx2 new_Px new_Sx
                clear newAssignment_Kb newAssignment_DoesProbeBindSite
                clear tempDoesProbeBindSite tempKb
            end
            end
            end 
            if (S_running>1)
                Num_of_Molecule_SitesTS(runningSeperateTS) = S_running;
                AN_ListTS{runningSeperateTS} = strcat('Nascent'," ",AssociatedTranscriptIDs_off{u});
                TS_SenseInfo(runningSeperateTS) = Nascent_Direction;
                NumPolymeraseOnTS(runningSeperateTS) = EquivalentNascentPolymerase;
                maxPolymeraseOnTS(runningSeperateTS) = maxNascentPolymerase;
                TS_ChrInfo{runningSeperateTS} = chr_num;
                runningSeperateTS = runningSeperateTS + 1;
            end
       end    
       clear Length_Nascent Nascent_Direction Associated_TranscriptIDs Associated_TranscriptIDs_off EquivalentNascentPolymerase
    end
end

OnTargetTS = zeros(1,runningSeperateTS);
for DNA_OnTargetIDs = Ix_DNA
    On_Events_ID = [];On_Events_Sites = [];On_Events_LocStart = [];
    for Sites = 1:length(Mol_ProbesAtEventsID{DNA_OnTargetIDs})
        if (~isempty(On_Events_ID))
            On_Events_LocStart = [On_Events_LocStart probes_location(MolProbesAtEvents{DNA_OnTargetIDs}{Sites}.')];
            On_Events_Probes = [On_Events_Probes MolProbesAtEvents{DNA_OnTargetIDs}{Sites}.'];
            On_Events_Sites = [On_Events_Sites Sites*ones(1,length(Mol_ProbesAtEventsID{DNA_OnTargetIDs}{Sites}))];
            On_Events_ID = [On_Events_ID length(On_Events_ID)+1:length(On_Events_ID)+length(Mol_ProbesAtEventsID{DNA_OnTargetIDs}{Sites})];
        else            
            On_Events_LocStart = probes_location(MolProbesAtEvents{DNA_OnTargetIDs}{Sites}.');
            On_Events_Probes = MolProbesAtEvents{DNA_OnTargetIDs}{Sites}.';
            On_Events_Sites = Sites*ones(1,length(Mol_ProbesAtEventsID{DNA_OnTargetIDs}{Sites}));
            On_Events_ID = 1:length(Mol_ProbesAtEventsID{DNA_OnTargetIDs}{Sites});    
        end                  
    end
    Name_Target = AN_List{DNA_OnTargetIDs};
    ChrIdx = find(strcmp(Name_Target,AN_Chr));
    if (strcmp(Organism,'Human'))
        switch ChrIdx 
            case 23
               chr_num = 'chrX';
            case 24
               chr_num = 'chrY';
            case 639
               chr_num = 'chrM';
            otherwise
                chr_num = strcat('chr',num2str(ChrIdx));
        end
    elseif strcmp(Organism,'Mouse') 
        switch ChrIdx 
            case 20
               chr_num = 'chrX';
            case 21
               chr_num = 'chrY';
            case 61
               chr_num = 'chrM';
            otherwise
                chr_num = strcat('chr',num2str(ChrIdx));
        end
    elseif strcmp(Organism,'Yeast')
       switch ChrIdx 
            case 17
               chr_num = 'chrM';
            otherwise
                chr_num = strcat('chr',num2str(ChrIdx));
        end
    end
    gene_strand = EMBLComp_db.strand{find(strcmp(chr_num,EMBLComp_db.chrom))};
    TSS = EMBLComp_db.txStart(find(strcmp(chr_num,EMBLComp_db.chrom)));
    TTS = EMBLComp_db.txEnd(find(strcmp(chr_num,EMBLComp_db.chrom)));
    AssociatedTranscriptIDs = EMBLComp_db.transcriptId(find(strcmp(chr_num,EMBLComp_db.chrom)));
    for i=1:length(onEMBL_transcriptID)
        Associated_OnTx{i} = find(cell2mat(arrayfun(@(x) find(strcmp(x,onEMBL_transcriptID{i})==1),AssociatedTranscriptIDs,'UniformOutput',false)));
    end
    try
        Associated_On = unique(cat(1,Associated_OnTx{:}));
    catch
        Associated_On = unique(cell2mat(Associated_OnTx{1}));
    end
    AssociatedTranscriptIDs_on = AssociatedTranscriptIDs(Associated_On);
    
    TSS = EMBLComp_db.txStart(find(strcmp(onEMBL_transcriptID,EMBLComp_db.transcriptId)));
    TTS = EMBLComp_db.txEnd(find(strcmp(onEMBL_transcriptID,EMBLComp_db.transcriptId)));
    gene_strand_on = EMBLComp_db.strand(find(strcmp(onEMBL_transcriptID,EMBLComp_db.transcriptId)));

    for u=1:length(onEMBL_transcriptID)
       Length_Nascent = abs(TTS(u) - TSS(u)) + 1;%Get directionality
       if (strcmp(gene_strand_on{u},'+'))%sense or anti-sense binding
           Nascent_Direction = 1;
       elseif (strcmp(gene_strand_on{u},'-'))
           Nascent_Direction = -1;
       end
       EquivalentNascentPolymerase = floor(NumberOfPolymerase*ExpOn(u)/max(ExpOn));   
       maxNascentPolymerase = floor(Length_Nascent/(PolymeraseOccupancy + PolymeraseSpacing));
       if (EquivalentNascentPolymerase>maxNascentPolymerase)
           EquivalentNascentPolymerase = maxNascentPolymerase;
       end
       if (EquivalentNascentPolymerase>0)
            S_running = 1;
            for n = 1:EquivalentNascentPolymerase
                BasePairRegionLength = floor(n*Length_Nascent/EquivalentNascentPolymerase);
                if (Nascent_Direction>0)%Sense   TTS-TSS
                    StartLoc = min([TSS(u) TTS(u)]);
                    Events_Sx1 = On_Events_ID(On_Events_LocStart>=StartLoc);
                    Events_Sx2 = On_Events_ID(On_Events_LocStart<=StartLoc + BasePairRegionLength); 
                else%Antisense
                    StartLoc = max([TSS(u) TTS(u)]);
                    Events_Sx1 = On_Events_ID(On_Events_LocStart<=StartLoc);  
                    Events_Sx2 = On_Events_ID(On_Events_LocStart>=StartLoc - BasePairRegionLength);
                end 
            Events_Sx = intersect(Events_Sx1,Events_Sx2); 
            Sx = On_Events_Sites(Events_Sx);
            Px = On_Events_Probes(Events_Sx);
            clear Events_Sx1 Events_Sx2 Events_Sx 
            Su = unique(Sx);
            Si = length(unique(Sx));
            if (Si>0)%If Probe Binds Sense or AntiSense since site map does not
                tempDoesProbeBindSite = squeeze(DoesProbeBindSite(:,DNA_OnTargetIDs,:));
                idx1 = sub2ind(size(tempDoesProbeBindSite),Px,Sx);
                newAssignment_DoesProbeBindSite = tempDoesProbeBindSite(idx1);
                if (size(DoesProbeBindSite_TS,2)<runningSeperateTS)
                    DoesProbeBindSite_TS(1,runningSeperateTS,1) = 0;
                end
                if (size(DoesProbeBindSite_TS,3)<max(new_Sx))
                    DoesProbeBindSite_TS(1,runningSeperateTS,max(new_Sx)) = 0;
                end
                DoesProbeBindSite_TS(sub2ind(size(DoesProbeBindSite_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx)) = newAssignment_DoesProbeBindSite;
                for m = 1:N_methods
                    tempKb = squeeze(Kb(:,DNA_OnTargetIDs,:,m));
                    tempHeq = squeeze(dHeq_mod(:,DNA_OnTargetIDs,:,m));
                    tempSeq = squeeze(dSeq_mod(:,DNA_OnTargetIDs,:,m));
                    tempCp = squeeze(dCp_mod(:,DNA_OnTargetIDs,:,m));
                    idx2 = sub2ind(size(tempKb),Px,Sx);
                    newAssignment_Kb = tempKb(idx2);
                    newAssignment_dHeq = tempHeq(idx2);
                    newAssignment_dSeq = tempSeq(idx2);
                    newAssignment_Cp = tempCp(idx2);
                    new_Px = Px;
                    new_Sx = cell2mat(arrayfun(@(x) find(x==Su)+S_running-1,Sx,'UniformOutput',false));
                    if (size(Kb_TS,2)<runningSeperateTS)
                        Kb_TS(1,runningSeperateTS,1,m) = 0;
                    end
                    if (size(dHeq_TS,2)<runningSeperateTS)
                        dHeq_TS(1,runningSeperateTS,1,m) = 0;
                    end
                    if (size(dSeq_TS,2)<runningSeperateTS)
                        dSeq_TS(1,runningSeperateTS,1,m) = 0;
                    end
                    if (size(dCp_TS,2)<runningSeperateTS)
                        dCp_TS(1,runningSeperateTS,1,m) = 0;
                    end
                    if (size(Kb_TS,3)<max(new_Sx))
                        Kb_TS(1,runningSeperateTS,max(new_Sx),m) = 0;
                    end
                    if (size(dHeq_TS,3)<max(new_Sx))
                        dHeq_TS(1,runningSeperateTS,max(new_Sx),m) = 0;
                    end
                    if (size(dSeq_TS,3)<max(new_Sx))
                        dSeq_TS(1,runningSeperateTS,max(new_Sx),m) = 0;
                    end
                    if (size(dCp_TS,3)<max(new_Sx))
                        dCp_TS(1,runningSeperateTS,max(new_Sx),m) = 0;
                    end
                    Kb_TS(sub2ind(size(Kb_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx,m*ones(1,length(new_Px)))) = newAssignment_Kb;     
                    dHeq_TS(sub2ind(size(dHeq_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx,m*ones(1,length(new_Px)))) = newAssignment_dHeq;  
                    dSeq_TS(sub2ind(size(dSeq_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx,m*ones(1,length(new_Px)))) = newAssignment_dSeq;    
                    dCp_TS(sub2ind(size(dCp_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx,m*ones(1,length(new_Px)))) = newAssignment_dCp;   
                end 
                for m = 1:N_methods2
                    tempHf = squeeze(dHf_mod(:,DNA_OnTargetIDs,:,m));
                    tempSf = squeeze(dSf_mod(:,DNA_OnTargetIDs,:,m));
                    tempHr = squeeze(dHr_mod(:,DNA_OnTargetIDs,:,m));
                    tempSr = squeeze(dSr_mod(:,DNA_OnTargetIDs,:,m));
                    idx2 = sub2ind(size(tempHf),Px,Sx);
                    newAssignment_Hf = tempHf(idx2);
                    newAssignment_Sf = tempSf(idx2);
                    newAssignment_Hr = tempHr(idx2);
                    newAssignment_Sr = tempSr(idx2); 
                    new_Px = Px;
                    new_Sx = cell2mat(arrayfun(@(x) find(x==Su)+S_running-1,Sx,'UniformOutput',false));
                    if (size(dHf_TS,2)<runningSeperateTS)
                        dHf_TS(1,runningSeperateTS,1,m) = 0;
                    end
                    if (size(dSf_TS,2)<runningSeperateTS)
                        dSf_TS(1,runningSeperateTS,1,m) = 0;
                    end
                    if (size(dHr_TS,2)<runningSeperateTS)
                        dHr_TS(1,runningSeperateTS,1,m) = 0;
                    end
                    if (size(dSr_TS,2)<runningSeperateTS)
                        dSr_TS(1,runningSeperateTS,1,m) = 0;
                    end
                    if (size(dHf_TS,3)<max(new_Sx))
                        dHf_TS(1,runningSeperateTS,max(new_Sx),m) = 0;
                    end
                    if (size(dSf_TS,3)<max(new_Sx))
                        dSf_TS(1,runningSeperateTS,max(new_Sx),m) = 0;
                    end
                    if (size(dHr_TS,3)<max(new_Sx))
                        dHr_TS(1,runningSeperateTS,max(new_Sx),m) = 0;
                    end
                    if (size(dSr_TS,3)<max(new_Sx))
                        dSr_TS(1,runningSeperateTS,max(new_Sx),m) = 0;
                    end
                    dHf_TS(sub2ind(size(dHf_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx,m*ones(1,length(new_Px)))) = newAssignment_Hf;     
                    dSf_TS(sub2ind(size(dSf_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx,m*ones(1,length(new_Px)))) = newAssignment_Sf;   
                    dHr_TS(sub2ind(size(dHr_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx,m*ones(1,length(new_Px)))) = newAssignment_Hr;   
                    dSr_TS(sub2ind(size(dSr_TS),new_Px,...
                        runningSeperateTS*ones(1,length(new_Px)),new_Sx,m*ones(1,length(new_Px)))) = newAssignment_Sr;   
                end  
%                 for ck = 1:length(Px)    
%                     DoesProbeBindSite_TS(Px(ck),runningSeperateTS,find(Sx(ck)==Su)+S_running-1) = DoesProbeBindSite(Px(ck),DNA_OffTargetIDs,Sx(ck));
%                     Kb_TS(Px(ck),runningSeperateTS,find(Sx(ck)==Su)+S_running-1) = Kb(Px(ck),DNA_OffTargetIDs,Sx(ck));
%                 end
                S_running = S_running + Si;
                clear BasePairRegionLength Px Sx Si Su idx1 idx2 new_Px new_Sx
                clear newAssignment_Kb newAssignment_DoesProbeBindSite
                clear tempDoesProbeBindSite tempKb
            end
            end  
            if (S_running>1)
                Num_of_Molecule_SitesTS(runningSeperateTS) = S_running;
                AN_ListTS{runningSeperateTS} = strcat('Nascent'," ",AssociatedTranscriptIDs_on{u});
                OnTargetTS(runningSeperateTS) = 1; 
                TS_SenseInfo(runningSeperateTS) = Nascent_Direction;
                NumPolymeraseOnTS(runningSeperateTS) = EquivalentNascentPolymerase;
                maxPolymeraseOnTS(runningSeperateTS) = maxNascentPolymerase;
                TS_ChrInfo{runningSeperateTS} = chr_num;
                runningSeperateTS = runningSeperateTS + 1;
            end
       end
       clear Length_Nascent Nascent_Direction Associated_TranscriptIDs Associated_TranscriptIDs_off EquivalentNascentPolymerase
    end
end
end


