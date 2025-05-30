function [IsSiteInLoop,Ks_TjLi_eq,dHs_TjLi_eq,dSs_TjLi_eq,dGs_TjLi_eq,dHs_TjLi_f,dSs_TjLi_f,dGs_TjLi_f,dHs_TjLi_r,dSs_TjLi_r,dGs_TjLi_r,dCp_TjLi_eq] = A_JH_GetRNASecondaryStructures_V2(settings,gene_table,DoesProbeBindSite,probes)
%Jason Hughes code to get/parse secondary structure of RNA in genome parsing 
%Into binding sites that probes could have blocked with how kinetic rates
%of individual loop formation for loops over 15bp long.
N_methods = 8;
N_methods2 = 3;
kb = 0.001987204259;%boltzman constant 
solveSeq = settings.SolveStructure;
RNASecondaryStructureRoot = settings.SecondaryStructureFileRoot;
clusterStatus = settings.clusterStatus;
RemoveMisMatches = settings.RemoveMisMatches;
SaltConcentration = settings.SaltConcentration;
T_hybrid = settings.HybridizationTemperature;
gene_table = sortrows(gene_table,[7 6],'ascend');
gene_table2 = gene_table(gene_table.Match>=15,:);
MinusStrandedHits = find(contains(gene_table2.Strand,'Minus'));
RNA_IDs_1 = find(contains(gene_table2.Name,'NM_'));
RNA_IDs_2 = find(contains(gene_table2.Name,'NR_'));
RNA_IDs_3 = find(contains(gene_table2.Name,'XM_'));
RNA_IDs_4 = find(contains(gene_table2.Name,'XR_'));
contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);
gene_table2 = gene_table2(setdiff(1:size(gene_table2,1),RNA_MissedFilteredHits),:);
gene_table2.Ax = min(gene_table2.SubjectIndices,[],2);
gene_table2.Bx = max(gene_table2.SubjectIndices,[],2);
gene_table3 = sortrows(gene_table2,[7 13],'ascend');
Names = unique(gene_table3.Name);
Names = convertCharsToStrings(Names);
uniNames = extractBefore(Names,'.');
NonDNA_IDs_1 = find(~contains(uniNames,'NC_'));%IDs
NonDNA_IDs_2 = find(~contains(uniNames,'NT_'));%IDs
NonDNA_IDs_3 = find(~contains(uniNames,'NW_'));%
NonDNA_IDs = intersect(intersect(NonDNA_IDs_1,NonDNA_IDs_2),NonDNA_IDs_3);
if strcmp(Organism,'Mouse') 
    switch clusterStatus
        case 0
            SEQdbRoot = settings.mLocalLocRoot;
        case 1
            SEQdbRoot = settings.mClusterLocRoot;
    end
elseif strcmp(Organism, 'Human')
    switch clusterStatus
        case 0
            SEQdbRoot = settings.hLocalLocRoot;
        case 1
            SEQdbRoot = settings.hClusterLocRoot;
    end
elseif strcmp(Organism, 'Yeast')
    switch clusterStatus
        case 0
            SEQdbRoot = settings.yLocalLocRoot;
        case 1
            SEQdbRoot = settings.yClusterLocRoot;
    end
else
    switch clusterStatus
       case 0
            SEQdbRoot = settings.oLocalLocRoot;
       case 1
            SEQdbRoot = settings.oClusterLocRoot;
    end
end
FNAFiles_FNA = dir([SEQdbRoot '**/*.fna']);
FNAFiles_FA = dir([SEQdbRoot '**/*.fa']);
if (~isempty(FNAFiles_FNA)&&~isempty(FNAFiles_FA))
    FNAFiles = cat(1,FNAFilesA,FNAFiles_FA);
elseif (~isempty(FNAFiles_FNA))
    FNAFiles = FNAFiles_FNA;
elseif (~isempty(FNAFiles_FA))
    FNAFiles = FNAFiles_FA;
end
for k=1:length(FNAFiles)
    FNAFileLoc{k} = strcat(FNAFiles(k).folder,'/',FNAFiles(k).name);
    BioIFobj{k} = BioIndexedFile('FASTA',FNAFileLoc{k});
    BioKeys{k}{:} = getKeys(BioIFobj{k});
end  
SeqBD = fastaread(FNAFileLoc{k});
SeqBD_Names = {SeqBD.Header};
RSSFiles = dir([RNASecondaryStructureRoot '**/*.dbn']);
RSSFiles_RFAM = find(contains({RSSFiles.name},'RFAM'));
RSSFileLoc = cell(length(RSSFiles_RFAM),1);
RSS_Name = cell(length(RSSFiles_RFAM),1);
RSS_Length = cell(length(RSSFiles_RFAM),1);
RSS_Seq = cell(length(RSSFiles_RFAM),1);
RSS_Bracket = cell(length(RSSFiles_RFAM),1);
for v=1:length(RSSFiles_RFAM)
    k = RSSFiles_RFAM(v);
    RSSFileLoc{k} = strcat(RSSFiles(k).folder,'/',RSSFiles(k).name);
    temp_struct = fileread(RSSFileLoc{k});
    temp_parsed = splitlines(temp_struct);
    RSS_Name{v} = temp_parsed{1};
    RSS_Length{v} = temp_parsed{2};
    RSS_Seq{v} =strrep(temp_parsed{4},'U','T');
    RSS_Bracket{v} = temp_parsed{5}(2:end-1);
end
Ks_TjLi = ndSparse.build([length(Names),1],0);
IsSiteInLoop = ndSparse.build([length(Names),1,size(probes,1)],0);
Lp_seq = cell(length(Names),1);
Lq_seq = cell(length(Names),1);
dCp_TjLi_eq = ndSparse.build([size(Names,1),1,N_methods],0);
dHs_TjLi_eq = ndSparse.build([length(Names),1,N_methods],0);
dSs_TjLi_eq = ndSparse.build([length(Names),1,N_methods],0);
dGs_TjLi_eq = ndSparse.build([length(Names),1,N_methods],0);
Ks_TjLi_eq = ndSparse.build([length(Names),1,N_methods],0);
dHs_TjLi_f = ndSparse.build([length(Names),1,N_methods2],0);
dSs_TjLi_f = ndSparse.build([length(Names),1,N_methods2],0);
dGs_TjLi_f = ndSparse.build([length(Names),1,N_methods2],0);
dHs_TjLi_r = ndSparse.build([length(Names),1,N_methods2],0);
dSs_TjLi_r = ndSparse.build([length(Names),1,N_methods2],0);
dGs_TjLi_r = ndSparse.build([length(Names),1,N_methods2],0);

           
for u=1:length(NonDNA_IDs)
    t = NonDNA_IDs(u);
    try
    tempName = Names(t);
    SeqBD_ID = find(strcmp(SeqBD_Names,tempName));
    if (~isempty(SeqBD_ID))
        checkSequence = SeqBD(strcmp(SeqBD_Names,tempName)).Sequence;
    else
        found = 0;
    end
    catch
        sdsd = 1;
    end
    if (solveSeq)
        found = 1;
        try
            seq1 = checkSequence;
            [RNAbrack,Energy] = rnafold(checkSequence);
        catch
            found = 0;
        end
    else
        found = find(strcmp(RSS_Seq,checkSequence)); 
        if (~isempty(found))
            seq1 = RSS_Seq{found};
            RNAbrack = RSS_Bracket{found};
        end
    end
    if (~isempty(found))  
        RNAnot = rnaconvert(RNAbrack);
        [rw] = find(RNAnot);
        %generate basepairs in loops
        [row,col] = ind2sub(size(RNAnot),rw);
        row = flip(row);col = flip(col);
        slope = diff(row);
        dir = abs(slope);
        q = [0 find(dir>1).'];
        Lp = cell(1,length(q));Lq =  cell(1,length(q));
        for v=1:length(q)
            if (v<length(q))
                Lp{v} = row(q(v)+1:q(v+1)).';
                Lq{v} = col(q(v)+1:q(v+1)).';
            else
                Lp{v} = row(q(v)+1:end).';   
                Lq{v} = col(q(v)+1:end).';
            end
        end
        temp_Events = find(strcmp(gene_table3.Name,tempName));
        temp_Ax = gene_table3.Ax(temp_Events);
        temp_Bx = gene_table3.Bx(temp_Events);
        temp_Probes = gene_table3.ProbeNum(temp_Events);
        %temp_Ranges = cell(1,length(temp_Ax));
        for m = 1:length(temp_Ax)
            temp_Ranges{m} = temp_Ax(m):temp_Bx(m);   
        end
        Lset = cell(1,length(q));
        for v = 1:length(q)
            Lset{v} = union(Lp{v},Lq{v});
            if (length(Lset{v})>=30)
            %IsEventInSite = cell(1,length(temp_Ax));
            IsEventInSite = [];
            for m = 1:length(temp_Ax)
               IsEventInSite(m) = double(sum(ismember(temp_Ranges{m},Lset{v}))>0);
            end
            Events_Wanted = find(IsEventInSite);
            if (~isempty(Events_Wanted))
               Probes_Wanted = temp_Probes(Events_Wanted);
               for p = Probes_Wanted
                    Potential_Sx = find(full(DoesProbeBindSite(p,t,:)));
                    %Mol_ProbesAtEventsID;%Figure out which potential Sites
                    if (~isempty(Potential_Sx))
                        IsSiteInLoop(t,v,Potential_Sx) = 1;
                    end
               end             
            end
            else
            end
            Lp_seq{t,v} = strrep(seq1(Lp{v}),'U','T');
            Lq_seq{t,v} = strrep(seq1(Lq{v}),'U','T');
            try
                [temp_dHeq, temp_dSeq, temp_dGeq, ...
                temp_dHf, temp_dSf, temp_dGf, ...
                temp_dHr, temp_dSr, temp_dGr,temp_dCpeq,~] = ...
                F_DeltaGibson_V3(strrep(Lp_seq{t,v},'-','N'),strrep(Lq_seq{t,v},'-','N'),SaltConcentration,T_hybrid);
                dCp_TjLi_eq(t,v,1:N_methods) = temp_dCpeq;    
                dHs_TjLi_eq(t,v,1:N_methods) = temp_dHeq;
                dSs_TjLi_eq(t,v,1:N_methods) = temp_dSeq;
                dGs_TjLi_eq(t,v,1:N_methods) = temp_dGeq;
                Ks_TjLi_eq(t,v,1:N_methods) = exp(-temp_dGeq/(kb*(T_hybrid+273.15)));
                dHs_TjLi_f(t,v,1:N_methods2) = temp_dHf;
                dSs_TjLi_f(t,v,1:N_methods2) = temp_dSf; 
                dGs_TjLi_f(t,v,1:N_methods2) = temp_dGf;
                dHs_TjLi_r(t,v,1:N_methods2) = temp_dHr;
                dSs_TjLi_r(t,v,1:N_methods2) = temp_dSr;
                dGs_TjLi_r(t,v,1:N_methods2) = temp_dGr;
            catch
                dCp_TjLi_eq(t,v,1:N_methods) = 0;    
                dHs_TjLi_eq(t,v,1:N_methods) = Inf;
                dSs_TjLi_eq(t,v,1:N_methods) = Inf;
                Ks_TjLi_eq(t,v,1:N_methods) = 0;
                dHs_TjLi_f(t,v,1:N_methods2) = Inf;
                dSs_TjLi_f(t,v,1:N_methods2) = Inf; 
                dGs_TjLi_f(t,v,1:N_methods2) = Inf;
                dHs_TjLi_r(t,v,1:N_methods2) = Inf;
                dSs_TjLi_r(t,v,1:N_methods2) = Inf;
                dGs_TjLi_r(t,v,1:N_methods2) = Inf;
            end   
        end
    end
end 
end