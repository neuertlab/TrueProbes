function [Kb_Match,Kon,Koff,dHeq_Match,dSeq_Match,dHf_Match,dSf_Match,dHr_Match,dSr_Match,dCpeq_Match,dHon_eq,dSon_eq,dHon_f,dSon_f,dHon_r,dSon_r,dCp_eq,Tm_on,Tm_Match] = A_JH_GenerateThermoInfo_V4(probes,gene_table,settings)
% This Function takes probe sequences for on/off-target duplexes in in blast table 
% and computes binding energy and rates associated with off-target and on-targets

N_methods = 8;
N_methods2 = 3;
kb = 0.001987204259;%boltzman constant
T_hybrid = settings.HybridizationTemperature;
SaltConcentration = settings.SaltConcentration;
RemoveMisMatches = settings.RemoveMisMatches;

gene_name1 = settings.GeneName;
transcript_ID = settings.transcript_IDs;
if (iscell(transcript_ID))
    transcript_IDon = cell(1,length(transcript_ID));
    for v = 1:length(transcript_ID)
       transcript_IDon{v} = transcript_ID{v};
    end
else 
    transcript_IDon{1} = transcript_ID; 
end
ChrNum = settings.ChrNum;
GeneChr = settings.GeneChr;
Organism = settings.Organism;
%% Identify location of on-target DNA and RNA molecules in list of DNA/RNA molecules
%generate ID vector to identify molecules that count as on-target hits
if (strcmp(Organism,'Human'))
    S8 = readtable(settings.HumanGenomeAssemblyReportFile);
    AN_Chr = S8.RefSeq_Accn(1:639);%%chromosomes 1-22, X (23), Y(24), patchs/scaffolds 25-638, chrM MT (639)
    if (strcmpi(ChrNum,'X'))
        AN_ON_Chr = AN_Chr{23};
    elseif (strcmpi(ChrNum,'Y'))
        AN_ON_Chr = AN_Chr{24};
    elseif (strcmpi(ChrNum,'M'))
        AN_ON_Chr = 'NC_012920.1';
    else
        AN_ON_Chr = AN_Chr{str2double(ChrNum)};
    end
elseif strcmp(Organism,'Mouse')  
    S16 = readtable(settings.MouseGenomeAssemblyReportFile); 
    AN_Chr = S16.RefSeq_Accn;%chromosomes 1-19, X (20), Y(21), patche scaffolds 22-60, chrM MT (61)
    if (strcmpi(ChrNum,'X'))
        AN_ON_Chr = AN_Chr{20};
    elseif (strcmpi(ChrNum,'Y'))
        AN_ON_Chr = AN_Chr{21};
    elseif (strcmpi(GeneChr,'MT'))
        AN_ON_Chr = 'NC_005089.1';
    else
        AN_ON_Chr = AN_Chr{str2double(ChrNum)};
    end
elseif strcmp(Organism,'Yeast')
    S21 = readtable(settings.YeastGenomeAssemblyReportFile); 
    AN_Chr = S21.RefSeq_Accn;%chromosome 1-16, chromosome M (Mitocondria) (17)
    if (strcmpi(GeneChr,'MT'))
        AN_ON_Chr = AN_Chr{17};
    else
        AN_ON_Chr = AN_Chr{roman2num(GeneChr)};
    end
else
    S22 = readtable(settings.CustomGenomeAssemblyReportFile); 
    AN_Chr = S22.RefSeq_Accn;%chromosome 1-16, chromosome M (Mitocondria) (17)
    if (strcmpi(ChrNum,'X'))
        AN_ON_Chr = AN_Chr{settings.Custom_X_ChromNumber};
    elseif (strcmpi(ChrNum,'Y'))
        AN_ON_Chr = AN_Chr{settings.Custom_Y_ChromNumber};
    elseif (strcmpi(GeneChr,'MT'))
        AN_ON_Chr = AN_Chr{settings.Custom_MT_ChromNumber};
    else
        AN_ON_Chr = AN_Chr{str2double(ChrNum)};
    end
end
AN_ON_Chr = extractBefore(AN_ON_Chr,'.');
transcript_IDon = extractBefore(transcript_IDon,'.');

RV = @(x) (seqrcomplement(x));
pi_seq = cell(1,size(probes,1));
for i=1:size(probes,1)
   pi_seq{i} = RV(probes{i,2});
end

% Extract Gene Taret Names
gene_table = sortrows(gene_table,[7 6],'ascend');
gene_table2 = gene_table(gene_table.Match>=settings.MinHomologySearchTargetSize,:);
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
clear gene_table gene_table2
Names = unique(gene_table3.Name);
Names = convertCharsToStrings(Names);
uniNames = extractBefore(Names,'.');

% Check if output file already exists
try
    if (settings.clusterStatus == 1)
        task_ID = getenv('SLURM_ARRAY_TASK_ID');
        if (isempty(task_ID))
            task_ID = 0;
        else
            task_ID = str2num(task_ID);
        end
        sleepRate = sum(1:task_ID);
        pause(5*sleepRate)
    end
    load([settings.FolderRootName '/' gene_name1 '_Tm' num2str(T_hybrid) '_RM' num2str(RemoveMisMatches) '_OnOffThermoInfo' settings.designerName '.mat'],'Kon','Koff','Kb_Match')    
    calcOnOff = 0;
catch ME
    disp(ME.message)
    calcOnOff = 1;
end        
if (calcOnOff)
    %gets sequences of probe and target in blast records
    targetMatch = arrayfun(@(x) strrep(gene_table3.Alignment{x}(3,:),'-','N'),1:size(gene_table3,1),'UniformOutput',false);%slow not so slow
    probeMatch = arrayfun(@(x) seqrcomplement(strrep(gene_table3.Alignment{x}(1,:),'-','N')),1:size(gene_table3,1),'UniformOutput',false);%slow 
    %Find Unique Pairs  (Do Mapping and Reverse Mapping)
    CalcDictionary = unique([targetMatch(:)' probeMatch(:)']);
    CalcDict2 = convertCharsToStrings(CalcDictionary);%1xN string
    indices = 1:length(CalcDict2);
    idMap = containers.Map(CalcDict2,indices);
    CalcPairs(:,1) = cell2mat(arrayfun(@(x) idMap(targetMatch{x}),1:size(gene_table3,1),'UniformOutput',false));
    CalcPairs(:,2) = cell2mat(arrayfun(@(x) idMap(probeMatch{x}),1:size(gene_table3,1),'UniformOutput',false));
    NamesZ = gene_table3.Name;
    ProbeNumZ = gene_table3.ProbeNum;
    MatchZ = gene_table3.Match;
    clear gene_table3
    uniqueCalcPairs = unique(CalcPairs,'rows');
    targetMatch_Unique = {CalcDictionary{uniqueCalcPairs(:,1)}};
    probeMatch_Unique = {CalcDictionary{uniqueCalcPairs(:,2)}};
    clear CalcDictionary CalcDict2 idMap indices
    uniquePair_BindingRate = zeros(size(uniqueCalcPairs,1),N_methods);
    %Re-Map Energies to non-unique sites
    %Unique_to_Original = @(x) find(strcmpi(targetMatch_Unique{x},targetMatch).*strcmpi(probeMatch_Unique{x},probeMatch)); 
    Kon = zeros(size(probes,1),N_methods);
    Tm_on = zeros(size(probes,1),N_methods+1);
    dHeq_Match = zeros(length(targetMatch),N_methods); 
    dSeq_Match = zeros(length(targetMatch),N_methods);
    dCpeq_Match = zeros(length(targetMatch),N_methods); 
    dHf_Match = zeros(length(targetMatch),3); 
    dSf_Match = zeros(length(targetMatch),3); 
    dHr_Match = zeros(length(targetMatch),3); 
    dSr_Match = zeros(length(targetMatch),3); 
    dHon_eq = zeros(size(probes,1),N_methods);
    dSon_eq = zeros(size(probes,1),N_methods);
    dCp_eq = zeros(size(probes,1),N_methods);
    dHon_f = zeros(size(probes,1),3);
    dSon_f = zeros(size(probes,1),3);
    dHon_r = zeros(size(probes,1),3);
    dSon_r = zeros(size(probes,1),3);
    uniquePair_Tm = zeros(size(probes,1),N_methods+1);
    uniquePair_dHeq = zeros(size(probes,1),N_methods);
    uniquePair_dSeq = zeros(size(probes,1),N_methods);
    uniquePair_dCpeq = zeros(size(uniqueCalcPairs,1),N_methods);
    uniquePair_dHf = zeros(size(uniqueCalcPairs,1),3);
    uniquePair_dSf = zeros(size(uniqueCalcPairs,1),3);
    uniquePair_dHr = zeros(size(uniqueCalcPairs,1),3);
    uniquePair_dSr = zeros(size(uniqueCalcPairs,1),3);
    Tm_Match = zeros(length(targetMatch),N_methods+1); 
    Kb_Match = zeros(length(targetMatch),N_methods); 
    try
        n_cores = str2num(getenv('SLURM_JOB_CPUS_PER_NODE')); 
        p = parpool(n_cores);
        %p.IdleTimeout = Inf;
        spmd
            warning('off','all')
        end
    catch
        %p = parpool;
        %p.IdleTimeout = Inf;
    end
    %compute binding thermodynamic/kinetic information for on-target hits parallelized
    parfor i=1:size(probes,1)
        [temp_dHeq{i}, temp_dSeq{i}, temp_dGeq{i}, ...
         temp_dHf{i}, temp_dSf{i}, ~, ...
         temp_dHr{i}, temp_dSr{i}, ~,temp_dCpeq{i},temp_Tm{i}] = ...
             F_DeltaGibson_V3(lower(pi_seq{i}),seqrcomplement(lower(pi_seq{i})),SaltConcentration,T_hybrid); 
        dHon_eq(i,:) = temp_dHeq{i};
        dSon_eq(i,:) = temp_dSeq{i};
        dHon_f(i,:) = temp_dHf{i};
        dSon_f(i,:) = temp_dSf{i};
        dHon_r(i,:) = temp_dHr{i};
        dSon_r(i,:) = temp_dSr{i};
        dCp_eq(i,:) = temp_dCpeq{i};
        Tm_on(i,:) = temp_Tm{i};
        Kon(i,:) = exp(-temp_dGeq{i}/(kb*(T_hybrid+273.15))); %Keq rates on-Target binding      
    end
    %compute binding thermodynamic/kinetic information for unique off-target hit sequences parallelized
    parfor v = 1:size(uniqueCalcPairs,1)
       [uniquePair_temp_dHeq{v}, uniquePair_temp_dSeq{v}, uniquePair_temp_dGeq{v}, ...
        uniquePair_temp_dHf{v}, uniquePair_temp_dSf{v}, ~, ...
        uniquePair_temp_dHr{v}, uniquePair_temp_dSr{v}, ~,uniquePair_temp_dCpeq{v},uniquePair_temp_Tm{v}] = ... 
           F_DeltaGibson_V3(targetMatch_Unique{v},probeMatch_Unique{v},SaltConcentration,T_hybrid);
       uniquePair_dHeq(v,:) = uniquePair_temp_dHeq{v}';
       uniquePair_dSeq(v,:) = uniquePair_temp_dSeq{v}';
       uniquePair_dHf(v,:) = uniquePair_temp_dHf{v}';
       uniquePair_dSf(v,:) = uniquePair_temp_dSf{v}';
       uniquePair_dHr(v,:) = uniquePair_temp_dHr{v}';
       uniquePair_dSr(v,:) = uniquePair_temp_dSr{v}';
       uniquePair_dCpeq(v,:) = uniquePair_temp_dCpeq{v}';
       uniquePair_Tm(v,:) = uniquePair_temp_Tm{v}';
       uniquePair_BindingRate(v,:) = exp(-uniquePair_temp_dGeq{v}/(kb*(T_hybrid+273.15)));
    end 
    parfor v = 1:size(uniqueCalcPairs,1)
       temp_indx{v}= find(strcmpi(targetMatch_Unique{v},targetMatch).*strcmpi(probeMatch_Unique{v},probeMatch));
       temp_val{v} = uniquePair_BindingRate(v,:);
       temp_valHeq{v} = uniquePair_dHeq(v,:);
       temp_valSeq{v} = uniquePair_dSeq(v,:);
       temp_valHf{v} = uniquePair_dHf(v,:);
       temp_valSf{v} = uniquePair_dSf(v,:);
       temp_valHr{v} = uniquePair_dHr(v,:);
       temp_valSr{v} = uniquePair_dSr(v,:);
       temp_valTm{v} = uniquePair_Tm(v,:);
       temp_valCpeq{v} = uniquePair_dCpeq(v,:); 
    end 
	
	%store results for on/off-target binding information matches to unique sequence pairs
    for v = 1:size(uniqueCalcPairs,1)
       Kb_Match(temp_indx{v},:) = repmat(temp_val{v},[length(temp_indx{v}) 1]);
       dHeq_Match(temp_indx{v},:) = repmat(temp_valHeq{v},[length(temp_indx{v}) 1]); 
       dSeq_Match(temp_indx{v},:) = repmat(temp_valSeq{v},[length(temp_indx{v}) 1]);
       dHf_Match(temp_indx{v},:) = repmat(temp_valHf{v},[length(temp_indx{v}) 1]);
       dSf_Match(temp_indx{v},:) = repmat(temp_valSf{v},[length(temp_indx{v}) 1]);
       dHr_Match(temp_indx{v},:) = repmat(temp_valHr{v},[length(temp_indx{v}) 1]);
       dSr_Match(temp_indx{v},:) = repmat(temp_valSr{v},[length(temp_indx{v}) 1]);
       dCpeq_Match(temp_indx{v},:) = repmat(temp_valCpeq{v},[length(temp_indx{v}) 1]);
       Tm_Match(temp_indx{v},:) = repmat(temp_valTm{v},[length(temp_indx{v}) 1]); 
    end   
        
Koff = ndSparse.build([size(probes,1),length(Names),N_methods],0);
indices = 1:length(uniNames);
idMap = containers.Map(uniNames,indices);
L = size(probes,1);
uniNamesZ = extractBefore(NamesZ,'.');  
%Might be better efficent to combine old Event way of getting Koff with new
%part for calculating energies minimally with deal

%build structure with on/off-target binding by each target in blast results.
for p=1:L
    try
    IDy = find(ProbeNumZ==p);%add step to remove Off-target from chromosome
    IDs = find(MatchZ<20);
    IDf = find(MatchZ==20);
    Indx_OFF = intersect(IDy,IDs);
    Indx_ON = intersect(IDy,IDf);%contains off hits
    Names_OFF = unique(extractBefore(NamesZ(Indx_OFF),'.'));
    Names_ON = unique(extractBefore(NamesZ(Indx_ON),'.'));
    Names_p = union(Names_OFF,Names_ON);
    OFF_uniNamesZ =  uniNamesZ(Indx_OFF);
    ON_uniNamesZ =  uniNamesZ(Indx_ON);
    IDon_temp = cell2mat(arrayfun(@(x) idMap(Names_ON{x}),1:size(Names_ON,1),'UniformOutput',false));
    IDoff_temp = cell2mat(arrayfun(@(x) idMap(Names_OFF{x}),1:size(Names_OFF,1),'UniformOutput',false));  
    actual_IDon_DNA = cell2mat(arrayfun(@(x) strcmp(Names_ON{x},AN_ON_Chr),1:size(Names_ON,1),'UniformOutput',false));
    actual_IDon_RNA = cell2mat(arrayfun(@(x) strcmp(Names_ON{x},transcript_IDon),1:size(Names_ON,1),'UniformOutput',false));
    actual_IDon = find(actual_IDon_DNA+actual_IDon_RNA>0);
    false_IDon = find(actual_IDon_DNA+actual_IDon_RNA==0);
    IDon = IDon_temp(actual_IDon);
    IDoff = union(IDoff_temp,IDon_temp(false_IDon)); 
    IDshared = intersect(IDon,IDoff);
    SharedName = {uniNames{IDshared}};
    falseNames = Names_ON(false_IDon);
    Indx_ON_False = Indx_ON(find(cell2mat(arrayfun(@(x) sum(strcmp(falseNames,ON_uniNamesZ(x))),1:length(ON_uniNamesZ),'Un',0))>0));
    for t = 1:length(Names_p)
       t1 =  find(strcmp(uniNames,Names_p{t}));
       Indx_OFF_1 = Indx_OFF(find(strcmp(OFF_uniNamesZ,Names_p{t})));
       Indx_OFF_2 = Indx_ON_False(find(strcmp(falseNames,Names_p{t}))); 
       if (isempty(IDshared))
           Koff(p,t1,:) = sum(Kb_Match(Indx_OFF_1,:),1);
       else     
           if (~strcmp(Names_p{t},uniNames{IDshared}))  
                Koff(p,t1,:) = sum(Kb_Match(Indx_OFF_1,:),1)+sum(Kb_Match(Indx_OFF_2),1);
           else% remove first ON event 
                Indx_OFF_3 = Indx_ON(find(cell2mat(arrayfun(@(x) sum(strcmp(ON_uniNamesZ,SharedName(x))),1:length(SharedName),'Un',0))>0));
                if (length(Indx_OFF_3)>1)
                    Indx_OFF_3 = Indx_OFF_3(2:end);
                    Koff(p,t1,:) = sum(Kb_Match(Indx_OFF_1,:),1)+sum(Kb_Match(Indx_OFF_2,:),1)+sum(Kb_Match(Indx_OFF_3,:),1);
                else
                    Koff(p,t1,:) = sum(Kb_Match(Indx_OFF_1,:),1)+sum(Kb_Match(Indx_OFF_2,:),1);
                end
            end
       end
    end
    catch
        swerwe = 1;
    end
end
end
end
