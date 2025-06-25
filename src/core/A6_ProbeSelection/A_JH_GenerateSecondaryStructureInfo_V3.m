function [Ks_eq,dHs_eq,dSs_eq,dHs_f,dSs_f,dHs_r,dSs_r,dCps_eq,Kd_eq,dHd_eq,dSd_eq,dHd_f,dSd_f,dHd_r,dSd_r,dCpd_eq,N_Self,N_Cross,ElapsedTimeForSecondaryStructureSeq,ElapsedEquilibriumRateEvaluation] = A_JH_GenerateSecondaryStructureInfo_V3(probes,FinalProbeSet,settings)
%Jason Hughes Software to get self binding and cross binding energy for probes used together in a probe set
N_methods = 8;
N_methods2 = 3;
kb = 0.001987204259;%boltzman constant
scr_mat = [-1,-1,-1,1;-1,-1,1,-1;-1,1,-1,-1;1,-1,-1,-1;];  
FinalProbeSet = sort(FinalProbeSet,'ascend');
RemoveMisMatches = settings.RemoveMisMatches;
SaltConcentration = settings.SaltConcentration;
T_hybrid = settings.HybridizationTemperature;
PrimerConcentration = settings.PrimerConcentration;
RV = @(x) (seqrcomplement(x));
pi_seq = cell(1,size(probes,1));
for i=1:size(probes,1)
   pi_seq{i} = RV(probes{i,2});
end
start = tic; 
warning('off','bioinfo:oligoprop:SeqLengthTooShort');
%Get Sequence for Hairpins and SelfDimers
HairSeq = cellfun(@(x) oligoprop(RV(x),'Salt',SaltConcentration,'Temp',T_hybrid).Hairpins,{probes{FinalProbeSet,2}},'UniformOutput',false);
locHair = cellfun(@(x) isstrprop(x,'upper'),{HairSeq{:}},'UniformOutput',false);
sublocHair = @(y) arrayfun(@(x) HairSeq{y}(x,(locHair{y}(x,:))),1:size(locHair{y},1),'UniformOutput',false);
HairSeqParsed = arrayfun(@(x) sublocHair(x),1:length(FinalProbeSet),'UniformOutput',false);
clear HairSeq locHair sublocHair
SelfDimerSeq = cellfun(@(x) oligoprop(RV(x),'Salt',SaltConcentration,'Temp',T_hybrid).Dimers,{probes{FinalProbeSet,2}},'UniformOutput',false);
locSelf = cellfun(@(x) isstrprop(x,'upper'),{SelfDimerSeq{:}},'UniformOutput',false);
sublocSelf = @(y) arrayfun(@(x) SelfDimerSeq{y}(x,(locSelf{y}(x,:))),1:size(locSelf{y},1),'UniformOutput',false);
SelfDimerSeqParsed = arrayfun(@(x) sublocSelf(x),1:length(FinalProbeSet),'UniformOutput',false);
clear SelfDimerSeq locSelf sublocSelf
%Find Unique Pairs of Hairpin and Self Hybridization Sequences
SelfSeqParsed = arrayfun(@(x) union(SelfDimerSeqParsed{x},HairSeqParsed{x}),1:length(FinalProbeSet),'UniformOutput',false);
%remove flips
for x = 1:length(SelfSeqParsed)   
    Self_Flip_Identity = cell2mat(arrayfun(@(y) strcmp(flip(SelfSeqParsed{x}{y}),SelfSeqParsed{x}).',1:length(SelfSeqParsed{x}),'UniformOutput',false));
    %Remove 
    Flip_Identity = triu(Self_Flip_Identity);
    [row,~] = find(Flip_Identity);
    if (~isempty(row))
        for sk = 1:length(row)  
            SelfSeqParsed{x}{row(sk)} = [];
        end
        SelfSeqParsed{x} = {SelfSeqParsed{x}{~cellfun(@isempty,SelfSeqParsed{x})}};
    end
end
selfNotExist = cellfun(@isempty,SelfSeqParsed);
selfId1 = find(selfNotExist==1);
for i = 1:length(selfId1)
   SelfSeqParsed{selfId1(i)}=cell(0,2); 
end
N_Self = cellfun(@(x) length(x), SelfSeqParsed);
clear row Flip_Identity Self_Flip_Identity
  
               
for i=FinalProbeSet
    for j=FinalProbeSet
        if (i ~= j)
            [~,tempAlign1]= swalign(pi_seq{i},pi_seq{j},'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
            TopCrossDimerSeqParsed1{i}{j} = tempAlign1(1,:);
            BotCrossDimerSeqParsed1{i}{j} = tempAlign1(3,:);
            TopCrossDimerSeqParsed1{j}{i} = tempAlign1(1,:);
            BotCrossDimerSeqParsed1{j}{i} = tempAlign1(3,:);
            [~,tempAlign2] = swalign(pi_seq{i},reverse(pi_seq{j}),'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
            TopCrossDimerSeqParsed2{i}{j} = tempAlign2(1,:);
            BotCrossDimerSeqParsed2{i}{j} = tempAlign2(3,:);
            TopCrossDimerSeqParsed2{j}{i} = tempAlign2(1,:);
            BotCrossDimerSeqParsed2{j}{i} = tempAlign2(3,:);
            [~,tempAlign3] = swalign(reverse(pi_seq{i}),pi_seq{j},'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
            TopCrossDimerSeqParsed3{i}{j} = tempAlign3(1,:);
            BotCrossDimerSeqParsed3{i}{j} = tempAlign3(3,:);
            TopCrossDimerSeqParsed3{j}{i} = tempAlign3(1,:);
            BotCrossDimerSeqParsed3{j}{i} = tempAlign3(3,:);
            [~,tempAlign4] = swalign(reverse(reverse(pi_seq{i})),reverse(pi_seq{j}),'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
            TopCrossDimerSeqParsed4{i}{j} = tempAlign4(1,:);
            BotCrossDimerSeqParsed4{i}{j} = tempAlign4(3,:);
            TopCrossDimerSeqParsed4{j}{i} = tempAlign4(1,:);
            BotCrossDimerSeqParsed4{j}{i} = tempAlign4(3,:); 
            clear tempAlign*
        else
            [~,tempAlign1] = swalign(pi_seq{i},reverse(pi_seq{j}),'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
            TopCrossDimerSeqParsed1{i}{j} = tempAlign1(1,:);
            BotCrossDimerSeqParsed1{i}{j} = tempAlign1(3,:);
            [~,tempAlign2] = swalign(pi_seq{i},pi_seq{j},'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
            TopCrossDimerSeqParsed2{i}{j} = tempAlign2(1,:);
            BotCrossDimerSeqParsed2{i}{j} = tempAlign2(3,:);
            clear tempAlign*
            TopCrossDimerSeqParsed3{i}{j} = '';
            BotCrossDimerSeqParsed3{i}{j} = '';
            TopCrossDimerSeqParsed4{i}{j} = '';
            BotCrossDimerSeqParsed4{i}{j} = '';        
        end
    end
end
ElapsedTimeForSecondaryStructureSeq = toc(start);   
%Find Unique Pairs of Hairpin and Self Hybridization Sequences

z0 = 1;z = 0;
for v = 1:length(FinalProbeSet)
count = length({TopCrossDimerSeqParsed1{FinalProbeSet(v)}{~cellfun(@isempty,TopCrossDimerSeqParsed1{FinalProbeSet(v)})}});
z = count + z;
CrossDictionaryTop1(z0:z) = {TopCrossDimerSeqParsed1{FinalProbeSet(v)}{~cellfun(@isempty,TopCrossDimerSeqParsed1{FinalProbeSet(v)})}};
z0 = z + 1;
end
z0 = 1;z = 0;
for v = 1:length(FinalProbeSet)
count = length({TopCrossDimerSeqParsed2{FinalProbeSet(v)}{~cellfun(@isempty,TopCrossDimerSeqParsed2{FinalProbeSet(v)})}});
z = count + z;
CrossDictionaryTop2(z0:z) = {TopCrossDimerSeqParsed2{FinalProbeSet(v)}{~cellfun(@isempty,TopCrossDimerSeqParsed2{FinalProbeSet(v)})}};
z0 = z + 1;
end
z0 = 1;z = 0;
for v = 1:length(FinalProbeSet)
count = length({TopCrossDimerSeqParsed3{FinalProbeSet(v)}{~cellfun(@isempty,TopCrossDimerSeqParsed3{FinalProbeSet(v)})}});
z = count + z;
CrossDictionaryTop3(z0:z) = {TopCrossDimerSeqParsed3{FinalProbeSet(v)}{~cellfun(@isempty,TopCrossDimerSeqParsed3{FinalProbeSet(v)})}};
z0 = z + 1;
end
z0 = 1;z = 0;
for v = 1:length(FinalProbeSet)
count = length({TopCrossDimerSeqParsed4{FinalProbeSet(v)}{~cellfun(@isempty,TopCrossDimerSeqParsed4{FinalProbeSet(v)})}});
z = count + z;
CrossDictionaryTop4(z0:z) = {TopCrossDimerSeqParsed4{FinalProbeSet(v)}{~cellfun(@isempty,TopCrossDimerSeqParsed4{FinalProbeSet(v)})}};
z0 = z + 1;
end
z0 = 1;z = 0;
for v = 1:length(FinalProbeSet)
count = length({BotCrossDimerSeqParsed1{FinalProbeSet(v)}{~cellfun(@isempty,BotCrossDimerSeqParsed1{FinalProbeSet(v)})}});
z = count + z;
CrossDictionaryBot1(z0:z) = {BotCrossDimerSeqParsed1{FinalProbeSet(v)}{~cellfun(@isempty,BotCrossDimerSeqParsed1{FinalProbeSet(v)})}};
z0 = z + 1;
end
z0 = 1;z = 0;
for v = 1:length(FinalProbeSet)
count = length({BotCrossDimerSeqParsed2{FinalProbeSet(v)}{~cellfun(@isempty,BotCrossDimerSeqParsed2{FinalProbeSet(v)})}});
z = count + z;
CrossDictionaryBot2(z0:z) = {BotCrossDimerSeqParsed2{FinalProbeSet(v)}{~cellfun(@isempty,BotCrossDimerSeqParsed2{FinalProbeSet(v)})}};
z0 = z + 1;
end
z0 = 1;z = 0;
for v = 1:length(FinalProbeSet)
count = length({BotCrossDimerSeqParsed3{FinalProbeSet(v)}{~cellfun(@isempty,BotCrossDimerSeqParsed3{FinalProbeSet(v)})}});
z = count + z;
CrossDictionaryBot3(z0:z) = {BotCrossDimerSeqParsed3{FinalProbeSet(v)}{~cellfun(@isempty,BotCrossDimerSeqParsed3{FinalProbeSet(v)})}};
z0 = z + 1;
end
z0 = 1;z = 0;
for v = 1:length(FinalProbeSet)
count = length({BotCrossDimerSeqParsed4{FinalProbeSet(v)}{~cellfun(@isempty,BotCrossDimerSeqParsed4{FinalProbeSet(v)})}});
z = count + z;
CrossDictionaryBot4(z0:z) = {BotCrossDimerSeqParsed4{FinalProbeSet(v)}{~cellfun(@isempty,BotCrossDimerSeqParsed4{FinalProbeSet(v)})}};
z0 = z + 1;
end
CrossDictionary = [CrossDictionaryTop1(:)' CrossDictionaryTop2(:)' CrossDictionaryTop3(:)' CrossDictionaryTop4(:)'...
    CrossDictionaryBot1(:)' CrossDictionaryBot2(:)' CrossDictionaryBot3(:)' CrossDictionaryBot4(:)'];
CrossDimerDictionary.Names = unique(CrossDictionary);
count = 1;
%Remove Pairs of Flips (Redundant matches, that look different but are a different pair both flipped)
for u = 1:length(FinalProbeSet)
    for v = 1:length(FinalProbeSet)
        try
            DictionaryPairs(count,1) = find(strcmp(CrossDimerDictionary.Names,TopCrossDimerSeqParsed1{FinalProbeSet(u)}{FinalProbeSet(v)}));
            DictionaryPairs(count,2) = find(strcmp(CrossDimerDictionary.Names,BotCrossDimerSeqParsed1{FinalProbeSet(u)}{FinalProbeSet(v)}));
            DictionaryPairs(count,3) = u;
            DictionaryPairs(count,4) = v;
            count = count + 1;
        catch
        end 
    end
end
for u = 1:length(FinalProbeSet)
    for v = 1:length(FinalProbeSet)
        try
            DictionaryPairs(count,1) = find(strcmp(CrossDimerDictionary.Names,TopCrossDimerSeqParsed2{FinalProbeSet(u)}{FinalProbeSet(v)}));
            DictionaryPairs(count,2) = find(strcmp(CrossDimerDictionary.Names,BotCrossDimerSeqParsed2{FinalProbeSet(u)}{FinalProbeSet(v)}));
            DictionaryPairs(count,3) = u;
            DictionaryPairs(count,4) = v;
            count = count + 1;
        catch
        end
    end
end
for u = 1:length(FinalProbeSet)
    for v = 1:length(FinalProbeSet)
        try
            DictionaryPairs(count,1) = find(strcmp(CrossDimerDictionary.Names,TopCrossDimerSeqParsed3{FinalProbeSet(u)}{FinalProbeSet(v)}));
            DictionaryPairs(count,2) = find(strcmp(CrossDimerDictionary.Names,BotCrossDimerSeqParsed3{FinalProbeSet(u)}{FinalProbeSet(v)}));
            DictionaryPairs(count,3) = u;
            DictionaryPairs(count,4) = v;
            count = count + 1;
        catch
        end
    end
end
for u = 1:length(FinalProbeSet)
    for v = 1:length(FinalProbeSet)
        try
            DictionaryPairs(count,1) = find(strcmp(CrossDimerDictionary.Names,TopCrossDimerSeqParsed4{FinalProbeSet(u)}{FinalProbeSet(v)}));
            DictionaryPairs(count,2) = find(strcmp(CrossDimerDictionary.Names,BotCrossDimerSeqParsed4{FinalProbeSet(u)}{FinalProbeSet(v)}));
            DictionaryPairs(count,3) = u;
            DictionaryPairs(count,4) = v;
            count = count + 1;
        catch
        end
    end
end
UniquePairs = unique(DictionaryPairs,'rows');
%Map Unique Pairs back to probes pairs
for u=1:length(FinalProbeSet)
    for v=1:length(FinalProbeSet)
       rows = find((UniquePairs(:,3)==u).*(UniquePairs(:,4)==v));
       for k = 1:length(rows)
            CrossDimerSeqParsed{u,v}{k,1} = CrossDimerDictionary.Names{UniquePairs(rows(k),1)};
            CrossDimerSeqParsed{u,v}{k,2} = CrossDimerDictionary.Names{UniquePairs(rows(k),2)};
       end
    end
end
%Remove Pairs of Flips (Redundant matches, that look different but are a different pair both flipped)
for u = 1:length(FinalProbeSet)
    for v = 1:length(FinalProbeSet)
        nl = size(CrossDimerSeqParsed{u,v},1);
        Cross_Flip_Identity1 = cell2mat(arrayfun(@(y) strcmp(flip(CrossDimerSeqParsed{u,v}{y,1}),{CrossDimerSeqParsed{u,v}{:,1}}).',1:nl,'UniformOutput',false));
        Cross_Flip_Identity2 = cell2mat(arrayfun(@(y) strcmp(flip(CrossDimerSeqParsed{u,v}{y,2}),{CrossDimerSeqParsed{u,v}{:,2}}).',1:nl,'UniformOutput',false));
        Flip_Identity1 = triu(Cross_Flip_Identity1);
        Flip_Identity2 = triu(Cross_Flip_Identity2);
        [row1,~] = find(Flip_Identity1);
        [row2,~] = find(Flip_Identity2);
        if (~isempty(row1)&&~isempty(row2))
            for w = 1:length(row1)
            CrossDimerSeqParsed{u,v}{row1(w),1} = [];
            end
            for w = 1:length(row2)
            CrossDimerSeqParsed{u,v}{row2(w),2} = [];
            end
        end
    end  
end
clear row Flip_Identity1 Flip_Identity2 Cross_Flip_Identity1 Cross_Flip_Identity2 row1 row2

crossNotExist = cellfun(@isempty,CrossDimerSeqParsed);
[crossId1,crossId2] = find(crossNotExist==1);
for i = 1:length(crossId1)
   CrossDimerSeqParsed{crossId1(i),crossId2(i)}=cell(0,2);
end
CrossDimerSeqParsed = cellfun(@(x) x(sum(cellfun(@isempty,x),2)'==0,1:2), CrossDimerSeqParsed,'Un',0);
N_Cross = cellfun(@(x) size(x,1), CrossDimerSeqParsed);

start2 = tic;
%Initialize matrix for storing binding affinity calculations
Ks_eq = ndSparse.build([size(probes,1),max([N_Self 1]),N_methods],0);
dHs_eq = ndSparse.build([size(probes,1),max([N_Self 1]),N_methods],0);
dSs_eq = ndSparse.build([size(probes,1),max([N_Self 1]),N_methods],0);
dHs_f = ndSparse.build([size(probes,1),max([N_Self 1]),N_methods2],0);
dSs_f = ndSparse.build([size(probes,1),max([N_Self 1]),N_methods2],0);
dHs_r = ndSparse.build([size(probes,1),max([N_Self 1]),N_methods2],0);
dSs_r = ndSparse.build([size(probes,1),max([N_Self 1]),N_methods2],0);
dCps_eq = ndSparse.build([size(probes,1),max([N_Self 1]),N_methods],0);

Kd_eq = ndSparse.build([size(probes,1),size(probes,1),max([max(N_Cross(:)) 1]),N_methods],0);
dHd_eq = ndSparse.build([size(probes,1),size(probes,1),max([max(N_Cross(:)) 1]),N_methods],0);
dSd_eq = ndSparse.build([size(probes,1),size(probes,1),max([max(N_Cross(:)) 1]),N_methods],0);
dCpd_eq = ndSparse.build([size(probes,1),size(probes,1),max([max(N_Cross(:)) 1]),N_methods],0);
dHd_f = ndSparse.build([size(probes,1),size(probes,1),max([max(N_Cross(:)) 1]),N_methods2],0);
dSd_f = ndSparse.build([size(probes,1),size(probes,1),max([max(N_Cross(:)) 1]),N_methods2],0);
dHd_r = ndSparse.build([size(probes,1),size(probes,1),max([max(N_Cross(:)) 1]),N_methods2],0);
dSd_r = ndSparse.build([size(probes,1),size(probes,1),max([max(N_Cross(:)) 1]),N_methods2],0);
%Compute binding affinity for self and cross-dimers (equilibrium, and transition-state forward/reverse rates)        
sequence_duplexes_thermo_generator_struct_Multi = struct();
sequence_duplexes_thermo_generator_struct_Multi.Model{1} = F_NearestNeighbors_Parser('Bres86','src/thirdparty/VarGibbs-4.1/P-BS86.par',[]);
sequence_duplexes_thermo_generator_struct_Multi.Model{2}  = F_NearestNeighbors_Parser('Sant96','src/thirdparty/VarGibbs-4.1/AOP-SL96.par',[]);
sequence_duplexes_thermo_generator_struct_Multi.Model{3}  = F_NearestNeighbors_Parser('Sant98','src/thirdparty/VarGibbs-4.1/AOP-SL98.par',[]);
sequence_duplexes_thermo_generator_struct_Multi.Model{4}   = F_NearestNeighbors_Parser('Sugi96','src/thirdparty/VarGibbs-4.1/P-SG96.par',[]);
sequence_duplexes_thermo_generator_struct_Multi.Model{5}   = F_NearestNeighbors_Parser('Sant04','src/thirdparty/VarGibbs-4.1/P-SL04.par',[]);
sequence_duplexes_thermo_generator_struct_Multi.Model{6}   = F_NearestNeighbors_Parser('Allawi97','src/thirdparty/VarGibbs-4.1/P-AL97.par',[]);
sequence_duplexes_thermo_generator_struct_Multi.Model{7}   = F_NearestNeighbors_Parser('Rejali21','src/thirdparty/VarGibbs-4.1/AOP-RJ21KE.par',[]);
sequence_duplexes_thermo_generator_struct_Multi.Model{8}   = F_NearestNeighbors_Parser('Martins24','src/thirdparty/VarGibbs-4.1/AOP-OW04-69.par','src/thirdparty/VarGibbs-4.1/AOP-MM-60.par');
sequence_duplexes_thermo_generator_structure = struct2table([sequence_duplexes_thermo_generator_struct_Multi.Model{:}]);

for i=FinalProbeSet
    v = find(FinalProbeSet==i);
    for j=1:N_Self(v)
         [temp_dHeq, temp_dSeq, temp_dGeq, ...
          temp_dHf, temp_dSf, ~, ...
          temp_dHr, temp_dSr, ~,temp_dCpeq,~] = ...
          F_DeltaGibson_V3(strrep(SelfSeqParsed{v}{j},'-','N'),strrep(SelfSeqParsed{v}{j},'-','N'),SaltConcentration,T_hybrid,PrimerConcentration,sequence_duplexes_thermo_generator_structure);
          dHs_eq(i,j,:) = temp_dHeq;
          dSs_eq(i,j,:) = temp_dSeq;
          Ks_eq(i,j,:) = exp(-temp_dGeq/(kb*(T_hybrid+273.15)));
          dHs_f(i,j,:) = temp_dHf;
          dSs_f(i,j,:) = temp_dSf; 
          dHs_r(i,j,:) = temp_dHr;
          dSs_r(i,j,:) = temp_dSr;
          dCps_eq(i,j,:) = temp_dCpeq;
    end
    for j=FinalProbeSet(FinalProbeSet>=i)
        w = find(FinalProbeSet==j);
        for k = 1:N_Cross(v,w)
            %try
            if (~isempty(CrossDimerSeqParsed{v,w}{k,1})&&~isempty(CrossDimerSeqParsed{v,w}{k,2}))
                [temp_dHeq, temp_dSeq, temp_dGeq, ...
                 temp_dHf, temp_dSf, ~, ...
                 temp_dHr, temp_dSr, ~,temp_dCpeq, ~] = ...
                 F_DeltaGibson_V3(strrep(CrossDimerSeqParsed{v,w}{k,1},'-','N'),reverse(strrep(CrossDimerSeqParsed{v,w}{k,2},'-','N')),SaltConcentration,T_hybrid,PrimerConcentration,sequence_duplexes_thermo_generator_structure);
                 dHd_eq(i,j,k,:) = temp_dHeq;
                 dSd_eq(i,j,k,:) = temp_dSeq;
                 Kd_eq(i,j,k,:) = exp(-temp_dGeq/(kb*(T_hybrid+273.15)));
                 dHd_f(i,j,k,:) = temp_dHf;
                 dSd_f(i,j,k,:) = temp_dSf; 
                 dHd_r(i,j,k,:) = temp_dHr;
                 dSd_r(i,j,k,:) = temp_dSr;
                 dHd_eq(j,i,k,:) = temp_dHeq;
                 dSd_eq(j,i,k,:) = temp_dSeq;
                 Kd_eq(j,i,k,:) = exp(-temp_dGeq/(kb*(T_hybrid+273.15)));
                 dHd_f(j,i,k,:) = temp_dHf;
                 dSd_f(j,i,k,:) = temp_dSf; 
                 dHd_r(j,i,k,:) = temp_dHr;
                 dSd_r(j,i,k,:) = temp_dSr;  
                 dCpd_eq(i,j,k,:) = temp_dCpeq;
            end
            %catch
             %   ss = 1;
            %end
        end
   end
end 
ElapsedEquilibriumRateEvaluation = toc(start2);
end