function [Ks_eq,dHs_eq,dSs_eq,dHs_f,dSs_f,dHs_r,dSs_r,dCps_eq,Kd_eq,dHd_eq,dSd_eq,dHd_f,dSd_f,dHd_r,dSd_r,dCpd_eq,N_Self,N_Cross,ElapsedTimeForSecondaryStructureSeq,ElapsedEquilibriumRateEvaluation] = A_JH_GenerateSecondaryStructureInfo_V3(probes,FinalProbeSet,settings)
%Jason Hughes Software to get self binding and cross binding energy for probes used together in a probe set
N_methods = 8;
N_methods2 = 3;
kb = 0.001987204259;%boltzman constant
scr_mat = [-1,-1,-1,1;-1,-1,1,-1;-1,1,-1,-1;1,-1,-1,-1;];  
FinalProbeSet = sort(FinalProbeSet,'ascend');
SaltConcentration = settings.SaltConcentration;
T_hybrid = settings.HybridizationTemperature;
PrimerConcentration = settings.PrimerConcentration;
RV = @(x) (seqrcomplement(x));
pi_seq = cell(1,size(probes,1));
for i=1:size(probes,1)
   pi_seq{i} = RV(probes{i,2});
end


%% Get Sequence for Hairpins and SelfDimers
start = tic; 
warning('off','bioinfo:oligoprop:SeqLengthTooShort');
HairSeq = cellfun(@(x) oligoprop(RV(x),'Salt',SaltConcentration,'Temp',T_hybrid).Hairpins,probes(FinalProbeSet,2)','UniformOutput',false);
locHair = cellfun(@(x) isstrprop(x,'upper'),HairSeq,'UniformOutput',false);
sublocHair = @(y) arrayfun(@(x) HairSeq{y}(x,(locHair{y}(x,:))),1:size(locHair{y},1),'UniformOutput',false);
HairSeqParsed = arrayfun(@(x) sublocHair(x),1:length(FinalProbeSet),'UniformOutput',false);
clear HairSeq locHair sublocHair
SelfDimerSeq = cellfun(@(x) oligoprop(RV(x),'Salt',SaltConcentration,'Temp',T_hybrid).Dimers,probes(FinalProbeSet,2)','UniformOutput',false);
locSelf = cellfun(@(x) isstrprop(x,'upper'),SelfDimerSeq,'UniformOutput',false);
sublocSelf = @(y) arrayfun(@(x) SelfDimerSeq{y}(x,(locSelf{y}(x,:))),1:size(locSelf{y},1),'UniformOutput',false);
SelfDimerSeqParsed = arrayfun(@(x) sublocSelf(x),1:length(FinalProbeSet),'UniformOutput',false);
clear SelfDimerSeq locSelf sublocSelf

%% Find Unique Pairs of Hairpin and Self Hybridization Sequences
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
        SelfSeqParsed{x} = SelfSeqParsed{x}(~cellfun(@isempty,SelfSeqParsed{x}));
    end
end
selfNotExist = cellfun(@isempty,SelfSeqParsed);
selfId1 = find(selfNotExist==1);
for i = 1:length(selfId1)
   SelfSeqParsed{selfId1(i)}=cell(0,2); 
end
N_Self = cellfun(@(x) length(x), SelfSeqParsed);
clear row Flip_Identity Self_Flip_Identity
  
TopCrossDimerSeqParsed1 = cell(1,length(FinalProbeSet));    
TopCrossDimerSeqParsed2 = cell(1,length(FinalProbeSet));  
TopCrossDimerSeqParsed3 = cell(1,length(FinalProbeSet));  
TopCrossDimerSeqParsed4 = cell(1,length(FinalProbeSet));  
BotCrossDimerSeqParsed1 = cell(1,length(FinalProbeSet));  
BotCrossDimerSeqParsed2 = cell(1,length(FinalProbeSet));  
BotCrossDimerSeqParsed3 = cell(1,length(FinalProbeSet));  
BotCrossDimerSeqParsed4 = cell(1,length(FinalProbeSet));  
for i=1:1:length(FinalProbeSet)
    for j=1:length(FinalProbeSet)
        if (i ~= j)
            [~,tempAlign1] = swalign(pi_seq{FinalProbeSet(i)},pi_seq{FinalProbeSet(j)},'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
            TopCrossDimerSeqParsed1{i}{j} = tempAlign1(1,:);
            BotCrossDimerSeqParsed1{i}{j} = tempAlign1(3,:);
            TopCrossDimerSeqParsed1{j}{i} = tempAlign1(1,:);
            BotCrossDimerSeqParsed1{j}{i} = tempAlign1(3,:);
            [~,tempAlign2] = swalign(pi_seq{FinalProbeSet(i)},reverse(pi_seq{FinalProbeSet(j)}),'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
            TopCrossDimerSeqParsed2{i}{j} = tempAlign2(1,:);
            BotCrossDimerSeqParsed2{i}{j} = tempAlign2(3,:);
            TopCrossDimerSeqParsed2{j}{i} = tempAlign2(1,:);
            BotCrossDimerSeqParsed2{j}{i} = tempAlign2(3,:);
            [~,tempAlign3] = swalign(reverse(pi_seq{FinalProbeSet(i)}),pi_seq{FinalProbeSet(j)},'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
            TopCrossDimerSeqParsed3{i}{j} = tempAlign3(1,:);
            BotCrossDimerSeqParsed3{i}{j} = tempAlign3(3,:);
            TopCrossDimerSeqParsed3{j}{i} = tempAlign3(1,:);
            BotCrossDimerSeqParsed3{j}{i} = tempAlign3(3,:);
            [~,tempAlign4] = swalign(reverse(reverse(pi_seq{FinalProbeSet(i)})),reverse(pi_seq{FinalProbeSet(j)}),'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
            TopCrossDimerSeqParsed4{i}{j} = tempAlign4(1,:);
            BotCrossDimerSeqParsed4{i}{j} = tempAlign4(3,:);
            TopCrossDimerSeqParsed4{j}{i} = tempAlign4(1,:);
            BotCrossDimerSeqParsed4{j}{i} = tempAlign4(3,:); 
            clear tempAlign*
        else
            [~,tempAlign1] = swalign(pi_seq{FinalProbeSet(i)},reverse(pi_seq{FinalProbeSet(j)}),'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
            TopCrossDimerSeqParsed1{i}{j} = tempAlign1(1,:);
            BotCrossDimerSeqParsed1{i}{j} = tempAlign1(3,:);
            [~,tempAlign2] = swalign(pi_seq{FinalProbeSet(i)},pi_seq{FinalProbeSet(j)},'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
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
%% Find Unique Pairs of Hairpin and Self Hybridization Sequences
z0 = 1;z = 0;
for v = 1:length(FinalProbeSet)
count = length(TopCrossDimerSeqParsed1{v}(~cellfun(@isempty,TopCrossDimerSeqParsed1{v})));
z = count + z;
CrossDictionaryTop1(z0:z) = TopCrossDimerSeqParsed1{v}(~cellfun(@isempty,TopCrossDimerSeqParsed1{v}));
z0 = z + 1;
end
z0 = 1;z = 0;
for v = 1:length(FinalProbeSet)
count = length(TopCrossDimerSeqParsed2{v}(~cellfun(@isempty,TopCrossDimerSeqParsed2{v})));
z = count + z;
CrossDictionaryTop2(z0:z) = TopCrossDimerSeqParsed2{v}(~cellfun(@isempty,TopCrossDimerSeqParsed2{v}));
z0 = z + 1;
end
z0 = 1;z = 0;
for v = 1:length(FinalProbeSet)
count = length(TopCrossDimerSeqParsed3{v}(~cellfun(@isempty,TopCrossDimerSeqParsed3{v})));
z = count + z;
CrossDictionaryTop3(z0:z) = TopCrossDimerSeqParsed3{v}(~cellfun(@isempty,TopCrossDimerSeqParsed3{v}));
z0 = z + 1;
end
z0 = 1;z = 0;
for v = 1:length(FinalProbeSet)
count = length(TopCrossDimerSeqParsed4{v}(~cellfun(@isempty,TopCrossDimerSeqParsed4{v})));
z = count + z;
CrossDictionaryTop4(z0:z) = TopCrossDimerSeqParsed4{v}(~cellfun(@isempty,TopCrossDimerSeqParsed4{v}));
z0 = z + 1;
end
z0 = 1;z = 0;
for v = 1:length(FinalProbeSet)
count = length(BotCrossDimerSeqParsed1{v}(~cellfun(@isempty,BotCrossDimerSeqParsed1{v})));
z = count + z;
CrossDictionaryBot1(z0:z) = BotCrossDimerSeqParsed1{v}(~cellfun(@isempty,BotCrossDimerSeqParsed1{v}));
z0 = z + 1;
end
z0 = 1;z = 0;
for v = 1:length(FinalProbeSet)
count = length(BotCrossDimerSeqParsed2{v}(~cellfun(@isempty,BotCrossDimerSeqParsed2{v})));
z = count + z;
CrossDictionaryBot2(z0:z) = BotCrossDimerSeqParsed2{v}(~cellfun(@isempty,BotCrossDimerSeqParsed2{v}));
z0 = z + 1;
end
z0 = 1;z = 0;
for v = 1:length(FinalProbeSet)
count = length(BotCrossDimerSeqParsed3{v}(~cellfun(@isempty,BotCrossDimerSeqParsed3{v})));
z = count + z;
CrossDictionaryBot3(z0:z) = BotCrossDimerSeqParsed3{v}(~cellfun(@isempty,BotCrossDimerSeqParsed3{v}));
z0 = z + 1;
end
z0 = 1;z = 0;
for v = 1:length(FinalProbeSet)
count = length(BotCrossDimerSeqParsed4{v}(~cellfun(@isempty,BotCrossDimerSeqParsed4{v})));
z = count + z;
CrossDictionaryBot4(z0:z) = BotCrossDimerSeqParsed4{v}(~cellfun(@isempty,BotCrossDimerSeqParsed4{v}));
z0 = z + 1;
end
CrossDictionary = [CrossDictionaryTop1(:)' CrossDictionaryTop2(:)' CrossDictionaryTop3(:)' CrossDictionaryTop4(:)'...
    CrossDictionaryBot1(:)' CrossDictionaryBot2(:)' CrossDictionaryBot3(:)' CrossDictionaryBot4(:)'];
CrossDimerDictionary.Names = unique(CrossDictionary);
count = 1;
%% Remove Pairs of Flips (Redundant matches, that look different but are a different pair both flipped)
for u = 1:length(FinalProbeSet)
    for v = 1:length(FinalProbeSet)
        try
            DictionaryPairs(count,1) = find(strcmp(CrossDimerDictionary.Names,TopCrossDimerSeqParsed1{u}{v}));
            DictionaryPairs(count,2) = find(strcmp(CrossDimerDictionary.Names,BotCrossDimerSeqParsed1{u}{v}));
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
            DictionaryPairs(count,1) = find(strcmp(CrossDimerDictionary.Names,TopCrossDimerSeqParsed2{u}{v}));
            DictionaryPairs(count,2) = find(strcmp(CrossDimerDictionary.Names,BotCrossDimerSeqParsed2{u}{v}));
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
            DictionaryPairs(count,1) = find(strcmp(CrossDimerDictionary.Names,TopCrossDimerSeqParsed3{u}{v}));
            DictionaryPairs(count,2) = find(strcmp(CrossDimerDictionary.Names,BotCrossDimerSeqParsed3{u}{v}));
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
            DictionaryPairs(count,1) = find(strcmp(CrossDimerDictionary.Names,TopCrossDimerSeqParsed4{u}{v}));
            DictionaryPairs(count,2) = find(strcmp(CrossDimerDictionary.Names,BotCrossDimerSeqParsed4{u}{v}));
            DictionaryPairs(count,3) = u;
            DictionaryPairs(count,4) = v;
            count = count + 1;
        catch
        end
    end
end
UniquePairs = unique(DictionaryPairs,'rows');
%% Map Unique Pairs back to probes pairs
CrossDimerSeqParsed = cell(length(FinalProbeSet),length(FinalProbeSet));
for u=1:length(FinalProbeSet)
    for v=1:length(FinalProbeSet)
       rows = find((UniquePairs(:,3)==u).*(UniquePairs(:,4)==v));
       for k = 1:length(rows)
            CrossDimerSeqParsed{u,v}{k,1} = CrossDimerDictionary.Names{UniquePairs(rows(k),1)};
            CrossDimerSeqParsed{u,v}{k,2} = CrossDimerDictionary.Names{UniquePairs(rows(k),2)};
       end
    end
end




%% Remove Pairs of Flips (Redundant matches, that look different but are a different pair both flipped)
for u = 1:length(FinalProbeSet)
    for v = 1:length(FinalProbeSet)
        nl = size(CrossDimerSeqParsed{u,v},1);
        Cross_Flip_Identity1 = cell2mat(arrayfun(@(y) strcmp(flip(CrossDimerSeqParsed{u,v}{y,1}),CrossDimerSeqParsed{u,v}(:,1)),1:nl,'UniformOutput',false));
        Cross_Flip_Identity2 = cell2mat(arrayfun(@(y) strcmp(flip(CrossDimerSeqParsed{u,v}{y,2}),CrossDimerSeqParsed{u,v}(:,2)),1:nl,'UniformOutput',false));
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

Self_V_List = cell2mat(cellfun(@(V,J) V*ones(1,J),num2cell(1:length(SelfSeqParsed)),num2cell(cellfun(@length,SelfSeqParsed)),'Un',0))';
Self_J_List = cell2mat(cellfun(@(V,J) 1:J,num2cell(1:length(SelfSeqParsed)),num2cell(cellfun(@length,SelfSeqParsed)),'Un',0))';
Self_SeqInput = arrayfun(@(x) strrep(SelfSeqParsed{Self_V_List(x)}{Self_J_List(x)},'-','N'),1:length(Self_V_List),'Un',0)';
Cross_V_List = cell2mat(reshape(cellfun(@(V,J) V*ones(1,J),num2cell(meshgrid(1:size(CrossDimerSeqParsed,1))),cellfun(@(z) size(z,1),CrossDimerSeqParsed,'Un',0),'Un',0),1,[]))';
Cross_W_List = cell2mat(reshape(cellfun(@(W,J) W*ones(1,J),num2cell(meshgrid(1:size(CrossDimerSeqParsed,1))'),cellfun(@(z) size(z,1),CrossDimerSeqParsed,'Un',0),'Un',0),1,[]))';
Cross_K_List = cell2mat(reshape(cellfun(@(V,J) 1:J,num2cell(meshgrid(1:size(CrossDimerSeqParsed,1))),cellfun(@(z) size(z,1),CrossDimerSeqParsed,'Un',0),'Un',0),1,[]))';
Cross_SeqInput1 = arrayfun(@(x) strrep(CrossDimerSeqParsed{Cross_V_List(x),Cross_W_List(x)}{Cross_K_List(x),1},'-','N'),1:length(Cross_V_List),'Un',0)';
Cross_SeqInput2 = arrayfun(@(x) reverse(strrep(CrossDimerSeqParsed{Cross_V_List(x),Cross_W_List(x)}{Cross_K_List(x),2},'-','N')),1:length(Cross_V_List),'Un',0)';

unique_secondary_structure_binding_seqs = unique([Self_SeqInput; Cross_SeqInput1; Cross_SeqInput2;...
    cellfun(@seqreverse,[Self_SeqInput; Cross_SeqInput1; Cross_SeqInput2],'Un',0)]);
binding_seqs_dictionary = dictionary(convertCharsToStrings(unique_secondary_structure_binding_seqs)',1:length(unique_secondary_structure_binding_seqs));
inverse_binding_seqs_dictionary = dictionary(1:length(unique_secondary_structure_binding_seqs),convertCharsToStrings(unique_secondary_structure_binding_seqs)');

Is_Self_Over_Cross = [1*ones(size(Self_V_List)); 0*ones(size(Cross_V_List))];
Secondary_Structure_VV_List = [Self_V_List; Cross_V_List];
Secondary_Structure_VW_List = [Self_V_List; Cross_W_List];
Secondary_Structure_JK_List = [Self_J_List; Cross_K_List];
Secondary_Structure_unique_Seq1_id_List = [binding_seqs_dictionary(convertCharsToStrings(Self_SeqInput)); binding_seqs_dictionary(convertCharsToStrings(Cross_SeqInput1))];
Secondary_Structure_unique_Seq2_id_List = [binding_seqs_dictionary(convertCharsToStrings(Self_SeqInput)); binding_seqs_dictionary(convertCharsToStrings(Cross_SeqInput2))];
seq_pair_ids_nonunique = [Secondary_Structure_unique_Seq1_id_List Secondary_Structure_unique_Seq2_id_List];
[unique_pair_ids,~,matched_unique_location] = unique(seq_pair_ids_nonunique,'rows');
unique_ordered_binding_paired_input_sequences = inverse_binding_seqs_dictionary(unique_pair_ids); 
unique_secondary_structure_pair_to_nonunique_entries = [Is_Self_Over_Cross Secondary_Structure_VV_List Secondary_Structure_VW_List ...
    Secondary_Structure_JK_List Secondary_Structure_unique_Seq1_id_List Secondary_Structure_unique_Seq2_id_List ...
    matched_unique_location];
SelfEQ_I_vector =struct('SelfEQ_I_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
SelfEQ_J_vector = struct('SelfEQ_J_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
SelfEQ_M_vector = struct('SelfEQ_M_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
SelfFR_I_vector =struct('SelfFR_I_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
SelfFR_J_vector = struct('SelfFR_J_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
SelfFR_M_vector = struct('SelfFR_M_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
Ks_eq_vector =struct('Ks_eq_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dHs_eq_vector = struct('dHs_eq_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dSs_eq_vector = struct('dSs_eq_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dHs_f_vector = struct('dHs_f_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dSs_f_vector = struct('dSs_f_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dHs_r_vector =struct('dHs_r_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dSs_r_vector = struct('dSs_r_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dCps_eq_vector = struct('dCps_eq_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
CrossEQ_I_vector =struct('CrossEQ_I_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
CrossEQ_J_vector =struct('CrossEQ_J_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
CrossEQ_K_vector =struct('CrossEQ_K_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
CrossEQ_M_vector =struct('CrossEQ_M_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
CrossFR_I_vector =struct('CrossFR_I_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
CrossFR_J_vector =struct('CrossFR_J_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
CrossFR_K_vector =struct('CrossFR_K_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
CrossFR_M_vector =struct('CrossFR_M_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
Kd_eq_vector =struct('Kd_eq_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dHd_eq_vector = struct('dHd_eq_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dSd_eq_vector = struct('dSd_eq_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dHd_f_vector = struct('dHd_f_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dSd_f_vector = struct('dSd_f_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dHd_r_vector =struct('dHd_r_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dSd_r_vector = struct('dSd_r_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));
dCpd_eq_vector = struct('dCpd_eq_vector',cell(1,size(unique_ordered_binding_paired_input_sequences,1)));

unique_ordered_binding_paired_input_sequences_List = parallel.pool.Constant(unique_ordered_binding_paired_input_sequences);
unique_secondary_structure_pair_to_nonunique_entries_List = parallel.pool.Constant(unique_secondary_structure_pair_to_nonunique_entries);
FinalProbeSet_List = parallel.pool.Constant(FinalProbeSet);
parfor unique_calc = 1:size(unique_ordered_binding_paired_input_sequences,1)
         [temp_dHeq, temp_dSeq, temp_dGeq, ...
          temp_dHf, temp_dSf, ~, ...
          temp_dHr, temp_dSr, ~,temp_dCpeq,~] = ...
          F_DeltaGibson_V3(char(unique_ordered_binding_paired_input_sequences_List.Value(unique_calc,1)),char(unique_ordered_binding_paired_input_sequences_List.Value(unique_calc,2)),SaltConcentration,T_hybrid,PrimerConcentration,sequence_duplexes_thermo_generator_structure);
        self_locs = find(double(unique_secondary_structure_pair_to_nonunique_entries_List.Value(:,7)==unique_calc).*...
            double(unique_secondary_structure_pair_to_nonunique_entries_List.Value(:,1)==1));
        cross_locs = find(double(unique_secondary_structure_pair_to_nonunique_entries_List.Value(:,7)==unique_calc).*...
            double(unique_secondary_structure_pair_to_nonunique_entries_List.Value(:,1)==0));
        if (~isempty(self_locs))
            V_vector = unique_secondary_structure_pair_to_nonunique_entries_List.Value(self_locs,2);
            J_vector = unique_secondary_structure_pair_to_nonunique_entries_List.Value(self_locs,4);
            SelfEQ_I_vector(unique_calc).SelfEQ_I_vector = repmat(FinalProbeSet_List.Value(V_vector),[1 N_methods])';
            SelfEQ_J_vector(unique_calc).SelfEQ_J_vector = repmat(reshape(J_vector,1,[]),[1 N_methods])';
            SelfEQ_M_vector(unique_calc).SelfEQ_M_vector = repelem((1:N_methods)',length(self_locs),1);
            SelfFR_I_vector(unique_calc).SelfFR_I_vector = repmat(FinalProbeSet_List.Value(V_vector),[1 N_methods2])';
            SelfFR_J_vector(unique_calc).SelfFR_J_vector = repmat(reshape(J_vector,1,[]),[1 N_methods2])';
            SelfFR_M_vector(unique_calc).SelfFR_M_vector = repelem((1:N_methods2)',length(self_locs),1);
            Ks_eq_vector(unique_calc).Ks_eq_vector = repelem(exp(-temp_dGeq/(kb*(T_hybrid+273.15))),length(self_locs),1);
            dHs_eq_vector(unique_calc).dHs_eq_vector = repelem(temp_dHeq,length(self_locs),1);
            dSs_eq_vector(unique_calc).dSs_eq_vector = repelem(temp_dSeq,length(self_locs),1);
            dHs_f_vector(unique_calc).dHs_f_vector = repelem(temp_dHf,length(self_locs),1);
            dSs_f_vector(unique_calc).dSs_f_vector = repelem(temp_dSf,length(self_locs),1);
            dHs_r_vector(unique_calc).dHs_r_vector = repelem(temp_dHr,length(self_locs),1);
            dSs_r_vector(unique_calc).dSs_r_vector = repelem(temp_dSr,length(self_locs),1);
            dCps_eq_vector(unique_calc).dCps_eq_vector = repelem(temp_dCpeq,length(self_locs),1);
        end
        if (~isempty(cross_locs))
            V_vector = unique_secondary_structure_pair_to_nonunique_entries_List.Value(cross_locs,2);
            W_vector = unique_secondary_structure_pair_to_nonunique_entries_List.Value(cross_locs,3);
            K_vector = unique_secondary_structure_pair_to_nonunique_entries_List.Value(cross_locs,4);
            CrossEQ_I_vector(unique_calc).CrossEQ_I_vector = repmat(FinalProbeSet_List.Value(V_vector),[1 N_methods])';
            CrossEQ_J_vector(unique_calc).CrossEQ_J_vector = repmat(FinalProbeSet_List.Value(W_vector),[1 N_methods])';
            CrossEQ_K_vector(unique_calc).CrossEQ_K_vector = repmat(reshape(K_vector,1,[]),[1 N_methods])';
            CrossEQ_M_vector(unique_calc).CrossEQ_M_vector = repelem((1:N_methods)',length(cross_locs),1);
            CrossFR_I_vector(unique_calc).CrossFR_I_vector = repmat(FinalProbeSet_List.Value(V_vector),[1 N_methods2])';
            CrossFR_J_vector(unique_calc).CrossFR_J_vector = repmat(FinalProbeSet_List.Value(W_vector),[1 N_methods2])';
            CrossFR_K_vector(unique_calc).CrossFR_K_vector = repmat(reshape(K_vector,1,[]),[1 N_methods2])';
            CrossFR_M_vector(unique_calc).CrossFR_M_vector = repelem((1:N_methods2)',length(cross_locs),1);
            Kd_eq_vector(unique_calc).Kd_eq_vector = repelem(exp(-temp_dGeq/(kb*(T_hybrid+273.15))),length(cross_locs),1);
            dHd_eq_vector(unique_calc).dHd_eq_vector = repelem(temp_dHeq,length(cross_locs),1);
            dSd_eq_vector(unique_calc).dSd_eq_vector = repelem(temp_dSeq,length(cross_locs),1);
            dHd_f_vector(unique_calc).dHd_f_vector = repelem(temp_dHf,length(cross_locs),1);
            dSd_f_vector(unique_calc).dSd_f_vector = repelem(temp_dSf,length(cross_locs),1);
            dHd_r_vector(unique_calc).dHd_r_vector = repelem(temp_dHr,length(cross_locs),1);
            dSd_r_vector(unique_calc).dSd_r_vector = repelem(temp_dSr,length(cross_locs),1);
            dCpd_eq_vector(unique_calc).dCpd_eq_vector = repelem(temp_dCpeq,length(cross_locs),1);
        end
end
SelfEQ_I_vector = vertcat(SelfEQ_I_vector(:).SelfEQ_I_vector);
SelfEQ_J_vector = vertcat(SelfEQ_J_vector(:).SelfEQ_J_vector);
SelfEQ_M_vector = vertcat(SelfEQ_M_vector(:).SelfEQ_M_vector);
SelfFR_I_vector = vertcat(SelfFR_I_vector(:).SelfFR_I_vector);
SelfFR_J_vector = vertcat(SelfFR_J_vector(:).SelfFR_J_vector);
SelfFR_M_vector = vertcat(SelfFR_M_vector(:).SelfFR_M_vector);
Ks_eq_vector = vertcat(Ks_eq_vector(:).Ks_eq_vector);
dHs_eq_vector = vertcat(dHs_eq_vector(:).dHs_eq_vector);
dSs_eq_vector = vertcat(dSs_eq_vector(:).dSs_eq_vector);
dCps_eq_vector = vertcat(dCps_eq_vector(:).dCps_eq_vector);
dHs_f_vector = vertcat(dHs_f_vector(:).dHs_f_vector);
dSs_f_vector = vertcat(dSs_f_vector(:).dSs_f_vector);
dHs_r_vector = vertcat(dHs_r_vector(:).dHs_r_vector);
dSs_r_vector = vertcat(dSs_r_vector(:).dSs_r_vector);
CrossEQ_I_vector = vertcat(CrossEQ_I_vector(:).CrossEQ_I_vector);
CrossEQ_J_vector = vertcat(CrossEQ_J_vector(:).CrossEQ_J_vector);
CrossEQ_K_vector = vertcat(CrossEQ_K_vector(:).CrossEQ_K_vector);
CrossEQ_M_vector = vertcat(CrossEQ_M_vector(:).CrossEQ_M_vector);
CrossFR_I_vector = vertcat(CrossFR_I_vector(:).CrossFR_I_vector);
CrossFR_J_vector = vertcat(CrossFR_J_vector(:).CrossFR_J_vector);
CrossFR_K_vector = vertcat(CrossFR_K_vector(:).CrossFR_K_vector);
CrossFR_M_vector = vertcat(CrossFR_M_vector(:).CrossFR_M_vector);
dKd_eq_vector = vertcat(Kd_eq_vector(:).Kd_eq_vector);
dHd_eq_vector = vertcat(dHd_eq_vector(:).dHd_eq_vector);
dSd_eq_vector = vertcat(dSd_eq_vector(:).dSd_eq_vector);
dCpd_eq_vector = vertcat(dCpd_eq_vector(:).dCpd_eq_vector);
dHd_f_vector = vertcat(dHd_f_vector(:).dHd_f_vector);
dSd_f_vector = vertcat(dSd_f_vector(:).dSd_f_vector);
dHd_r_vector = vertcat(dHd_r_vector(:).dHd_r_vector);
dSd_r_vector = vertcat(dSd_r_vector(:).dSd_r_vector);
SelfEQ_IJM_vector = [SelfEQ_I_vector SelfEQ_J_vector SelfEQ_M_vector];
SelfFR_IJM_vector = [SelfFR_I_vector SelfFR_J_vector SelfFR_M_vector];
CrossEQ_IJKM_vector = [CrossEQ_I_vector CrossEQ_J_vector CrossEQ_K_vector CrossEQ_M_vector];
CrossFR_IJKM_vector = [CrossFR_I_vector CrossFR_J_vector CrossFR_K_vector CrossFR_M_vector];

Ks_eq = ndSparse.build(SelfEQ_IJM_vector,Ks_eq_vector,[size(probes,1),max([N_Self 1]),N_methods]);
dHs_eq = ndSparse.build(SelfEQ_IJM_vector,dHs_eq_vector,[size(probes,1),max([N_Self 1]),N_methods]);
dSs_eq = ndSparse.build(SelfEQ_IJM_vector,dSs_eq_vector,[size(probes,1),max([N_Self 1]),N_methods]);
dCps_eq = ndSparse.build(SelfEQ_IJM_vector,dCps_eq_vector,[size(probes,1),max([N_Self 1]),N_methods]);
dHs_f = ndSparse.build(SelfFR_IJM_vector,dHs_f_vector,[size(probes,1),max([N_Self 1]),N_methods2]);
dSs_f = ndSparse.build(SelfFR_IJM_vector,dSs_f_vector,[size(probes,1),max([N_Self 1]),N_methods2]);
dHs_r = ndSparse.build(SelfFR_IJM_vector,dHs_r_vector,[size(probes,1),max([N_Self 1]),N_methods2]);
dSs_r = ndSparse.build(SelfFR_IJM_vector,dSs_r_vector,[size(probes,1),max([N_Self 1]),N_methods2]);
Kd_eq = ndSparse.build(CrossEQ_IJKM_vector,dKd_eq_vector,[size(probes,1),size(probes,1),max([max(N_Cross(:)) 1]),N_methods]);
dHd_eq = ndSparse.build(CrossEQ_IJKM_vector,dHd_eq_vector,[size(probes,1),size(probes,1),max([max(N_Cross(:)) 1]),N_methods]);
dSd_eq = ndSparse.build(CrossEQ_IJKM_vector,dSd_eq_vector,[size(probes,1),size(probes,1),max([max(N_Cross(:)) 1]),N_methods]);
dCpd_eq = ndSparse.build(CrossEQ_IJKM_vector,dCpd_eq_vector,[size(probes,1),size(probes,1),max([max(N_Cross(:)) 1]),N_methods]);
dHd_f = ndSparse.build(CrossFR_IJKM_vector,dHd_f_vector,[size(probes,1),size(probes,1),max([max(N_Cross(:)) 1]),N_methods2]);
dSd_f = ndSparse.build(CrossFR_IJKM_vector,dSd_f_vector,[size(probes,1),size(probes,1),max([max(N_Cross(:)) 1]),N_methods2]);
dHd_r = ndSparse.build(CrossFR_IJKM_vector,dHd_r_vector,[size(probes,1),size(probes,1),max([max(N_Cross(:)) 1]),N_methods2]);
dSd_r = ndSparse.build(CrossFR_IJKM_vector,dSd_r_vector,[size(probes,1),size(probes,1),max([max(N_Cross(:)) 1]),N_methods2]);

ElapsedEquilibriumRateEvaluation = toc(start2);
end