function [Ks_eq,dHs_eq,dSs_eq,dHs_f,dSs_f,dHs_r,dSs_r,dCps_eq,Kd_eq,dHd_eq,dSd_eq,dHd_f,dSd_f,dHd_r,dSd_r,dCpd_eq,N_Self,N_Cross,ElapsedTimeForSecondaryStructureSeq,ElapsedEquilibriumRateEvaluation] = A_JH_GenerateSecondaryStructureInfo_V3(probes,FinalProbeSet,settings,overrides)
%Jason Hughes Software to get self binding and cross binding energy for probes used together in a probe set

if nargin < 4
    overrides = struct();
    overrides.includeSelf = true;
    overrides.includeCross = true;
end
gapOpen = 5;
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
if overrides.includeSelf
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
else
    N_Self = 0;
end
if overrides.includeCross
    TopCrossDimerSeqParsed1 = cell(1,length(FinalProbeSet));
    TopCrossDimerSeqParsed2 = cell(1,length(FinalProbeSet));
    TopCrossDimerSeqParsed3 = cell(1,length(FinalProbeSet));
    TopCrossDimerSeqParsed4 = cell(1,length(FinalProbeSet));
    BotCrossDimerSeqParsed1 = cell(1,length(FinalProbeSet));
    BotCrossDimerSeqParsed2 = cell(1,length(FinalProbeSet));
    BotCrossDimerSeqParsed3 = cell(1,length(FinalProbeSet));
    BotCrossDimerSeqParsed4 = cell(1,length(FinalProbeSet));
    cross_pair_combos = [nchoosek(1:length(FinalProbeSet),2);repmat(1:length(FinalProbeSet),[2 1])'];
    cross_pair_combos_with_flip =  [nchoosek(1:length(FinalProbeSet),2);fliplr(nchoosek(1:length(FinalProbeSet),2));repmat(1:length(FinalProbeSet),[2 1])'];
    cross_pair_combos_List = parallel.pool.Constant(cross_pair_combos);
    tempAlign1 = cell(1,size(cross_pair_combos,1));
    tempAlign2 = cell(1,size(cross_pair_combos,1));
    tempAlign3 = cell(1,size(cross_pair_combos,1));
    tempAlign4 = cell(1,size(cross_pair_combos,1));
    FinalProbeSet_List = parallel.pool.Constant(FinalProbeSet);
    pi_seq_List = parallel.pool.Constant(pi_seq);
    parfor nn = 1:size(cross_pair_combos,1)
        [~,tempAlign1{nn}] = swalign(pi_seq_List.Value{FinalProbeSet_List.Value(cross_pair_combos_List.Value(nn,1))},pi_seq_List.Value{FinalProbeSet_List.Value(cross_pair_combos_List.Value(nn,2))},'SCORINGMATRIX',scr_mat,'GAPOPEN',gapOpen,'ALPHA','NT');
        [~,tempAlign2{nn}] = swalign(pi_seq_List.Value{FinalProbeSet_List.Value(cross_pair_combos(nn,1))},reverse(pi_seq_List.Value{FinalProbeSet_List.Value(cross_pair_combos_List.Value(nn,2))}),'SCORINGMATRIX',scr_mat,'GAPOPEN',gapOpen,'ALPHA','NT');
        [~,tempAlign3{nn}] = swalign(reverse(pi_seq_List.Value{FinalProbeSet_List.Value(cross_pair_combos_List.Value(nn,1))}),pi_seq_List.Value{FinalProbeSet_List.Value(cross_pair_combos_List.Value(nn,2))},'SCORINGMATRIX',scr_mat,'GAPOPEN',gapOpen,'ALPHA','NT');
        [~,tempAlign4{nn}] = swalign(reverse(pi_seq_List.Value{FinalProbeSet_List.Value(cross_pair_combos_List.Value(nn,1))}),reverse(pi_seq_List.Value{FinalProbeSet_List.Value(cross_pair_combos_List.Value(nn,2))}),'SCORINGMATRIX',scr_mat,'GAPOPEN',gapOpen,'ALPHA','NT');
    end
    parfor_tempAlignIntermediateC = cell(size(cross_pair_combos,1), 1);
    parfor ii = 1:size(cross_pair_combos,1)
        if (~isequal(cross_pair_combos_List.Value(ii,1),cross_pair_combos_List.Value(ii,2)))
            Current_Z = {tempAlign1{ii}(1,:);tempAlign1{ii}(3,:); tempAlign2{ii}(1,:);tempAlign2{ii}(3,:); tempAlign3{ii}(1,:);tempAlign3{ii}(3,:);tempAlign4{ii}(1,:);tempAlign4{ii}(3,:)};
            parfor_tempAlignIntermediateC{ii} = {[cross_pair_combos_List.Value(ii,1); cross_pair_combos_List.Value(ii,2)],[cross_pair_combos_List.Value(ii,2); cross_pair_combos_List.Value(ii,1)],[Current_Z Current_Z]'};
        else
            Current_Z = {tempAlign2{ii}(1,:);tempAlign2{ii}(3,:);tempAlign1{ii}(1,:);tempAlign1{ii}(3,:);'';'';'';''};
             parfor_tempAlignIntermediateC{ii} = {cross_pair_combos_List.Value(ii,1),cross_pair_combos_List.Value(ii,2),(Current_Z)'};
        end
    end
    clear tempAlign*
    positionX = cell2mat(cellfun(@(x) x{1}, parfor_tempAlignIntermediateC, 'Un', 0));
    positionY = cell2mat(cellfun(@(x) x{2}, parfor_tempAlignIntermediateC, 'Un', 0));
    positionZ = CATnWrapper(cellfun(@(x) x{3}, parfor_tempAlignIntermediateC, 'Un', 0),1);
    clear parfor_tempAlignIntermediate*
    positionY_List = parallel.pool.Constant(positionY);
    positionZ_List = parallel.pool.Constant(positionZ);
    parfor ii = 1:length(FinalProbeSet)
        TopCrossDimerSeqParsed1{ii}(positionY_List.Value(positionX==ii)) = positionZ_List.Value((positionX==ii),1);
        BotCrossDimerSeqParsed1{ii}(positionY_List.Value(positionX==ii)) = positionZ_List.Value((positionX==ii),2);
        TopCrossDimerSeqParsed2{ii}(positionY_List.Value(positionX==ii)) = positionZ_List.Value((positionX==ii),3);
        BotCrossDimerSeqParsed2{ii}(positionY_List.Value(positionX==ii)) = positionZ_List.Value((positionX==ii),4);
        TopCrossDimerSeqParsed3{ii}(positionY_List.Value(positionX==ii)) = positionZ_List.Value((positionX==ii),5);
        BotCrossDimerSeqParsed3{ii}(positionY_List.Value(positionX==ii)) = positionZ_List.Value((positionX==ii),6);
        TopCrossDimerSeqParsed4{ii}(positionY_List.Value(positionX==ii)) = positionZ_List.Value((positionX==ii),7);
        BotCrossDimerSeqParsed4{ii}(positionY_List.Value(positionX==ii)) = positionZ_List.Value((positionX==ii),8);
    end
    CrossDictionary = ...
        [CATnWrapper(arrayfun(@(nn) TopCrossDimerSeqParsed1{nn}(~cellfun(@isempty,TopCrossDimerSeqParsed1{nn})), 1:length(FinalProbeSet),'Un',0),2) ...
        CATnWrapper(arrayfun(@(nn) TopCrossDimerSeqParsed2{nn}(~cellfun(@isempty,TopCrossDimerSeqParsed2{nn})), 1:length(FinalProbeSet),'Un',0),2) ...
        CATnWrapper(arrayfun(@(nn) TopCrossDimerSeqParsed3{nn}(~cellfun(@isempty,TopCrossDimerSeqParsed3{nn})), 1:length(FinalProbeSet),'Un',0),2) ...
        CATnWrapper(arrayfun(@(nn) TopCrossDimerSeqParsed4{nn}(~cellfun(@isempty,TopCrossDimerSeqParsed4{nn})), 1:length(FinalProbeSet),'Un',0),2) ...
        CATnWrapper(arrayfun(@(nn) BotCrossDimerSeqParsed1{nn}(~cellfun(@isempty,BotCrossDimerSeqParsed1{nn})), 1:length(FinalProbeSet),'Un',0),2) ...
        CATnWrapper(arrayfun(@(nn) BotCrossDimerSeqParsed2{nn}(~cellfun(@isempty,BotCrossDimerSeqParsed2{nn})), 1:length(FinalProbeSet),'Un',0),2) ...
        CATnWrapper(arrayfun(@(nn) BotCrossDimerSeqParsed3{nn}(~cellfun(@isempty,BotCrossDimerSeqParsed3{nn})), 1:length(FinalProbeSet),'Un',0),2) ...
        CATnWrapper(arrayfun(@(nn) BotCrossDimerSeqParsed4{nn}(~cellfun(@isempty,BotCrossDimerSeqParsed4{nn})), 1:length(FinalProbeSet),'Un',0),2)];
    CrossDimerDictionary.Names = unique(CrossDictionary);
    %% Find Unique Pairs of Hairpin and Self Hybridization Sequences
    tempPair1_Exists = arrayfun(@(nn) sum(strcmp(CrossDimerDictionary.Names,TopCrossDimerSeqParsed1{cross_pair_combos_with_flip(nn,1)}{cross_pair_combos_with_flip(nn,2)})),1:size(cross_pair_combos_with_flip,1));
    tempPair2_Exists = arrayfun(@(nn) sum(strcmp(CrossDimerDictionary.Names,TopCrossDimerSeqParsed2{cross_pair_combos_with_flip(nn,1)}{cross_pair_combos_with_flip(nn,2)})),1:size(cross_pair_combos_with_flip,1));
    tempPair3_Exists = arrayfun(@(nn) sum(strcmp(CrossDimerDictionary.Names,TopCrossDimerSeqParsed3{cross_pair_combos_with_flip(nn,1)}{cross_pair_combos_with_flip(nn,2)})),1:size(cross_pair_combos_with_flip,1));
    tempPair4_Exists = arrayfun(@(nn) sum(strcmp(CrossDimerDictionary.Names,TopCrossDimerSeqParsed4{cross_pair_combos_with_flip(nn,1)}{cross_pair_combos_with_flip(nn,2)})),1:size(cross_pair_combos_with_flip,1));
    DictionaryPairs = ...
        [CATnWrapper(arrayfun(@(nn) [find(strcmp(CrossDimerDictionary.Names,TopCrossDimerSeqParsed1{cross_pair_combos_with_flip(nn,1)}{cross_pair_combos_with_flip(nn,2)})) ...
        find(strcmp(CrossDimerDictionary.Names,BotCrossDimerSeqParsed1{cross_pair_combos_with_flip(nn,1)}{cross_pair_combos_with_flip(nn,2)})) cross_pair_combos_with_flip(nn,1) cross_pair_combos_with_flip(nn,2)],find(tempPair1_Exists),'Un',0),1);...
        CATnWrapper(arrayfun(@(nn) [find(strcmp(CrossDimerDictionary.Names,TopCrossDimerSeqParsed2{cross_pair_combos_with_flip(nn,1)}{cross_pair_combos_with_flip(nn,2)})) ...
        find(strcmp(CrossDimerDictionary.Names,BotCrossDimerSeqParsed2{cross_pair_combos_with_flip(nn,1)}{cross_pair_combos_with_flip(nn,2)})) cross_pair_combos_with_flip(nn,1) cross_pair_combos_with_flip(nn,2)],find(tempPair2_Exists),'Un',0),1);...
        CATnWrapper(arrayfun(@(nn) [find(strcmp(CrossDimerDictionary.Names,TopCrossDimerSeqParsed3{cross_pair_combos_with_flip(nn,1)}{cross_pair_combos_with_flip(nn,2)})) ...
        find(strcmp(CrossDimerDictionary.Names,BotCrossDimerSeqParsed3{cross_pair_combos_with_flip(nn,1)}{cross_pair_combos_with_flip(nn,2)})) cross_pair_combos_with_flip(nn,1) cross_pair_combos_with_flip(nn,2)],find(tempPair3_Exists),'Un',0),1);...
        CATnWrapper(arrayfun(@(nn) [find(strcmp(CrossDimerDictionary.Names,TopCrossDimerSeqParsed4{cross_pair_combos_with_flip(nn,1)}{cross_pair_combos_with_flip(nn,2)})) ...
        find(strcmp(CrossDimerDictionary.Names,BotCrossDimerSeqParsed4{cross_pair_combos_with_flip(nn,1)}{cross_pair_combos_with_flip(nn,2)})) cross_pair_combos_with_flip(nn,1) cross_pair_combos_with_flip(nn,2)],find(tempPair4_Exists),'Un',0),1)];
    clear tempPair*
    UniquePairs = unique(DictionaryPairs,'rows');
    clear DictionaryPairs
    %% Map Unique Pairs back to probes pairs
    CrossDimerSeqParsed = cell(length(FinalProbeSet),length(FinalProbeSet));
    Row_Vector =  CATnWrapper(arrayfun(@(nn)  find((UniquePairs(:,3)==cross_pair_combos_with_flip(nn,1)).*(UniquePairs(:,4)==cross_pair_combos_with_flip(nn,2))),1:size(cross_pair_combos_with_flip,1),'Un',0),1);
    K_Vector =  CATnWrapper(arrayfun(@(nn)  reshape(1:size(find((UniquePairs(:,3)==cross_pair_combos_with_flip(nn,1)).*(UniquePairs(:,4)==cross_pair_combos_with_flip(nn,2))),1),[],1),1:size(cross_pair_combos_with_flip,1),'Un',0),1);
    U_Vector =  CATnWrapper(arrayfun(@(nn)  cross_pair_combos_with_flip(nn,1)*ones(size(find((UniquePairs(:,3)==cross_pair_combos_with_flip(nn,1)).*(UniquePairs(:,4)==cross_pair_combos_with_flip(nn,2))))),1:size(cross_pair_combos_with_flip,1),'Un',0),1);
    V_Vector = CATnWrapper(arrayfun(@(nn)  cross_pair_combos_with_flip(nn,2)*ones(size(find((UniquePairs(:,3)==cross_pair_combos_with_flip(nn,1)).*(UniquePairs(:,4)==cross_pair_combos_with_flip(nn,2))))),1:size(cross_pair_combos_with_flip,1),'Un',0),1);
    K_Vector_List = parallel.pool.Constant(K_Vector);
    Row_Vector_List = parallel.pool.Constant(Row_Vector);
    UniquePairs_List = parallel.pool.Constant(UniquePairs);
    CrossDimerDictionaryNames_List = parallel.pool.Constant(CrossDimerDictionary.Names);
    L = length(FinalProbeSet);
    parfor ii = 1:length(FinalProbeSet)
        for ij = 1:L
            CrossDimerSeqParsed{ii,ij}(K_Vector_List.Value((U_Vector==ii).*(V_Vector==ij)==1),1) = CrossDimerDictionaryNames_List.Value(UniquePairs_List.Value(Row_Vector_List.Value((U_Vector==ii).*(V_Vector==ij)==1),1));
            CrossDimerSeqParsed{ii,ij}(K_Vector_List.Value((U_Vector==ii).*(V_Vector==ij)==1),2) = CrossDimerDictionaryNames_List.Value(UniquePairs_List.Value(Row_Vector_List.Value((U_Vector==ii).*(V_Vector==ij)==1),2));
        end
    end
    clear Row_Vector* K_Vector* U_Vector* V_Vector* CrossDimerDictionary*
    %% Remove Pairs of Flips (Redundant matches, that look different but are a different pair both flipped)
    parfor_tempPar = cell(size(cross_pair_combos,1), 1);
    cross_pair_combos_with_flip_List = parallel.pool.Constant(cross_pair_combos_with_flip);
    CrossDimerSeqParsed_List = parallel.pool.Constant(CrossDimerSeqParsed);
    parfor nn = 1:size(cross_pair_combos_with_flip,1)
        u = cross_pair_combos_with_flip_List.Value(nn,1);
        v = cross_pair_combos_with_flip_List.Value(nn,2);
        [row1,~] = find(triu(cell2mat(arrayfun(@(y) strcmp(flip(CrossDimerSeqParsed_List.Value{u,v}{y,1}),CrossDimerSeqParsed_List.Value{u,v}(:,1)),1:size(CrossDimerSeqParsed_List.Value{u,v},1),'UniformOutput',false))));
        [row2,~] = find(triu(cell2mat(arrayfun(@(y) strcmp(flip(CrossDimerSeqParsed_List.Value{u,v}{y,2}),CrossDimerSeqParsed_List.Value{u,v}(:,2)),1:size(CrossDimerSeqParsed_List.Value{u,v},1),'UniformOutput',false))));
        parfor_tempPar{nn} = {u,v, row1, row2};
    end
    clear CrossDimerSeqParsed_List cross_pair_combos_with_flip_List
    U_Vector = cross_pair_combos_with_flip(:,1);
    V_Vector = cross_pair_combos_with_flip(:,2);
    parfor_tempPar_List = parallel.pool.Constant(parfor_tempPar);
    clear parfor_tempPar
    parfor ii = 1:length(FinalProbeSet)
        for ij = 1:L
            CrossDimerSeqParsed{ii,ij}( cell2mat(arrayfun(@(nn) parfor_tempPar_List.Value{nn}{3},find((U_Vector == ii).*(V_Vector==ij)),'Un',0))    ,1) = arrayfun(@(x) [],find((U_Vector == ii).*(V_Vector==ij)),'Un',0);
            CrossDimerSeqParsed{ii,ij}( cell2mat(arrayfun(@(nn) parfor_tempPar_List.Value{nn}{4},find((U_Vector == ii).*(V_Vector==ij)),'Un',0))  ,  2) = arrayfun(@(x) [],find((U_Vector == ii).*(V_Vector==ij)),'Un',0);
        end
    end
    clear row Flip_Identity1 Flip_Identity2 Cross_Flip_Identity1 Cross_Flip_Identity2 row1 row2 parfor_tempPar_List
    crossNotExist = cellfun(@isempty,CrossDimerSeqParsed);
    [crossId1,crossId2] = find(crossNotExist==1);
    for i = 1:length(crossId1)
        CrossDimerSeqParsed{crossId1(i),crossId2(i)}=cell(0,2);
    end
    CrossDimerSeqParsed = cellfun(@(x) x(sum(cellfun(@isempty,x),2)'==0,1:2), CrossDimerSeqParsed,'Un',0);
    N_Cross = cellfun(@(x) size(x,1), CrossDimerSeqParsed);
else
    N_Cross = 0;
end
ElapsedTimeForSecondaryStructureSeq = toc(start);

start2 = tic;
%% Initialize matrix for storing binding affinity calculations
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
Self_V_List = [];
Self_J_List = [];
Cross_V_List  = [];
Cross_W_List  = [];
Cross_K_List = [];
if overrides.includeSelf
    Self_V_List = cell2mat(cellfun(@(V,J) V*ones(1,J),num2cell(1:length(SelfSeqParsed)),num2cell(cellfun(@length,SelfSeqParsed)),'Un',0))';
    Self_J_List = cell2mat(cellfun(@(V,J) 1:J,num2cell(1:length(SelfSeqParsed)),num2cell(cellfun(@length,SelfSeqParsed)),'Un',0))';
    Self_SeqInput = arrayfun(@(x) strrep(SelfSeqParsed{Self_V_List(x)}{Self_J_List(x)},'-','N'),1:length(Self_V_List),'Un',0)';
end
if overrides.includeCross
    Cross_V_List = cell2mat(reshape(cellfun(@(V,J) V*ones(1,J),num2cell(meshgrid(1:size(CrossDimerSeqParsed,1))),cellfun(@(z) size(z,1),CrossDimerSeqParsed,'Un',0),'Un',0),1,[]))';
    Cross_W_List = cell2mat(reshape(cellfun(@(W,J) W*ones(1,J),num2cell(meshgrid(1:size(CrossDimerSeqParsed,1))'),cellfun(@(z) size(z,1),CrossDimerSeqParsed,'Un',0),'Un',0),1,[]))';
    Cross_K_List = cell2mat(reshape(cellfun(@(V,J) 1:J,num2cell(meshgrid(1:size(CrossDimerSeqParsed,1))),cellfun(@(z) size(z,1),CrossDimerSeqParsed,'Un',0),'Un',0),1,[]))';
    Cross_SeqInput1 = arrayfun(@(x) strrep(CrossDimerSeqParsed{Cross_V_List(x),Cross_W_List(x)}{Cross_K_List(x),1},'-','N'),1:length(Cross_V_List),'Un',0)';
    Cross_SeqInput2 = arrayfun(@(x) reverse(strrep(CrossDimerSeqParsed{Cross_V_List(x),Cross_W_List(x)}{Cross_K_List(x),2},'-','N')),1:length(Cross_V_List),'Un',0)';
end
if overrides.includeSelf && overrides.includeCross
    unique_secondary_structure_binding_seqs = unique([Self_SeqInput; Cross_SeqInput1; Cross_SeqInput2;...
        cellfun(@seqreverse,[Self_SeqInput; Cross_SeqInput1; Cross_SeqInput2],'Un',0)]);
elseif overrides.includeSelf && ~overrides.includeCross
    unique_secondary_structure_binding_seqs = unique([Self_SeqInput;cellfun(@seqreverse,Self_SeqInput,'Un',0)]);
elseif ~overrides.includeSelf && overrides.includeCross
    unique_secondary_structure_binding_seqs = unique([Cross_SeqInput1; Cross_SeqInput2;...
        cellfun(@seqreverse,[Cross_SeqInput1; Cross_SeqInput2],'Un',0)]);
end
if overrides.includeSelf || overrides.includeCross
    binding_seqs_dictionary = dictionary(convertCharsToStrings(unique_secondary_structure_binding_seqs)',1:length(unique_secondary_structure_binding_seqs));
    inverse_binding_seqs_dictionary = dictionary(1:length(unique_secondary_structure_binding_seqs),convertCharsToStrings(unique_secondary_structure_binding_seqs)');
end

Is_Self_Over_Cross = [1*ones(size(Self_V_List)); 0*ones(size(Cross_V_List))];
Secondary_Structure_VV_List = [Self_V_List; Cross_V_List];
Secondary_Structure_VW_List = [Self_V_List; Cross_W_List];
Secondary_Structure_JK_List = [Self_J_List; Cross_K_List];

if overrides.includeSelf && overrides.includeCross
    Secondary_Structure_unique_Seq1_id_List = [binding_seqs_dictionary(convertCharsToStrings(Self_SeqInput)); binding_seqs_dictionary(convertCharsToStrings(Cross_SeqInput1))];
    Secondary_Structure_unique_Seq2_id_List = [binding_seqs_dictionary(convertCharsToStrings(Self_SeqInput)); binding_seqs_dictionary(convertCharsToStrings(Cross_SeqInput2))];
elseif overrides.includeSelf && ~overrides.includeCross
    Secondary_Structure_unique_Seq1_id_List = [binding_seqs_dictionary(convertCharsToStrings(Self_SeqInput))];
    Secondary_Structure_unique_Seq2_id_List = [binding_seqs_dictionary(convertCharsToStrings(Self_SeqInput))];
elseif ~overrides.includeSelf && overrides.includeCross
    Secondary_Structure_unique_Seq1_id_List = [binding_seqs_dictionary(convertCharsToStrings(Cross_SeqInput1))];
    Secondary_Structure_unique_Seq2_id_List = [binding_seqs_dictionary(convertCharsToStrings(Cross_SeqInput2))];
else
    Secondary_Structure_unique_Seq1_id_List = [];
    Secondary_Structure_unique_Seq2_id_List = [];
end

seq_pair_ids_nonunique = [Secondary_Structure_unique_Seq1_id_List Secondary_Structure_unique_Seq2_id_List];
[unique_pair_ids,~,matched_unique_location] = unique(seq_pair_ids_nonunique,'rows');
unique_secondary_structure_pair_to_nonunique_entries = [Is_Self_Over_Cross Secondary_Structure_VV_List Secondary_Structure_VW_List ...
    Secondary_Structure_JK_List Secondary_Structure_unique_Seq1_id_List Secondary_Structure_unique_Seq2_id_List matched_unique_location];
unique_ordered_binding_paired_input_sequences = [];
if overrides.includeSelf || overrides.includeCross
    unique_ordered_binding_paired_input_sequences = inverse_binding_seqs_dictionary(unique_pair_ids);
end
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