function [Ks,Kd,ElapsedTimeForSecondaryStructureSeq,ElapsedEquilibriumRateEvaluation] = A_JH_GenerateSecondaryStructureInfo(probes,FinalProbeSet,settings)
%Jason Hughes Software to get self binding and cross binding energy for
%probes used together in a probe set

kb = 0.001987204259;%boltzman constant
scr_mat = [-1,-1,-1,1;-1,-1,1,-1;-1,1,-1,-1;1,-1,-1,-1;];  
RemoveMisMatches = settings.RemoveMisMatches;
SaltConcentration = settings.SaltConcentration;
T_hybrid = settings.HybridizationTemperature;
RV = @(x) (seqrcomplement(x));
for i=1:size(probes,1)
   pi_seq{i} = RV(probes{i,2});
end
start = tic; 

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
N_Self = cellfun(@(x) length(x), SelfSeqParsed);
clear row Flip_Identity Self_Flip_Identity
    %this is slow 
	
%Get Cross-Dimerization by finding sequence alignment between probe, and its complement                    
%Generates Cross-Dimerization for different configurations of two probes binding
for i=FinalProbeSet
    for j=FinalProbeSet
        if (i ~= j)
            tempAlign1 = swalignMod(pi_seq{i},pi_seq{j},'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
            TopCrossDimerSeqParsed1{i}{j} = tempAlign1(1,:);
            BotCrossDimerSeqParsed1{i}{j} = tempAlign1(3,:);
            TopCrossDimerSeqParsed1{j}{i} = tempAlign1(1,:);
            BotCrossDimerSeqParsed1{j}{i} = tempAlign1(3,:);
            tempAlign2 = swalignMod(pi_seq{i},reverse(pi_seq{j}),'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
            TopCrossDimerSeqParsed2{i}{j} = tempAlign2(1,:);
            BotCrossDimerSeqParsed2{i}{j} = tempAlign2(3,:);
            TopCrossDimerSeqParsed2{j}{i} = tempAlign2(1,:);
            BotCrossDimerSeqParsed2{j}{i} = tempAlign2(3,:);
            tempAlign3 = swalignMod(reverse(pi_seq{i}),pi_seq{j},'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
            TopCrossDimerSeqParsed3{i}{j} = tempAlign3(1,:);
            BotCrossDimerSeqParsed3{i}{j} = tempAlign3(3,:);
            TopCrossDimerSeqParsed3{j}{i} = tempAlign3(1,:);
            BotCrossDimerSeqParsed3{j}{i} = tempAlign3(3,:);
            tempAlign4 = swalignMod(reverse(reverse(pi_seq{i})),reverse(pi_seq{j}),'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
            TopCrossDimerSeqParsed4{i}{j} = tempAlign4(1,:);
            BotCrossDimerSeqParsed4{i}{j} = tempAlign4(3,:);
            TopCrossDimerSeqParsed4{j}{i} = tempAlign4(1,:);
            BotCrossDimerSeqParsed4{j}{i} = tempAlign4(3,:); 
            clear tempAlign*
        else
            tempAlign1 = swalignMod(pi_seq{i},reverse(pi_seq{j}),'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
            TopCrossDimerSeqParsed1{i}{j} = tempAlign1(1,:);
            BotCrossDimerSeqParsed1{i}{j} = tempAlign1(3,:);
            tempAlign2 = swalignMod(pi_seq{i},pi_seq{j},'SCORINGMATRIX',scr_mat,'GAPOPEN',5,'ALPHA','NT');
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
%Store CrossDimerization Pairs into dictionary
count = 1;
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
N_Cross = cellfun(@(x) length(x), CrossDimerSeqParsed);

start2 = tic;
SG = cell(1,size(probes,1));
CG = cell(size(probes,1),size(probes,1));
Ks = ndSparse.build([size(probes,1),1],0);
Kd = ndSparse.build([size(probes,1),size(probes,1)],0);
%Compute binding affinity for self and cross-dimers
for i=FinalProbeSet
    v = find(FinalProbeSet==i);
    for j=1:N_Self(v)
        try
            SG{i}(j) = F_DeltaGibson(strrep(SelfSeqParsed{v}{j},'-','N'),strrep(SelfSeqParsed{v}{j},'-','N'),SaltConcentration,T_hybrid,RemoveMisMatches);
        catch
            SG{i}(j) = Inf;
        end
    end
    if (~isempty(SG{i}))
        Ks(i) = sum(exp(-SG{i}/(kb*(T_hybrid+273.15))));%Keq rates Self-Dimerization
    else
        Ks(i) = 0;
    end
    for j=FinalProbeSet
        w = find(FinalProbeSet==j);
        for k = 1:N_Cross(v,w)
            try
                if (~isempty(CrossDimerSeqParsed{v,w}{k,1})&&~isempty(CrossDimerSeqParsed{v,w}{k,2}))
                    CG{i,j}{k} = F_DeltaGibson(strrep(CrossDimerSeqParsed{v,w}{k,1},'-','N'),reverse(strrep(CrossDimerSeqParsed{v,w}{k,2},'-','N')),SaltConcentration,T_hybrid,RemoveMisMatches);
                else
                    CG{i,j}{k} = Inf;    
                end
            catch
                CG{i,j}{k} = Inf;
            end
        end
        if (~isempty(CG{i,j}))
            Kd(i,j) = sum(exp(-cell2mat(CG{i,j})/(kb*(T_hybrid+273.15)))); %K rates Cross-Dimerization
            Kd(j,i) = Kd(i,j);
        else
            Kd(i,j) = 0;
            Kd(j,i) = 0;
        end
   end
end 
ElapsedEquilibriumRateEvaluation = toc(start2);
end