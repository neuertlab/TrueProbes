function [dHeq, dSeq, dGeq, dHf, dSf, dGf, dHr, dSr, dGr, dCpeq, Tm] = F_DeltaGibson_V3(seq1,seq2,C,T,P)
% This function computes terms for DeltaGibson Free energy changes for any two pairs of DNA sequences binding.
% This function compute DeltaGibson for both binding equilibrium and forward/reverse transition-state rates.
% Additionally this compute Tm for DNA sequence binding.
% Additionally includeded are meso-scale mechanistic model which explicity includes terms for mismatches.


N_models = 8;
dGeq = zeros(N_models,1);
dHeq = zeros(N_models,1);
dSeq = zeros(N_models,1);
Tm = zeros(N_models-1,1);
dCpeq = zeros(N_models,1);
dGf = zeros(3,1);
dHf = zeros(3,1);
dSf = zeros(3,1);
dGr = zeros(3,1);
dHr = zeros(3,1);
dSr = zeros(3,1);
format('long');
%upper strand is RNA with 5->3 and lower strand DNA (3'->5')
%Seq1 5->3
%Seq2 5->3
%DNA/DNA;
%RNA/DNA;
%RNA/RNA;
Np = 400;
T_vector = 25:5:75;
kb = 8.617333262145*10^-5;%eV/K
b = 1/(kb*(273.15));%T
%Tm = dH/dS
%Keq = exp(-dH/RT(1-T/Tm))
%http://matlab.cheme.cmu.edu/2011/12/20/the-gibbs-free-energy-of-a-reacting-mixture-and-the-equilibrium-composition/#13
%dG = dH-TdS
%at melting temp dG=0
%dH = TmdS
%dG = (Tm-T)dS
%fraction bound
%theta = 1/(1+exp(dG/RT))
%theta = 1/(1+exp(-dTdS/RT))
RV = @(x) (seqrcomplement(x));
SaltConcentration = C;
Temperature = T;
PrimerConc = P;
if (strcmp(seq1,seq2))%for when they are the same
    seq1 = strip(seq1,'N');
    seq2 = strip(seq2,'N');
else
    seq2i = seq2;
    seq2 = strip(seq2,'N');
    st = strfind(seq2i,seq2);
    seq1 = seq1(st:st+length(seq2)-1);
end
seq2c = reverse(seq2);

%Removes mismatch/Drop out effect of mismatches
%Cherepinsky, V., Hashmi, G., & Mishra, B. (2010). Competitive hybridization models. Physical Review E, 82(5), 051914.
%Find Mismatches
IsMisMatched = [];
for i=1:length(seq1)
    IsMisMatched(i) = ~strcmpi(Match(seq1(i)),seq2c(i));
end
allLocs = 1:length(seq1);
MisMatchLoc = find(IsMisMatched);
MatchedLoc = setdiff(allLocs,MisMatchLoc);
seq1Matched = seq1(MatchedLoc);
seq2Matched = seq2c(MatchedLoc);
%delta Sd wrong?  correct sometimes also incorrect double
%wrong for more gcs?


if (~isempty(MatchedLoc))
    seq1_Data = F_NearestNeighborWrapper_V2(RV(seq1Matched),SaltConcentration,Temperature,PrimerConc);
    seq2_Data = F_NearestNeighborWrapper_V2(RV(seq2Matched),SaltConcentration,Temperature,PrimerConc);
    dHseq1 = seq1_Data.Thermo(:,1);
    dHseq2 = seq2_Data.Thermo(:,1);
    dSseq1 = seq1_Data.Thermo(:,2)./1000;
    dSseq2 = seq2_Data.Thermo(:,2)./1000;
    dGseq1 = seq1_Data.Thermo(:,3);
    dGseq2 = seq2_Data.Thermo(:,3);
    seq1_TData = F_NearestNeighborTransitionWrapper_V2(RV(seq1Matched),SaltConcentration,Temperature,PrimerConc);
    seq2_TData = F_NearestNeighborTransitionWrapper_V2(RV(seq2Matched),SaltConcentration,Temperature,PrimerConc); 
    dHseq1_T = seq1_TData.Thermo(:,1);
    dHseq2_T = seq2_TData.Thermo(:,1);
    dSseq1_T = seq1_TData.Thermo(:,2)./1000;
    dSseq2_T = seq2_TData.Thermo(:,2)./1000;
    dGseq1_T = seq1_TData.Thermo(:,3);
    dGseq2_T = seq2_TData.Thermo(:,3);
    dHeq(1:N_models-1) = 0.5*(dHseq1 + dHseq2);   
    dSeq(1:N_models-1) = 0.5*(dSseq1 + dSseq2);   
    dGeq(1:N_models-1) = 0.5*(dGseq1 + dGseq2);   
    dHt = 0.5*(dHseq1_T + dHseq2_T);   
    dSt = 0.5*(dSseq1_T + dSseq2_T);   
    dGt = 0.5*(dGseq1_T + dGseq2_T);   
    % association units has moles diff molecules
    %fL fM fH rL rM rH eqL eqM eqH
    dHf = dHt(1:3);
    dHr = dHt(4:6);
    dSf = dSt(1:3);
    dSr = dSt(4:6);
    dGf = dGt(1:3);
    dGr = dGt(4:6);    
    Tm = 0.5*(seq1_Data.Tm+seq2_Data.Tm);
else
    dHeq(1:N_models-1) = Inf*ones(N_models-1,1);
    dSeq(1:N_models-1) = Inf*ones(N_models-1,1);
    dGeq(1:N_models-1) = Inf*ones(N_models-1,1);
    Tm(1:N_models-1) = NaN*ones(N_models-1,1);
    dHf = Inf*ones(3,1);
    dHr = Inf*ones(3,1);
    dSf = Inf*ones(3,1);
    dSr = Inf*ones(3,1);
    dGf = Inf*ones(3,1);
    dGr = Inf*ones(3,1); 
end  
if (isempty(seq1))
    dHeq = Inf*ones(N_models,1);
    dSeq = Inf*ones(N_models,1);
    dGeq = Inf*ones(N_models,1);
    Tm = NaN*ones(N_models-1,1);
    dHf = Inf*ones(3,1);
    dHr = Inf*ones(3,1);
    dSf = Inf*ones(3,1);
    dSr = Inf*ones(3,1);
    dGf = Inf*ones(3,1);
    dGr = Inf*ones(3,1); 
end
% 
% if (~isempty(seq1))  
%     %Calculate energy given mismatch
%     %Triplicate MisMatch Model (estimate with mismatch delta G)
%     %Mesoscopic model
%     %Oliveira, L. M., Long, A. S., Brown, T., Fox, K. R., & Weber, G. (2020). Melting temperature measurement and mesoscopic evaluation of single, double and triple DNA mismatches. Chemical science, 11(31), 8273-8287.
%     %Hydrogen bonds 
%     %V(y_i) = D_i*(e^(-y_i/l_i)-1)^2, y_i is distance between bases.
%     %D is potential depth, l is potential width
%     %stacking interaction
%     %W(y_i,y_i+1) = ki,i+1/2*(y_i^2-2y_i*y_i+1*cos(theta)+y_i+1^2)
%     %could add term in W to get k/2(1+p*exp(-alpha(yn+yn+1))(yn+1-yn)^2
%     %term for anharmonic potential to describe stacking energy,
%     %quantitative description of melting. 
%     %k is elastic constant, theta is smal angle 0.01 radians
%     %U(y_i,y_i+1) = W(y_i,y_i+1) + V(y_i)
%     %Z is partion function 
%     %K = exp(-BU),  Z = integral all y, K(y1,y2)*K(y2,y3)...K(yN,y1)
%     %average <y> of each base pair
%     %A = angstrum
%     %Hy for y<d, max value 30A, p=0.5, alpha = 0.35A^-1, 
%     %average yi between 1 and 4.5A
%     optsFull = detectImportOptions('DatabaseData/ContextDependent_MisMatchStackingPotentials.txt');
%     optsFull.DataLines=[1,Inf];
%     optsFull.VariableNamesLine = 0;
%     k_charMisMatch = readtable('DatabaseData/ContextDependent_MisMatchStackingPotentials.txt',optsFull);%mismatch and correct matches
%     k_charMisMatchUnparsedList = strjoin(table2cell(k_charMisMatch)," ");
%     extraSpaces = strfind(k_charMisMatchUnparsedList," -");
%     k_charMisMatchUnparsedList(extraSpaces)=[];
%     k_MisMatchParsedList = split(k_charMisMatchUnparsedList," ");
%     emptyParsedListIndx = find(cellfun(@isempty,k_MisMatchParsedList));
%     k_MisMatchParsedList(emptyParsedListIndx) = [];
%     N_Cases = length(k_MisMatchParsedList)/3;
%     KDependent_CaseName = {k_MisMatchParsedList{3*[1:N_Cases]-2}};
%     KDependent_NearestNeighbours = {k_MisMatchParsedList{3*[1:N_Cases]-1}};
%     KDependent_StackingPotentials = {k_MisMatchParsedList{3*[1:N_Cases]}};
%     dependentSeqLength = 1/2*(cellfun(@length,KDependent_NearestNeighbours)-3);
%     dependentTrimerIDs = find(dependentSeqLength==3);
%     dependentTetramerIDs = find(dependentSeqLength==4);
%     KDependentTrimer_NearestNeighbours =  KDependent_NearestNeighbours(dependentTrimerIDs);
%     KDependentTrimer_StackingPotentials = KDependent_StackingPotentials(dependentTrimerIDs);
%     KDependentTetramer_NearestNeighbours =  KDependent_NearestNeighbours(dependentTetramerIDs);
%     KDependentTetramer_StackingPotentials = KDependent_StackingPotentials(dependentTetramerIDs);
%     k_tableMisMatch = readtable('DatabaseData/ContextIndependent_MisMatchStackingPotentials.txt');%mismatch and correct matches
%     k_structMisMatch = table2struct(k_tableMisMatch);
%     KIndependent_NearestNeighbours = {k_structMisMatch.Var1 k_structMisMatch.Var3 k_structMisMatch.Var5 k_structMisMatch.Var7 ...
%         k_structMisMatch.Var9 k_structMisMatch.Var11 k_structMisMatch.Var13 k_structMisMatch.Var15};
%     KIndependent_StackingPotentials = {k_structMisMatch.Var2 k_structMisMatch.Var4 k_structMisMatch.Var6 k_structMisMatch.Var8 ...
%         k_structMisMatch.Var10 k_structMisMatch.Var12 k_structMisMatch.Var14 k_structMisMatch.Var16};
%     D_tableCanonical = readtable('DatabaseData/CanonicalMorsePotentials.txt');%first two base pairs match
%     D_structCanonical = table2struct(D_tableCanonical);
%     DCanonical_NearestNeighbours = {D_structCanonical.Var1 D_structCanonical.Var3 D_structCanonical.Var5 D_structCanonical.Var7};
%     DCanonical_NearestNeighbours1 = cellfun(@(x) strcat(x,'-',upper(seqcomplement(x))),DCanonical_NearestNeighbours,'UniformOutput',false);
%     DCanonical_NearestNeighbours2 = cellfun(@(x) strcat(reverse(x),'-',upper(seqcomplement(reverse(x)))),DCanonical_NearestNeighbours([3 4 5 7 8 10]),'UniformOutput',false);
%     
%     DCanonical_HBPotentials1 = {D_structCanonical.Var2 D_structCanonical.Var4 D_structCanonical.Var6 D_structCanonical.Var8};
%     DCanonical_HBPotentials2 = DCanonical_HBPotentials1([3 4 5 7 8 10]);
%     DCanonical_NearestNeighbours = [DCanonical_NearestNeighbours1 DCanonical_NearestNeighbours2];
%     DCanonical_HBPotentials = [DCanonical_HBPotentials1 DCanonical_HBPotentials2];
%     D_tableContext = readtable('DatabaseData/ContextDependentMorsePotentials.txt');%no correct matches, or correct adjacent
%     D_structContext = table2struct(D_tableContext);
%     DContext_NearestNeighbours = {D_structContext.Var2 D_structContext.Var5 D_structContext.Var8 D_structContext.Var11};
%     DContext_HBPotentials = {D_structContext.Var3 D_structContext.Var6 D_structContext.Var9 D_structContext.Var12};
%     DContext_NearestNeighbours = cellfun(@(x) upper(x),DContext_NearestNeighbours,'UniformOutput',false);
%     KDependentTrimer_NearestNeighbours = cellfun(@(x) upper(x),KDependentTrimer_NearestNeighbours,'UniformOutput',false);
%     KDependentTetramer_NearestNeighbours = cellfun(@(x) upper(x),KDependentTetramer_NearestNeighbours,'UniformOutput',false);
%     %1eV = 1.60217734*10^-22kJ
%     %1eV * (1.60217734*10^-22)*(6.0223*10^23)*0.23900573614kcal/kJ
%     %1eV = 23.0612 kcal/mol
%     %V equation, y relative displacement between base pairs
%     %sum(m/2(dudt^2+dvdt^2)+k/2*((un-un-1)^2+(vn-vn-1)^2)+D(e^-a(un-vn)-1)^2, m
%     %average mass
%     %displacement of nth base of first strand un, of second strand vn
%     %displacement at position at equilibrium along direction of hydrongen
%     %bond
%     %un and vn have different signs
%     %diameter of dna helix 20A, chain full turn at 34A
%     %for cytosine to guanine
%     %   two h-bond distance 2.9A and 3rd 3A.
%     %  distance across 10.8A
%     % energy hybrogen bond 4-29kJ/mol
%     %distance between nucleic acid on same strand  3.4A
%     %y is distance from equilibrium seperation
%     eVtoG = 23.0612;L = 0.03;%nm  y range -.1 to 20nm 
%     theta = 0.01; p=0;alpha=0.035;Ymin = -0.1;Ymax = 20;
%     UCanonical_HydrogenBonding = @(sq2,y) (0.001*getD(sq2,DCanonical_NearestNeighbours,DCanonical_HBPotentials)*(exp(-y/L)-1)^2);
%     UContext_HydrogenBonding = @(sq3,y) (0.001*getD(sq3,DContext_NearestNeighbours,DContext_HBPotentials)*(exp(-y/L)-1)^2);
%     UIndependent_Stacking = @(sq2,y1,y2) (getK(sq2,KIndependent_NearestNeighbours,KIndependent_StackingPotentials)/2*(y1-y2)^2);
%     UContext_Stacking = @(sq3,sq4,y1,y2) ((1+p*exp(-alpha*(y1+y2)))*...
%         (getK(sq3,KDependentTrimer_NearestNeighbours,KDependentTrimer_StackingPotentials)+...
%                getK(sq4,KDependentTetramer_NearestNeighbours,KDependentTetramer_StackingPotentials))/2*(y1-y2)^2);
%     UCI_Hamiltonian = @(sq2,y1,y2) (UCanonical_HydrogenBonding(sq2,y1) + UIndependent_Stacking(sq2,y1,y2));
%     UCC_Hamiltonian = @(sq2,sq3,sq4,y1,y2) (UCanonical_HydrogenBonding(sq2,y1) + UContext_Stacking(sq3,sq4,y1,y2));
%     UI_Hamiltonian = @(sq2,sq3,y1,y2) (UContext_HydrogenBonding(sq3,y1) + UIndependent_Stacking(sq2,y1,y2));
%     UC_Hamiltonian = @(sq3,sq4,y1,y2) (UContext_HydrogenBonding(sq3,y1) + UContext_Stacking(sq3,sq4,y1,y2));
%     %D is in meV (10^-3 eV) and k is in eV/nm^2, 
%     %could compute over all possible displacement y, or pick specific
%     %values for different base pair combinations, or average
%     for i=1:length(seq1)-1
%        %get each D abd k for each base pair and store in matrix
%            sq2 = strcat(upper(seq1(i:i+1)),'-',upper(seq2c(i:i+1)));
%        if (i<length(seq1)-2)
%            sq3  = strcat('{',upper(seq1(i:i+2)),'/',upper(seq2c(i:i+2)),'}');
%            sq3D = strcat('[',upper(seq1(i:i+2)),'/',upper(seq2c(i:i+2)),']');
%        else
%            sq3  = [];%add wildcards for edge cases
%            sq3D = [];%add wildcards for edge cases
%        end
%        if (i<length(seq1)-3)
%            sq4D = strcat('[',upper(seq1(i:i+3)),'/',upper(seq2c(i:i+3)),']'); 
%        else
%            sq4D = [];%add wildcards for edge cases
%        end
%        isCanon = sum(strcmp(sq2,DCanonical_NearestNeighbours));%is first two perfect match
%        %UC context dependent then independent
%        if (isCanon)
%            D_C(i) = getD(sq2,DCanonical_NearestNeighbours,DCanonical_HBPotentials);
%            K_C(i,i+1) = getK(sq2,KIndependent_NearestNeighbours,KIndependent_StackingPotentials);           
%        else
%            D_C(i) = getD(sq3,DContext_NearestNeighbours,DContext_HBPotentials);
%            K_C(i,i+1) = getK(sq3D,KDependentTrimer_NearestNeighbours,KDependentTrimer_StackingPotentials)+...
%                getK(sq4D,KDependentTetramer_NearestNeighbours,KDependentTetramer_StackingPotentials);
%        end
%     end
%     VoltageText = '$V(y_{i}) = D_{i}(e^{-y_{i}/\lambda_{i}}-1)^2$';
%     StackingText = '$W(y_{i},y_{i+1}) = \frac{k_{i,i+1}}{2}(1+\rho e^{-\alpha(y_{i}+y_{i+1})})(y_{i}^2-2y_{i}y_{i+1}\cos\theta+y_{i+1}^2)$';
%     HamiltonianText = '$U(y_{i},y_{i+1}) = V(y_{i})+W(y_{i},y_{i+1})$';
%     PartitionFuncText = '$Z_{y}(T) = \prod_{i=1}^{N}e^{\frac{-U(y_{i},y_{i+1})}{k_{B}T}}$';
%     EquilibriumDisplacementText = '$\langle y_{m}(T) \rangle = \int_{1}^{1} \prod_{i=1}^{N}e^{\frac{-U(y_{i},y_{i+1})}{k_{B}T}} dy$';
%     EnergyText = '$\Delta G = \Delta G_{init} + \sum_{i=1}^{N} U(\langle y_{i} \rangle,\langle y_{i+1}\rangle)$';
%       
%     Y = sym('y_%d',[1 length(seq1)]);
%     for i=1:length(seq1)-1
%        V_C(i) = D_C(i)/1000*(exp(-Y(i)/L)-1)^2;
%     end
%        V_C(length(seq1)) = 0;
%     for i=1:length(seq1)-1
%        W_C(i,i+1) = K_C(i,i+1)/2*(1+p*exp(-alpha*(Y(i)+Y(i+1))))*(Y(i)^2-2*Y(i)*Y(i+1)*cos(theta)+Y(i+1)^2);
%     end
%     for i=1:length(seq1)-1
%        U_C(i,i+1) = W_C(i,i+1) + V_C(i);
%     end
%     UC2 = sum(U_C);
%     UC3 = sum(UC2);
%     Gibson_Handle = matlabFunction(UC3,'vars',{Y});
%     UC4 = exp(-b*UC3);
%     Kernal_C = UC4;
%        %Kernal_CD = exp(-b*sum(U_CD,'all'));%integral prod(KFunc) from ymin to ymax for each y combined.
%     KS_C = Kernal_C;
% 
%     N_VecMat = eye(length(seq1));
%     T_hybrid = sym('T_hybrid%d',[1 1]);nK = sym('nK%d',[1 length(seq1)]);
%     PartialKernel_KC = cell(1,length(seq1));
%     PartialKernel_Handle = cell(1,length(seq1));
%     for i=1:length(seq1)
%         if (i<length(seq1))
%             PartialKernel_KC{i} = Y(i)^nK(i)*exp(-U_C(i,i+1)/(kb*(T_hybrid(1)+273.15)));   
%             PartialKernel_Handle{i} = matlabFunction(PartialKernel_KC{i},'vars',{Y(i:i+1) nK(i) T_hybrid});
%         else
%             PartialKernel_KC{i} = Y(i)^nK(i);
%             PartialKernel_Handle{i} = matlabFunction(PartialKernel_KC{i},'vars',{Y(i) nK(i)});  
%         end    
%     end
%     [Xi,Yi] = meshgrid(linspace(Ymin,Ymax,Np),linspace(Ymin,Ymax,Np));
%     Ck0 = cell(1,length(seq1));
%     Ck_N = cell(1,length(seq1));
%     Zy = zeros(1,length(seq1));
%     MZ = zeros(length(seq1),length(T_vector));
%     for v = 1:length(seq1)
%         if (v==1)
%             Ck0{v} = arrayfun(@(y) arrayfun(@(x) trapz(Xi(x,:),PartialKernel_Handle{v}([Xi(x,:)' Yi(x,:)'],0,T_vector(y))),1:size(Xi,1)),1:length(T_vector),'Un',0);
%             for j = 1:length(seq1)
%                Ck_N{j}{v} = arrayfun(@(y) arrayfun(@(x) trapz(Xi(x,:),PartialKernel_Handle{v}([Xi(x,:)' Yi(x,:)'],N_VecMat(j,v),T_vector(y))),1:size(Xi,1)),1:length(T_vector),'Un',0);
%             end
%         elseif (v<length(seq1))
%             Ck0{v} = arrayfun(@(y) arrayfun(@(x) trapz(Xi(x,:),Ck0{v-1}{y}'.*PartialKernel_Handle{v}([Xi(x,:)' Yi(x,:)'],0,T_vector(y))),1:size(Xi,1)),1:length(T_vector),'Un',0);
%             for j = 1:length(seq1)
%                Ck_N{j}{v} = arrayfun(@(y) arrayfun(@(x) trapz(Xi(x,:),Ck_N{j}{v-1}{y}'.*PartialKernel_Handle{v}([Xi(x,:)' Yi(x,:)'],N_VecMat(j,v),T_vector(y))),1:size(Xi,1)),1:length(T_vector),'Un',0);
%             end
%         else
%             Zy = arrayfun(@(y) trapz(Xi(1,:),Ck0{v-1}{y}'.*PartialKernel_Handle{v}(Xi(1,:)',0)),1:length(T_vector));
%             for j = 1:length(seq1)
%                MZ(j,:) = arrayfun(@(y) trapz(Xi(1,:),Ck_N{j}{v-1}{y}'.*PartialKernel_Handle{v}(Xi(1,:)',N_VecMat(j,v))),1:length(T_vector));
%             end
%          end       
%     end
%     Y_EQ = MZ./Zy;
%     Gibson_Energy = arrayfun(@(x) Gibson_Handle(Y_EQ(:,x)'),1:length(T_vector));
%     dGinit = -4.5;
%     dG_PDB = -Gibson_Energy*eVtoG +dGinit;  
%     
%     [phi,~,mu] = polyfit((T_vector+273.15),dG_PDB,1); 
%     pS = fit_reverseNormalization(phi,mu);%dG = pS(1)*T+pS(2) = dH-TdS 
%     dH_PDB = pS(2);
%     dS_PDB = -pS(1);
%      
%     dG_PDB0 = dH_PDB-T*dS_PDB;  
%     dHeq(N_models) = dH_PDB;
%     dSeq(N_models) = dS_PDB;
%     dGeq(N_models) = dG_PDB0;
% end
       %https://www.mathworks.com/matlabcentral/answers/56324-nested-numerical-integral-in-matlab
       %https://www.mathworks.com/matlabcentral/answers/373345-nested-integration-with-multiple-variables
       %https://www.mathworks.com/matlabcentral/answers/168330-how-to-do-numerical-multiple-integral-more-than-triple-by-using-matlab
       %https://www.mathworks.com/matlabcentral/answers/380906-is-it-possible-to-calculate-a-10-dimensional-integral
       %https://www.mathworks.com/matlabcentral/answers/419419-numerical-integration-for-nested-integral
       %average displacement
       
%     U_Total = 0;
%     for i=1:length(seq1)-1%Stacking sq2, first 2
%        sq2 = strcat(upper(seq1(i:i+1)),'-',upper(seq2c(i:i+1)));
%        if (i<length(seq1)-2)
%            sq3 = strcat('{',upper(seq1(i:i+2)),'/',upper(seq2c(i:i+2)),'}'); %if at end 
%            sq3D = strcat('[',upper(seq1(i:i+2)),'/',upper(seq2c(i:i+2)),']');
%        else
%            sq3 = [];%add wild cards for edge cases
%            sq3D = [];%add wild cards for edge cases
%        end
%        if (i<length(seq1)-3)
%           sq4D = strcat('[',upper(seq1(i:i+3)),'/',upper(seq2c(i:i+3)),']');
%        else
%           sq4D = [];%add wild cards for edge cases
%        end
%        isCanon = sum(strcmp(sq2,DCanonical_NearestNeighbours));%is first two perfect match
%        if (isCanon==1)
%            U_Total = U_Total + UCI_Hamiltonian(sq2,Ydis_C(i),Ydis_C(i+1));
%        else
%            U_Total = U_Total + UC_Hamiltonian(sq3D,sq4D,Ydis_C(i),Ydis_C(i+1));
%        end
%     end
%     %Trapp et al. BMC Biophysics 2011, 4:20
%     dGinit = -4.5;
%     dG = -U_Total*eVtoG +dGinit;
    %Triplicate NN model SNP
    %Li, G., Quan, Y., Wang, X., Liu, R., Bie, L., Gao, J., & Zhang, H. Y. (2019). Trinucleotide base pair stacking free energy for understanding TF-DNA recognition and the functions of SNPs. Frontiers in chemistry, 6, 666.
    %E
    %Hopfinger, M. C., Kirkpatrick, C. C., & Znosko, B. M. (2020). Predictions and analyses of RNA nearest neighbor parameters for modified nucleotides. Nucleic acids research, 48(16), 8901-8913.
    %Prog. Theor. Exp. Phys. 2020, 063J02
end
function m = Match(s)
if (strcmpi(s,'a')==1)
    m = 't';
elseif (strcmpi(s,'t')==1)
    m = 'a';
elseif (strcmpi(s,'g')==1)
    m = 'c';
elseif (strcmpi(s,'c')==1)
    m = 'g';
elseif (strcmpi(s,'n')==1)
    m = 'x';
else 
    m = [];
end

end
function K = getK(seq,K_Context,K_values)
% This function parses meso-scale model triple nucleotide binding sequence 
% cases to get parameter values for harmonic coupling.
   Idx = find(strcmp(seq,K_Context));
   if (length(seq)==5&&isempty(Idx))
      Idx2 = find(strcmp(reverse(seq),K_Context));
      Idx = union(Idx,Idx2);
   end
   if (length(seq)==9&&isempty(Idx))
      seqr = strcat(seq(1),reverse(seq(6:8)),seq(5),reverse(seq(2:4)),seq(9));
      Idx2 = find(strcmp(seqr,K_Context));
      Idx = union(Idx,Idx2);
   end
   if (length(seq)==11&&isempty(Idx))
      seqr = strcat(seq(1),reverse(seq(7:10)),seq(6),reverse(seq(2:5)),seq(11));
      Idx2 = find(strcmp(seqr,K_Context));
      Idx = union(Idx,Idx2);
   end
   if (~isempty(Idx))
       K = mean(cell2mat(K_values(Idx)));
   else
       K = 0;
   end
end
function D = getD(seq,K_Context,K_values)
% This function parses meso-scale model triple nucleotide binding sequence 
% cases to get parameter values for Morse Potential Depth.
   Idx = find(strcmp(seq,K_Context));
   if (length(seq)==5&&isempty(Idx))
      Idx2 = find(strcmp(reverse(seq),K_Context));
      Idx = union(Idx,Idx2);
   end
   if (length(seq)==9&&isempty(Idx))
      seqr = strcat(seq(1),reverse(seq(6:8)),seq(5),reverse(seq(2:4)),seq(9));
      Idx2 = find(strcmp(seqr,K_Context));
      Idx = union(Idx,Idx2);
   end
   if (length(seq)==11&&isempty(Idx))
      seqr = strcat(seq(1),reverse(seq(7:10)),seq(6),reverse(seq(2:5)),seq(11));
      Idx2 = find(strcmp(seqr,K_Context));
      Idx = union(Idx,Idx2);
   end
   if (~isempty(Idx))
       D = mean(cell2mat(K_values(Idx)));
   else
       D = 0;
   end
end

