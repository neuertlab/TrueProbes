function oligo = F_NearestNeighborTransitionWrapper_V2(sequence,SaltConcentration,Temperature,PrimerConc)
% This function is modified from matlab function OLIGOPROP to calculates properties of DNA oligonucleotide.
% This adds additional Transition-state Gibson Free Energy Calculations from different NN models 
sequence = convertStringsToChars(sequence);
seq = upper(sequence);
if(isnumeric(seq))
    seq = upper(int2nt(seq));
else
    seq = upper(seq);
end
if(any(seq=='N'))
    nFlag=1;
else
    nFlag=0;
end



T = str2sym(sprintf('T%d(t)',1)); % temperature in Celsius
S = sym('salt');          % salt concentration in moles per liter (M)
P= sym('primerConc');    % concentration of primers in mole per liter (M)

% compute sequence properties
numSeq = double(nt2int(seq));
[tm, tmdelta, NN, NNdelta]  = get_tm_NN(numSeq, S, P, nFlag);
NN(:,3) = NN(:,1) - ((Temperature+273.15)  .* (NN(:,2)./1000)); % DeltaG
if (~nFlag) % no ambiguous symbols 'N'
    NNdelta(:,3) = zeros(9,1);
else % occurrences of ambiguous N symbols
    NNdelta(:,3) = NNdelta(:,1) - ((Temperature+273.15)  .* (NNdelta(:,2) ./1000));
end
%salt correction
% build output structure

tm0 = [arrayfun(@(x) double(subs(tm(x),'salt',SaltConcentration)),1:3) arrayfun(@(x) double(subs(subs(tm(x),'salt',SaltConcentration),'primerConc',PrimerConc)),4:9)];
dGF  = NN(1:3,1) - ((T+273.15)  .* (NN(1:3,2)./1000)); % DeltaG
dGR  = NN(4:6,1) - ((T+273.15)  .* (NN(4:6,2)./1000)); % DeltaG
dGE  = NN(7:9,1) - ((T+273.15)  .* (NN(7:9,2)./1000)); % DeltaG

oligo.Tm = tm0;
oligo.Tmdelta = tmdelta;
oligo.Thermo = NN;
oligo.Thermodelta = NNdelta;
oligo.Tm_Function = tm;
oligo.dGf = dGF;
oligo.dGr = dGR;
oligo.dGeq = dGE;

end
% calculate melting temperature and thermo values (37 degrees C)
function [tm, tmdelta, NN, NNdelta]  = get_tm_NN(numSeq, salt, primerConc, nFlag)
% melting temperatures and thermodynamic values are returned as average +/- delta level.
% If no ambiguous symbols are present, the delta level is zero.

selfCompFlag = all(5-numSeq == numSeq(end:-1:1));

if(selfCompFlag) % self complementary sequence
    b = 1; % correction value for nearest neighbor melting temperature calculation
else
    b = 4;
end

tmlo_delta = zeros(3,1);
tmmd_delta = zeros(3,1);
tmhi_delta = zeros(3,1);

if (~nFlag) % no ambiguous symbols

    [NNl, NNldelta] = near_neigh(numSeq, length(numSeq), selfCompFlag, nFlag,-1);
    [NNm, NNmdelta] = near_neigh(numSeq, length(numSeq), selfCompFlag, nFlag,0);
    [NNu, NNudelta] = near_neigh(numSeq, length(numSeq), selfCompFlag, nFlag,1);
    
    tm_lo = (NNl(:,1) * 1000 ./ (NNl(:,2) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15; %TM NEAREST NEIGHBOR
    tm_md = (NNm(:,1) * 1000 ./ (NNm(:,2) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15; %TM NEAREST NEIGHBOR
    tm_hi = (NNu(:,1) * 1000 ./ (NNu(:,2) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15; %TM NEAREST NEIGHBOR
    
else % occurrences of 'N'
    [NNl, NNldelta] = near_neigh(numSeq, length(numSeq), selfCompFlag, nFlag,-1);
    [NNm, NNmdelta] = near_neigh(numSeq, length(numSeq), selfCompFlag, nFlag,0);
    [NNu, NNudelta] = near_neigh(numSeq, length(numSeq), selfCompFlag, nFlag,1);
    
    tm_lo = (((NNl(:,1)+ NNldelta(:,1)) * 1000 ./ ((NNl(:,2)+ NNldelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15 + ...
        ((NNl(:,1)- NNldelta(:,1)) * 1000 ./ ((NNl(:,2)- NNldelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15)* 1/2  ; % NEAREST NEIGHBOR
    tm_md = (((NNm(:,1)+ NNmdelta(:,1)) * 1000 ./ ((NNm(:,2)+ NNmdelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15 + ...
        ((NNm(:,1)- NNmdelta(:,1)) * 1000 ./ ((NNm(:,2)- NNmdelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15)* 1/2  ; % NEAREST NEIGHBOR
    tm_hi = (((NNu(:,1)+ NNudelta(:,1)) * 1000 ./ ((NNu(:,2)+ NNudelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15 + ...
        ((NNu(:,1)- NNudelta(:,1)) * 1000 ./ ((NNu(:,2)- NNudelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15)* 1/2  ; % NEAREST NEIGHBOR
    tmlo_delta =(((NNl(:,1)+ NNldelta(:,1)) * 1000 ./ ((NNl(:,2)+ NNldelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15 - ...
        (((NNl(:,1)- NNldelta(:,1)) * 1000 ./ ((NNl(:,2)- NNldelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15)) * 1/2 ;
    tmmd_delta =(((NNm(:,1)+ NNmdelta(:,1)) * 1000 ./ ((NNm(:,2)+ NNmdelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15 - ...
        (((NNm(:,1)- NNmdelta(:,1)) * 1000 ./ ((NNm(:,2)- NNmdelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15)) * 1/2 ;
    tmhi_delta =(((NNu(:,1)+ NNudelta(:,1)) * 1000 ./ ((NNu(:,2)+ NNudelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15 - ...
        (((NNu(:,1)- NNudelta(:,1)) * 1000 ./ ((NNu(:,2)- NNudelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15)) * 1/2 ;
end
tm = [tm_lo; tm_md; tm_hi]';
tmdelta = [tmlo_delta; tmmd_delta; tmhi_delta]';
NN = [NNl; NNm; NNu];
NNdelta = [NNldelta; NNmdelta; NNudelta];

tm = tm([1 4 7 2 5 8 3 6 9]);
tmdelta = tmdelta([1 4 7 2 5 8 3 6 9]);
NN  = NN([1 4 7 2 5 8 3 6 9],:);
NNdelta = NNdelta([1 4 7 2 5 8 3 6 9],:);
end
% compute thermo values using Nearest Neighbor methods
%one typo in energy eq. for GG/CC
function [NN, NNdelta] = near_neigh(seq, seq_length, selfCompFlag, nFlag,bounded)
%AA/TT AT/TA TA/AT CA/GT GT/CA CT/GA GA/CT CG/GC GC/CG GG/CC Initiation Terminal_AT    
DeltaFm_H = zeros(1,12);%kcal/mol
DeltaFm_HErr = zeros(1,12);%kcal/mol
DeltaFm_S = [-0.08 -0.16 -0.41 -0.35 -0.02 -0.29 -0.22 0.08 -0.02 0.47 -25.1 0.32];%cal/molK
DeltaFm_SErr = [0.11 0.23 0.23 0.21 0.20 0.27 0.20 0.29 0.24 0.16 0.7 0.15];%cal/molK
DeltaRm_H = [9.2 8.6 5.6 11.9 9.6 10.2 8.1 14.5 11.2 9.8 -14.8 -1.0];%kcal/mol
DeltaRm_HErr = [0.6 1.7 1.4 1.4 1.6 1.2 1.3 1.5 1.7 1.3 4.2 1.0];
DeltaRm_S = [26.5 24.2 15.8 33.5 25.9 28.9 21.7 40.7 29.2 26.6 -66.8 -2.3];%cal/molK
DeltaRm_SErr = [1.8 5.3 4.3 4.4 5 3.6 4.0 4.9 5.5 4.2 13.3 3.1];%cal/molK
DeltaSymmetry = [0 -1.4];

DeltaFm_Init = [0 -25.1];
DeltaFm_InitErr = [0 0.7];
DeltaRm_Init = [-14.8 -66.8];
DeltaRm_InitErr = [4.2 13.3];
DeltaFm_TermAT = [0 0.32];
DeltaFm_TermATErr = [0 0.15];
DeltaRm_TermAT = [-1 -2.3];
DeltaRm_TermATErr = [1 3.1];

DeltaF_H = DeltaFm_H + bounded*DeltaFm_HErr;
DeltaF_S = DeltaFm_S + bounded*DeltaFm_SErr;
DeltaR_H = DeltaRm_H + bounded*DeltaRm_HErr;
DeltaR_S = DeltaRm_S + bounded*DeltaRm_SErr;
DeltaF_Init = DeltaFm_Init + bounded*DeltaFm_InitErr;
DeltaR_Init = DeltaRm_Init + bounded*DeltaRm_InitErr;
DeltaF_TermAT = DeltaFm_TermAT + bounded*DeltaFm_TermATErr;
DeltaR_TermAT = DeltaRm_TermAT + bounded*DeltaRm_TermATErr;


DeltaE_H = DeltaF_H-DeltaR_H;
DeltaE_S = DeltaF_S-DeltaR_S;
DeltaE_Init = DeltaF_Init - DeltaR_Init;
DeltaE_TermAT = DeltaF_TermAT - DeltaR_TermAT;
corrs = [1 5 7 2;4 10 8 6;7 9 10 5;3 7 4 1];   
%salt sequence dependent vs sequence independent correction
%concentration dependent correction to energy and Tm
%     A   C  G  T
    %(5'->3')/(3'-5)
    %opposite bonds
    %AA/TT (1)
    %AT/TA (2)   TA/AT (3)
    %CA/GT (4)   GT/CA (5)
    %CT/GA (6)   GA/CT (7)
    %CG/GC (8)   GC/CG (9)
    %GG/CC (10)
	%first check first part in front of slash, if not flip and check second
    %A   AA(1)      AC(5)   AG(7)  AT(2)    
    %C   CA(4)      CC(10)  CG(8)  CT(6)
    %G   GA(7)      GC(9)   GG(10) GT(5)
    %T   TA(3)      TC(7)   TG(4)  TT(1)
    
MethodF_H = DeltaF_H(corrs);MethodF_S = DeltaF_S(corrs);    
MethodR_H = DeltaR_H(corrs);MethodR_S = DeltaR_S(corrs); 
MethodE_H = DeltaE_H(corrs);MethodE_S = DeltaE_S(corrs); 
if (~nFlag)
    % nearest neighbor parameters from: Panjkovich and Melo, Bioinformatics  Vol 21 no 6 pp 711-722 2004 [1]
    % rows corresponds to A,C,G,T respectively; columns correspond to A,C,G,T respectively
    ind = sub2ind([4 4],seq(1:seq_length-1),seq(2:seq_length));
else
    % nearest neighbor parameters as in [1] with added average values for
    % all possible combinations involving 'N'
    % rows corresponds to A,C,G,T,N respectively; columns correspond to
    % A,C,G,T,N respectively
    
    % MethodF_H(1,5) = mean(DeltaF_H(corrs(1,:)));
    % MethodF_H(2,5) = mean(DeltaF_H(corrs(2,:)));
    % MethodF_H(3,5) = mean(DeltaF_H(corrs(3,:)));
    % MethodF_H(4,5) = mean(DeltaF_H(corrs(4,:)));
    % MethodF_H(5,1) = mean(DeltaF_H(corrs(:,1)));
    % MethodF_H(5,2) = mean(DeltaF_H(corrs(:,2)));    
    % MethodF_H(5,3) = mean(DeltaF_H(corrs(:,3)));
    % MethodF_H(5,4) = mean(DeltaF_H(corrs(:,4)));   
    % MethodF_H(5,5) = mean(DeltaF_H(corrs(:)));
    % MethodF_S(1,5) = mean(DeltaF_S(corrs(1,:)));
    % MethodF_S(2,5) = mean(DeltaF_S(corrs(2,:)));
    % MethodF_S(3,5) = mean(DeltaF_S(corrs(3,:)));
    % MethodF_S(4,5) = mean(DeltaF_S(corrs(4,:)));
    % MethodF_S(5,1) = mean(DeltaF_S(corrs(:,1)));
    % MethodF_S(5,2) = mean(DeltaF_S(corrs(:,2)));    
    % MethodF_S(5,3) = mean(DeltaF_S(corrs(:,3)));
    % MethodF_S(5,4) = mean(DeltaF_S(corrs(:,4)));   
    % MethodF_S(5,5) = mean(DeltaF_S(corrs(:)));
    % 
    % MethodR_H(1,5) = mean(DeltaR_H(corrs(1,:)));
    % MethodR_H(2,5) = mean(DeltaR_H(corrs(2,:)));
    % MethodR_H(3,5) = mean(DeltaR_H(corrs(3,:)));
    % MethodR_H(4,5) = mean(DeltaR_H(corrs(4,:)));
    % MethodR_H(5,1) = mean(DeltaR_H(corrs(:,1)));
    % MethodR_H(5,2) = mean(DeltaR_H(corrs(:,2)));    
    % MethodR_H(5,3) = mean(DeltaR_H(corrs(:,3)));
    % MethodR_H(5,4) = mean(DeltaR_H(corrs(:,4)));   
    % MethodR_H(5,5) = mean(DeltaR_H(corrs(:)));
    % MethodR_S(1,5) = mean(DeltaR_S(corrs(1,:)));
    % MethodR_S(2,5) = mean(DeltaR_S(corrs(2,:)));
    % MethodR_S(3,5) = mean(DeltaR_S(corrs(3,:)));
    % MethodR_S(4,5) = mean(DeltaR_S(corrs(4,:)));
    % MethodR_S(5,1) = mean(DeltaR_S(corrs(:,1)));
    % MethodR_S(5,2) = mean(DeltaR_S(corrs(:,2)));    
    % MethodR_S(5,3) = mean(DeltaR_S(corrs(:,3)));
    % MethodR_S(5,4) = mean(DeltaR_S(corrs(:,4)));   
    % MethodR_S(5,5) = mean(DeltaR_S(corrs(:)));
    
    % MethodE_H(1,5) = mean(DeltaE_H(corrs(1,:)));
    % MethodE_H(2,5) = mean(DeltaE_H(corrs(2,:)));
    % MethodE_H(3,5) = mean(DeltaE_H(corrs(3,:)));
    % MethodE_H(4,5) = mean(DeltaE_H(corrs(4,:)));
    % MethodE_H(5,1) = mean(DeltaE_H(corrs(:,1)));
    % MethodE_H(5,2) = mean(DeltaE_H(corrs(:,2)));    
    % MethodE_H(5,3) = mean(DeltaE_H(corrs(:,3)));
    % MethodE_H(5,4) = mean(DeltaE_H(corrs(:,4)));   
    % MethodE_H(5,5) = mean(DeltaE_H(corrs(:)));
    % MethodE_S(1,5) = mean(DeltaE_S(corrs(1,:)));
    % MethodE_S(2,5) = mean(DeltaE_S(corrs(2,:)));
    % MethodE_S(3,5) = mean(DeltaE_S(corrs(3,:)));
    % MethodE_S(4,5) = mean(DeltaE_S(corrs(4,:)));
    % MethodE_S(5,1) = mean(DeltaE_S(corrs(:,1)));
    % MethodE_S(5,2) = mean(DeltaE_S(corrs(:,2)));    
    % MethodE_S(5,3) = mean(DeltaE_S(corrs(:,3)));
    % MethodE_S(5,4) = mean(DeltaE_S(corrs(:,4)));   
    % MethodE_S(5,5) = mean(DeltaE_S(corrs(:)));

    MethodF_H(1:4,5) = mean(DeltaF_H(corrs(1:4,:)),2);
    MethodF_H(5,1:4) = mean(DeltaF_H(corrs(:,1:4)),1);
    MethodF_H(5,5) = mean(DeltaF_H(corrs(:)));
    MethodF_S(1:4,5) = mean(DeltaF_S(corrs(1:4,:)),2);
    MethodF_S(5,1:4) = mean(DeltaF_S(corrs(:,1:4)),1);
    MethodF_S(5,5) = mean(DeltaF_S(corrs(:)));
    MethodR_H(1:4,5) = mean(DeltaR_H(corrs(1:4,:)),2);
    MethodR_H(5,1:4) = mean(DeltaR_H(corrs(:,1:4)),1);
    MethodR_H(5,5) = mean(DeltaR_H(corrs(:)));
    MethodR_S(1:4,5) = mean(DeltaR_S(corrs(1:4,:)),2);
    MethodR_S(5,1:4) = mean(DeltaR_S(corrs(:,1:4)),1);
    MethodR_S(5,5) = mean(DeltaR_S(corrs(:)));
    MethodE_H(1:4,5) = mean(DeltaE_H(corrs(1:4,:)),2);
    MethodE_H(5,1:4) = mean(DeltaE_H(corrs(:,1:4)),1);
    MethodE_H(5,5) = mean(DeltaE_H(corrs(:)));
    MethodE_S(1:4,5) = mean(DeltaE_S(corrs(1:4,:)),2);
    MethodE_S(5,1:4) = mean(DeltaE_S(corrs(:,1:4)),1);
    MethodE_S(5,5) = mean(DeltaE_S(corrs(:)));
    seq(seq==15)=5; % substitute numeric value of 'N' with 5
    ind = sub2ind([5 5],seq(1:seq_length-1),seq(2:seq_length));
end
%Tm^-1 = R/dHln(Ct)+dS/dH, 
%Tm = DH/(DS+R*ln(Ct)),  Ct/4 for self complementart
% NN is 4x2 matrix. Columns are DeltaH and DeltaS. Rows correspond to
% methods by Bres86, SantaLucia96, SantaLucia98 and Sugimoto96.
NN = [sum(MethodF_H(ind)),sum(MethodF_S(ind));...
    sum(MethodR_H(ind)),sum(MethodR_S(ind));...
    sum(MethodE_H(ind)),sum(MethodE_S(ind))];

% Corrections: all AT pairs, any GC pairs, symmetry, initiation
if(~nFlag) % ambiguous symbols 'N' not present

    %86 initiation 5kcal for GC and 6kcal for A-T
    %self symetry 0.4, comp of 0
    % only AT pairs or any GC pairs?
    %specify that heat is zero  intination coverted into entropy units
    NN =  NN + [DeltaF_Init; DeltaR_Init; DeltaE_Init];%initiation

%for lucia 96  symetry and AT/GC AT asumes at least once
    % symmetry
    if(selfCompFlag)
        NN = NN + [DeltaSymmetry; DeltaSymmetry; DeltaSymmetry];
    end

    % initiation with terminal  5'
    if(seq(1) == 2 || seq(1) == 3)
        NN(1,:) = NN(1,:) + [0 0]; %GC
        NN(2,:) = NN(2,:) + [0 0]; %GC
        NN(3,:) = NN(3,:) + [0 0]; %GC
    elseif(seq(1) == 1 || seq(1) == 4)
        NN(1,:) = NN(1,:) + [DeltaF_TermAT]; %AT
        NN(2,:) = NN(2,:) + [DeltaR_TermAT]; %AT
        NN(3,:) = NN(3,:) + [DeltaE_TermAT]; %AT
    end

    % initiation with terminal  3'
    if(seq(end) == 2 || seq(end) == 3)
        NN(1,:) = NN(1,:) + [0 0]; %GC
        NN(2,:) = NN(2,:) + [0 0]; %GC
        NN(3,:) = NN(3,:) + [0 0]; %GC
    elseif(seq(end) == 1 || seq(end) == 4)
        NN(1,:) = NN(1,:) + [DeltaF_TermAT]; %AT
        NN(2,:) = NN(2,:) + [DeltaR_TermAT]; %AT
        NN(3,:) = NN(3,:) + [DeltaE_TermAT]; %AT
    end

    NNdelta=zeros(3,2);

else % 'N' symbols are present

    % only AT pairs or any GC pairs?
    NN1 = NN;  % case when all Ns are G/C
    NN2 = NN; % case when all Ns are A/T
    NN1 =  NN1 + [DeltaF_Init; DeltaR_Init; DeltaE_Init];%initiation
    NN2 =  NN2 + [DeltaF_Init; DeltaR_Init; DeltaE_Init];%initiation

    % symmetry
    if(selfCompFlag)
         NN1 = NN1 + [DeltaSymmetry; DeltaSymmetry; DeltaSymmetry];
         NN2 = NN2 + [DeltaSymmetry; DeltaSymmetry; DeltaSymmetry];
    end

    % initiation with terminal 5'(only Sant98)
    if(seq(1) == 2 || seq(1) == 3 || seq(1) == 5)
        NN2(1,:) = NN2(1,:) + [0 0]; %GC
        NN2(2,:) = NN2(2,:) + [0 0]; %GC
        NN2(3,:) = NN2(3,:) + [0 0]; %GC
    elseif(seq(1) == 1 || seq(1) == 4 || seq(1) == 5)
        NN1(1,:) = NN1(1,:) + [DeltaF_TermAT]; %AT
        NN1(2,:) = NN1(2,:) + [DeltaR_TermAT]; %AT
        NN1(3,:) = NN1(3,:) + [DeltaE_TermAT]; %AT
    end

    % initiation with terminal 3'(only Sant98)
    if(seq(end) == 2 || seq(end) == 3 || seq(end) == 5)
        NN2(1,:) = NN2(1,:) + [0 0]; %GC
        NN2(2,:) = NN2(2,:) + [0 0]; %GC
        NN2(3,:) = NN2(3,:) + [0 0]; %GC
    elseif(seq(end) == 1 || seq(end) == 4 || seq(end) == 5)
        NN1(1,:) = NN1(1,:) + [DeltaF_TermAT]; %AT
        NN1(2,:) = NN1(2,:) + [DeltaR_TermAT]; %AT
        NN1(3,:) = NN1(3,:) + [DeltaE_TermAT]; %AT
    end

    NN = (NN1+NN2)/2; % avg
    NNdelta = (max(NN1,NN2)- min(NN1,NN2))/2; % delta level
end
end
