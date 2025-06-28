function oligo = F_NearestNeighborWrapper_V2(sequence,SaltConcentration,Temperature,PrimerConc)
%MERFISH uses 2004 SantaLucia Hicks
% PaintSHOP defaults to BioPython default='DNA_NN3',
%AA/TT AT/TA TA/AT CA/GT GT/CA CT/GA GA/CT CG/GC GC/CG GG/CC Initiation Terminal_AT  


%  Different Paper different parameters for binding
        %Probe Design Methods Used different models 
        %    Model 1, Model 2, Model 3
        %    Model 1, Model 1, Model 1
        %    Compare using same model as them
        %    PaintSHOP (Model 1)
        %    MERFISH   (Model 2)
        %    OligoStan (Model 3) -> Probes Design  
        %    PaintSHOP/MERFISH/Oligostan/Ours (Model X)
        %    PaintSHOP OFF-Target dG -> 0-100
        %    MERFISH off-targets? 
        
        %    Model -> dG,  predict Tm
        %                         ->  
        %    Perform Better If Use Same Model
        %            Expression, Off-target, interactions
        %    Potential Models Performce Similarly
        %            Model 1 = Model 2 + Value
        
        %    Best Model (Newer Papers)
        
        %    Most Recent Model from other software 2004
        
        %    New Paper 2017-2023.
        
        %  Statements Even when using same models predicted to work better
        %  Use more updated model so more accurate
        
        % How are newer model validated
        
        
        % dG    dH  dS  dCp   -> dG/K(T)
        
        %Keq(Probe,Target,Site)  Model
        %Keq{Model}(Probe,Target,Site)
        
        %Keq(Probe,Target,Site,Model)  
        
        %Keq(Probe,Target,Site)
        
        %Method1.Keq vs Method2.Keq
        % Keq(:,:,:,[1 2])   speed vs memory usage.
        % Keq(:,:,y,:)
        % Salt Concentration.
       
        
        % Dimensions 
%BioPython 
%DNA/DNA
    %DNA_NN1 Breslauer 86
    %DNA_NN2 Sugimoto 96
    %DNA_NN3 Allawi& SantaLucia 1997
    %DNA_NN4 SantaLucia 2004
%RNA/RNA
    %RNA_NN1 Freier 1996
    %RNA_NN2 Xia 1998
    %RNA_NN3 Chen 2012
%RNA/DNA
    %R_DNA_NN1 Sugimoto 1995
%Terminal Mismatches
    %DNA_TMM1 SantaLucia & Peyret 2001
%Internal Mismatches (including isonine mismatches)
    %DNA_IMM1 Allawi & SantraLucia 1997-1998
    %DNA_IMM1 Peyret 1999
    %DNA_IMM1 Watkins & SantaLucia 2005
%Dangling Ends
    %DNA_DE1 Bommarito et al 2000
    %RNA_DE1 Tuner & Mathews 2010
    
    
%Biochemistry
%. 1995 Sep 5;34(35):11211-6. doi: 10.1021/bi00035a029.
%Thermodynamic parameters to predict stability of RNA/DNA hybrid duplexes    
%N Sugimoto 1, S Nakano, M Katoh, A Matsumura, H Nakamuta, T Ohmichi, M Yoneyama, M Sasaki

%Nucleic Acids Res
% 2020 Dec 2;48(21):12042-12054. doi: 10.1093/nar/gkaa572.
%Improved nearest-neighbor parameters for the stability of RNA/DNA hybrids under a physiological condition    
%Banerjee D, Tateishi-Karimata H, Ohyama T, Ghosh S, Endoh T, Takahashi S, Sugimoto N.
    
%Biophys Rep (N Y)
%. 2023 Mar 2;3(2):100101. doi: 10.1016/j.bpr.2023.100101. eCollection 2023 Jun 14.
%MeltR software provides facile determination of nucleic acid thermodynamics
%Jacob P Sieg 1 2, Sebastian J Arteaga 3, Brent M Znosko 3, Philip C Bevilacqua 1 2 4

%Nucleic Acids Res
%. 2023 May 22;51(9):4101-4111. doi: 10.1093/nar/gkad020.
%Nearest-neighbor parameters for the prediction of RNA duplex stability in diverse in vitro and cellular-like crowding conditions
%Saptarshi Ghosh 1, Shuntaro Takahashi 1, Dipanwita Banerjee 1, Tatsuya Ohyama 1, Tamaki Endoh 1, Hisae Tateishi-Karimata 1, Naoki Sugimoto


%Bres86
%SantaLucia96
%SantaLucia98
%Sugimoto96
%SantaLucia04
%Allawi97
%Rejali21
warning('off','bioinfo:oligoprop:SeqLengthTooShort');
oligo = oligoprop(sequence,'Salt',SaltConcentration,'Temp',Temperature,'PrimerConc',PrimerConc);
TmBlock = oligo.Tm;
TmDeltaBlock = oligo.Tmdelta;
ThermoBlock = oligo.Thermo;
ThermoDeltaBlock = oligo.Thermodelta;

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
% compute sequence properties
numSeq = double(nt2int(seq));
baseNum = [sum(numSeq == 1) sum(numSeq == 2) sum(numSeq == 3) sum(numSeq == 4) sum(numSeq == 15)];
if (~nFlag) % no ambiguous symbols 'N'
    [tm, tmdelta, NN, NNdelta]  = get_tm_NN(numSeq, baseNum, SaltConcentration, PrimerConc, nFlag);
    NN(:,3) = NN(:,1) - ((Temperature+273.15)  .* (NN(:,2)./1000)); % DeltaG
    NNdelta(:,3) = zeros(3,1);
else % occurrences of ambiguous N symbols
    % warning('bioinfo:oligoprop:ambiguousCheck', 'Ambiguous symbols N in input sequence.');
    % average between case when all Ns are C/G and case when all Ns are A/T
    [tm, tmdelta, NN, NNdelta]  = get_tm_NN(numSeq, baseNum, SaltConcentration, PrimerConc, nFlag);
    NN(:,3) = NN(:,1) - ((Temperature+273.15)  .* (NN(:,2)./1000)); % DeltaG
    NNdelta(:,3) = NNdelta(:,1) - ((Temperature+273.15)  .* (NNdelta(:,2)./1000));
end
        
% This adds additional Gibson Free Energy Calculations from different NN models 
Tm = [TmBlock tm];
TmDelta = [TmDeltaBlock tmdelta];
Thermo = [ThermoBlock; NN];
ThermoDelta = [ThermoDeltaBlock; NNdelta];
oligo.Tm = Tm;
oligo.Tmdelta = TmDelta;
oligo.Thermo = Thermo;
oligo.Thermodelta = ThermoDelta;
oligo = rmfield(oligo,{'Hairpins','Dimers'});
%NN(:,2) = NN(:,2) + 0.368*(length(seq))*log(salt);
%NN(:,3) = NN(:,3) - 0.114*(length(seq))*log(salt);
%dG =dH - T*dS
%  dS(1M) + 0.368*N/2*log(Na+)
%length(seq)
%salt adjustment sait lucia 2004
% add physiological conditions, different types of salt conc.
% add energy from free loops and secondary/tertirary structure nearest neigbor
%chemical modifications and RNA/DNA, RNA/RNA interactions'
%new 2023 papers

%./vargibbs -ct=400 -o=mytest -par=data/P-SL98.par -seq=ACGT -calc=prediction


end

% FUNCTIONS

% calculate melting temperature and thermo values (37 degrees C)
function [tm, tmdelta, NN, NNdelta]  = get_tm_NN(numSeq, baseNum, salt, primerConc, nFlag)
% melting temperatures and thermodynamic values are returned as average +/- delta level.
% If no ambiguous symbols are present, the delta level is zero.

selfCompFlag = all(5-numSeq == numSeq(end:-1:1));

if(selfCompFlag) % self complementary sequence
    b = 1; % correction value for nearest neighbor melting temperature calculation
else
    b = 4;
end

tmdelta = zeros(1,3);

if (~nFlag) % no ambiguous symbols
    [NN, NNdelta] = near_neigh(numSeq, length(numSeq), selfCompFlag, nFlag);
    tm = (NN(:,1) * 1000 ./ (NN(:,2) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15; %TM NEAREST NEIGHBOR

else % occurrences of 'N'
    [NN, NNdelta] = near_neigh(numSeq, length(numSeq), selfCompFlag, nFlag);
    tm = (((NN(:,1)+ NNdelta(:,1)) * 1000 ./ ((NN(:,2)+ NNdelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15 + ...
        ((NN(:,1)- NNdelta(:,1)) * 1000 ./ ((NN(:,2)- NNdelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15)* 1/2  ; % NEAREST NEIGHBOR
    tmdelta(1:3)=(((NN(:,1)+ NNdelta(:,1)) * 1000 ./ ((NN(:,2)+ NNdelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15 - ...
        (((NN(:,1)- NNdelta(:,1)) * 1000 ./ ((NN(:,2)- NNdelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15)) * 1/2 ;
end
tm = (tm)';


end
% compute thermo values using Nearest Neighbor methods
function [NN, NNdelta] = near_neigh(seq, seq_length, selfCompFlag, nFlag)
Sant04H_Table = [-7.6 -7.2 -7.2 -8.5 -8.4 -7.8 -8.2 -10.6 -9.8 -8.0];%kcal/mol
Sant04S_Table = [-21.3 -20.4 -21.3 -22.7 -22.4 -21.0 -22.2 -27.2 -24.4 -19.9];%cal/molK
Sant04InitwTermGC = [0.2 -5.7];%H %S
Sant04InitwTermAT = [2.2 6.9];%H %SS
Sant04InitTermATwoTerminalAnnotation = [0 0];%H %S
Sant04InitTermGCwoTerminalAnnotation = [0 0];%H %S
Sant04SymmetryCorrection = [0 -1.4];%H %S



Allawi97H_Table = [-7.9 -7.2 -7.2 -8.5 -8.4 -7.8 -8.2 -10.6 -9.8 -8.0];%kcal/mol
Allawi97S_Table = [-22.2 -20.4 -21.3 -22.7 -22.4 -21.0 -22.2 -27.2 -24.4 -19.9];%cal/molK
Allawi97Herr_Table = [0.2 0.7 0.9 0.6 0.5 0.6 0.6 0.6 0.4 0.9];%kcal/mol
Allawi97Serr_Table = [0.8 2.4 2.4 2.0 2.0 2.0 1.7 2.6 2.0 1.8];%cal/molK
Allawi97InitwTermGC = [0.1 -2.8];%H %S
Allawi97InitGCErr = [1.1 0.2];%H %S
Allawi97InitwTermAT = [2.3 4.1];%H %S
Allawi97InitTermATErr = [1.3 0.2];%H %S
Allawi97SymmetryCorrection = [0 -1.4];%H %S
Allawi97InitTermATwoTerminalAnnotation = [0 0];%H %S
Allawi97InitTermGCwoTerminalAnnotation = [0 0];%H %S


Rejali21H_Table = [-9.2 -8.6 -5.6 -11.9 -10.2 -8.1 -14.5 -11.2 -9.8 14.8 1.0];
Rejali21S_Table = [-26.6 -24.4 -16.2 -33.9 -25.9 -29.2 -21.9 -40.6 -29.2 -26.2 41.7 2.6];
Rejali21Herr_Table = [0.6 1.7 1.4 1.4 1.6 1.2 1.3 1.5 1.7 1.3 4.2 1.0];
Rejali21Serr_Table = [1.8 5.3 4.3 4.4 5.0 3.6 4.0 4.9 5.5 4.2 13.3 3.1];
Rejali21InitwTermGC = [14.8 41.7];
Rejali21InitwTermAT = [1.0 2.6];
Rejali21SymmetryCorrection = [0 -1.4];%H %S
Rejali21InitTermATwoTerminalAnnotation = [0 0];%H %S
Rejali21InitTermGCwoTerminalAnnotation = [0 0];%H %S



 
corrs =[1 5  7  2;4 10 8  6;7 9  10 5;3 7  4  1];  
Sant04_H = Sant04H_Table(corrs);Sant04_S = Sant04S_Table(corrs);  
Allawi97_H = Allawi97H_Table(corrs);Allawi97_S = Allawi97S_Table(corrs); 
Rejali21_H = Rejali21H_Table(corrs);Rejali21_S = Rejali21S_Table(corrs);  

if (~nFlag)
    % nearest neighbor parameters from: Panjkovich and Melo, Bioinformatics  Vol 21 no 6 pp 711-722 2004 [1]
    % rows corresponds to A,C,G,T respectively; columns correspond to A,C,G,T respectively
    ind = sub2ind([4 4],seq(1:seq_length-1),seq(2:seq_length));
else
    % nearest neighbor parameters as in [1] with added average values for
    % all possible combinations involving 'N'
    % rows corresponds to A,C,G,T,N respectively; columns correspond to
    % A,C,G,T,N respectively
    Sant04_H(1:4,5) = mean(Sant04H_Table(corrs(1:4,:)),2);
    Sant04_H(5,1:4) = mean(Sant04H_Table(corrs(:,1:4)),1);
    Sant04_H(5,5) = mean(Sant04H_Table(corrs(:)));
    Sant04_S(1:4,5) = mean(Sant04S_Table(corrs(1:4,:)),2);
    Sant04_S(5,1:4) = mean(Sant04S_Table(corrs(:,1:4)),1);
    Sant04_S(5,5) = mean(Sant04S_Table(corrs(:)));
    Allawi97_H(1:4,5) = mean(Allawi97H_Table(corrs(1:4,:)),2);
    Allawi97_H(5,1:4) = mean(Allawi97H_Table(corrs(:,1:4)),1);
    Allawi97_H(5,5) = mean(Allawi97H_Table(corrs(:)));
    Allawi97_S(1:4,5) = mean(Allawi97S_Table(corrs(1:4,:)),2);
    Allawi97_S(5,1:4) = mean(Allawi97S_Table(corrs(:,1:4)),1);
    Allawi97_S(5,5) = mean(Allawi97S_Table(corrs(:)));
    Rejali21_H(1:4,5) = mean(Rejali21H_Table(corrs(1:4,:)),2);
    Rejali21_H(5,1:4) = mean(Rejali21H_Table(corrs(:,1:4)),1);
    Rejali21_H(5,5) = mean(Rejali21H_Table(corrs(:)));
    Rejali21_S(1:4,5) = mean(Rejali21S_Table(corrs(1:4,:)),2);
    Rejali21_S(5,1:4) = mean(Rejali21S_Table(corrs(:,1:4)),1);
    Rejali21_S(5,5) = mean(Rejali21S_Table(corrs(:)));
    % Sant04_H(1,5) = mean(Sant04H_Table(corrs(1,:)));
    % Sant04_H(2,5) = mean(Sant04H_Table(corrs(2,:)));
    % Sant04_H(3,5) = mean(Sant04H_Table(corrs(3,:)));
    % Sant04_H(4,5) = mean(Sant04H_Table(corrs(4,:)));
    % Sant04_H(5,1) = mean(Sant04H_Table(corrs(:,1)));
    % Sant04_H(5,2) = mean(Sant04H_Table(corrs(:,2)));    
    % Sant04_H(5,3) = mean(Sant04H_Table(corrs(:,3)));
    % Sant04_H(5,4) = mean(Sant04H_Table(corrs(:,4)));   
    % Sant04_S(1,5) = mean(Sant04S_Table(corrs(1,:)));
    % Sant04_S(2,5) = mean(Sant04S_Table(corrs(2,:)));
    % Sant04_S(3,5) = mean(Sant04S_Table(corrs(3,:)));
    % Sant04_S(4,5) = mean(Sant04S_Table(corrs(4,:)));
    % Sant04_S(5,1) = mean(Sant04S_Table(corrs(:,1)));
    % Sant04_S(5,2) = mean(Sant04S_Table(corrs(:,2)));    
    % Sant04_S(5,3) = mean(Sant04S_Table(corrs(:,3)));
    % Sant04_S(5,4) = mean(Sant04S_Table(corrs(:,4)));   
    %Sant04_S(5,5) = mean(Sant04S_Table(corrs(:)));
    % Allawi97_H(1,5) = mean(Allawi97H_Table(corrs(1,:)));
    % Allawi97_H(2,5) = mean(Allawi97H_Table(corrs(2,:)));
    % Allawi97_H(3,5) = mean(Allawi97H_Table(corrs(3,:)));
    % Allawi97_H(4,5) = mean(Allawi97H_Table(corrs(4,:)));
    % Allawi97_H(5,1) = mean(Allawi97H_Table(corrs(:,1)));
    % Allawi97_H(5,2) = mean(Allawi97H_Table(corrs(:,2)));    
    % Allawi97_H(5,3) = mean(Allawi97H_Table(corrs(:,3)));
    % Allawi97_H(5,4) = mean(Allawi97H_Table(corrs(:,4)));   
    % Allawi97_H(5,5) = mean(Allawi97H_Table(corrs(:)));
    % Allawi97_S(1,5) = mean(Allawi97S_Table(corrs(1,:)));
    % Allawi97_S(2,5) = mean(Allawi97S_Table(corrs(2,:)));
    % Allawi97_S(3,5) = mean(Allawi97S_Table(corrs(3,:)));
    % Allawi97_S(4,5) = mean(Allawi97S_Table(corrs(4,:)));
    % Allawi97_S(5,1) = mean(Allawi97S_Table(corrs(:,1)));
    % Allawi97_S(5,2) = mean(Allawi97S_Table(corrs(:,2)));    
    % Allawi97_S(5,3) = mean(Allawi97S_Table(corrs(:,3)));
    % Allawi97_S(5,4) = mean(Allawi97S_Table(corrs(:,4)));   
    % Allawi97_S(5,5) = mean(Allawi97S_Table(corrs(:)));
    % Rejali21_H(1,5) = mean(Rejali21H_Table(corrs(1,:)));
    % Rejali21_H(2,5) = mean(Rejali21H_Table(corrs(2,:)));
    % Rejali21_H(3,5) = mean(Rejali21H_Table(corrs(3,:)));
    % Rejali21_H(4,5) = mean(Rejali21H_Table(corrs(4,:)));
    % Rejali21_H(5,1) = mean(Rejali21H_Table(corrs(:,1)));
    % Rejali21_H(5,2) = mean(Rejali21H_Table(corrs(:,2)));    
    % Rejali21_H(5,3) = mean(Rejali21H_Table(corrs(:,3)));
    % Rejali21_H(5,4) = mean(Rejali21H_Table(corrs(:,4)));   
    % Rejali21_H(5,5) = mean(Rejali21H_Table(corrs(:)));
    % Rejali21_S(1,5) = mean(Rejali21S_Table(corrs(1,:)));
    % Rejali21_S(2,5) = mean(Rejali21S_Table(corrs(2,:)));
    % Rejali21_S(3,5) = mean(Rejali21S_Table(corrs(3,:)));
    % Rejali21_S(4,5) = mean(Rejali21S_Table(corrs(4,:)));
    % Rejali21_S(5,1) = mean(Rejali21S_Table(corrs(:,1)));
    % Rejali21_S(5,2) = mean(Rejali21S_Table(corrs(:,2)));    
    % Rejali21_S(5,3) = mean(Rejali21S_Table(corrs(:,3)));
    % Rejali21_S(5,4) = mean(Rejali21S_Table(corrs(:,4)));   
    % Rejali21_S(5,5) = mean(Rejali21S_Table(corrs(:)));
    
    seq(seq==15)=5; % substitute numeric value of 'N' with 5
    ind = sub2ind([5 5],seq(1:seq_length-1),seq(2:seq_length));
end

% NN is 3x2 matrix. Columns are DeltaH and DeltaS. 
NN = [sum(Sant04_H(ind)),sum(Sant04_S(ind)); ...
    sum(Allawi97_H(ind)),sum(Allawi97_S(ind)); ...
    sum(Rejali21_H(ind)),sum(Rejali21_S(ind));];

% Corrections: all AT pairs, any GC pairs, symmetry, initiation
if(~nFlag) % ambiguous symbols 'N' not present
    % only AT pairs or any GC pairs?
    if(all((seq ==  1)|(seq == 4)))
        NN =  NN + [Sant04InitTermATwoTerminalAnnotation; Allawi97InitTermATwoTerminalAnnotation; Rejali21InitTermATwoTerminalAnnotation];
    else
        NN = NN +  [Sant04InitTermGCwoTerminalAnnotation; Allawi97InitTermGCwoTerminalAnnotation; Rejali21InitTermGCwoTerminalAnnotation];
    end
    % symmetry
    if(selfCompFlag)
        NN = NN + [Sant04SymmetryCorrection; Allawi97SymmetryCorrection; Rejali21SymmetryCorrection];
    end
    % initiation with terminal  5'
    if(seq(1) == 2 || seq(1) == 3)
        NN(1,:) = NN(1,:) + Sant04InitwTermGC;%GC
        NN(2,:) = NN(2,:) + Allawi97InitwTermGC;%GC
        NN(3,:) = NN(3,:) + Rejali21InitwTermGC;%GC
    elseif(seq(1) == 1 || seq(1) == 4)
        NN(1,:) = NN(1,:) + Sant04InitwTermAT;%AT
        NN(2,:) = NN(2,:) + Allawi97InitwTermAT;%AT
        NN(3,:) = NN(3,:) + Rejali21InitwTermAT;%AT
    end
    
    % initiation with terminal  3'
    if(seq(end) == 2 || seq(end) == 3)
        NN(1,:) = NN(1,:) + Sant04InitwTermGC;%GC
        NN(2,:) = NN(2,:) + Allawi97InitwTermGC;%GC
        NN(3,:) = NN(3,:) + Rejali21InitwTermGC;%GC
    elseif(seq(end) == 1 || seq(end) == 4)
        NN(1,:) = NN(1,:) + Sant04InitwTermAT;
        NN(2,:) = NN(2,:) + Allawi97InitwTermAT;%AT
        NN(3,:) = NN(3,:) + Rejali21InitwTermAT;%AT
    end

    NNdelta=zeros(3,2);

else % 'N' symbols are present

    % only AT pairs or any GC pairs?
    NN1 = NN;  % case when all Ns are G/C
    NN2 = NN; % case when all Ns are A/T
    if(all((seq == 1)|(seq == 4)|(seq == 5)))
        NN1 =  NN1 + [Sant04InitTermATwoTerminalAnnotation; Allawi97InitTermATwoTerminalAnnotation; Rejali21InitTermATwoTerminalAnnotation];

    else
        NN2 = NN2 +  [Sant04InitTermGCwoTerminalAnnotation; Allawi97InitTermGCwoTerminalAnnotation; Rejali21InitTermGCwoTerminalAnnotation];
    end

    % symmetry
    if(selfCompFlag)
        NN1 = NN1 + [Sant04SymmetryCorrection; Allawi97SymmetryCorrection; Rejali21SymmetryCorrection];
        NN2 = NN2 + [Sant04SymmetryCorrection; Allawi97SymmetryCorrection; Rejali21SymmetryCorrection];
    end

    % initiation with terminal 5'(only Sant98/Sant04/Allawi97)
    if(seq(1) == 2 || seq(1) == 3 || seq(1) == 5)
        NN2(1,:) = NN2(1,:) + Sant04InitwTermGC;
        NN2(2,:) = NN2(2,:) + Allawi97InitwTermGC;
        NN2(3,:) = NN2(3,:) + Rejali21InitwTermGC;
    elseif(seq(1) == 1 || seq(1) == 4 || seq(1) == 5)
        NN1(1,:) = NN1(1,:) + Sant04InitwTermAT;
        NN1(2,:) = NN1(2,:) + Allawi97InitwTermAT;
        NN1(3,:) = NN1(3,:) + Rejali21InitwTermAT;
    end

    % initiation with terminal 3'(only Sant98/Sant04/Allawi97)
    if(seq(end) == 2 || seq(end) == 3 || seq(end) == 5)
        NN2(1,:) = NN2(1,:) + Sant04InitwTermGC;
        NN2(2,:) = NN2(2,:) + Allawi97InitwTermGC;
        NN2(3,:) = NN2(3,:) + Rejali21InitwTermGC;
    elseif(seq(end) == 1 || seq(end) == 4 || seq(end) == 5)
        NN1(1,:) = NN1(1,:) + Sant04InitwTermAT;
        NN1(2,:) = NN1(2,:) + Allawi97InitwTermAT;
        NN1(3,:) = NN1(3,:) + Rejali21InitwTermAT;
    end

    NN = (NN1+NN2)/2; % avg
    NNdelta = (max(NN1,NN2)- min(NN1,NN2))/2; % delta level
end
end