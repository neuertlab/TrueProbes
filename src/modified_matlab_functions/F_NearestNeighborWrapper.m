function oligo = F_NearestNeighborWrapper(sequence,varargin)
% This function is modified from matlab function OLIGOPROP to calculates properties of DNA oligonucleotide.
% This adds additional Gibson Free Energy Calculations from different NN models 


% add physiological conditions, different types of salt conc.
% add energy from free loops and secondary/tertirary structure nearest neigbor
%chemical modifications and RNA/DNA, RNA/RNA interactions'
%new 2023 papers

% also proprogate error estimates into solution using dH and dS +- errors.

% also how to add dCp for larger temperature deviations from reference

% salt conc changes of K, Mg, Na, etc on Tm and dG, dH and dS

%model based on GC% and probe length

%model with mismatches internal loops, dangaling ends, and internal
%dangling loops

% kon(T) using uni binding sequential steps

%   OLIGOPROP(SeqNT) returns oligonucleotide properties for a DNA sequence in
%   a structure with the following fields:
%
%       GC: percent GC content for the DNA oligonucleotide. If SEQ
%       contains ambiguous N characters, GC field represents the midpoint
%       value, and its uncertainty is expressed by GCdelta. Ambiguous N
%       characters are considered to potentially be any nucleotide.
%
%       GCdelta: difference between GC (midpoint value) and either the
%       maximum or minimum value GC could assume. The maximum and minimum
%       values are calculated by assuming all N characters are G/C or not G/C,
%       respectively. Therefore, GCdelta defines the possible range of GC content.
%
%       Hairpins: H-by-length(SeqNT) matrix of characters displaying all
%       potential hairpin structures for the sequence SeqNT. Each row is a
%       potential hairpin structure of the sequence, with the hairpin
%       forming nucleotides designated by capital letters. H is the number
%       of potential hairpin structures for the sequence. Ambiguous N
%       characters are considered to potentially complement any nucleotide.
%
%       Dimers: D-by-length(SeqNT) matrix of characters displaying all
%       potential dimers for the sequence SeqNT. Each row is a potential
%       dimer of the sequence, with the self-dimerizing nucleotides
%       designated by capital letters. D is the number of potential dimers
%       for the sequence. Ambiguous N characters are considered to
%       potentially complement any nucleotide.
%
%       MolWeight: Molecular weight of the DNA oligonucleotide. If SeqNT
%       contains ambiguous N characters, MolWeight is the midpoint value,
%       and its uncertainty is expressed by MolWeightdelta.  Ambiguous N
%       characters in SeqNT are considered to potentially be any
%       nucleotide.
%
%       MolWeightdelta: difference between MolWeight (midpoint value) and
%       either the maximum or minimum value MolWeight could assume. The
%       maximum and minimum values are calculated by assuming all N
%       characters are G or C, respectively. Therefore, MolWeightdelta
%       defines the possible range of molecular weight for SeqNT.
%
%       Tm: vector with melting temperature values, in degrees Celsius,
%       calculated by six different methods, listed in the following order:
%            Basic (Marmur et al., 1962)
%            Salt adjusted (Howley et al., 1979)
%            Nearest-neighbor (Breslauer et al., 1986)
%            Nearest-neighbor (SantaLucia Jr. et al., 1996)
%            Nearest-neighbor (SantaLucia Jr., 1998)
%            Nearest-neighbor (Sugimoto et al., 1996)
%       Ambiguous N characters in SeqNT are considered to potentially be
%       any nucleotide. If SeqNT contains ambiguous N characters, Tm is the
%       midpoint value, and its uncertainty is expressed by Tmdelta.
%
%       Tmdelta: vector containing the differences between Tm (midpoint
%       value) and either the maximum or minimum value Tm could assume for
%       each of the six methods. Therefore, Tmdelta defines the possible
%       range of melting temperatures for SeqNT.
%
%       Thermo: 4-by-3 matrix of thermodynamic calculations. The rows
%       correspond to nearest-neighbor parameters from:
%            Breslauer et al., 1986
%            SantaLucia Jr. et al., 1996
%            SantaLucia Jr., 1998
%            Sugimoto et al., 1996
%       The columns correspond to:
%            delta H - enthalpy in kilocalories per mole, kcal/mol
%            delta S - entropy in calories per mole-degrees Kelvin, cal/(K*mol)
%            delta G - free energy in kilocalories per mole, kcal/mol
%       Ambiguous N characters in SeqNT are considered to potentially be
%       any nucleotide. If SeqNT contains ambiguous N characters, Thermo is
%       the midpoint value, and its uncertainty is expressed by Thermodelta.
%
%       Thermodelta: 4-by-3 matrix containing the differences between
%       Thermo (midpoint value) and either the maximum or minimum value
%       Thermo could assume for each calculation and method. Therefore,
%       Thermodelta defines the possible range of thermodynamic
%       values for SeqNT.
%
%   OLIGOPROP(...,'SALT',SALT) specifies a salt concentration in moles per
%   liter for melting temperature calculations. The default is 0.05M.
%
%   OLIGOPROP(...,'TEMP',TEMP) specifies the temperature for nearest
%   neighbor calculations of free energy. The default is 25 degrees Celsius.
%
%   OLIGOPROP(...,'PRIMERCONC',PRIMERCONC) specifies what concentration
%   will be used for melting temperature calculations. The default is 50e-6M.
%
%   OLIGOPROP(...,'HPBASE',HAIRPINBASE) specifies the minimum number of
%   paired bases forming the neck of the hairpin. The default is 4.
%
%   OLIGOPROP(...,'HPLOOP',HAIRPINLOOP) specifies the minimum number of
%   bases forming a hairpin. The default is 2.
%
%   OLIGOPROP(...,'DIMERLENGTH',LENGTH) the minimum number of aligned bases
%   between the sequence and its reverse. The default is 4.
%
%   Examples:
%
%       S1 = oligoprop(randseq(25))
%       S2 = oligoprop('ACGTAGAGGACGTN')
%
%   See also NTDENSITY, PALINDROMES, RANDSEQ, ISOELECTRIC,
%   MOLWEIGHT.

%   Copyright 2005-2012 The MathWorks, Inc.

%   References:
%   [1] Panjkovich A., Melo F., "Comparison of different melting
%       temperature calculation methods for short DNA sequences",
%       Bioinformatics 21(6):711-722 (2004).
%   [2] Breslauer K.J., Frank R., Blocker H., Marky L.A., "Predicting DNA
%       duplex stability from the base sequence" PNAS 83:3746-3750 (1986).
%   [3] Sugimoto N., Nakano S., Yoneyama M., Honda K., "Improved
%       thermodynamic parameters and helix initiation factor to predict
%       stability of DNA duplexes", Nucleic Acids Research
%       24(22):4501-4505 (1996).
%   [4] SantaLucia J., "A unified view of polymer, dumbbell, and
%       oligonucleotide DNA nearest-neighbor thermodynamics", PNAS
%       95:1460-1465 (1998).
%   [5] SantaLucia, J., Allawi, H.T., Seneviratne P.A., "Improved
%       Nearest-Neighbor Parameters for Predicting DNA Duplex Stability",
%       Biochemistry 35:3555-3562 (1996).
%   [6] Howley P.M., Israel M.A., Law M.F., Martin M.A., "A rapid method
%       for detecting and mapping homology between heterologous DNAs.
%       Evaluation of polyomavirus genomes", Journal of Biological
%       Chemistry, 254(11):4876-4883 (1979).
%   [7] Marmur J., Doty P., "Determination of the base composition of
%       deoxyribonucleic acid from its thermal denaturation temperature",
%       JMB 5:109-118 (1962).
%   [8] Chen S.H., Lin C.Y., Lo C.Z., Cho C.S., Hsiung C.A., "Primer Design
%       Assistant (PDA): a Web-based Primer Design Tool" Nucleic Acids
%       Research 31:3751-3754 (2003).
%   [9] http://www.basic.northwestern.edu/biotools/oligocalc.html for
%       weight calculations

% determine the format of the sequence
if nargin > 0
    sequence = convertStringsToChars(sequence);
end

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if isstruct(sequence)
    sequence = bioinfoprivate.seqfromstruct(sequence);
end
seq = upper(sequence);


if(isnumeric(seq))
    seq = upper(int2nt(seq));
else
    seq = upper(seq);
end

if (~all(seq=='A'|seq=='C'|seq=='G'|seq=='T'|seq=='N'))
    error(message('bioinfo:oligoprop:IncorrectSequenceType'));
elseif(any(seq=='N'))
    nFlag=1;
else
    nFlag=0;
end

% defaults parameters
temp = 25;            % temperature in Celsius
salt = 0.05;          % salt concentration in moles per liter (M)
primerConc= 50e-6;    % concentration of primers in mole per liter (M)
weight = [313.21 289.18 329.21 304.2]; % molecular weight A,C,G,T respectively

% check arguments
if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:oligoprop:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'SALT','PRIMERCONC','TEMP'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strcmpi(upper(pname), okargs)); %#ok
        if isempty(k)
            error(message('bioinfo:oligoprop:UnknownParameterName', pname));
        else
            switch(k)
                case 1  % SALT
                    if (pval > -1 && isnumeric(pval) && isreal(pval))
                        salt =  pval;
                    else
                       % error(message('bioinfo:oligoprop:BadSaltConcentrationParam'));
                    end
                case 2  % PRIMERCONC
                    if (pval > -1 && isnumeric(pval) && isreal(pval))
                        primerConc =  pval;
                    else
                       % error(message('bioinfo:oligoprop:BadConcentrationParam'));
                    end
                case 3  % Temp
                    if (isnumeric(pval) && isreal(pval))
                        temp =  pval;
                    else
                    %    error(message('bioinfo:oligoprop:BadtempParam'));
                    end
            end
        end
    end
end

% compute sequence properties
numSeq = double(nt2int(seq));
baseNum = [sum(numSeq == 1) sum(numSeq == 2) sum(numSeq == 3) sum(numSeq == 4) sum(numSeq == 15)];

if(numel(numSeq)<8)
    %warning(message('bioinfo:oligoprop:SeqLengthTooShort'));
end

if (~nFlag) % no ambiguous symbols 'N'

    GC = 100 * ((baseNum(2) + baseNum(3)) / length(numSeq));
    GCdelta = 0;

    w = sum(baseNum(1:4) .* weight) - 61.96;
    wdelta = 0;

    [tm tmdelta NN NNdelta]  = get_tm_NN(numSeq, baseNum, salt, primerConc, nFlag);
    NN(:,3) = NN(:,1) - ((temp+273.15)  .* (NN(:,2)./1000)); % DeltaG
    NNdelta(:,3) = zeros(7,1);

else % occurrences of ambiguous N symbols
    % warning('bioinfo:oligoprop:ambiguousCheck', 'Ambiguous symbols N in input sequence.');

    % average between case when all Ns are C/G and case when all Ns are A/T
    GC = 100 * ((baseNum(2) + baseNum(3) + baseNum(5)/2)/ length(numSeq));
    GCdelta = 100 * 1/2 * (baseNum(5)/length(numSeq));

    % avg between case when all 'N' are G and case when all 'N' are C
    w = sum(baseNum(1:4) .* weight) - 61.96 + baseNum(5) * (weight(3) + weight(2)) / 2;
    wdelta = baseNum(5) * (weight(3) - weight(2)) / 2;

    [tm tmdelta NN NNdelta]  = get_tm_NN(numSeq, baseNum, salt, primerConc, nFlag);
    NN(:,3) = NN(:,1) - ((temp+273.15)  .* (NN(:,2)./1000)); % DeltaG
    NNdelta(:,3) = NNdelta(:,1) - ((temp+273.15)  .* (NNdelta(:,2) ./1000));
end

% build output structure
oligo.GC = GC;
oligo.GCdelta = GCdelta;
oligo.MolWeight = w;
oligo.MolWeightdelta = wdelta;
oligo.Tm = tm;
oligo.Tmdelta = tmdelta;
oligo.Thermo = NN;
oligo.Thermodelta = NNdelta;
end

% FUNCTIONS

% calculate melting temperature and thermo values (37 degrees C)
function [tm tmdelta NN NNdelta]  = get_tm_NN(numSeq, baseNum, salt, primerConc, nFlag)
% melting temperatures and thermodynamic values are returned as average +/- delta level.
% If no ambiguous symbols are present, the delta level is zero.

selfCompFlag = all(5-numSeq == numSeq(end:-1:1));

if(selfCompFlag) % self complementary sequence
    b = 1; % correction value for nearest neighbor melting temperature calculation
else
    b = 4;
end

tmdelta = zeros(1,9);

if (~nFlag) % no ambiguous symbols
    if (sum(baseNum)<14)
        basic = 2 * (baseNum(1) + baseNum(4)) + 4 * (baseNum(2) + baseNum(3)); % TM BASIC [9]
        saltadj = basic - 16.6 * log10(0.05) + 16.6 * log10(salt); % TM SALT ADJUSTED [9]
    else
        basic = 64.9 + (41 * ((baseNum(3) + baseNum(2) - 16.4) / sum(baseNum))); %TM BASIC [1],[9]
        saltadj = 100.5 + (41 * ((baseNum(3) + baseNum(2))/ sum(baseNum))) - (820/sum(baseNum)) + (16.6 * log10(salt)); %TM SALT ADJUSTED [9]
    end

    [NN, NNdelta] = near_neigh(numSeq, length(numSeq), selfCompFlag, nFlag);
    tm = (NN(:,1) * 1000 ./ (NN(:,2) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15; %TM NEAREST NEIGHBOR

else % occurrences of 'N'
    if(sum(baseNum)<14)
        basic = 2 * (baseNum(1) + baseNum(4)) + 4 * (baseNum(2) + baseNum(3)) + 3 * baseNum(5); % TM BASIC [9]
        tmdelta(1) = baseNum(5);
        saltadj = basic - 16.6 * log10(0.05) + 16.6 * log10(salt); % TM SALT ADJUSTED [9]
        tmdelta(2) = 3 * baseNum(5);
    else
        basic = 64.9 + (41 * (baseNum(3) + baseNum(2) - 16.4 + baseNum(5)/2)/sum(baseNum)); % avg TM BASIC
        tmdelta(1)= 1/2 * 41 * (baseNum(5)/sum(baseNum));
        saltadj = 100.5 - (820/sum(baseNum)) + (16.6 * log10(salt)) +  41 * ((baseNum(2) + baseNum(3)+ baseNum(5)/2)/sum(baseNum));% avg TM SALT ADJUSTED
        tmdelta(2)= 41 * 1/2 * (baseNum(5) / sum(baseNum));
    end
    [NN, NNdelta] = near_neigh(numSeq, length(numSeq), selfCompFlag, nFlag);
    tm = (((NN(:,1)+ NNdelta(:,1)) * 1000 ./ ((NN(:,2)+ NNdelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15 + ...
        ((NN(:,1)- NNdelta(:,1)) * 1000 ./ ((NN(:,2)- NNdelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15)* 1/2  ; % NEAREST NEIGHBOR
    tmdelta(3:9)=(((NN(:,1)+ NNdelta(:,1)) * 1000 ./ ((NN(:,2)+ NNdelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15 - ...
        (((NN(:,1)- NNdelta(:,1)) * 1000 ./ ((NN(:,2)- NNdelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15)) * 1/2 ;
end
tm = [basic; saltadj; tm]';
end
% compute thermo values using Nearest Neighbor methods
function [NN NNdelta] = near_neigh(seq, seq_length, selfCompFlag, nFlag)
Sant04H_Table = [-7.6 -7.2 -7.2 -8.5 -8.4 -7.8 -8.2 -10.6 -9.8 -8.0];%kcal/mol
Sant04S_Table = [-21.3 -20.4 -21.3 -22.7 -22.4 -21.0 -22.2 -27.2 -24.4 -19.9];%cal/molK
Allawi97H_Table = [-7.9 -7.2 -7.2 -8.5 -8.4 -7.8 -8.2 -10.6 -9.8 -8.0];%kcal/mol
Allawi97S_Table = [-22.2 -20.4 -21.3 -22.7 -22.4 -21.0 -22.2 -27.2 -24.4 -19.9];%cal/molK
Allawi97Herr_Table = [0.2 0.7 0.9 0.6 0.5 0.6 0.6 0.6 0.4 0.9];%kcal/mol
Allawi97Serr_Table = [0.8 2.4 2.4 2.0 2.0 2.0 1.7 2.6 2.0 1.8];%cal/molK
Allawi97InitTermGC = [0.1 -2.8];%H %S
Allawi97InitGCErr = [1.1 0.2];%H %S
Allawi97InitTermAT = [2.3 4.1];%H %S
Allawi97InitTermATErr = [1.3 0.2];%H %S
Rejali21H_Table = [-9.2 -8.6 -5.6 -11.9 -10.2 -8.1 -14.5 -11.2 -9.8 14.8 1.0];
Rejali21S_Table = [-26.6 -24.4 -16.2 -33.9 -25.9 -29.2 -21.9 -40.6 -29.2 -26.2 41.7 2.6];
Rejali21Herr_Table = [0.6 1.7 1.4 1.4 1.6 1.2 1.3 1.5 1.7 1.3 4.2 1.0];
Rejali21Serr_Table = [1.8 5.3 4.3 4.4 5.0 3.6 4.0 4.9 5.5 4.2 13.3 3.1];

Rejali21InitTermGC = [14.8 41.7];
Rejali21InitTermAT = [1.0 2.6];



corrs =[1 5  7  2;4 10 8  6;7 9  10 5;3 7  4  1];  
Sant04_H = Sant04H_Table(corrs);Sant04_S = Sant04S_Table(corrs);  
Allawi97_H = Allawi97H_Table(corrs);Allawi97_S = Allawi97S_Table(corrs); 
Rejali21_H = Rejali21H_Table(corrs);Rejali21_S = Rejali21S_Table(corrs);  
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
if (~nFlag)
    % nearest neighbor parameters from: Panjkovich and Melo, Bioinformatics  Vol 21 no 6 pp 711-722 2004 [1]
    % rows corresponds to A,C,G,T respectively; columns correspond to A,C,G,T respectively
    Bres86_H = [-9.1,-6.5,-7.8,-8.6,;-5.8,-11,-11.9,-7.8,;-5.6,-11.1,-11,-6.5,;-6,-5.6,-5.8,-9.1,];
    Bres86_S = [-24,-17.3,-20.8,-23.9,;-12.9,-26.6,-27.8,-20.8,;-13.5,-26.7,-26.6,-17.3,;-16.9,-13.5,-12.9,-24,];
    Sant96_H = [-8.4,-8.6,-6.1,-6.5,;-7.4,-6.7,-10.1,-6.1,;-7.7,-11.1,-6.7,-8.6,;-6.3,-7.7,-7.4,-8.4,];
    Sant96_S = [-23.6,-23,-16.1,-18.8,;-19.3,-15.6,-25.5,-16.1,;-20.3,-28.4,-15.6,-23,;-18.5,-20.3,-19.3,-23.6,];
    Sant98_H = [-7.9,-8.4,-7.8,-7.2,;-8.5,-8,-10.6,-7.8,;-8.2,-9.8,-8,-8.4,;-7.2,-8.2,-8.5,-7.9,];
    Sant98_S = [-22.2,-22.4,-21,-20.4,;-22.7,-19.9,-27.2,-21,;-22.2,-24.4,-19.9,-22.4,;-21.3,-22.2,-22.7,-22.2,];
    Sugi96_H = [-8,-9.4,-6.6,-5.6,;-8.2,-10.9,-11.8,-6.6,;-8.8,-10.5,-10.9,-9.4,;-6.6,-8.8,-8.2,-8,];
    Sugi96_S = [-21.9,-25.5,-16.4,-15.2,;-21,-28.4,-29,-16.4,;-23.5,-26.4,-28.4,-25.5,;-18.4,-23.5,-21,-21.9,];

    ind = sub2ind([4 4],seq(1:seq_length-1),seq(2:seq_length));
else
    % nearest neighbor parameters as in [1] with added average values for
    % all possible combinations involving 'N'
    % rows corresponds to A,C,G,T,N respectively; columns correspond to
    % A,C,G,T,N respectively
    Bres86_H = [-9.1,-6.5,-7.8,-8.6,-8;-5.8,-11,-11.9,-7.8,-9.125; -5.6,-11.1,-11,-6.5,-8.55;-6,-5.6,-5.8,-9.1,-6.625;-6.625,-8.55,-9.125,-8,-8.075];
    Bres86_S = [-24,-17.3,-20.8,-23.9,-21.5;-12.9,-26.6,-27.8,-20.8,-22.025;-13.5,-26.7,-26.6,-17.3,-21.025;-16.9,-13.5,-12.9,-24,-16.825;-16.825,-21.025,-22.025,-21.5,-20.3438];
    Sant96_H = [-8.4,-8.6,-6.1,-6.5,-7.4;-7.4,-6.7,-10.1,-6.1,-7.575;-7.7,-11.1,-6.7,-8.6,-8.525;-6.3,-7.7,-7.4,-8.4,-7.45;-7.45,-8.525,-7.575,-7.4,-7.7375];
    Sant96_S = [-23.6,-23,-16.1,-18.8,-20.375;-19.3,-15.6,-25.5,-16.1,-19.125;-20.3,-28.4,-15.6,-23,-21.825;-18.5,-20.3,-19.3,-23.6,-20.425;-20.425,-21.825,-19.125,-20.375,-20.4375];
    Sant98_H = [-7.9,-8.4,-7.8,-7.2,-7.825;-8.5,-8,-10.6,-7.8,-8.725;-8.2,-9.8,-8,-8.4,-8.6;-7.2,-8.2,-8.5,-7.9,-7.95;-7.95,-8.6,-8.725,-7.825,-8.275];
    Sant98_S = [-22.2,-22.4,-21,-20.4,-21.5;-22.7,-19.9,-27.2,-21,-22.7;-22.2,-24.4,-19.9,-22.4,-22.225;-21.3,-22.2,-22.7,-22.2,-22.1;-22.1,-22.225,-22.7,-21.5,-22.1313];
    Sugi96_H = [-8,-9.4,-6.6,-5.6,-7.4;-8.2,-10.9,-11.8,-6.6,-9.375;-8.8,-10.5,-10.9,-9.4,-9.9;-6.6,-8.8,-8.2,-8,-7.9;-7.9,-9.9,-9.375,-7.4,-8.6437];
    Sugi96_S = [-21.9,-25.5,-16.4,-15.2,-19.75;-21,-28.4,-29,-16.4,-23.7;-23.5,-26.4,-28.4,-25.5,-25.95;-18.4,-23.5,-21,-21.9,-21.2;-21.2,-25.95,-23.7,-19.75,-22.65];
    Sant04_H(1,5) = mean(Sant04H_Table(corrs(1,:)));
    Sant04_H(2,5) = mean(Sant04H_Table(corrs(2,:)));
    Sant04_H(3,5) = mean(Sant04H_Table(corrs(3,:)));
    Sant04_H(4,5) = mean(Sant04H_Table(corrs(4,:)));
    Sant04_H(5,1) = mean(Sant04H_Table(corrs(:,1)));
    Sant04_H(5,2) = mean(Sant04H_Table(corrs(:,2)));    
    Sant04_H(5,3) = mean(Sant04H_Table(corrs(:,3)));
    Sant04_H(5,4) = mean(Sant04H_Table(corrs(:,4)));   
    Sant04_H(5,5) = mean(Sant04H_Table(corrs(:)));
    Sant04_S(1,5) = mean(Sant04S_Table(corrs(1,:)));
    Sant04_S(2,5) = mean(Sant04S_Table(corrs(2,:)));
    Sant04_S(3,5) = mean(Sant04S_Table(corrs(3,:)));
    Sant04_S(4,5) = mean(Sant04S_Table(corrs(4,:)));
    Sant04_S(5,1) = mean(Sant04S_Table(corrs(:,1)));
    Sant04_S(5,2) = mean(Sant04S_Table(corrs(:,2)));    
    Sant04_S(5,3) = mean(Sant04S_Table(corrs(:,3)));
    Sant04_S(5,4) = mean(Sant04S_Table(corrs(:,4)));   
    Sant04_S(5,5) = mean(Sant04S_Table(corrs(:)));
    
    Allawi97_H(1,5) = mean(Allawi97H_Table(corrs(1,:)));
    Allawi97_H(2,5) = mean(Allawi97H_Table(corrs(2,:)));
    Allawi97_H(3,5) = mean(Allawi97H_Table(corrs(3,:)));
    Allawi97_H(4,5) = mean(Allawi97H_Table(corrs(4,:)));
    Allawi97_H(5,1) = mean(Allawi97H_Table(corrs(:,1)));
    Allawi97_H(5,2) = mean(Allawi97H_Table(corrs(:,2)));    
    Allawi97_H(5,3) = mean(Allawi97H_Table(corrs(:,3)));
    Allawi97_H(5,4) = mean(Allawi97H_Table(corrs(:,4)));   
    Allawi97_H(5,5) = mean(Allawi97H_Table(corrs(:)));
    Allawi97_S(1,5) = mean(Allawi97S_Table(corrs(1,:)));
    Allawi97_S(2,5) = mean(Allawi97S_Table(corrs(2,:)));
    Allawi97_S(3,5) = mean(Allawi97S_Table(corrs(3,:)));
    Allawi97_S(4,5) = mean(Allawi97S_Table(corrs(4,:)));
    Allawi97_S(5,1) = mean(Allawi97S_Table(corrs(:,1)));
    Allawi97_S(5,2) = mean(Allawi97S_Table(corrs(:,2)));    
    Allawi97_S(5,3) = mean(Allawi97S_Table(corrs(:,3)));
    Allawi97_S(5,4) = mean(Allawi97S_Table(corrs(:,4)));   
    Allawi97_S(5,5) = mean(Allawi97S_Table(corrs(:)));
    
    Rejali21_H(1,5) = mean(Rejali21H_Table(corrs(1,:)));
    Rejali21_H(2,5) = mean(Rejali21H_Table(corrs(2,:)));
    Rejali21_H(3,5) = mean(Rejali21H_Table(corrs(3,:)));
    Rejali21_H(4,5) = mean(Rejali21H_Table(corrs(4,:)));
    Rejali21_H(5,1) = mean(Rejali21H_Table(corrs(:,1)));
    Rejali21_H(5,2) = mean(Rejali21H_Table(corrs(:,2)));    
    Rejali21_H(5,3) = mean(Rejali21H_Table(corrs(:,3)));
    Rejali21_H(5,4) = mean(Rejali21H_Table(corrs(:,4)));   
    Rejali21_H(5,5) = mean(Rejali21H_Table(corrs(:)));
    Rejali21_S(1,5) = mean(Rejali21S_Table(corrs(1,:)));
    Rejali21_S(2,5) = mean(Rejali21S_Table(corrs(2,:)));
    Rejali21_S(3,5) = mean(Rejali21S_Table(corrs(3,:)));
    Rejali21_S(4,5) = mean(Rejali21S_Table(corrs(4,:)));
    Rejali21_S(5,1) = mean(Rejali21S_Table(corrs(:,1)));
    Rejali21_S(5,2) = mean(Rejali21S_Table(corrs(:,2)));    
    Rejali21_S(5,3) = mean(Rejali21S_Table(corrs(:,3)));
    Rejali21_S(5,4) = mean(Rejali21S_Table(corrs(:,4)));   
    Rejali21_S(5,5) = mean(Rejali21S_Table(corrs(:)));
    
    seq(seq==15)=5; % substitute numeric value of 'N' with 5
    ind = sub2ind([5 5],seq(1:seq_length-1),seq(2:seq_length));
end

% NN is 5x2 matrix. Columns are DeltaH and DeltaS. Rows correspond to
% methods by Bres86, SantaLucia96, SantaLucia98, Sugimoto96, and
% SantaLucia04
NN = [sum(Bres86_H(ind)),sum(Bres86_S(ind)); ...
    sum(Sant96_H(ind)),sum(Sant96_S(ind)); ...
    sum(Sant98_H(ind)),sum(Sant98_S(ind)); ...
    sum(Sugi96_H(ind)),sum(Sugi96_S(ind)); ...
    sum(Sant04_H(ind)),sum(Sant04_S(ind)); ...
    sum(Allawi97_H(ind)),sum(Allawi97_S(ind)); ...
    sum(Rejali21_H(ind)),sum(Rejali21_S(ind));];

% Corrections: all AT pairs, any GC pairs, symmetry, initiation
if(~nFlag) % ambiguous symbols 'N' not present

    % only AT pairs or any GC pairs?
    if(all((seq == 1)|(seq == 4)))
        NN =  NN + [0 -20.13; 0 -9; 0 0; 0.6 -9; 0 0; 0 0; 0 0];
    else
        NN = NN +  [0 -16.77; 0 -5.9; 0 0; 0.6 -9; 0 0; 0 0; 0 0];
    end

    % symmetry
    if(selfCompFlag)
        NN = NN + [0 -1.34 ; 0 -1.4 ;0 -1.4; 0 -1.4 ; 0 -1.4; 0 -1.4; 0 -1.4];
    end

    % initiation with terminal  5'
    if(seq(1) == 2 || seq(1) == 3)
        NN(3,:) = NN(3,:) + [0.1 -2.8];%GC
        NN(5,:) = NN(5,:) + [0.2 -5.7];%GC
        NN(6,:) = NN(6,:) + Allawi97InitTermGC;%GC
        NN(7,:) = NN(7,:) + Rejali21InitTermGC;%GC
    elseif(seq(1) == 1 || seq(1) == 4)
        NN(3,:) = NN(3,:) + [2.3 4.1];%AT
        NN(5,:) = NN(5,:) + [2.2 6.9];%AT
        NN(6,:) = NN(6,:) + Allawi97InitTermAT;%AT
        NN(7,:) = NN(7,:) + Rejali21InitTermAT;%AT
    end
    
    % initiation with terminal  3'
    if(seq(end) == 2 || seq(end) == 3)
        NN(3,:) = NN(3,:) + [0.1 -2.8];
        NN(5,:) = NN(5,:) + [0.2 -5.7];
        NN(6,:) = NN(6,:) + Allawi97InitTermGC;%GC
        NN(7,:) = NN(7,:) + Rejali21InitTermGC;%GC
    elseif(seq(end) == 1 || seq(end) == 4)
        NN(3,:) = NN(3,:) + [2.3 4.1];
        NN(5,:) = NN(5,:) + [2.2 6.9];
        NN(6,:) = NN(6,:) + Allawi97InitTermAT;%AT
        NN(7,:) = NN(7,:) + Rejali21InitTermAT;%AT
    end

    NNdelta=zeros(7,2);

else % 'N' symbols are present

    % only AT pairs or any GC pairs?
    NN1 = NN;  % case when all Ns are G/C
    NN2 = NN; % case when all Ns are A/T
    if(all((seq == 1)|(seq == 4)|(seq == 5)))
        NN1 = NN1 + [0 -20.13; 0 -9; 0 0; 0.6 -9; 0 0; 0 0; 0 0];
    else
        NN2 = NN2 +  [0 -16.77; 0 -5.9; 0 0; 0.6 -9; 0 0; 0 0; 0 0];
    end

    % symmetry
    if(selfCompFlag)
        NN1 = NN1 + [0 -1.34 ; 0 -1.4 ;0 -1.4; 0 -1.4 ; 0 -1.4; 0 -1.4; 0 -1.4];
        NN2 = NN2 + [0 -1.34 ; 0 -1.4 ;0 -1.4; 0 -1.4 ; 0 -1.4; 0 -1.4; 0 -1.4];
    end

    % initiation with terminal 5'(only Sant98/Sant04/Allawi97)
    if(seq(1) == 2 || seq(1) == 3 || seq(1) == 5)
        NN2(3,:) = NN2(3,:) + [0.1 -2.8];
        NN2(5,:) = NN2(5,:) + [0.2 -5.7];
        NN2(6,:) = NN2(6,:) + Allawi97InitTermGC;
        NN2(7,:) = NN2(7,:) + Rejali21InitTermGC;
    elseif(seq(1) == 1 || seq(1) == 4 || seq(1) == 5)
        NN1(3,:) = NN1(3,:) + [2.3 4.1];
        NN1(5,:) = NN1(5,:) + [2.2 6.9];
        NN1(6,:) = NN1(6,:) + Allawi97InitTermAT;
        NN1(7,:) = NN1(7,:) + Rejali21InitTermAT;
    end

    % initiation with terminal 3'(only Sant98/Sant04/Allawi97)
    if(seq(end) == 2 || seq(end) == 3 || seq(end) == 5)
        NN2(3,:) = NN2(3,:) + [0.1 -2.8];
        NN2(5,:) = NN2(5,:) + [0.2 -5.7];
        NN2(6,:) = NN2(6,:) + Allawi97InitTermGC;
        NN2(7,:) = NN2(7,:) + Rejali21InitTermGC;
    elseif(seq(end) == 1 || seq(end) == 4 || seq(end) == 5)
        NN1(3,:) = NN1(3,:) + [2.3 4.1];
        NN1(5,:) = NN1(5,:) + [2.2 6.9];
        NN1(6,:) = NN1(6,:) + Allawi97InitTermAT;
        NN1(7,:) = NN1(7,:) + Rejali21InitTermAT;
    end

    NN = (NN1+NN2)/2; % avg
    NNdelta = (max(NN1,NN2)- min(NN1,NN2))/2; % delta level
end
end