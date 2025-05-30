function dG = F_TransitionStateDeltaGibson_V2(seq1,seq2,C,T_vector)
Nmodels = 5;
Ntemps = length(T_vector);
dG = zeros(Nmodels,Ntemps);
format('long');
%upper strand is RNA with 5->3 and lower strand DNA (3'->5')
%Seq1 5->3
%Seq2 5->3
%DNA/DNA;
%RNA/DNA;
%RNA/RNA;
Np = 400;
R = 0.001987204259;%gas constant [Energy Kcal/mol K
kb = 1.380649*10^-23;%bolzman constant J/K
h = 6.62607015*10^-34;%planks constant J/Hz
approx = 4;
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

%NN Model
%Length and GC% approximation with salt conc.

%k = k*kbT/h*exp(-dGtrans/RT)*(c0)^`-m
%dG_trans = dH_trans-T*dS_trans
%dG = dG_transA-dG_transB
%dG_total_trans = dG_trans_intiation +

%bNNs are reported (5'->3')/(3'-5) at 37C
%Params  salt is Na
%AA/TT AT/TA TA/AT CA/GT GT/CA CT/GA GA/CT CG/GC GC/CG GG/CC Initiation Terminal_AT    
DeltaH_Transition_Dissociation_1M_Mean = [9.2 8.6 5.6 11.9 9.6 10.2 8.1 14.5 11.2 9.8 -14.8 -1.0];%kcal/mol
DeltaH_Transition_Dissociation_1M_Err = [0.6 1.7 1.4 1.4 1.6 1.2 1.3 1.5 1.7 1.3 4.2 1.0];
DeltaS_Transition_Dissociation_1M_Mean = [26.5 24.2 15.8 33.5 25.9 28.9 21.7 40.7 29.2 266.6 -66.8 -2.3];%cal/molK
DeltaS_Transition_Dissociation_1M_Err = [1.8 5.3 4.3 4.4 5 3.6 4.0 4.9 5.5 4.2 13.3 3.1];%cal/molK     
DeltaG_Transition_Dissociation_1M_Mean = [1.03 1.09 0.71 1.51 1.56 1.26 1.37 1.87 2.13 1.57 5.87 -0.29];%kcal/mol
DeltaG_Transition_Dissociation_1M_Err = [0.04 0.09 0.10 0.10 0.10 0.09 0.08 0.12 0.08 0.09 0.32 0.05];%kcal/mol
DeltaS_Transition_Association_1M_Mean = [-0.08 -0.16 -0.41 -0.35 -0.02 -0.29 -0.22 0.08 -0.02 0.47 -25.1 0.32];%cal/molK
DeltaS_Transition_Association_1M_Err = [0.11 0.23 0.23 0.21 0.20 0.27 0.20 0.29 0.24 0.16 0.7 0.15];%cal/molK     
DeltaG_Transition_Association_1M_Mean = [0.03 0.05 0.13 0.11 0.01 0.09 0.07 -0.02 0.01 -0.15 7.79 -0.10];%kcal/mol
DeltaG_Transition_Association_1M_Err = [0.03 0.07 0.07 0.07 0.06 0.08 0.06 0.09 0.07 0.05 0.21 0.05];%kcal/mol
%2.2mM MgCl2 2.2*10^-3M  2.2*10^3 microM  
DeltaH_Transition_Dissociation_4400microM_Mean = [9.6 11.0 3.7 10.0 12.0 9.8 7.4 13.2 11.6 11.3 -10.5 -3.0];%kcal/mol
DeltaH_Transition_Dissociation_4400microM_Err = [0.5 1.5 1.6 1.2 1.3 1.3 1.1 1.5 1.3 1.2 3.9 0.7];
DeltaS_Transition_Dissociation_4400microM_Mean = [28.0 32.8 10.0 27.6 34.0 27.2 20.1 36.0 31.3 31.2 -54.0 -8.7];%cal/molK
DeltaS_Transition_Dissociation_4400microM_Err = [1.5 4.8 5.0 3.8 4 4.2 3.7 5.0 4.2 3.5 12.3 2.3];%cal/molK     
DeltaG_Transition_Dissociation_4400microM_Mean = [0.93 0.82 0.63 1.46 1.46 1.33 1.20 2.00 1.88 1.59 6.21 -0.28];%kcal/mol
DeltaG_Transition_Dissociation_4400microM_Err = [0.05 0.12 0.13 0.09 0.10 0.12 0.09 0.14 0.10 0.11 0.37 0.06];%kcal/mol
DeltaS_Transition_Association_4400microM_Mean = [-0.37 -0.24 -0.58 -0.24 -0.24 -0.52 -0.22 0.11 0.02 0.49 -27.5 0.31];%cal/molK
DeltaS_Transition_Association_4400microM_Err = [0.12 0.34 0.34 0.27 0.26 0.37 0.28 0.38 0.27 0.24 0.9 0.22];%cal/molK     
DeltaG_Transition_Association_4400microM_Mean = [0.11 0.08 0.18 0.08 0.07 0.16 0.07 -0.03 -0.01 -0.15 8.51 -0.10];%kcal/mol
DeltaG_Transition_Association_4400microM_Err = [0.04 0.11 0.10 0.90 0.08 0.12 0.09 0.13 0.08 0.07 0.28 0.07];%kcal/mol

DeltaH_Transition_Dissociation_1M_Upper = DeltaH_Transition_Dissociation_1M_Mean + DeltaH_Transition_Dissociation_1M_Err;
DeltaS_Transition_Dissociation_1M_Upper = DeltaS_Transition_Dissociation_1M_Mean + DeltaS_Transition_Dissociation_1M_Err;   
DeltaG_Transition_Dissociation_1M_Upper = DeltaG_Transition_Dissociation_1M_Mean + DeltaG_Transition_Dissociation_1M_Err;
DeltaS_Transition_Association_1M_Upper = DeltaS_Transition_Association_1M_Mean + DeltaS_Transition_Association_1M_Err;
DeltaG_Transition_Association_1M_Upper = DeltaG_Transition_Association_1M_Mean + DeltaG_Transition_Association_1M_Err;
DeltaH_Transition_Dissociation_1M_Lower = DeltaH_Transition_Dissociation_1M_Mean - DeltaH_Transition_Dissociation_1M_Err;
DeltaS_Transition_Dissociation_1M_Lower = DeltaS_Transition_Dissociation_1M_Mean - DeltaS_Transition_Dissociation_1M_Err;   
DeltaG_Transition_Dissociation_1M_Lower = DeltaG_Transition_Dissociation_1M_Mean - DeltaG_Transition_Dissociation_1M_Err;
DeltaS_Transition_Association_1M_Lower = DeltaS_Transition_Association_1M_Mean - DeltaS_Transition_Association_1M_Err;
DeltaG_Transition_Association_1M_Lower = DeltaG_Transition_Association_1M_Mean - DeltaG_Transition_Association_1M_Err;

DeltaH_Transition_Dissociation_4400microM_Upper = DeltaH_Transition_Dissociation_4400microM_Mean + DeltaH_Transition_Dissociation_4400microM_Err;
DeltaS_Transition_Dissociation_4400microM_Upper = DeltaS_Transition_Dissociation_4400microM_Mean + DeltaS_Transition_Dissociation_4400microM_Err;   
DeltaG_Transition_Dissociation_4400microM_Upper = DeltaG_Transition_Dissociation_4400microM_Mean + DeltaG_Transition_Dissociation_4400microM_Err;
DeltaS_Transition_Association_4400microM_Upper = DeltaS_Transition_Association_4400microM_Mean + DeltaS_Transition_Association_4400microM_Err;
DeltaG_Transition_Association_4400microM_Upper = DeltaG_Transition_Association_4400microM_Mean + DeltaG_Transition_Association_4400microM_Err;
DeltaH_Transition_Dissociation_4400microM_Lower = DeltaH_Transition_Dissociation_4400microM_Mean - DeltaH_Transition_Dissociation_4400microM_Err;
DeltaS_Transition_Dissociation_4400microM_Lower = DeltaS_Transition_Dissociation_4400microM_Mean - DeltaS_Transition_Dissociation_4400microM_Err;   
DeltaG_Transition_Dissociation_4400microM_Lower = DeltaG_Transition_Dissociation_4400microM_Mean - DeltaG_Transition_Dissociation_4400microM_Err;
DeltaS_Transition_Association_4400microM_Lower = DeltaS_Transition_Association_4400microM_Mean - DeltaS_Transition_Association_4400microM_Err;
DeltaG_Transition_Association_4400microM_Lower = DeltaG_Transition_Association_4400microM_Mean - DeltaG_Transition_Association_4400microM_Err;

A = kb/h;
%dG = dH-TdS;
%dG = dGa-dGd;
%k = kbT/h * exp(-dGTT/RT)  dGTT =dHTT -T*dSTT
% dG/RT = dH/RT - dS/R
%TdS/RT = dS/R
R = 0.001987204259;%gas constant [Energy Kcal/mol K
T = str2sym(sprintf('T%d(t)',1)); % temperature in Celsius
S = sym('salt');          % salt concentration in moles per liter (M)
P= sym('primerConc');    % concentration of primers in mole per liter (M)

if (~isempty(MatchedLoc))
    dGseq1 = F_NearestNeighborTransitionWrapper(RV(seq1Matched));
    dGseq2 = F_NearestNeighborTransitionWrapper(RV(seq2Matched));
    
    dGf = 0.5*(dGseq1.dGf + dGseq2.dGf);   
    dGr = 0.5*(dGseq1.dGr + dGseq2.dGr);  
    dGeq = 0.5*(dGseq1.dGeq + dGseq2.dGeq);  
else
    dGf = Inf;
    dGr = Inf;
    dGeq = Inf;
end  
if (isempty(seq1))
   dGf = Inf;
   dGr = Inf;
   dGeq = Inf;
end

    Keq = exp(-dGeq/(R*(T+273.15)));
    kf_arr = exp(-dGf/(R*(T+273.15)));
    kr_arr = exp(-dGr/(R*(T+273.15)));
    kf_eyr = kb*(T+273.15)/h*exp(-dGf/(R*(T+273.15)));
    kr_eyr = kb*(T+273.15)/h*exp(-dGr/(R*(T+273.15)));
    %   dH(m)  dS(m)
    %   dG = dH(m)-T*dS(m)
    %   dG/RT = dH(m)/RT -dS(m)/R;
    %   dG/RT = 1/R*[dH(m)/T-dS(m)]
%     dGF  = NN(1:3,1) - ((T+273.15)  .* (NN(1:3,2)./1000)); % DeltaG
% dGR  = NN(4:6,1) - ((T+273.15)  .* (NN(4:6,2)./1000)); % DeltaG
% dGE  = NN(7:9,1) - ((T+273.15)  .* (NN(7:9,2)./1000)); % DeltaG
%     %dH(m)
mean(double(subs(kf_arr,T,25)))
%unimolecular kon calculation

    %kr = kf/Keq(T)  

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

