function [Sol] = F_TemperatureColocalizationStats(Pset,probes,gene_table,settings,DoesProbeBindSite,ExpressionMatrix,C_var)
% Given Probe Set Get Temperature To Use for optimal performance
GuessConc = 10^-20;
ConcRange = 3; DilutionFactor = 0.5;    
MaxProbeConc = 1;   
ProbeConcSeries = 10.^(DilutionFactor*[0:-1:-ConcRange+1])';
ProbeConc =  10^5;P_1FLAPs = [0 1];
T_low_celsius = 4;T_high_celsius = 75;
N_methods = 8;
R = 0.001987204259;%gas constant [Energy Kcal/mol K
kb = 1.380649*10^-23;%bolzman constant J/K
h = 6.62607015*10^-34;%planks constant J/Hz
Tref = 37+273.15;

LocMax = max(cell2mat(cellfun(@(x) x,{probes{:,3}},'UniformOutput',false)));  
Lpmin = min(cell2mat(cellfun(@length,{probes{:,2}},'UniformOutput',false)));
TargetLength = LocMax + Lpmin - 1;
theoryMaxProbes = floor(TargetLength/(Lpmin+settings.ProbeSpacing));
maxProbes = floor(TargetLength/(Lpmin+settings.ProbeSpacing));

GeneName = settings.GeneName;
GeneTarget = settings.transcript_IDs;
gene_table = sortrows(gene_table,[7 6],'ascend');
gene_table2 = gene_table(gene_table.Match>=15,:);
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
Names = unique(gene_table3.Name);
Names = convertCharsToStrings(Names);
uniNames = extractBefore(Names,'.');

Isoforms = find(contains(Names,GeneName));
Desired_Isoforms = find(contains(uniNames,extractBefore(GeneTarget{:},'.')));
UnDesired_Isoforms = setdiff(Isoforms,Desired_Isoforms);
offRNA_minusIso = setdiff(offRNA,UnDesired_Isoforms);

ON_IDs = Desired_Isoforms;
OFF_IDs = offRNA_minusIso;

Tset = find(sum(squeeze(sum(DoesProbeBindSite(Pset,:,:),1)),2)>0);
Smax = max(full(sum(sum(DoesProbeBindSite(Pset,Tset,:),3),1)));

Corr_ONDNA = arrayfun(@(x) find(Tset==offRNA(x)),1:length(onDNA),'Un',0);
Corr_ONRNA_IsoSpecific = arrayfun(@(x) find(Tset==onRNA(x)),1:length(onRNA),'Un',0);
Corr_ONRNA_IsoAgnostic = arrayfun(@(x) find(Tset==Isoforms(x)),1:length(Isoforms),'Un',0);
Corr_OFFDNA = arrayfun(@(x) find(Tset==offDNA(x)),1:length(offDNA),'Un',0);
Corr_OFFRNA_withIsos = arrayfun(@(x) find(Tset==offRNA(x)),1:length(offRNA),'Un',0);
Corr_OFFRNA_minusIso = arrayfun(@(x) find(Tset==offRNA_minusIso(x)),1:length(offRNA_minusIso),'Un',0);
Corr_ONDNA = [Corr_ONDNA{:}];
Corr_ONRNA_IsoSpecific = [Corr_ONRNA_IsoSpecific{:}];
Corr_ONRNA_IsoAgnostic = [Corr_ONRNA_IsoAgnostic{:}];
Corr_OFFDNA = [Corr_OFFDNA{:}];
Corr_OFFRNA_withIsos = [Corr_OFFRNA_withIsos{:}];
Corr_OFFRNA_minusIso = [Corr_OFFRNA_minusIso{:}];
Corr_OFFRNA_otherIso = setdiff(Corr_ONRNA_IsoAgnostic,Corr_ONRNA_IsoSpecific);

PxF_Count_G = @(G,PxF_ExprMatrix) sum(PxF_ExprMatrix(G,2:end),1);
PxF_Dist_G = @(G,PxF_ExprMatrix) sum(PxF_ExprMatrix(G,2:end),1)/sum(sum(PxF_ExprMatrix(G,2:end),1),2);
PxF_Mean_G = @(G,PxF_ExprMatrix) dot(1:Smax,sum(PxF_ExprMatrix(G,2:end),1))/sum(sum(PxF_ExprMatrix(G,2:end),1),2);
PxF_Var_G = @(G,PxF_ExprMatrix) dot((linspace(1,Smax,Smax)-PxF_Mean_G(G,PxF_ExprMatrix)).^2,sum(PxF_ExprMatrix(G,2:end),1))/sum(sum(PxF_ExprMatrix(G,2:end),1),2);
PxF_Fano_G = @(G,PxF_ExprMatrix) PxF_Var_G(G,PxF_ExprMatrix)/PxF_Mean_G(G,PxF_ExprMatrix);
PxF_Moment_G = @(G,PxF_ExprMatrix,n) dot((linspace(1,Smax,Smax)-PxF_Mean_G(G,PxF_ExprMatrix)).^n,sum(PxF_ExprMatrix(G,2:end),1))/(sum(sum(PxF_ExprMatrix(G,2:end),1),2)*PxF_Var_G(G,PxF_ExprMatrix)^(n/2));   

Js = @(x) find(sum(squeeze(sum(DoesProbeBindSite(x,:,:),1)),2)>0);
Js_Sites = @(x) find(sum(sum(DoesProbeBindSite(x,Js(x),:),1),2)>0);
Js_RNA = @(x) NonDNA_IDs(ismember(NonDNA_IDs,Js(x)));
Js_DNA = @(x) DNA_IDs(ismember(DNA_IDs,Js(x)));
AddDim3 = @(x,C) permute(repmat(C,[1 1 1 x]),[1 4 2 3]);
AddDim2 = @(x,C) permute(repmat(C,[1 1 1 x]),[4 1 2 3]);
AddDim1 = @(x,C) permute(repmat(C,[1 1 x]),[3 1 2]);

dCp_mod = C_var{1}; 
dH0eq_mod = C_var{2}; 
dS0eq_mod = C_var{3}; 
dH0f_mod = C_var{4}; 
dS0f_mod = C_var{5}; 
dH0r_mod = C_var{6}; 
dS0r_mod = C_var{7}; 
dH0eq_Complement = C_var{8}; 
dCp_Complement = C_var{9}; 
dS0eq_Complement = C_var{10}; 
dH0f_Complement = C_var{11}; 
dS0f_Complement = C_var{12}; 
dH0r_Complement = C_var{13}; 
dS0r_Complement = C_var{14}; 

[SelfSeqParsed,CrossDimerSeqParsed] = A_JH_GenerateSecondaryStructureSeq(probes,Pset,settings);
CrossDimerSeqParsed = cellfun(@(x) x(sum(cellfun(@isempty,x),2)'==0,1:2), CrossDimerSeqParsed,'Un',0);
Ns = cellfun(@(x) length(x), SelfSeqParsed);
Nc = cellfun(@(x) size(x,1), CrossDimerSeqParsed);
c_alpha = ones(size(probes,1),size(probes,1))+diag(ones(size(probes,1),1));
Ns_val = zeros(size(probes,1),max(Ns(:)));
Nc_val = zeros(size(probes,1),size(probes,1),max(Nc(:)));
for i=1:length(Pset)
    if (Ns(i)>0)
        Ns_val(Pset(i),1:Ns(i)) = 1; 
    end
    for j = 1:length(Pset)
        if (j<=i)
            if (Nc(i,j)>0)
                Nc_val(Pset(i),Pset(j),1:Nc(i,j)) = 1; 
            end
        end
    end
end

[~,dH0s_eq,dS0s_eq,...
 dH0s_f,dS0s_f,dH0s_r,dS0s_r,dCps,...
 ~,dH0d_eq,dS0d_eq,...
 dH0d_f,dS0d_f,dH0d_r,dS0d_r,dCpd,~,~] = A_JH_GenerateSecondaryStructureInfo_V2(probes,Pset,settings);

dHs_eq = @(T,Tref,p,n,m) full(dH0s_eq(p,n,m)) + full(dCps(p,n,m))*(T-Tref);
dSs_eq = @(T,Tref,p,n,m) full(dS0s_eq(p,n,m)) + full(dCps(p,n,m))*log(T/Tref);
dGs_eq = @(T,Tref,p,n,m) dHs_eq(T,Tref,p,n,m) - T*dSs_eq(T,Tref,p,n,m);  
Ks_eq = @(T,Tref,p,n,m) Ns_val(p,n)*exp(-dGs_eq(T,Tref,p,n,m)/(R*T));
dHd_eq = @(T,Tref,p,q,n,m) full(dH0d_eq(p,q,n,m)) + full(dCpd(p,q,n,m))*(T-Tref);
dSd_eq = @(T,Tref,p,q,n,m) full(dS0d_eq(p,q,n,m)) + full(dCpd(p,q,n,m))*log(T/Tref);
dGd_eq = @(T,Tref,p,q,n,m) dHd_eq(T,Tref,p,q,n,m) - T*dSd_eq(T,Tref,p,q,n,m); 
Kd_eq = @(T,Tref,p,q,n,m) Nc_val(p,q,n)*repmat(c_alpha(p,q),[1 1 length(n)]).*exp(-dGd_eq(T,Tref,p,q,n,m)/(R*T));
dHeq_mod = @(T,Tref,p,q,r,m) full(dH0eq_mod(p,q,r,m)) + full(dCp_mod(p,q,r,m))*(T-Tref);
dSeq_mod = @(T,Tref,p,q,r,m) full(dS0eq_mod(p,q,r,m)) + full(dCp_mod(p,q,r,m))*log(T/Tref);
dGeq_mod = @(T,Tref,p,q,r,m) dHeq_mod(T,Tref,p,q,r,m) - T*dSeq_mod(T,Tref,p,q,r,m); 
Keq_mod = @(T,Tref,p,q,r,m) full(DoesProbeBindSite(p,q,r))*exp(-dGeq_mod(T,Tref,p,q,r,m)/(R*T));
dHeq_Complement = @(T,Tref,p,r,m) full(dH0eq_Complement(p,r,m)) + full(dCp_Complement(p,r,m))*(T-Tref);
dSeq_Complement = @(T,Tref,p,r,m) full(dS0eq_Complement(p,r,m)) + full(dCp_Complement(p,r,m))*log(T/Tref);
dGeq_Complement = @(T,Tref,p,r,m) dHeq_Complement(T,Tref,p,r,m) - T*dSeq_Complement(T,Tref,p,r,m); 
Keq_Complement = @(T,Tref,p,r,m) double(full(sum(DoesProbeBindSite(p,q,r),1))>0)*exp(-dGeq_Complement(T,Tref,p,r,m)/(R*T));

%Two-Color Tests
PColor2_Split = [1*ones(1,floor(length(Pset)/2)) 2*ones(1,ceil(length(Pset)/2))];
if (mod(length(Pset),2)==0)
    PColor2_Weave = repmat([1 2],[1 floor(length(Pset)/2)]);
else
    PColor2_Weave = [repmat([1 2],[1 floor(length(Pset)/2)]) 1];
end
N_Channels = length(unique(PColor2_Split));
Ix_Split = cell(1,N_Channels);
Ix_Weave = cell(1,N_Channels);
vsi_Split = cell(1,N_Channels);vsi_Weave = cell(1,N_Channels);
for k = 1:N_Channels
    Ix_Split{k} = Pset(PColor2_Split==k);    
    Ix_Weave{k} = Pset(PColor2_Weave==k);   
    vsi_Split{k} = @(Tx) find(sum(squeeze(DoesProbeBindSite(Ix_Split{k},Tset(Tx),:)),1)>0);
    vsi_Weave{k} = @(Tx) find(sum(squeeze(DoesProbeBindSite(Ix_Weave{k},Tset(Tx),:)),1)>0);
end
vsi_Split_Uni{1} = @(Tx) setdiff(vsi_Split{1}(Tx),vsi_Split{2}(Tx));
vsi_Weave_Uni{1} = @(Tx) setdiff(vsi_Weave{1}(Tx),vsi_Weave{2}(Tx));
vsi_Split_Uni{2} = @(Tx) setdiff(vsi_Split{2}(Tx),vsi_Split{1}(Tx));
vsi_Weave_Uni{2} = @(Tx) setdiff(vsi_Weave{2}(Tx),vsi_Weave{1}(Tx));
overlap_Split = @(Tx) intersect(vsi_Split{1}(Tx),vsi_Split{2}(Tx));
overlap_Weave = @(Tx) intersect(vsi_Weave{1}(Tx),vsi_Weave{2}(Tx)); 
PxFD_Split = cell(1,N_Channels);PxFD_Weave = cell(1,N_Channels);



%% Time Varying Base Case 
if (solution_case ==1)
    CProbes_Free = arrayfun(@(N) str2sym(sprintf('CProbes_Free%d(t)',N)), 1:length(Pset)).';   
    ProbeConc = sym('ProbeConc',[length(Pset) 1]);
    T = str2sym('T');
    M = str2sym('Model');
    C = str2sym('CellType');
    CTargets_Total_TargetSite = repmat(ExpressionMatrix,[1 size(DoesProbeBindSite,3)]);
    h_sub_RNA_base = @(T,Tref,x,pf,m) 1+squeeze(multiprod(permute(squeeze(Keq_mod(T,Tref,x,Js_RNA(x),Js_Sites(x),m)).*DoesProbeBindSite(x,Js_RNA(x),Js_Sites(x)),[2 1 3]),pf));
    h_sub_DNA_base = @(T,Tref,x,pf,m) 1+squeeze(multiprod(permute(squeeze(Keq_mod(T,Tref,x,Js_DNA(x),Js_Sites(x),m)).*DoesProbeBindSite(x,Js_DNA(x),Js_Sites(x)),[2 1 3]),pf));
    h_sub_RNA = @(T,Tref,x,pf,m,cell) CTargets_Total_TargetSite(Js_RNA(x),Js_Sites(x),cell)./h_sub_RNA_base(T,Tref,x,pf,m);
    h_sub_DNA = @(T,Tref,x,pf,m,cell) (sqrt(1+4*CTargets_Total_TargetSite(Js_DNA(x),Js_Sites(x),cell)./h_sub_DNA_base(T,Tref,x,pf,m))-1)./(2*Keq_Complement(T,Tref,Js_DNA(x),Js_Sites(x),m));    
    h_sub_DNA_complement = @(T,Tref,x,pf,m,cell) CTargets_Total_TargetSite(Js_DNA(x),Js_Sites(x),cell)./(1+h_sub_DNA(T,Tref,x,pf,m,cell)./Keq_Complement(T,Tref,Js_DNA(x),Js_Sites(x),m));   
    h_RNA = @(T,Tref,x,k,C,m) full(sum(sum(squeeze(Keq_mod(T,Tref,x(k),Js_RNA(x),Js_Sites(x),m)).*C,2,'omitnan'),1,'omitnan'));
    h_DNA = @(T,Tref,x,k,C,m) full(sum(sum(squeeze(Keq_mod(T,Tref,x(k),Js_DNA(x),Js_Sites(x),m)).*C,2,'omitnan'),1,'omitnan')); 
    h_sub_RNA1_Func = @(T,CProbes_Free,Pset) removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(cellType) sum(CATnWrapper(arrayfun(@(m) h_sub_RNA(T,Tref,Pset,CProbes_Free,m,cellType)*kroneckerDelta(M,m)*kroneckerDelta(C,cellType),1:N_methods,'Un',0),3),3), 1:size(ExpressionMatrix,2),'Un',0),3),3));
    h_sub_DNA1_Func = @(T,CProbes_Free,Pset) removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(cellType) sum(CATnWrapper(arrayfun(@(m) h_sub_DNA(T,Tref,Pset,CProbes_Free,m,cellType)*kroneckerDelta(M,m)*kroneckerDelta(C,cellType),1:N_methods,'Un',0),3),3), 1:size(ExpressionMatrix,2),'Un',0),3),3));
    h_sub_DNA1_complement_Func = @(T,CProbes_Free,Pset) removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(cellType) sum(CATnWrapper(arrayfun(@(m) h_sub_DNA_complement(T,Tref,Pset,CProbes_Free,m,cellType)*kroneckerDelta(M,m)*kroneckerDelta(C,cellType),1:N_methods,'Un',0),3),3), 1:size(ExpressionMatrix,2),'Un',0),3),3));
    h_RNAx_Func = @(T,CProbes_Free,h_sub_RNA1,Pset) CATnWrapper(arrayfun(@(J) sum(CATnWrapper(arrayfun(@(m) h_RNA(T,Tref,Pset,J,h_sub_RNA1,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),3),3),1:length(Pset),'Un',0),1);
    h_DNAx_Func = @(T,CProbes_Free,h_sub_DNA1,Pset) CATnWrapper(arrayfun(@(J) sum(CATnWrapper(arrayfun(@(m) h_DNA(T,Tref,Pset,J,h_sub_DNA1,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),3),3),1:length(Pset),'Un',0),1);
    K_S_Func = @(T,P) sum(sum(CATnWrapper(arrayfun(@(m) Ks_eq(T,Tref,P,1:max(Ns),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),3),3),2);
    K_CD_Func = @(T,P) sum(sum(CATnWrapper(arrayfun(@(m) Kd_eq(T,Tref,P,P,1:max(Nc(:)),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4),3);
    K_S = K_S_Func(T,Pset);
    K_CD = K_CD_Func(T,Pset);
    h_sub_RNA1 = h_sub_RNA1_Func(T,CProbes_Free,Pset);
    h_sub_DNA1 = h_sub_DNA1_Func(T,CProbes_Free,Pset);
    h_sub_DNA1_complement = h_sub_DNA1_complement_Func(T,CProbes_Free,Pset);
    h_RNAx = h_RNAx_Func(T,CProbes_Free,h_sub_RNA1,Pset);
    h_DNAx = h_DNAx_Func(T,CProbes_Free,h_sub_DNA1,Pset);
    ODE_basic = diff(CProbes_Free) == ProbeConc./(1+K_S+K_CD*CProbes_Free+h_RNAx + h_DNAx)-CProbes_Free;
    [M0_basic,F0_basic] = massMatrixForm(ODE_basic,CProbes_Free);
    M_basic = odeFunction(M0_basic,CProbes_Free,'Sparse',true);
    F_basic = odeFunction(F0_basic,CProbes_Free,ProbeConc,T,M,C);
pRNA_TargetSites_Bound = @(T,Pset,CProbes_Free,m) (squeeze(multiprod(permute(Keq_mod(T,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),m).*DoesProbeBindSite(Pset,Js_RNA(Pset),Js_Sites(Pset)),[2 1 3]),CProbes_Free)).*h_sub_RNA1)./...
    (squeeze(multiprod(permute(Keq_mod(T,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),m).*DoesProbeBindSite(Pset,Js_RNA(Pset),Js_Sites(Pset)),[2 1 3]),CProbes_Free)).*h_sub_RNA1+h_sub_RNA1);
pDNA_TargetSites_Bound = @(T,Pset,CProbes_Free,m) (squeeze(multiprod(permute(Keq_mod(T,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),m).*DoesProbeBindSite(Pset,Js_DNA(Pset),Js_Sites(Pset)),[2 1 3]),CProbes_Free)).*h_sub_DNA1)./...
 (squeeze(multiprod(permute(Keq_mod(T,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),m).*DoesProbeBindSite(Pset,Js_DNA(Pset),Js_Sites(Pset)),[2 1 3]),CProbes_Free)).*h_sub_DNA1+h_sub_DNA1+Keq_Complement(T,Tref,Js_DNA(Pset),Js_Sites(Pset),m).*h_sub_DNA1.*h_sub_DNA1_complement);
p_TargetSites_Bound(Js_RNA(Pset),Js_Sites(Pset)) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pRNA_TargetSites_Bound(T,Pset,CProbes_Free,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),3),3));
p_TargetSites_Bound(Js_DNA(Pset),Js_Sites(Pset)) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pDNA_TargetSites_Bound(T,Pset,CProbes_Free,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),3),3));
p_TargetSites_Unbound = 1 - p_TargetSites_Bound;   
p_TargetSites_Bound_Func = matlabFunction(p_TargetSites_Bound,'vars',{CProbes_Free,T,M,C});
p_TargetSites_Unbound_Func = matlabFunction(p_TargetSites_Unbound,'vars',{CProbes_Free,T,M,C});
p_TargetSites_Bound_Weave = cell(1,N_Channels);
p_TargetSites_Unbound_Weave = cell(1,N_Channels);
p_TargetSites_Bound_Weave_Func = cell(1,N_Channels);
p_TargetSites_Unbound_Weave_Func = cell(1,N_Channels);
p_TargetSites_Bound_Split = cell(1,N_Channels);
p_TargetSites_Unbound_Split = cell(1,N_Channels);
p_TargetSites_Bound_Split_Func = cell(1,N_Channels);
p_TargetSites_Unbound_Split_Func = cell(1,N_Channels);
for k = 1:N_Channels
    p_TargetSites_Bound_Weave{k}(Js_RNA(Ix_Weave{k}),Js_Sites(Ix_Weave{k})) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pRNA_TargetSites_Bound(T,Ix_Weave{k},CProbes_Free(ismember(Pset,Ix_Weave{k})),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),3),3));
    p_TargetSites_Bound_Weave{k}(Js_DNA(Ix_Weave{k}),Js_Sites(Ix_Weave{k})) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pDNA_TargetSites_Bound(T,Ix_Weave{k},CProbes_Free(ismember(Pset,Ix_Weave{k})),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),3),3));
    p_TargetSites_Unbound_Weave{k} = 1 - p_TargetSites_Bound_Weave{k};   
    p_TargetSites_Bound_Weave_Func{k} = matlabFunction(p_TargetSites_Bound_Weave{k},'vars',{CProbes_Free,T,M,C});
    p_TargetSites_Unbound_Weave_Func{k} = matlabFunction(p_TargetSites_Unbound_Weave{k},'vars',{CProbes_Free,T,M,C});
    p_TargetSites_Bound_Split{k}(Js_RNA(Ix_Split{k}),Js_Sites(Ix_Split{k})) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pRNA_TargetSites_Bound(T,Ix_Split{k},CProbes_Free(ismember(Pset,Ix_Split{k})),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),3),3));
    p_TargetSites_Bound_Split{k}(Js_DNA(Ix_Split{k}),Js_Sites(Ix_Split{k})) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pDNA_TargetSites_Bound(T,Ix_Split{k},CProbes_Free(ismember(Pset,Ix_Split{k})),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),3),3));
    p_TargetSites_Unbound_Split{k} = 1 - p_TargetSites_Bound_Split{k};   
    p_TargetSites_Bound_Split_Func{k} = matlabFunction(p_TargetSites_Bound_Split{k},'vars',{CProbes_Free,T,M,C});
    p_TargetSites_Unbound_Split_Func{k} = matlabFunction(p_TargetSites_Unbound_Split{k},'vars',{CProbes_Free,T,M,C});
end
Mean_G0_2D = @(G,c,mean_probes_Target) sum(ExpressionMatrix(G,c).*mean_probes_Target(G))/sum(ExpressionMatrix(G,c));
 Var_G0_2D = @(G,c,var_probes_Target) sum(ExpressionMatrix(G,c).*var_probes_Target(G))/sum(ExpressionMatrix(G,c));
 Mean_G_2D = @(G,c,mean_probes_Target,p_TargetSites_Zero) sum(ExpressionMatrix(G,c).*mean_probes_Target(G))/sum(ExpressionMatrix(G,c).*(1-p_TargetSites_Zero(G)));
Var_G_2D = @(G,c,var_probes_Target,p_TargetSites_Zero) sum(ExpressionMatrix(G,c).*var_probes_Target(G))/sum(ExpressionMatrix(G,c).*(1-p_TargetSites_Zero(G)));%edge cases
sConc_G_2D = @(G,c,p_TargetSites_Zero) sum(ExpressionMatrix(G,c).*(1-p_TargetSites_Zero(G))); %count spots ~
pConc_G_2D = @(G,c,mean_probes_Target) sum(ExpressionMatrix(G,c).*mean_probes_Target(G)); %count probes  ~
sFrac_G_2D = @(G,c,p_TargetSites_Zero) 100*sConc_G_2D(G,c,p_TargetSites_Zero)/sConc_G_2D(1:length(Tset),c,p_TargetSites_Zero);
pFrac_G_2D = @(G,c,mean_probes_Target) 100*pConc_G_2D(G,c,mean_probes_Target)/pConc_G_2D(1:length(Tset),c,mean_probes_Target);
elseif (solution_case == 2)%% Simul MultiTissue
    CProbes_Free_MultiTissue = CATnWrapper(arrayfun(@(M) arrayfun(@(N) str2sym(sprintf('CProbes_Free%d_%d(t)',M,N)), 1:length(Pset)),1:size(ExpressionMatrix,2),'Un',0),1);
    ProbeConc_MultiTissue = sym('ProbeConc',[size(ExpressionMatrix,2) length(Pset)]);
    T = str2sym('T');
    M = str2sym('Model');
    EKernel_MultiTissue = permute(ExpressionMatrix,[2 1]);
    CTargets_Total_TargetSite_MultiTissue = repmat(permute(ExpressionMatrix,[2 1]),[1 1 size(DoesProbeBindSite,3)]);  
    h_sub_RNA_base_MultiTissue = @(T,Tref,x,pf,m) 1+permute(multiprod(permute(Keq_mod(T,Tref,x,Js_RNA(x),Js_Sites(x),m).*DoesProbeBindSite(x,Js_RNA(x),Js_Sites(x)),[2 1 3]),pf.'),[2 1 3]);
    h_sub_DNA_base_MultiTissue = @(T,Tref,x,pf,m) 1+permute(multiprod(permute(Keq_mod(T,Tref,x,Js_DNA(x),Js_Sites(x),m).*DoesProbeBindSite(x,Js_DNA(x),Js_Sites(x)),[2 1 3]),pf.'),[2 1 3]); 
    h_sub_RNA_MultiTissue = @(T,Tref,x,pf,m) CTargets_Total_TargetSite_MultiTissue(1:size(EKernel_MultiTissue,1),Js_RNA(x),Js_Sites(x))./h_sub_RNA_base_MultiTissue(T,Tref,x,pf,m);
    h_sub_DNA_MultiTissue = @(T,Tref,x,pf,Kernel,m) (sqrt(1+4*CTargets_Total_TargetSite_MultiTissue(1:size(EKernel_MultiTissue,1),Js_DNA(x),Js_Sites(x))./h_sub_DNA_base_MultiTissue(T,Tref,x,pf,m))-1)./(2*AddDim1(size(Kernel,1),Keq_Complement(T,Tref,Js_DNA(x),Js_Sites(x),m)));
    h_sub_DNA_complement_MultiTissue = @(T,Tref,x,pf,Kernel,m) CTargets_Total_TargetSite_MultiTissue(1:size(EKernel_MultiTissue,1),Js_DNA(x),Js_Sites(x))./(1+h_sub_DNA_MultiTissue(T,Tref,x,pf,Kernel,m)./AddDim1(size(Kernel,1),Keq_Complement(T,Tref,Js_DNA(x),Js_Sites(x),m)));   
    h_RNA_MultiTissue = @(T,Tref,x,k,C,Kernel,m) full(squeeze(sum(squeeze(sum(AddDim2(size(Kernel,1),squeeze(Keq_mod(T,Tref,x(k),Js_RNA(x),Js_Sites(x),m))).*C,3,'omitnan')),2,'omitnan')));
    h_DNA_MultiTissue = @(T,Tref,x,k,C,Kernel,m) full(squeeze(sum(squeeze(sum(AddDim2(size(Kernel,1),squeeze(Keq_mod(T,Tref,x(k),Js_DNA(x),Js_Sites(x),m))).*C,3,'omitnan')),2,'omitnan')));      
    h_sub_RNA1_MultiTissue_Func = @(T,CProbes_Free,Pset) removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) h_sub_RNA_MultiTissue(T,Tref,Pset,CProbes_Free,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
    h_sub_DNA1_MultiTissue_Func = @(T,CProbes_Free,Pset,Kernel) removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) h_sub_DNA_MultiTissue(T,Tref,Pset,CProbes_Free,Kernel,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
    h_sub_DNA1_complement_MultiTissue_Func = @(T,CProbes_Free,Pset,Kernel) removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) h_sub_DNA_complement_MultiTissue(T,Tref,Pset,CProbes_Free,Kernel,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
    h_RNAx_MultiTissue_Func = @(T,CProbes_Free,h_sub_RNA1,Pset,Kernel) CATnWrapper(arrayfun(@(J) sum(CATnWrapper(arrayfun(@(m) h_RNA_MultiTissue(T,Tref,Pset,J,h_sub_RNA1,Kernel,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4),1:length(Pset),'Un',0),2);
    h_DNAx_MultiTissue_Func = @(T,CProbes_Free,h_sub_DNA1,Pset,Kernel) CATnWrapper(arrayfun(@(J) sum(CATnWrapper(arrayfun(@(m) h_DNA_MultiTissue(T,Tref,Pset,J,h_sub_DNA1,Kernel,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4),1:length(Pset),'Un',0),2);
    K_S_Func = @(T,P) sum(sum(CATnWrapper(arrayfun(@(m) Ks_eq(T,Tref,P,1:max(Ns),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),3),3),2);
    K_CD_Func = @(T,P) sum(sum(CATnWrapper(arrayfun(@(m) Kd_eq(T,Tref,P,P,1:max(Nc(:)),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4),3);    
    K_S = K_S_Func(T,Pset);
    K_CD = K_CD_Func(T,Pset);
    h_sub_RNA1_MultiTissue = h_sub_RNA1_MultiTissue_Func(T,CProbes_Free_MultiTissue,Pset);
    h_sub_DNA1_MultiTissue = h_sub_DNA1_MultiTissue_Func(T,CProbes_Free_MultiTissue,Pset,EKernel_MultiTissue);
    h_sub_DNA1_complement_MultiTissue = h_sub_DNA1_complement_MultiTissue_Func(T,CProbes_Free_MultiTissue,Pset,EKernel_MultiTissue);
    h_RNAx_MultiTissue = h_RNAx_MultiTissue_Func(T,CProbes_Free_MultiTissue,h_sub_RNA1_MultiTissue,Pset,EKernel_MultiTissue);
    h_DNAx_MultiTissue = h_DNAx_MultiTissue_Func(T,CProbes_Free_MultiTissue,h_sub_DNA1_MultiTissue,Pset,EKernel_MultiTissue);
    ODE_MultiTissue = diff(CProbes_Free_MultiTissue) == ProbeConc_MultiTissue./(1+AddDim1(size(EKernel_MultiTissue,1),K_S)+multiprod(K_CD,CProbes_Free_MultiTissue')'+ h_RNAx_MultiTissue + h_DNAx_MultiTissue)-CProbes_Free_MultiTissue;  
    [M0_MultiTissue,F0_MultiTissue] = massMatrixForm(ODE_MultiTissue,CProbes_Free_MultiTissue);
    M_MultiTissue = odeFunction(M0_MultiTissue,CProbes_Free_MultiTissue,'Sparse',true);
    F_MultiTissue = odeFunction(F0_MultiTissue,CProbes_Free_MultiTissue,ProbeConc_MultiTissue,T,M);    
pRNA_TargetSites_Bound_MultiTissue = @(T,Pset,CProbes_Free,m) ... 
    (permute(multiprod(permute(Keq_mod(T,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),m).*DoesProbeBindSite(Pset,Js_RNA(Pset),Js_Sites(Pset)),[2 1 3]),CProbes_Free'),[2 1 3]).*h_sub_RNA1_MultiTissue)./(...
     permute(multiprod(permute(Keq_mod(T,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),m).*DoesProbeBindSite(Pset,Js_RNA(Pset),Js_Sites(Pset)),[2 1 3]),CProbes_Free'),[2 1 3]).*h_sub_RNA1_MultiTissue + h_sub_RNA1_MultiTissue);
pDNA_TargetSites_Bound_MultiTissue = @(T,Pset,CProbes_Free,m) ...
    (permute(multiprod(permute(Keq_mod(T,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),m).*DoesProbeBindSite(Pset,Js_DNA(Pset),Js_Sites(Pset)),[2 1 3]),CProbes_Free'),[2 1 3]).*h_sub_DNA1_MultiTissue)./(...
     permute(multiprod(permute(Keq_mod(T,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),m).*DoesProbeBindSite(Pset,Js_DNA(Pset),Js_Sites(Pset)),[2 1 3]),CProbes_Free'),[2 1 3]).*h_sub_DNA1_MultiTissue + ...
    AddDim1(size(EKernel_MultiTissue,1),Keq_Complement(T,Tref,Js_DNA(Pset),Js_Sites(Pset),m)).*h_sub_DNA1_MultiTissue.*h_sub_DNA1_complement_MultiTissue + h_sub_DNA1_MultiTissue);
p_TargetSites_Bound_MultiTissue(1:size(EKernel_MultiTissue,1),Js_RNA(Pset),Js_Sites(Pset)) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pRNA_TargetSites_Bound_MultiTissue(T,Pset,CProbes_Free_MultiTissue,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
p_TargetSites_Bound_MultiTissue(1:size(EKernel_MultiTissue,1),Js_DNA(Pset),Js_Sites(Pset)) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pDNA_TargetSites_Bound_MultiTissue(T,Pset,CProbes_Free_MultiTissue,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
p_TargetSites_Unbound_MultiTissue = 1 - p_TargetSites_Bound_MultiTissue;
p_TargetSites_Bound_MultiTissue_Func = matlabFunction(p_TargetSites_Bound_MultiTissue,'vars',{CProbes_Free_MultiTissue,T,M});
p_TargetSites_Unbound_MultiTissue_Func = matlabFunction(p_TargetSites_Unbound_MultiTissue,'vars',{CProbes_Free_MultiTissue,T,M});
p_TargetSites_Bound_MultiTissue_Weave = cell(1,N_Channels);
p_TargetSites_Unbound_MultiTissue_Weave = cell(1,N_Channels);
p_TargetSites_Bound_MultiTissue_Weave_Func = cell(1,N_Channels);
p_TargetSites_Unbound_MultiTissue_Weave_Func = cell(1,N_Channels);
p_TargetSites_Bound_MultiTissue_Split = cell(1,N_Channels);
p_TargetSites_Unbound_MultiTissue_Split = cell(1,N_Channels);
p_TargetSites_Bound_MultiTissue_Split_Func = cell(1,N_Channels);
p_TargetSites_Unbound_MultiTissue_Split_Func = cell(1,N_Channels);
for k = 1:N_Channels
    p_TargetSites_Bound_MultiTissue_Weave{k}(1:size(EKernel_MultiTissue,1),Js_RNA(Ix_Weave{k}),Js_Sites(Ix_Weave{k})) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pRNA_TargetSites_Bound_MultiTissue(T,Ix_Weave{k},CProbes_Free_MultiTissue(:,ismember(Pset,Ix_Weave{k})),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
    p_TargetSites_Bound_MultiTissue_Weave{k}(1:size(EKernel_MultiTissue,1),Js_DNA(Ix_Weave{k}),Js_Sites(Ix_Weave{k})) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pDNA_TargetSites_Bound_MultiTissue(T,Ix_Weave{k},CProbes_Free_MultiTissue(:,ismember(Pset,Ix_Weave{k})),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
    p_TargetSites_Unbound_MultiTissue_Weave{k} = 1 - p_TargetSites_Bound_MultiTissue_Weave{k};
    p_TargetSites_Bound_MultiTissue_Weave_Func{k} = matlabFunction(p_TargetSites_Bound_MultiTissue_Weave{k},'vars',{CProbes_Free_MultiTissue,T,M});
    p_TargetSites_Unbound_MultiTissue_Weave_Func{k} = matlabFunction(p_TargetSites_Unbound_MultiTissue_Weave{k},'vars',{CProbes_Free_MultiTissue,T,M});
    p_TargetSites_Bound_MultiTissue_Split{k}(1:size(EKernel_MultiTissue,1),Js_RNA(Ix_Split{k}),Js_Sites(Ix_Split{k})) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pRNA_TargetSites_Bound_MultiTissue(T,Ix_Split{k},CProbes_Free_MultiTissue(:,ismember(Pset,Ix_Split{k})),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
    p_TargetSites_Bound_MultiTissue_Split{k}(1:size(EKernel_MultiTissue,1),Js_DNA(Ix_Split{k}),Js_Sites(Ix_Split{k})) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pDNA_TargetSites_Bound_MultiTissue(T,Ix_Split{k},CProbes_Free_MultiTissue(:,ismember(Pset,Ix_Split{k})),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
    p_TargetSites_Unbound_MultiTissue_Split{k} = 1 - p_TargetSites_Bound_MultiTissue_Split{k};
    p_TargetSites_Bound_MultiTissue_Split_Func{k} = matlabFunction(p_TargetSites_Bound_MultiTissue_Split{k},'vars',{CProbes_Free_MultiTissue,T,M});
    p_TargetSites_Unbound_MultiTissue_Split_Func{k} = matlabFunction(p_TargetSites_Unbound_MultiTissue_Split{k},'vars',{CProbes_Free_MultiTissue,T,M});
end
Mean_G0_3D_MultiTissue = @(G,mean_probes_Target) sum(EKernel_MultiTissue(:,G).*mean_probes_Target(:,G),2)./sum(EKernel_MultiTissue(:,G),2);
 Var_G0_3D_MultiTissue = @(G,var_probes_Target) sum(EKernel_MultiTissue(:,G).*var_probes_Target(:,G),2)./sum(EKernel_MultiTissue(:,G),2);
 Mean_G_3D_MultiTissue = @(G,mean_probes_Target,p_TargetSites_Zero) sum(EKernel_MultiTissue(:,G).*mean_probes_Target(:,G),2)./sum(EKernel_MultiTissue(:,G).*(1-p_TargetSites_Zero(:,G)),2);
Var_G_3D_MultiTissue = @(G,var_probes_Target,p_TargetSites_Zero) sum(EKernel_MultiTissue(:,G).*var_probes_Target(:,G),2)./sum(EKernel_MultiTissue(:,G).*(1-p_TargetSites_Zero(:,G)),2);%edge cases
sConc_G_3D_MultiTissue = @(G,p_TargetSites_Zero) sum(EKernel_MultiTissue(:,G).*(1-p_TargetSites_Zero(:,G)),2); %count spots ~
pConc_G_3D_MultiTissue = @(G,mean_probes_Target) sum(EKernel_MultiTissue(:,G).*mean_probes_Target(:,G),2); %count probes  ~
sFrac_G_3D_MultiTissue = @(G,p_TargetSites_Zero) 100*sConc_G_3D_MultiTissue(G,p_TargetSites_Zero)./sConc_G_3D_MultiTissue(1:length(Tset),p_TargetSites_Zero);
pFrac_G_3D_MultiTissue = @(G,mean_probes_Target) 100*pConc_G_3D_MultiTissue(G,mean_probes_Target)./pConc_G_3D_MultiTissue(1:length(Tset),mean_probes_Target);
elseif (solution_case ==3) %% Simul Multi Probe Conc in Multi Tissue
    CProbes_Free_MultiProbeConc = CATnWrapper(arrayfun(@(M) arrayfun(@(N) str2sym(sprintf('CProbes_Free%d_%d(t)',M,N)), 1:length(Pset)),1:size(ExpressionMatrix,2)*ConcRange,'Un',0),1);
    ProbeConcMatrix_MultiProbeConc = sym('ProbeConc',[size(ExpressionMatrix,2)*ConcRange length(Pset)]);
    T = str2sym('T');
    M = str2sym('Model');
    EKernel_MultiProbeConc = permute(repmat(ExpressionMatrix,[1 ConcRange]),[2 1]); 
    CTargets_Total_TargetSite_MultiProbeConc = permute(repmat(ExpressionMatrix,[1 ConcRange size(DoesProbeBindSite,3)]),[2 1 3]);
    h_sub_RNA_base_MultiProbeConc = @(T,Tref,x,pf,m) 1+permute(multiprod(permute(Keq_mod(T,Tref,x,Js_RNA(x),Js_Sites(x),m).*DoesProbeBindSite(x,Js_RNA(x),Js_Sites(x)),[2 1 3]),pf.'),[2 1 3]);
    h_sub_DNA_base_MultiProbeConc = @(T,Tref,x,pf,m) 1+permute(multiprod(permute(Keq_mod(T,Tref,x,Js_DNA(x),Js_Sites(x),m).*DoesProbeBindSite(x,Js_DNA(x),Js_Sites(x)),[2 1 3]),pf.'),[2 1 3]);
    h_sub_RNA_MultiProbeConc = @(T,Tref,x,pf,m) CTargets_Total_TargetSite_MultiProbeConc(1:size(EKernel_MultiProbeConc,1),Js_RNA(x),Js_Sites(x))./h_sub_RNA_base_MultiProbeConc(T,Tref,x,pf,m);
    h_sub_DNA_MultiProbeConc = @(T,Tref,x,pf,Kernel,m) (sqrt(1+4*CTargets_Total_TargetSite_MultiProbeConc(1:size(EKernel_MultiProbeConc,1),Js_DNA(x),Js_Sites(x))./h_sub_DNA_base_MultiProbeConc(T,Tref,x,pf,m))-1)./(2*AddDim1(size(Kernel,1),Keq_Complement(T,Tref,Js_DNA(x),Js_Sites(x),m)));
    h_sub_DNA_complement_MultiProbeConc = @(T,Tref,x,pf,Kernel,m) CTargets_Total_TargetSite_MultiProbeConc(1:size(EKernel_MultiProbeConc,1),Js_DNA(x),Js_Sites(x))./(1+h_sub_DNA_MultiProbeConc(T,Tref,x,pf,Kernel,m)./AddDim1(size(Kernel,1),Keq_Complement(T,Tref,Js_DNA(x),Js_Sites(x),m)));        
    h_RNA_MultiProbeConc = @(T,Tref,x,k,C,Kernel,m) full(squeeze(sum(squeeze(sum(AddDim2(size(Kernel,1),squeeze(Keq_mod(T,Tref,x(k),Js_RNA(x),Js_Sites(x),m))).*C,3,'omitnan')),2,'omitnan')));
    h_DNA_MultiProbeConc = @(T,Tref,x,k,C,Kernel,m) full(squeeze(sum(squeeze(sum(AddDim2(size(Kernel,1),squeeze(Keq_mod(T,Tref,x(k),Js_DNA(x),Js_Sites(x),m))).*C,3,'omitnan')),2,'omitnan')));      
    h_sub_RNA1_MultiProbeConc_Func = @(T,CProbes_Free,Pset) removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) h_sub_RNA_MultiProbeConc(T,Tref,Pset,CProbes_Free,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
    h_sub_DNA1_MultiProbeConc_Func = @(T,CProbes_Free,Pset,Kernel) removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) h_sub_DNA_MultiProbeConc(T,Tref,Pset,CProbes_Free,Kernel,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
    h_sub_DNA1_complement_MultiProbeConc_Func = @(T,CProbes_Free,Pset,Kernel) removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) h_sub_DNA_complement_MultiProbeConc(T,Tref,Pset,CProbes_Free,Kernel,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
    h_RNAx_MultiProbeConc_Func = @(T,CProbes_Free,h_sub_RNA1,Pset,Kernel) CATnWrapper(arrayfun(@(J) sum(CATnWrapper(arrayfun(@(m) h_RNA_MultiProbeConc(T,Tref,Pset,J,h_sub_RNA1,Kernel,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4),1:length(Pset),'Un',0),2);
    h_DNAx_MultiProbeConc_Func = @(T,CProbes_Free,h_sub_DNA1,Pset,Kernel) CATnWrapper(arrayfun(@(J) sum(CATnWrapper(arrayfun(@(m) h_DNA_MultiProbeConc(T,Tref,Pset,J,h_sub_DNA1,Kernel,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4),1:length(Pset),'Un',0),2);
    K_S_Func = @(T,P) sum(sum(CATnWrapper(arrayfun(@(m) Ks_eq(T,Tref,P,1:max(Ns),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),3),3),2);
    K_CD_Func = @(T,P) sum(sum(CATnWrapper(arrayfun(@(m) Kd_eq(T,Tref,P,P,1:max(Nc(:)),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4),3);
    K_S = K_S_Func(T,Pset);
    K_CD = K_CD_Func(T,Pset);
    h_sub_RNA1_MultiProbeConc = h_sub_RNA1_MultiProbeConc_Func(T,CProbes_Free_MultiProbeConc,Pset);
    h_sub_DNA1_MultiProbeConc = h_sub_DNA1_MultiProbeConc_Func(T,CProbes_Free_MultiProbeConc,Pset,EKernel_MultiProbeConc);
    h_sub_DNA1_complement_MultiProbeConc = h_sub_DNA1_complement_MultiProbeConc_Func(T,CProbes_Free_MultiProbeConc,Pset,EKernel_MultiProbeConc);
    h_RNAx_MultiProbeConc = h_RNAx_MultiProbeConc_Func(T,CProbes_Free_MultiProbeConc,h_sub_RNA1_MultiProbeConc,Pset,EKernel_MultiProbeConc);
    h_DNAx_MultiProbeConc = h_DNAx_MultiProbeConc_Func(T,CProbes_Free_MultiProbeConc,h_sub_DNA1_MultiProbeConc,Pset,EKernel_MultiProbeConc);  
    ODE_MultiProbeConc = diff(CProbes_Free_MultiProbeConc) == MaxProbeConc*ProbeConcMatrix_MultiProbeConc./(1+AddDim1(size(EKernel_MultiProbeConc,1),K_S)+multiprod(K_CD,CProbes_Free_MultiProbeConc')'+ h_RNAx_MultiProbeConc + h_DNAx_MultiProbeConc)-CProbes_Free_MultiProbeConc;
    [M0_MultiProbeConc,F0_MultiProbeConc] = massMatrixForm(ODE_MultiProbeConc,CProbes_Free_MultiProbeConc);
    M_MultiProbeConc = odeFunction(M0_MultiProbeConc,CProbes_Free_MultiProbeConc,'Sparse',true);
    F_MultiProbeConc = odeFunction(F0_MultiProbeConc,CProbes_Free_MultiProbeConc,ProbeConcMatrix_MultiProbeConc,T,M);   
pRNA_TargetSites_Bound_MultiProbeConc = @(T,Pset,CProbes_Free,m) ...
    (permute(multiprod(permute(Keq_mod(T,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),m).*DoesProbeBindSite(Pset,Js_RNA(Pset),Js_Sites(Pset)),[2 1 3]),CProbes_Free'),[2 1 3]).*h_sub_RNA1_MultiProbeConc)./(...
     permute(multiprod(permute(Keq_mod(T,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),m).*DoesProbeBindSite(Pset,Js_RNA(Pset),Js_Sites(Pset)),[2 1 3]),CProbes_Free'),[2 1 3]).*h_sub_RNA1_MultiProbeConc + h_sub_RNA1_MultiProbeConc);
pDNA_TargetSites_Bound_MultiProbeConc = @(T,Pset,CProbes_Free,m) ...
    (permute(multiprod(permute(Keq_mod(T,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),m).*DoesProbeBindSite(Pset,Js_DNA(Pset),Js_Sites(Pset)),[2 1 3]),CProbes_Free'),[2 1 3]).*h_sub_DNA1_MultiProbeConc)./(...
     permute(multiprod(permute(Keq_mod(T,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),m).*DoesProbeBindSite(Pset,Js_DNA(Pset),Js_Sites(Pset)),[2 1 3]),CProbes_Free'),[2 1 3]).*h_sub_DNA1_MultiProbeConc + ...
    AddDim1(size(EKernel_MultiProbeConc,1),Keq_Complement(T,Tref,Js_DNA(Pset),Js_Sites(Pset),m)).*h_sub_DNA1_MultiProbeConc.*h_sub_DNA1_complement_MultiProbeConc + h_sub_DNA1_MultiProbeConc);
p_TargetSites_Bound_MultiProbeConc(1:size(EKernel_MultiProbeConc,1),Js_RNA(Pset),Js_Sites(Pset)) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pRNA_TargetSites_Bound_MultiProbeConc(T,Pset,CProbes_Free_MultiProbeConc,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
p_TargetSites_Bound_MultiProbeConc(1:size(EKernel_MultiProbeConc,1),Js_DNA(Pset),Js_Sites(Pset)) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pDNA_TargetSites_Bound_MultiProbeConc(T,Pset,CProbes_Free_MultiProbeConc,m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
p_TargetSites_Unbound_MultiProbeConc = 1 - p_TargetSites_Bound_MultiProbeConc;
p_TargetSites_Bound_MultiProbeConc_Func = matlabFunction(p_TargetSites_Bound_MultiProbeConc,'vars',{CProbes_Free_MultiProbeConc,T,M});
p_TargetSites_Unbound_MultiProbeConc_Func = matlabFunction(p_TargetSites_Unbound_MultiProbeConc,'vars',{CProbes_Free_MultiProbeConc,T,M});
p_TargetSites_Bound_MultiProbeConc_Weave = cell(1,N_Channels);
p_TargetSites_Unbound_MultiProbeConc_Weave = cell(1,N_Channels);
p_TargetSites_Bound_MultiProbeConc_Weave_Func = cell(1,N_Channels);
p_TargetSites_Unbound_MultiProbeConc_Weave_Func = cell(1,N_Channels);
p_TargetSites_Bound_MultiProbeConc_Split = cell(1,N_Channels);
p_TargetSites_Unbound_MultiProbeConc_Split = cell(1,N_Channels);
p_TargetSites_Bound_MultiProbeConc_Split_Func = cell(1,N_Channels);
p_TargetSites_Unbound_MultiProbeConc_Split_Func = cell(1,N_Channels);
for k = 1:N_Channels
    p_TargetSites_Bound_MultiProbeConc_Weave{k}(1:size(EKernel_MultiProbeConc,1),Js_RNA(Ix_Weave{k}),Js_Sites(Ix_Weave{k})) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pRNA_TargetSites_Bound_MultiProbeConc(T,Ix_Weave{k},CProbes_Free_MultiProbeConc(:,ismember(Pset,Ix_Weave{k})),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
    p_TargetSites_Bound_MultiProbeConc_Weave{k}(1:size(EKernel_MultiProbeConc,1),Js_DNA(Ix_Weave{k}),Js_Sites(Ix_Weave{k})) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pDNA_TargetSites_Bound_MultiProbeConc(T,Ix_Weave{k},CProbes_Free_MultiProbeConc(:,ismember(Pset,Ix_Weave{k})),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
    p_TargetSites_Unbound_MultiProbeConc_Weave{k} = 1 - p_TargetSites_Bound_MultiProbeConc_Weave{k};
    p_TargetSites_Bound_MultiProbeConc_Weave_Func{k} = matlabFunction(p_TargetSites_Bound_MultiProbeConc_Weave{k},'vars',{CProbes_Free_MultiProbeConc,T,M});
    p_TargetSites_Unbound_MultiProbeConc_Weave_Func{k} = matlabFunction(p_TargetSites_Unbound_MultiProbeConc_Weave{k},'vars',{CProbes_Free_MultiProbeConc,T,M});
    p_TargetSites_Bound_MultiProbeConc_Split{k}(1:size(EKernel_MultiProbeConc,1),Js_RNA(Ix_Split{k}),Js_Sites(Ix_Split{k})) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pRNA_TargetSites_Bound_MultiProbeConc(T,Ix_Split{k},CProbes_Free_MultiProbeConc(:,ismember(Pset,Ix_Split{k})),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
    p_TargetSites_Bound_MultiProbeConc_Split{k}(1:size(EKernel_MultiProbeConc,1),Js_DNA(Ix_Split{k}),Js_Sites(Ix_Split{k})) = removeNaNvInfWrapper(sum(CATnWrapper(arrayfun(@(m) pDNA_TargetSites_Bound_MultiProbeConc(T,Ix_Split{k},CProbes_Free_MultiProbeConc(:,ismember(Pset,Ix_Split{k})),m)*kroneckerDelta(M,m),1:N_methods,'Un',0),4),4));
    p_TargetSites_Unbound_MultiProbeConc_Split{k} = 1 - p_TargetSites_Bound_MultiProbeConc_Split{k};
    p_TargetSites_Bound_MultiProbeConc_Split_Func{k} = matlabFunction(p_TargetSites_Bound_MultiProbeConc_Split{k},'vars',{CProbes_Free_MultiProbeConc,T,M});
    p_TargetSites_Unbound_MultiProbeConc_Split_Func{k} = matlabFunction(p_TargetSites_Unbound_MultiProbeConc_Split{k},'vars',{CProbes_Free_MultiProbeConc,T,M});
end
Mean_G0_3D_MultiProbeConc = @(G,mean_probes_Target) sum(EKernel_MultiProbeConc(:,G).*mean_probes_Target(:,G),2)./sum(EKernel_MultiProbeConc(:,G),2);
Var_G0_3D_MultiProbeConc = @(G,var_probes_Target) sum(EKernel_MultiProbeConc(:,G).*var_probes_Target(:,G),2)./sum(EKernel_MultiProbeConc(:,G),2);
Mean_G_3D_MultiProbeConc = @(G,mean_probes_Target,p_TargetSites_Zero) sum(EKernel_MultiProbeConc(:,G).*mean_probes_Target(:,G),2)./sum(EKernel_MultiProbeConc(:,G).*(1-p_TargetSites_Zero(:,G)),2);
Var_G_3D_MultiProbeConc = @(G,var_probes_Target,p_TargetSites_Zero) sum(EKernel_MultiProbeConc(:,G).*var_probes_Target(:,G),2)./sum(EKernel_MultiProbeConc(:,G).*(1-p_TargetSites_Zero(:,G)),2);%edge cases
sConc_G_3D_MultiProbeConc = @(G,p_TargetSites_Zero) sum(EKernel_MultiProbeConc(:,G).*(1-p_TargetSites_Zero(:,G)),2); %count spots ~
pConc_G_3D_MultiProbeConc = @(G,mean_probes_Target) sum(EKernel_MultiProbeConc(:,G).*mean_probes_Target(:,G),2); %count probes  ~
sFrac_G_3D_MultiProbeConc = @(G,p_TargetSites_Zero) 100*sConc_G_3D_MultiProbeConc(G,p_TargetSites_Zero)./sConc_G_3D_MultiProbeConc(1:length(Tset),p_TargetSites_Zero);
pFrac_G_3D_MultiProbeConc = @(G,mean_probes_Target) 100*pConc_G_3D_MultiProbeConc(G,mean_probes_Target)./pConc_G_3D_MultiProbeConc(1:length(Tset),mean_probes_Target);
end
tspan = [0 10^6];
specific_2D = @(A,k) A(k,:);
specific_3D = @(A,k) A(:,k,:);
y_MomGen1 = @(bound) 1 - 2*bound;
y_MomGen2 = @(bound,unbound) 1 - 6*bound.*unbound;
y_Zero_2D = @(unbound) prod(unbound,2);
y_Zero_3D = @(unbound) prod(unbound,3);
y_Bound_ON_2D = @(ON_IDs,bound) specific_2D(bound,ON_IDs);
y_Bound_OFF_2D = @(OFF_IDs,bound) specific_2D(bound,OFF_IDs);
y_Unbound_ON_2D = @(ON_IDs,unbound) specific_2D(unbound,ON_IDs);
y_Unbound_OFF_2D = @(OFF_IDs,unbound) specific_2D(unbound,OFF_IDs); 
y_Bound_ON_3D = @(ON_IDs,bound) specific_3D(bound,ON_IDs);
y_Bound_OFF_3D = @(OFF_IDs,bound) specific_3D(bound,OFF_IDs);
y_Unbound_ON_3D = @(ON_IDs,unbound) specific_3D(unbound,ON_IDs);
y_Unbound_OFF_3D = @(OFF_IDs,unbound) specific_3D(unbound,OFF_IDs);
mean_probes_Target_2D = @(p_TargetSites_Bound) sum(p_TargetSites_Bound,2,'omitnan');
var_probes_Target_2D = @(p_TargetSites_Bound,p_TargetSites_Unbound) sum(p_TargetSites_Unbound.*p_TargetSites_Bound,2,'omitnan');
fano_probes_Target_2D = @(p_TargetSites_Bound,p_TargetSites_Unbound) sum(p_TargetSites_Unbound.*p_TargetSites_Bound,2,'omitnan')./sum(p_TargetSites_Bound,2,'omitnan');
skew_probes_Target_2D = @(p_TargetSites_Bound,p_TargetSites_Unbound,p_TargetSites_MomGen1) sum(p_TargetSites_MomGen1.*p_TargetSites_Unbound.*p_TargetSites_Bound,2,'omitnan');
kurt_probes_Target_2D = @(p_TargetSites_Bound,p_TargetSites_Unbound,p_TargetSites_MomGen2) sum(p_TargetSites_MomGen2.*p_TargetSites_Unbound.*p_TargetSites_Bound,2,'omitnan');
mean_probes_Target_3D = @(p_TargetSites_Bound) sum(p_TargetSites_Bound,3,'omitnan');
var_probes_Target_3D = @(p_TargetSites_Bound,p_TargetSites_Unbound) sum(p_TargetSites_Unbound.*p_TargetSites_Bound,3,'omitnan');
fano_probes_Target_3D = @(p_TargetSites_Bound,p_TargetSites_Unbound) sum(p_TargetSites_Unbound.*p_TargetSites_Bound,3,'omitnan')./sum(p_TargetSites_Bound,3,'omitnan');
skew_probes_Target_3D = @(p_TargetSites_Bound,p_TargetSites_Unbound,p_TargetSites_MomGen1) sum(p_TargetSites_MomGen1.*p_TargetSites_Unbound.*p_TargetSites_Bound,3,'omitnan');
kurt_probes_Target_3D = @(p_TargetSites_Bound,p_TargetSites_Unbound,p_TargetSites_MomGen2) sum(p_TargetSites_MomGen2.*p_TargetSites_Unbound.*p_TargetSites_Bound,3,'omitnan');

if (solution_case == 1)
    PC0 = repmat(ProbeConc_Init,[1 length(Pset)]);  
    CPF0 = GuessConc*ones(length(Pset),1);  
    opts_basic = @(T,M,C) odeset('mass', M_basic, 'RelTol', 10^(-6),'AbsTol', 10^(-6), 'InitialSlope', F_basic(0,CPF0,PC0,T,M,C));
    y_basic = @(T,M,C) ode45(@(t,y) F_basic(t,y,PC0,T,M,C),tspan,CPF0,opts_basic(T,M,C)).y(:,end);
    y_Bound = @(T,M,C) p_TargetSites_Bound_Func(y_basic(T,M,C),T,M,C); 
    y_Unbound = @(T,M,C) p_TargetSites_Unbound_Func(y_basic(T,M,C),T,M,C);
    y_Bound_Weave = @(T,M,C,k) p_TargetSites_Bound_Weave_Func{k}(y_basic(T,M,C),T,M,C); 
    y_Unbound_Weave = @(T,M,C,k) p_TargetSites_Unbound_Weave_Func{k}(y_basic(T,M,C),T,M,C);
    y_Bound_Split = @(T,M,C,k) p_TargetSites_Bound_Split_Func{k}(y_basic(T,M,C),T,M,C); 
    y_Unbound_Split = @(T,M,C,k) p_TargetSites_Unbound_Split_Func{k}(y_basic(T,M,C),T,M,C);
    for Ti = T_low_celsius:T_high_celsius
       for m = 1:N_methods
          for c = 1:size(ExpressionMatrix,2) 
             for k = 1:N_Channels  
                Sol(Ti,m,c).Bound_Weave{k} = y_Bound_Weave(Ti+273.15,m,c,k);
                Sol(Ti,m,c).Unbound_Weave{k} = y_Unbound_Weave(Ti+273.15,m,c,k);
                Sol(Ti,m,c).Bound_Split{k} = y_Bound_Split(Ti+273.15,m,c,k);
                Sol(Ti,m,c).Unbound_Split{k} = y_Unbound_Split(Ti+273.15,m,c,k);  
                PxFD_Weave{k} = arrayfun(@(J) F_DiscretePossionBinomial(P_1FLAPs(2)*Sol(Ti,m,c).Bound_Weave{k}(Tset(J),vsi_Weave_Uni{k}(J))),1:length(Tset),'Un',0);% unique 1 each 
                PxFD_Split{k} = arrayfun(@(J) F_DiscretePossionBinomial(P_1FLAPs(2)*Sol(Ti,m,c).Bound_Split{k}(Tset(J),vsi_Split_Uni{k}(J))),1:length(Tset),'Un',0);% unique 1 each   
             end 
             Sol(Ti,m,c).Bound = y_Bound(Ti+273.15,m,c);
             Sol(Ti,m,c).Unbound = y_Unbound(Ti+273.15,m,c);
             Sol(Ti,m,c).MomGen1 = y_MomGen1(Sol(Ti,m,c).Bound);
             Sol(Ti,m,c).MomGen2 = y_MomGen2(Sol(Ti,m,c).Bound,Sol(Ti,m,c).Unbound);
             Sol(Ti,m,c).Zero = y_Zero_2D(Sol(Ti,m,c).Unbound);
             Sol(Ti,m,c).Bound_ON = y_Bound_ON_2D(ON_IDs,Sol(Ti,m,c).Bound);
             Sol(Ti,m,c).Bound_OFF = y_Bound_OFF_2D(OFF_IDs,Sol(Ti,m,c).Bound);
             Sol(Ti,m,c).Unbound_ON = y_Unbound_ON_2D(ON_IDs,Sol(Ti,m,c).Unbound);
             Sol(Ti,m,c).Unbound_OFF = y_Unbound_OFF_2D(OFF_IDs,Sol(Ti,m,c).Unbound);
             Sol(Ti,m,c).Mean = mean_probes_Target_2D(Sol(Ti,m,c).Bound);
             Sol(Ti,m,c).Var  = var_probes_Target_2D(Sol(Ti,m,c).Bound,Sol(Ti,m,c).Unbound);
             Sol(Ti,m,c).Fano = fano_probes_Target_2D(Sol(Ti,m,c).Bound,Sol(Ti,m,c).Unbound);
             Sol(Ti,m,c).Skew = skew_probes_Target_2D(Sol(Ti,m,c).Bound,Sol(Ti,m,c).Unbound,Sol(Ti,m,c).MomGen1);
             Sol(Ti,m,c).Kurt = kurt_probes_Target_2D(Sol(Ti,m,c).Bound,Sol(Ti,m,c).Unbound,Sol(Ti,m,c).MomGen2);
             Sol(Ti,m,c).Mean0_ON = Mean_G0_2D(ON_IDs,c,Sol(Ti,m,c).Mean); 
             Sol(Ti,m,c).STD0_ON = Var_G0_2D(ON_IDs,c,Sol(Ti,m,c).Var);
             Sol(Ti,m,c).Mean_ON = Mean_G_2D(ON_IDs,c,Sol(Ti,m,c).Mean,Sol(Ti,m,c).Zero); 
             Sol(Ti,m,c).STD_ON = Var_G_2D(ON_IDs,c,Sol(Ti,m,c).Var,Sol(Ti,m,c).Zero); 
             Sol(Ti,m,c).sConc_ON = sConc_G_2D(ON_IDs,c,Sol(Ti,m,c).Zero);
             Sol(Ti,m,c).pConc_ON = pConc_G_2D(ON_IDs,c,Sol(Ti,m,c).Mean); 
             Sol(Ti,m,c).sFrac_ON = sFrac_G_2D(ON_IDs,c,Sol(Ti,m,c).Zero); 
             Sol(Ti,m,c).pFrac_ON = pFrac_G_2D(ON_IDs,c,Sol(Ti,m,c).Mean); 
             Sol(Ti,m,c).Mean0_OFF =  Mean_G0_2D(OFF_IDs,c,Sol(Ti,m,c).Mean); 
             Sol(Ti,m,c).STD0_OFF = Var_G0_2D(OFF_IDs,c,Sol(Ti,m,c).Var);
             Sol(Ti,m,c).Mean_OFF = Mean_G_2D(OFF_IDs,c,Sol(Ti,m,c).Mean,Sol(Ti,m,c).Zero); 
             Sol(Ti,m,c).STD_OFF = Var_G_2D(OFF_IDs,c,Sol(Ti,m,c).Var,Sol(Ti,m,c).Zero); 
             Sol(Ti,m,c).sConc_OFF = sConc_G_2D(OFF_IDs,c,Sol(Ti,m,c).Zero);
             Sol(Ti,m,c).pConc_OFF = pConc_G_2D(OFF_IDs,c,Sol(Ti,m,c).Mean); 
             Sol(Ti,m,c).sFrac_OFF = sFrac_G_2D(OFF_IDs,c,Sol(Ti,m,c).Zero); 
             Sol(Ti,m,c).pFrac_OFF = pFrac_G_2D(OFF_IDs,c,Sol(Ti,m,c).Mean);
             Sol(Ti,m,c).PxF = arrayfun(@(J) F_DiscretePossionBinomial(P_1FLAPs(2)*Sol(Ti,m,c).Bound(Tset(J),find(sum(squeeze(DoesProbeBindSite(Pset,Tset(J),:)),1)>0))),1:length(Tset),'Un',0);    
             Sol(Ti,m,c).PxF_Matrix = vertcat(Sol(Ti,m,c).PxF{:});
             Sol(Ti,m,c).PxF_ExprMatrix = ExpressionMatrix(Tset,c).*Sol(Ti,m,c).PxF_Matrix;    
             PxF2D_Weave = arrayfun(@(J) PxFD_Weave{2}{J}.*PxFD_Weave{1}{J}',1:length(Tset),'Un',0);
             PxF2D_Split = arrayfun(@(J) PxFD_Split{2}{J}.*PxFD_Split{1}{J}',1:length(Tset),'Un',0);
             temp_Weave2 = arrayfun(@(J) PxF2D_Weave{J}./sum(sum(PxF2D_Weave{J},1),2),1:length(Tset),'Un',0);
             temp_Split2 = arrayfun(@(J) PxF2D_Split{J}./sum(sum(PxF2D_Split{J},1),2),1:length(Tset),'Un',0);
             PxF2D_Weave = temp_Weave2;
             PxF2D_Split = temp_Split2;
             p1_Split = arrayfun(@(J) Sol(Ti,m,c).Bound_Split{1}(Tset(J), overlap_Split(J)),1:length(Tset),'Un',0);%1st dimension k = 1
             p2_Split = arrayfun(@(J) Sol(Ti,m,c).Bound_Split{2}(Tset(J), overlap_Split(J)),1:length(Tset),'Un',0);%2nd dimension k = 2
             p1_Weave = arrayfun(@(J) Sol(Ti,m,c).Bound_Weave{1}(Tset(J), overlap_Weave(J)),1:length(Tset),'Un',0);%1st dimension k = 1
             p2_Weave = arrayfun(@(J) Sol(Ti,m,c).Bound_Weave{2}(Tset(J), overlap_Weave(J)),1:length(Tset),'Un',0);%2nd dimension k = 2
             Split_Basis_Matrix = arrayfun(@(J) [p1_Split{J};p2_Split{J}],1:length(Tset),'Un',0);
             Weave_Basis_Matrix = arrayfun(@(J) [p1_Weave{J};p2_Weave{J}],1:length(Tset),'Un',0);
             PxF2D_Split = arrayfun(@(J) F_DiscretePossionBinomial_Base2(PxF2D_Split{J},Split_Basis_Matrix{J}),1:length(Tset),'Un',0);
             PxF2D_Weave = arrayfun(@(J) F_DiscretePossionBinomial_Base2(PxF2D_Weave{J},Weave_Basis_Matrix{J}),1:length(Tset),'Un',0);
             Dim_Split = cell2mat(arrayfun(@(J) max(size(PxF2D_Split{J}))+1,1:length(Tset),'Un',0));
             Dim_Weave = cell2mat(arrayfun(@(J) max(size(PxF2D_Weave{J}))+1,1:length(Tset),'Un',0));
             DimL = max([Dim_Split Dim_Weave]);
             PxF2D_Split = arrayfun(@(J) padarray(PxF2D_Split{J},[max(size(PxF2D_Split{J}))+1-size(PxF2D_Split{J},1) max(size(PxF2D_Split{J}))+1-size(PxF2D_Split{J},2)],0,'post'),1:length(Tset),'Un',0);
             PxF2D_Weave = arrayfun(@(J) padarray(PxF2D_Weave{J},[max(size(PxF2D_Weave{J}))+1-size(PxF2D_Weave{J},1) max(size(PxF2D_Weave{J}))+1-size(PxF2D_Weave{J},2)],0,'post'),1:length(Tset),'Un',0);
             PxF2D_Expr_Split = arrayfun(@(J) ExpressionMatrix(Tset(J),c)*PxF2D_Split{J},1:length(Tset),'Un',0);
             PxF2D_Expr_Weave = arrayfun(@(J) ExpressionMatrix(Tset(J),c)*PxF2D_Weave{J},1:length(Tset),'Un',0); 
             Sol(Ti,m,c).PxF2D_Split = arrayfun(@(J) padarray(PxF2D_Split{J},[DimL-size(PxF2D_Split{J},1) DimL-size(PxF2D_Split{J},2)],0,'post'),1:length(Tset),'Un',0);
             Sol(Ti,m,c).PxF2D_Expr_Split = arrayfun(@(J) padarray(PxF2D_Expr_Split{J},[DimL-size(PxF2D_Expr_Split{J},1) DimL-size(PxF2D_Expr_Split{J},2)],0,'post'),1:length(Tset),'Un',0);
             Sol(Ti,m,c).PxF2D_Weave = arrayfun(@(J) padarray(PxF2D_Weave{J},[DimL-size(PxF2D_Weave{J},1) DimL-size(PxF2D_Weave{J},2)],0,'post'),1:length(Tset),'Un',0);
             Sol(Ti,m,c).PxF2D_Expr_Weave = arrayfun(@(J) padarray(PxF2D_Expr_Weave{J},[DimL-size(PxF2D_Expr_Weave{J},1) DimL-size(PxF2D_Expr_Weave{J},2)],0,'post'),1:length(Tset),'Un',0);
             %check ability to distiguish isoforms (P(x) similarity)
             if (length(Isoforms)>1)
                 for u = 1:length(Isoforms)
                    for v = 1:length(Isoforms)
                    Sol(Ti,m,c).D_Isoforms(u,v) = 1/2*trapz(0:Smax,abs(Sol(Ti,m,c).PxF{Corr_ONRNA_IsoAgnostic(u)}-Sol(Ti,m,c).PxF{Corr_ONRNA_IsoAgnostic(v)}));
                    Sol(Ti,m,c).D_Isoforms_Split(u,v) = 1/2*trapz(0:DimL-1,trapz(0:DimL-1,abs(Sol(Ti,m,c).PxF2D_Split{Corr_ONRNA_IsoAgnostic(u)}-Sol(Ti,m,c).PxF2D_Split{Corr_ONRNA_IsoAgnostic(v)})));
                    Sol(Ti,m,c).D_Isoforms_Weave(u,v) = 1/2*trapz(0:DimL-1,trapz(0:DimL-1,abs(Sol(Ti,m,c).PxF2D_Weave{Corr_ONRNA_IsoAgnostic(u)}-Sol(Ti,m,c).PxF2D_Weave{Corr_ONRNA_IsoAgnostic(v)})));
                    end
                 end
             else
                 Sol(Ti,m,c).D_Isoforms = NaN;
                 Sol(Ti,m,c).D_Isoforms_Split = NaN;
                 Sol(Ti,m,c).D_Isoforms_Weave = NaN;
             end
             Sol(Ti,m,c).MeanIsoIdentifiability = mean(Sol(Ti,m,c).D_Isoforms(triu(true(size(Sol(Ti,m,c).D_Isoforms)),1)));%Optimize Isoform identifiability (mean,var)
             Sol(Ti,m,c).SplitMeanIsoIdentifiability = mean(Sol(Ti,m,c).D_Isoforms_Split(triu(true(size(Sol(Ti,m,c).D_Isoforms_Split)),1)));%Optimize Isoform identifiability (mean,var)    
             Sol(Ti,m,c).WeaveMeanIsoIdentifiability = mean(Sol(Ti,m,c).D_Isoforms_Weave(triu(true(size(Sol(Ti,m,c).D_Isoforms_Weave)),1)));%Optimize Isoform identifiability (mean,var)    
             %Gibbs Sampling
             W_Norm = sum(Sol(Ti,m,c).PxF_ExprMatrix(:,2:end),1);%count dist
             W_ONDNA = sum(Sol(Ti,m,c).PxF_ExprMatrix(Corr_ONDNA,2:end),1);%count dist
             W_ONRNA_IsoSpecific = sum(Sol(Ti,m,c).PxF_ExprMatrix(Corr_ONRNA_IsoSpecific,2:end),1);%count dist
             W_ONRNA_IsoAgnostic = sum(Sol(Ti,m,c).PxF_ExprMatrix(Corr_ONRNA_IsoAgnostic,2:end),1);%count dist
             W_OFFDNA = sum(Sol(Ti,m,c).PxF_ExprMatrix(Corr_OFFDNA,2:end),1);%count dist
             W_OFFRNA_withIsos = sum(Sol(Ti,m,c).PxF_ExprMatrix(Corr_OFFRNA_withIsos,2:end),1);%count dist
             W_OFFRNA_minusIso = sum(Sol(Ti,m,c).PxF_ExprMatrix(Corr_OFFRNA_minusIso,2:end),1);%count dist
             W_OFFRNA_otherIso = sum(Sol(Ti,m,c).PxF_ExprMatrix(Corr_OFFRNA_otherIso,2:end),1);%count dist 
             Spot_Detect_Max = sum(W_Norm);
             Spot_Detect_ON_IsoSpecific = sum(W_ONDNA) + sum(W_ONRNA_IsoSpecific);
             Spot_Detect_ON_IsoAgnostic = sum(W_ONDNA) + sum(W_ONRNA_IsoAgnostic);
             Spot_Detect_OFF_withIsos = sum(W_OFFDNA) + sum(W_OFFRNA_withIsos);
             Spot_Detect_OFF_minusIso = sum(W_OFFDNA) + sum(W_OFFRNA_minusIso);
             Spot_Detect_OFF_otherIso = sum(W_OFFDNA) + sum(W_OFFRNA_otherIso);
             Probe_Detect_Max = dot(1:length(W_Norm),W_Norm); 
             Probe_Detect_ON_IsoSpecific = dot(1:length(W_ONDNA),W_ONDNA) + dot(1:length(W_ONRNA_IsoSpecific),W_ONRNA_IsoSpecific);
             Probe_Detect_ON_IsoAgnostic = dot(1:length(W_ONDNA),W_ONDNA) + dot(1:length(W_ONRNA_IsoAgnostic),W_ONRNA_IsoAgnostic);
             Probe_Detect_OFF_withIsos = dot(1:length(W_OFFDNA),W_OFFDNA) + dot(1:length(W_OFFRNA_withIsos),W_OFFRNA_withIsos);
             Probe_Detect_OFF_minusIso = dot(1:length(W_OFFDNA),W_OFFDNA) + dot(1:length(W_OFFRNA_minusIso),W_OFFRNA_minusIso);
             Probe_Detect_OFF_otherIso = dot(1:length(W_OFFRNA_otherIso),W_OFFRNA_otherIso);
             Sol(Ti,m,c).Spot.Detection.Max = Spot_Detect_Max;
             Sol(Ti,m,c).Spot.Detection.ON_IsoSpecific = Spot_Detect_ON_IsoSpecific;
             Sol(Ti,m,c).Spot.Detection.ON_IsoAgnostic = Spot_Detect_ON_IsoAgnostic;
             Sol(Ti,m,c).Spot.Detection.OFF_withIsos = Spot_Detect_OFF_withIsos;
             Sol(Ti,m,c).Spot.Detection.OFF_minusIso = Spot_Detect_OFF_minusIso;
             Sol(Ti,m,c).Spot.Detection.OFF_otherIso = Spot_Detect_OFF_otherIso;
             Sol(Ti,m,c).Spot.Percentage.ON_IsoSpecific = Spot_Detect_ON_IsoSpecific/Spot_Detect_Max;
             Sol(Ti,m,c).Spot.Percentage.ON_IsoAgnostic = Spot_Detect_ON_IsoAgnostic/Spot_Detect_Max;
             Sol(Ti,m,c).Spot.Percentage.OFF_withIsos = Spot_Detect_OFF_withIsos/Spot_Detect_Max;
             Sol(Ti,m,c).Spot.Percentage.OFF_minusIso = Spot_Detect_OFF_minusIso/Spot_Detect_Max;
             Sol(Ti,m,c).Spot.Percentage.OFF_otherIso = Spot_Detect_OFF_otherIso/Spot_Detect_Max;
             Sol(Ti,m,c).Probe.Detection.Max = Probe_Detect_Max;
             Sol(Ti,m,c).Probe.Detection.ON_IsoSpecific = Probe_Detect_ON_IsoSpecific;
             Sol(Ti,m,c).Probe.Detection.ON_IsoAgnostic = Probe_Detect_ON_IsoAgnostic;
             Sol(Ti,m,c).Probe.Detection.OFF_withIsos = Probe_Detect_OFF_withIsos;
             Sol(Ti,m,c).Probe.Detection.OFF_minusIso = Probe_Detect_OFF_minusIso;
             Sol(Ti,m,c).Probe.Detection.OFF_otherIso = Probe_Detect_OFF_otherIso;
             Sol(Ti,m,c).Probe.Percentage.ON_IsoSpecific = Probe_Detect_ON_IsoSpecific/Probe_Detect_Max;
             Sol(Ti,m,c).Probe.Percentage.ON_IsoAgnostic = Probe_Detect_ON_IsoAgnostic/Probe_Detect_Max;
             Sol(Ti,m,c).Probe.Percentage.OFF_withIsos = Probe_Detect_OFF_withIsos/Probe_Detect_Max;
             Sol(Ti,m,c).Probe.Percentage.OFF_minusIso = Probe_Detect_OFF_withIsos/Probe_Detect_Max;
             Sol(Ti,m,c).Probe.Percentage.OFF_otherIso = Probe_Detect_OFF_otherIso/Probe_Detect_Max;
             for groups = 1:5
                 switch groups 
                    case 1%RNA_ON_Target
                    Gx = arrayfun(@(x) find(Tset==onRNA(x)),1:length(onRNA),'Un',0);
                    withZero = 1;
                    case 2%RNA_ON_Target_IsoAgnostic
                    Gx = arrayfun(@(x) find(Tset==Isoforms(x)),1:length(Isoforms),'Un',0);
                    withZero = 1;
                    case 3%RNA_OFF_Target_withIso
                    Gx = arrayfun(@(x) find(Tset==offRNA(x)),1:length(offRNA),'Un',0);
                    withZero = 0;
                    case 4%RNA_OFF_Target_withoutIso
                    Gx = arrayfun(@(x) find(Tset==offRNA_minusIso(x)),1:length(offRNA_minusIso),'Un',0);
                    withZero = 0;
                    case 5%RNA_ON_EachIso
                    Gx = arrayfun(@(x) find(Tset==Isoforms(x)),1:length(Isoforms),'Un',0);
                    withZero = 1;                    
                    case 6%DNA_ON_Target
                    Gx = arrayfun(@(x) find(Tset==onDNA(x)),1:length(onDNA),'Un',0);
                    withZero = 1;
                    case 7%DNA_OFF_Target
                    Gx = arrayfun(@(x) find(Tset==offDNA(x)),1:length(offDNA),'Un',0);
                    withZero = 0;
                 end
                 Gx = [Gx{:}];%ON-DNA, ON-RNA, ON-IsoSpec, ON-IsoAgnostic, OFF-RNAwoIso, OFF-RNAwIso, OFF-DNA    
                 if (groups~=5)
                 Sol(Ti,m,c).CountY(groups,:) = PxF_Count_G(Gx,Sol(Ti,m,c).PxF_ExprMatrix);
                 Sol(Ti,m,c).CountCY(groups,:) = cumsum(PxF_Count_G(Gx,Sol(Ti,m,c).PxF_ExprMatrix));
                 Sol(Ti,m,c).DistY(groups,:) = PxF_Dist_G(Gx,Sol(Ti,m,c).PxF_ExprMatrix);
                 Sol(Ti,m,c).DistCY(groups,:) = cumsum(PxF_Dist_G(Gx,Sol(Ti,m,c).PxF_ExprMatrix));
                 Sol(Ti,m,c).DistX(groups,:) = 0:Smax;
                 Sol(Ti,m,c).Mean(groups) = PxF_Mean_G(Gx,Sol(Ti,m,c).PxF_ExprMatrix);
                 Sol(Ti,m,c).Fano(groups) = PxF_Fano_G(Gx,Sol(Ti,m,c).PxF_ExprMatrix);
                 Sol(Ti,m,c).Skew(groups) = PxF_Moment_G(Gx,Sol(Ti,m,c).PxF_ExprMatrix,3);
                 Sol(Ti,m,c).Kurt(groups) = PxF_Moment_G(Gx,Sol(Ti,m,c).PxF_ExprMatrix,4);           
                 else
                     for vs = 1:length(Gx)
                     Sol(Ti,m,c).CountY(groups-1+vs,:) = PxF_Count_G(Gx(vs),Sol(Ti,m,c).PxF_ExprMatrix);
                     Sol(Ti,m,c).CountCY(groups-1+vs,:) = cumsum(PxF_Count_G(Gx(vs),Sol(Ti,m,c).PxF_ExprMatrix));
                     Sol(Ti,m,c).DistY(groups-1+vs,:) = PxF_Dist_G(Gx(vs),Sol(Ti,m,c).PxF_ExprMatrix); 
                     Sol(Ti,m,c).DistCY(groups-1+vs,:) = cumsum(PxF_Dist_G(Gx(vs),Sol(Ti,m,c).PxF_ExprMatrix));
                     Sol(Ti,m,c).DistX(groups-1+vs,:) = 0:Smax;
                     Sol(Ti,m,c).Mean(groups-1+vs) = PxF_Mean_G(Gx(vs),Sol(Ti,m,c).PxF_ExprMatrix);
                     Sol(Ti,m,c).Fano(groups-1+vs) = PxF_Fano_G(Gx(vs),Sol(Ti,m,c).PxF_ExprMatrix);
                     Sol(Ti,m,c).Skew(groups-1+vs) = PxF_Moment_G(Gx(vs),Sol(Ti,m,c).PxF_ExprMatrix,3);
                     Sol(Ti,m,c).Kurt(groups-1+vs) = PxF_Moment_G(Gx(vs),Sol(Ti,m,c).PxF_ExprMatrix,4);
                     end
                 end  
                 if (groups~=5)
                 Mat_Split = sum(CATnWrapper(Sol(Ti,m,c).PxF2D_Expr_Split{Gx},3),3);
                 Mat_Weave = sum(CATnWrapper(Sol(Ti,m,c).PxF2D_Expr_Weave{Gx},3),3);
                 Sol(Ti,m,c) = SubStats_2D(groups,Mat_Split,Sol(Ti,m,c),1,withZero,maxProbes);
                 Sol(Ti,m,c) = SubStats_2D(groups,Mat_Weave,Sol(Ti,m,c),2,withZero,maxProbes);
                 else
                     for vs = 1:length(Gx)
                     Sol(Ti,m,c) = SubStats_2D(groups-1+vs,Sol(Ti,m,c).PxF2D_Expr_Split{Gx(vs)},Sol(Ti,m,c),1,withZero,maxProbes);
                     Sol(Ti,m,c) = SubStats_2D(groups-1+vs,Sol(Ti,m,c).PxF2D_Expr_Weave{Gx(vs)},Sol(Ti,m,c),2,withZero,maxProbes);
                     end        
                 end
             end   
          end
       end  
    end      
elseif (solution_case == 2)
    PC0 = repmat(ProbeConc_Init,[size(ExpressionMatrix,2) length(Pset)]);  
    CPF0 = GuessConc*ones(size(EKernel_MultiTissue,1),length(Pset));
    opts_MultiTissue = @(T,M) odeset('mass', M_MultiTissue, 'RelTol', 10^(-6),'AbsTol', 10^(-6), 'InitialSlope', F_MultiTissue(0,reshape(CPF0,[numel(CPF0) 1]),PC0,T,M));    
    y_MultiTissue = @(T,M) reshape(ode45(@(t,y) F_MultiTissue(t,y,PC0,T,M),tspan,reshape(CPF0,[numel(CPF0) 1]),opts_MultiTissue(T,M)).y(:,end),size(CPF0));
    y_Bound_MultiTissue = @(T,M) p_TargetSites_Bound_MultiTissue_Func(y_MultiTissue(T,M),T,M);
    y_Unbound_MultiTissue = @(T,M) p_TargetSites_Unbound_MultiTissue_Func(y_MultiTissue(T,M),T,M);
    y_Bound_MultiTissue_Weave = @(T,M,k) p_TargetSites_Bound_MultiTissue_Weave_Func{k}(y_MultiTissue(T,M),T,M);
    y_Unbound_MultiTissue_Weave = @(T,M,k) p_TargetSites_Unbound_MultiTissue_Weave_Func{k}(y_MultiTissue(T,M),T,M);
    y_Bound_MultiTissue_Split = @(T,M,k) p_TargetSites_Bound_MultiTissue_Split_Func{k}(y_MultiTissue(T,M),T,M);
    y_Unbound_MultiTissue_Split = @(T,M,k) p_TargetSites_Unbound_MultiTissue_Split_Func{k}(y_MultiTissue(T,M),T,M);
    for Ti = T_low_celsius:T_high_celsius
       for m = 1:N_methods
          for k = 1:N_Channels
             Sol(Ti,m).Bound_Weave{k} = y_Bound_MultiTissue_Weave(Ti+273.15,m,k);
             Sol(Ti,m).Unbound_Weave{k} = y_Unbound_MultiTissue_Weave(Ti+273.15,m,k);
             Sol(Ti,m).Bound_Split{k} = y_Bound_MultiTissue_Split(Ti+273.15,m,k);
             Sol(Ti,m).Unbound_Split{k} = y_Unbound_MultiTissue_Split(Ti+273.15,m,k);
             PxFD_Weave{k} = arrayfun(@(c) arrayfun(@(J) F_DiscretePossionBinomial(P_1FLAPs(2)*Sol(Ti,m).Bound_Weave{k}(c,Tset(J),vsi_Weave_Uni{k}(J))),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
             PxFD_Split{k} = arrayfun(@(c) arrayfun(@(J) F_DiscretePossionBinomial(P_1FLAPs(2)*Sol(Ti,m).Bound_Split{k}(c,Tset(J),vsi_Split_Uni{k}(J))),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          end
          Sol(Ti,m).Bound = y_Bound_MultiTissue(Ti+273.15,m);
          Sol(Ti,m).Unbound = y_Unbound_MultiTissue(Ti+273.15,m); 
          Sol(Ti,m).MomGen1 = y_MomGen1(Sol(Ti,m).Bound);
          Sol(Ti,m).MomGen2 = y_MomGen2(Sol(Ti,m).Bound,Sol(Ti,m).Unbound);
          Sol(Ti,m).Zero = y_Zero_3D(Sol(Ti,m).Unbound);
          Sol(Ti,m).Bound_ON = y_Bound_ON_3D(ON_IDs,Sol(Ti,m).Bound);
          Sol(Ti,m).Bound_OFF = y_Bound_OFF_3D(OFF_IDs,Sol(Ti,m).Bound);
          Sol(Ti,m).Unbound_ON = y_Unbound_ON_3D(ON_IDs,Sol(Ti,m).Unbound);
          Sol(Ti,m).Unbound_OFF = y_Unbound_OFF_3D(OFF_IDs,Sol(Ti,m).Unbound);
          Sol(Ti,m).Mean = mean_probes_Target_3D(Sol(Ti,m).Bound);
          Sol(Ti,m).Var  = var_probes_Target_3D(Sol(Ti,m).Bound,Sol(Ti,m).Unbound);
          Sol(Ti,m).Fano = fano_probes_Target_3D(Sol(Ti,m).Bound,Sol(Ti,m).Unbound);
          Sol(Ti,m).Skew = skew_probes_Target_3D(Sol(Ti,m).Bound,Sol(Ti,m).Unbound,Sol(Ti,m).MomGen1);
          Sol(Ti,m).Kurt = kurt_probes_Target_3D(Sol(Ti,m).Bound,Sol(Ti,m).Unbound,Sol(Ti,m).MomGen2);
          Sol(Ti,m).Mean0_ON = Mean_G0_3D_MultiTissue(ON_IDs,Sol(Ti,m).Mean); 
          Sol(Ti,m).STD0_ON = Var_G0_3D_MultiTissue(ON_IDs,Sol(Ti,m).Var);
          Sol(Ti,m).Mean_ON = Mean_G_3D_MultiTissue(ON_IDs,Sol(Ti,m).Mean,Sol(Ti,m).Zero); 
          Sol(Ti,m).STD_ON = Var_G_3D_MultiTissue(ON_IDs,Sol(Ti,m).Var,Sol(Ti,m).Zero); 
          Sol(Ti,m).sConc_ON = sConc_G_3D_MultiTissue(ON_IDs,Sol(Ti,m).Zero);
          Sol(Ti,m).pConc_ON = pConc_G_3D_MultiTissue(ON_IDs,Sol(Ti,m).Mean); 
          Sol(Ti,m).sFrac_ON = sFrac_G_3D_MultiTissue(ON_IDs,Sol(Ti,m).Zero); 
          Sol(Ti,m).pFrac_ON = pFrac_G_3D_MultiTissue(ON_IDs,Sol(Ti,m).Mean); 
          Sol(Ti,m).Mean0_OFF =  Mean_G0_3D_MultiTissue(OFF_IDs,Sol(Ti,m).Mean); 
          Sol(Ti,m).STD0_OFF = Var_G0_3D_MultiTissue(OFF_IDs,Sol(Ti,m).Var);
          Sol(Ti,m).Mean_OFF = Mean_G_3D_MultiTissue(OFF_IDs,Sol(Ti,m).Mean,Sol(Ti,m).Zero); 
          Sol(Ti,m).STD_OFF = Var_G_3D_MultiTissue(OFF_IDs,Sol(Ti,m).Var,Sol(Ti,m).Zero); 
          Sol(Ti,m).sConc_OFF = sConc_G_3D_MultiTissue(OFF_IDs,Sol(Ti,m).Zero);
          Sol(Ti,m).pConc_OFF = pConc_G_3D_MultiTissue(OFF_IDs,Sol(Ti,m).Mean); 
          Sol(Ti,m).sFrac_OFF = sFrac_G_3D_MultiTissue(OFF_IDs,Sol(Ti,m).Zero); 
          Sol(Ti,m).pFrac_OFF = pFrac_G_3D_MultiTissue(OFF_IDs,Sol(Ti,m).Mean); 
          Sol(Ti,m).PxF = arrayfun(@(c) arrayfun(@(J) F_DiscretePossionBinomial(P_1FLAPs(2)*Sol(Ti,m).Bound(c,Tset(J),find(sum(squeeze(DoesProbeBindSite(Pset,Tset(J),:)),1)>0))),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);       
          Sol(Ti,m).PxF_Matrix = CATnWrapper(arrayfun(@(c) vertcat(Sol(Ti,m).PxF{c}{:}),1:size(EKernel_MultiTissue,1),'Un',0),3); 
          Sol(Ti,m).PxF_ExprMatrix = CATnWrapper(arrayfun(@(c) ExpressionMatrix(Tset,c).*squeeze(Sol(Ti,m).PxF_Matrix(:,:,c)),1:size(EKernel_MultiTissue,1),'Un',0),3);    
          PxF2D_Weave = arrayfun(@(c) arrayfun(@(J) PxFD_Weave{2}{c}{J}.*PxFD_Weave{1}{c}{J}',1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          PxF2D_Split = arrayfun(@(c) arrayfun(@(J) PxFD_Split{2}{c}{J}.*PxFD_Split{1}{c}{J}',1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          temp_Weave2 = arrayfun(@(c) arrayfun(@(J) PxF2D_Weave{c}{J}./sum(sum(PxF2D_Weave{c}{J},1),2),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          temp_Split2 = arrayfun(@(c) arrayfun(@(J) PxF2D_Split{c}{J}./sum(sum(PxF2D_Split{c}{J},1),2),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          PxF2D_Weave = temp_Weave2;
          PxF2D_Split = temp_Split2;  
          p1_Split = arrayfun(@(c) arrayfun(@(J) Sol(Ti,m).Bound_Split{1}(c,Tset(J),overlap_Split(J)),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          p2_Split = arrayfun(@(c) arrayfun(@(J) Sol(Ti,m).Bound_Split{2}(c,Tset(J),overlap_Split(J)),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          p1_Weave = arrayfun(@(c) arrayfun(@(J) Sol(Ti,m).Bound_Weave{1}(c,Tset(J),overlap_Weave(J)),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          p2_Weave = arrayfun(@(c) arrayfun(@(J) Sol(Ti,m).Bound_Weave{2}(c,Tset(J),overlap_Weave(J)),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          Split_Basis_Matrix = arrayfun(@(c) arrayfun(@(J) [p1_Split{c}{J};p2_Split{c}{J}],1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          Weave_Basis_Matrix = arrayfun(@(c) arrayfun(@(J) [p1_Weave{c}{J};p2_Weave{c}{J}],1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);  
          PxF2D_Split = arrayfun(@(c) arrayfun(@(J) F_DiscretePossionBinomial_Base2(PxF2D_Split{c}{J},Split_Basis_Matrix{c}{J}),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          PxF2D_Weave = arrayfun(@(c) arrayfun(@(J) F_DiscretePossionBinomial_Base2(PxF2D_Weave{c}{J},Weave_Basis_Matrix{c}{J}),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);  
          Dim_Split = arrayfun(@(c) cell2mat(arrayfun(@(J) max(size(PxF2D_Split{c}{J}))+1,1:length(Tset),'Un',0)),1:size(EKernel_MultiTissue,1),'Un',0);
          Dim_Weave = arrayfun(@(c) cell2mat(arrayfun(@(J) max(size(PxF2D_Weave{c}{J}))+1,1:length(Tset),'Un',0)),1:size(EKernel_MultiTissue,1),'Un',0);
          DimL = max(cell2mat(arrayfun(@(c) max([Dim_Split{c} Dim_Weave{c}]),1:size(EKernel_MultiTissue,1),'Un',0)));
          PxF2D_Split = arrayfun(@(c) arrayfun(@(J) padarray(PxF2D_Split{c}{J},[max(size(PxF2D_Split{c}{J}))+1-size(PxF2D_Split{c}{J},1) max(size(PxF2D_Split{c}{J}))+1-size(PxF2D_Split{c}{J},2)],0,'post'),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          PxF2D_Weave = arrayfun(@(c) arrayfun(@(J) padarray(PxF2D_Weave{c}{J},[max(size(PxF2D_Weave{c}{J}))+1-size(PxF2D_Weave{c}{J},1) max(size(PxF2D_Weave{c}{J}))+1-size(PxF2D_Weave{c}{J},2)],0,'post'),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          PxF2D_Expr_Split = arrayfun(@(c) arrayfun(@(J) ExpressionMatrix(Tset(J),c)*PxF2D_Split{c}{J},1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          PxF2D_Expr_Weave = arrayfun(@(c) arrayfun(@(J) ExpressionMatrix(Tset(J),c)*PxF2D_Weave{c}{J},1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          Sol(Ti,m).PxF2D_Split = arrayfun(@(c) arrayfun(@(J) padarray(PxF2D_Split{c}{J},[DimL-size(PxF2D_Split{c}{J},1) DimL-size(PxF2D_Split{c}{J},2)],0,'post'),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          Sol(Ti,m).PxF2D_Expr_Split = arrayfun(@(c) arrayfun(@(J) padarray(PxF2D_Expr_Split{c}{J},[DimL-size(PxF2D_Expr_Split{c}{J},1) DimL-size(PxF2D_Expr_Split{c}{J},2)],0,'post'),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          Sol(Ti,m).PxF2D_Weave = arrayfun(@(c) arrayfun(@(J) padarray(PxF2D_Weave{c}{J},[DimL-size(PxF2D_Weave{c}{J},1) DimL-size(PxF2D_Weave{c}{J},2)],0,'post'),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          Sol(Ti,m).PxF2D_Expr_Weave = arrayfun(@(c) arrayfun(@(J) padarray(PxF2D_Expr_Weave{c}{J},[DimL-size(PxF2D_Expr_Weave{c}{J},1) DimL-size(PxF2D_Expr_Weave{c}{J},2)],0,'post'),1:length(Tset),'Un',0),1:size(EKernel_MultiTissue,1),'Un',0);
          %check ability to distiguish isoforms (P(x) similarity)
          if (length(Isoforms)>1)
              for u = 1:length(Isoforms)
                 for v = 1:length(Isoforms)
                     for c = 1:size(EKernel_MultiTissue,1)
                     Sol(Ti,m).D_Isoforms{c}(u,v) = 1/2*trapz(0:Smax,abs(Sol(Ti,m).PxF{c}{Corr_ONRNA_IsoAgnostic(u)}-Sol(Ti,m).PxF{c}{Corr_ONRNA_IsoAgnostic(v)}));
                     Sol(Ti,m).D_Isoforms_Split{c}(u,v) = 1/2*trapz(0:DimL-1,trapz(0:DimL-1,abs(Sol(Ti,m).PxF2D_Split{c}{Corr_ONRNA_IsoAgnostic(u)}-Sol(Ti,m).PxF2D_Split{c}{Corr_ONRNA_IsoAgnostic(v)})));
                     Sol(Ti,m).D_Isoforms_Weave{c}(u,v) = 1/2*trapz(0:DimL-1,trapz(0:DimL-1,abs(Sol(Ti,m).PxF2D_Weave{c}{Corr_ONRNA_IsoAgnostic(u)}-Sol(Ti,m).PxF2D_Weave{c}{Corr_ONRNA_IsoAgnostic(v)})));
                     end
                 end
              end
          else
              for c = 1:size(EKernel_MultiTissue,1)
              Sol(Ti,m).D_Isoforms{c} = NaN;
              Sol(Ti,m).D_Isoforms_Split{c} = NaN;
              Sol(Ti,m).D_Isoforms_Weave{c} = NaN;
              end
          end
          Sol(Ti,m).MeanIsoIdentifiability = cell2mat(arrayfun(@(c) mean(Sol(Ti,m).D_Isoforms{c}(triu(true(size(Sol(Ti,m).D_Isoforms{c})),1))),1:size(EKernel_MultiTissue,1),'Un',0));
          Sol(Ti,m).SplitMeanIsoIdentifiability = cell2mat(arrayfun(@(c) mean(Sol(Ti,m).D_Isoforms_Split{c}(triu(true(size(Sol(Ti,m).D_Isoforms_Split{c})),1))),1:size(EKernel_MultiTissue,1),'Un',0));
          Sol(Ti,m).WeaveMeanIsoIdentifiability = cell2mat(arrayfun(@(c) mean(Sol(Ti,m).D_Isoforms_Weave{c}(triu(true(size(Sol(Ti,m).D_Isoforms_Weave{c})),1))),1:size(EKernel_MultiTissue,1),'Un',0));  
          %Gibbs Sampling
          W_Norm = squeeze(sum(Sol(Ti,m).PxF_ExprMatrix(:,2:end,:),1));%count dist
          W_ONDNA = squeeze(sum(Sol(Ti,m).PxF_ExprMatrix(Corr_ONDNA,2:end,:),1));%count dist
          W_ONRNA_IsoSpecific = squeeze(sum(Sol(Ti,m).PxF_ExprMatrix(Corr_ONRNA_IsoSpecific,2:end,:),1));%count dist
          W_ONRNA_IsoAgnostic = squeeze(sum(Sol(Ti,m).PxF_ExprMatrix(Corr_ONRNA_IsoAgnostic,2:end,:),1));%count dist
          W_OFFDNA = squeeze(sum(Sol(Ti,m).PxF_ExprMatrix(Corr_OFFDNA,2:end,:),1));%count dist
          W_OFFRNA_withIsos = squeeze(sum(Sol(Ti,m).PxF_ExprMatrix(Corr_OFFRNA_withIsos,2:end,:),1));%count dist
          W_OFFRNA_minusIso = squeeze(sum(Sol(Ti,m).PxF_ExprMatrix(Corr_OFFRNA_minusIso,2:end,:),1));%count dist
          W_OFFRNA_otherIso = squeeze(sum(Sol(Ti,m).PxF_ExprMatrix(Corr_OFFRNA_otherIso,2:end,:),1));%count dist    
          Spot_Detect_Max = sum(W_Norm,1);
          Spot_Detect_ON_IsoSpecific = sum(W_ONDNA,1) + sum(W_ONRNA_IsoSpecific,1);
          Spot_Detect_ON_IsoAgnostic = sum(W_ONDNA,1) + sum(W_ONRNA_IsoAgnostic,1);
          Spot_Detect_OFF_withIsos = sum(W_OFFDNA,1) + sum(W_OFFRNA_withIsos,1);
          Spot_Detect_OFF_minusIso = sum(W_OFFDNA,1) + sum(W_OFFRNA_minusIso,1);
          Spot_Detect_OFF_otherIso = sum(W_OFFDNA,1) + sum(W_OFFRNA_otherIso,1);
          Probe_Detect_Max = cell2mat(arrayfun(@(c) dot(1:size(W_Norm,1),W_Norm(:,c)),1:size(EKernel_MultiTissue,1),'Un',0));
          Probe_Detect_ON_IsoSpecific = cell2mat(arrayfun(@(c) dot(1:size(W_ONDNA,1),W_ONDNA(:,c)) + dot(1:size(W_ONRNA_IsoSpecific,1),W_ONRNA_IsoSpecific(:,c)),1:size(EKernel_MultiTissue,1),'Un',0));
          Probe_Detect_ON_IsoAgnostic = cell2mat(arrayfun(@(c) dot(1:size(W_ONDNA,1),W_ONDNA(:,c)) + dot(1:size(W_ONRNA_IsoAgnostic,1),W_ONRNA_IsoAgnostic(:,c)),1:size(EKernel_MultiTissue,1),'Un',0));
          Probe_Detect_OFF_withIsos = cell2mat(arrayfun(@(c) dot(1:size(W_OFFDNA,1),W_OFFDNA(:,c)) + dot(1:size(W_OFFRNA_withIsos,1),W_OFFRNA_withIsos(:,c)),1:size(EKernel_MultiTissue,1),'Un',0));
          Probe_Detect_OFF_minusIso = cell2mat(arrayfun(@(c) dot(1:size(W_OFFDNA,1),W_OFFDNA(:,c)) + dot(1:size(W_OFFRNA_minusIso,1),W_OFFRNA_minusIso(:,c)),1:size(EKernel_MultiTissue,1),'Un',0));
          Probe_Detect_OFF_otherIso = cell2mat(arrayfun(@(c) dot(1:size(W_OFFRNA_otherIso,1),W_OFFRNA_otherIso(:,c)),1:size(EKernel_MultiTissue,1),'Un',0));  
          Sol(Ti,m).Spot.Detection.Max = Spot_Detect_Max;
          Sol(Ti,m).Spot.Detection.ON_IsoSpecific = Spot_Detect_ON_IsoSpecific;
          Sol(Ti,m).Spot.Detection.ON_IsoAgnostic = Spot_Detect_ON_IsoAgnostic;
          Sol(Ti,m).Spot.Detection.OFF_withIsos = Spot_Detect_OFF_withIsos;
          Sol(Ti,m).Spot.Detection.OFF_minusIso = Spot_Detect_OFF_minusIso;
          Sol(Ti,m).Spot.Detection.OFF_otherIso = Spot_Detect_OFF_otherIso;
          Sol(Ti,m).Spot.Percentage.ON_IsoSpecific = Spot_Detect_ON_IsoSpecific./Spot_Detect_Max;
          Sol(Ti,m).Spot.Percentage.ON_IsoAgnostic = Spot_Detect_ON_IsoAgnostic./Spot_Detect_Max;
          Sol(Ti,m).Spot.Percentage.OFF_withIsos = Spot_Detect_OFF_withIsos./Spot_Detect_Max;
          Sol(Ti,m).Spot.Percentage.OFF_minusIso = Spot_Detect_OFF_minusIso./Spot_Detect_Max;
          Sol(Ti,m).Spot.Percentage.OFF_otherIso = Spot_Detect_OFF_otherIso./Spot_Detect_Max;
          Sol(Ti,m).Probe.Detection.Max = Probe_Detect_Max;
          Sol(Ti,m).Probe.Detection.ON_IsoSpecific = Probe_Detect_ON_IsoSpecific;
          Sol(Ti,m).Probe.Detection.ON_IsoAgnostic = Probe_Detect_ON_IsoAgnostic;
          Sol(Ti,m).Probe.Detection.OFF_withIsos = Probe_Detect_OFF_withIsos;
          Sol(Ti,m).Probe.Detection.OFF_minusIso = Probe_Detect_OFF_minusIso;
          Sol(Ti,m).Probe.Detection.OFF_otherIso = Probe_Detect_OFF_otherIso;
          Sol(Ti,m).Probe.Percentage.ON_IsoSpecific = Probe_Detect_ON_IsoSpecific./Probe_Detect_Max;
          Sol(Ti,m).Probe.Percentage.ON_IsoAgnostic = Probe_Detect_ON_IsoAgnostic./Probe_Detect_Max;
          Sol(Ti,m).Probe.Percentage.OFF_withIsos = Probe_Detect_OFF_withIsos./Probe_Detect_Max;
          Sol(Ti,m).Probe.Percentage.OFF_minusIso = Probe_Detect_OFF_withIsos./Probe_Detect_Max;
          Sol(Ti,m).Probe.Percentage.OFF_otherIso = Probe_Detect_OFF_otherIso./Probe_Detect_Max;
          for groups = 1:5
              switch groups 
                 case 1%RNA_ON_Target
                    Gx = arrayfun(@(x) find(Tset==onRNA(x)),1:length(onRNA),'Un',0);
                    withZero = 1;
                    case 2%RNA_ON_Target_IsoAgnostic
                    Gx = arrayfun(@(x) find(Tset==Isoforms(x)),1:length(Isoforms),'Un',0);
                    withZero = 1;
                    case 3%RNA_OFF_Target_withIso
                    Gx = arrayfun(@(x) find(Tset==offRNA(x)),1:length(offRNA),'Un',0);
                    withZero = 0;
                    case 4%RNA_OFF_Target_withoutIso
                    Gx = arrayfun(@(x) find(Tset==offRNA_minusIso(x)),1:length(offRNA_minusIso),'Un',0);
                    withZero = 0;
                    case 5%RNA_ON_EachIso
                    Gx = arrayfun(@(x) find(Tset==Isoforms(x)),1:length(Isoforms),'Un',0);
                    withZero = 1;                    
                    case 6%DNA_ON_Target
                    Gx = arrayfun(@(x) find(Tset==onDNA(x)),1:length(onDNA),'Un',0);
                    withZero = 1;
                    case 7%DNA_OFF_Target
                    Gx = arrayfun(@(x) find(Tset==offDNA(x)),1:length(offDNA),'Un',0);
                    withZero = 0;
              end
              Gx = [Gx{:}];%ON-DNA, ON-RNA, ON-IsoSpec, ON-IsoAgnostic, OFF-RNAwoIso, OFF-RNAwIso, OFF-DNA    
              if (groups~=5)
                 Sol(Ti,m).CountY(groups,1:Smax+1,1:size(EKernel_MultiTissue,1)) = CATnWrapper(arrayfun(@(c)PxF_Count_G(Gx,squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c)))',1:size(EKernel_MultiTissue,1),'Un',0),2);
                 Sol(Ti,m).CountCY(groups,1:Smax+1,1:size(EKernel_MultiTissue,1)) = CATnWrapper(arrayfun(@(c)cumsum(PxF_Count_G(Gx,squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c))))',1:size(EKernel_MultiTissue,1),'Un',0),2);
                 Sol(Ti,m).DistY(groups,1:Smax+1,1:size(EKernel_MultiTissue,1)) = CATnWrapper(arrayfun(@(c)PxF_Dist_G(Gx,squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c)))',1:size(EKernel_MultiTissue,1),'Un',0),2);
                 Sol(Ti,m).DistCY(groups,1:Smax+1,1:size(EKernel_MultiTissue,1)) = CATnWrapper(arrayfun(@(c)cumsum(PxF_Dist_G(Gx,squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c))))',1:size(EKernel_MultiTissue,1),'Un',0),2);
                 Sol(Ti,m).DistX(groups,1:Smax+1,1:size(EKernel_MultiTissue,1)) = repmat(0:Smax,[size(EKernel_MultiTissue,1) 1])';
                 Sol(Ti,m).Mean(groups,1:size(EKernel_MultiTissue,1)) = cell2mat(arrayfun(@(c)PxF_Mean_G(Gx,squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c))),1:size(EKernel_MultiTissue,1),'Un',0));
                 Sol(Ti,m).Fano(groups,1:size(EKernel_MultiTissue,1)) = cell2mat(arrayfun(@(c)PxF_Fano_G(Gx,squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c))),1:size(EKernel_MultiTissue,1),'Un',0));
                 Sol(Ti,m).Skew(groups,1:size(EKernel_MultiTissue,1)) = cell2mat(arrayfun(@(c)PxF_Moment_G(Gx,squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c)),3),1:size(EKernel_MultiTissue,1),'Un',0));
                 Sol(Ti,m).Kurt(groups,1:size(EKernel_MultiTissue,1)) = cell2mat(arrayfun(@(c)PxF_Moment_G(Gx,squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c)),4),1:size(EKernel_MultiTissue,1),'Un',0));        
             else
                 for vs = 1:length(Gx)
                     Sol(Ti,m).CountY(groups-1+vs,1:Smax+1,1:size(EKernel_MultiTissue,1)) = CATnWrapper(arrayfun(@(c) PxF_Count_G(Gx(vs),squeeze(Sol(Ti,m,c).PxF_ExprMatrix(:,:,c)))',1:size(EKernel_MultiTissue,1),'Un',0),2);
                     Sol(Ti,m).CountCY(groups-1+vs,1:Smax+1,1:size(EKernel_MultiTissue,1)) = CATnWrapper(arrayfun(@(c) cumsum(PxF_Count_G(Gx(vs),squeeze(Sol(Ti,m,c).PxF_ExprMatrix(:,:,c))))',1:size(EKernel_MultiTissue,1),'Un',0),2);
                     Sol(Ti,m).DistY(groups-1+vs,1:Smax+1,1:size(EKernel_MultiTissue,1)) = CATnWrapper(arrayfun(@(c) PxF_Dist_G(Gx(vs),squeeze(Sol(Ti,m,c).PxF_ExprMatrix))',1:size(EKernel_MultiTissue,1),'Un',0),2);
                     Sol(Ti,m).DistCY(groups-1+vs,1:Smax+1,1:size(EKernel_MultiTissue,1)) = CATnWrapper(arrayfun(@(c) cumsum(PxF_Dist_G(Gx(vs),squeeze(Sol(Ti,m,c).PxF_ExprMatrix(:,:,c))))',1:size(EKernel_MultiTissue,1),'Un',0),2);
                     Sol(Ti,m).DistX(groups-1+vs,1:Smax+1,1:size(EKernel_MultiTissue,1)) = repmat(0:Smax,[size(EKernel_MultiTissue,1) 1])';
                     Sol(Ti,m).Mean(groups-1+vs,1:size(EKernel_MultiTissue,1)) = cell2mat(arrayfun(@(c) PxF_Mean_G(Gx(vs),squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c))),1:size(EKernel_MultiTissue,1),'Un',0));   
                     Sol(Ti,m).Fano(groups-1+vs,1:size(EKernel_MultiTissue,1)) = cell2mat(arrayfun(@(c) PxF_Fano_G(Gx(vs),squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c))),1:size(EKernel_MultiTissue,1),'Un',0));   
                     Sol(Ti,m).Skew(groups-1+vs,1:size(EKernel_MultiTissue,1)) = cell2mat(arrayfun(@(c) PxF_Moment_G(Gx(vs),squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c)),3),1:size(EKernel_MultiTissue,1),'Un',0));   
                     Sol(Ti,m).Kurt(groups-1+vs,1:size(EKernel_MultiTissue,1)) = cell2mat(arrayfun(@(c) PxF_Moment_G(Gx(vs),squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c)),4),1:size(EKernel_MultiTissue,1),'Un',0));   
                 end
             end  
             if (groups~=5)
                 Mat_Split = arrayfun(@(c) sum(CATnWrapper(Sol(Ti,m).PxF2D_Expr_Split{c}{Gx},3),3),1:size(EKernel_MultiTissue,1),'Un',0); 
                 Mat_Weave = arrayfun(@(c) sum(CATnWrapper(Sol(Ti,m).PxF2D_Expr_Weave{c}{Gx},3),3),1:size(EKernel_MultiTissue,1),'Un',0); 
                 Sol(Ti,m) = SubStats_3D(groups,Mat_Split,Sol(Ti,m),1,withZero,maxProbes);
                 Sol(Ti,m) = SubStats_3D(groups,Mat_Weave,Sol(Ti,m),2,withZero,maxProbes);
             else
                for vs = 1:length(Gx)
                     Mat_Split2 = arrayfun(@(c) Sol(Ti,m).PxF2D_Expr_Split{c}{Gx(vs)},1:size(EKernel_MultiTissue,1),'Un',0);
                     Mat_Weave2 = arrayfun(@(c) Sol(Ti,m).PxF2D_Expr_Weave{c}{Gx(vs)},1:size(EKernel_MultiTissue,1),'Un',0);
                     Sol(Ti,m,c) = SubStats_3D(groups-1+vs,Mat_Split2,Sol(Ti,m),1,withZero,maxProbes);
                     Sol(Ti,m,c) = SubStats_3D(groups-1+vs,Mat_Weave2,Sol(Ti,m),2,withZero,maxProbes);
                end        
             end
          end   
       end       
    end
elseif (solution_case == 3)
    PC0 = repmat(ProbeConcSeries,[size(ExpressionMatrix,2) length(Pset)]);   
    CPF0 = GuessConc*ones(size(EKernel_MultiProbeConc,1),length(Pset)); 
    opts_MultiTissue = @(T,M) odeset('mass', M_MultiProbeConc, 'RelTol', 10^(-6),'AbsTol', 10^(-6), 'InitialSlope', F_MultiProbeConc(0,reshape(CPF0,[numel(CPF0) 1]),PC0,T,M));   
    y_MultiProbeConc = @(T,M) reshape(ode45(@(t,y) F_MultiTissue(t,y,PC0,T,M),tspan,reshape(CPF0,[numel(CPF0) 1]),opts_MultiTissue(T,M)).y(:,end),size(CPF0)); 
    y_Bound_MultiProbeConc = @(T,M) p_TargetSites_Bound_MultiProbeConc_Func(y_MultiProbeConc(T,M),T,M);
    y_Unbound_MultiProbeConc = @(T,M) p_TargetSites_Unbound_MultiProbeConc_Func(y_MultiProbeConc(T,M),T,M);
    y_Bound_MultiProbeConc_Weave = @(T,M,k) p_TargetSites_Bound_MultiProbeConc_Weave_Func{k}(y_MultiProbeConc(T,M),T,M);
    y_Unbound_MultiProbeConc_Weave = @(T,M,k) p_TargetSites_Unbound_MultiProbeConc_Weave_Func{k}(y_MultiProbeConc(T,M),T,M);
    y_Bound_MultiProbeConc_Split = @(T,M,k) p_TargetSites_Bound_MultiProbeConc_Split_Func{k}(y_MultiProbeConc(T,M),T,M);
    y_Unbound_MultiProbeConc_Split = @(T,M,k) p_TargetSites_Unbound_MultiProbeConc_Split_Func{k}(y_MultiProbeConc(T,M),T,M);
    for Ti = T_low_celsius:T_high_celsius
       for m = 1:N_methods
          for k = 1:N_Channels
             Sol(Ti,m).Bound_Weave{k} = y_Bound_MultiProbeConc_Weave(Ti+273.15,m,k);
             Sol(Ti,m).Unbound_Weave{k} = y_Unbound_MultiProbeConc_Weave(Ti+273.15,m,k);
             Sol(Ti,m).Bound_Split{k} = y_Bound_MultiProbeConc_Split(Ti+273.15,m,k);
             Sol(Ti,m).Unbound_Split{k} = y_Unbound_MultiProbeConc_Split(Ti+273.15,m,k); 
             PxFD_Weave{k} = arrayfun(@(c) arrayfun(@(J) F_DiscretePossionBinomial(P_1FLAPs(2)*Sol(Ti,m).Bound_Weave{k}(c,Tset(J),vsi_Weave_Uni{k}(J))),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
             PxFD_Split{k} = arrayfun(@(c) arrayfun(@(J) F_DiscretePossionBinomial(P_1FLAPs(2)*Sol(Ti,m).Bound_Split{k}(c,Tset(J),vsi_Split_Uni{k}(J))),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          end
          Sol(Ti,m).Bound = y_Bound_MultiProbeConc(Ti+273.15,m);
          Sol(Ti,m).Unbound = y_Unbound_MultiProbeConc(Ti+273.15,m); 
          Sol(Ti,m).MomGen1 = y_MomGen1(Sol(Ti,m).Bound);
          Sol(Ti,m).MomGen2 = y_MomGen2(Sol(Ti,m).Bound,Sol(Ti,m).Unbound);
          Sol(Ti,m).Zero = y_Zero_3D(Sol(Ti,m).Unbound);
          Sol(Ti,m).Bound_ON = y_Bound_ON_3D(ON_IDs,Sol(Ti,m).Bound);
          Sol(Ti,m).Bound_OFF = y_Bound_OFF_3D(OFF_IDs,Sol(Ti,m).Bound);
          Sol(Ti,m).Unbound_ON = y_Unbound_ON_3D(ON_IDs,Sol(Ti,m).Unbound);
          Sol(Ti,m).Unbound_OFF = y_Unbound_OFF_3D(OFF_IDs,Sol(Ti,m).Unbound);
          Sol(Ti,m).Mean = mean_probes_Target_3D(Sol(Ti,m).Bound);
          Sol(Ti,m).Var  = var_probes_Target_3D(Sol(Ti,m).Bound,Sol(Ti,m).Unbound);
          Sol(Ti,m).Fano = fano_probes_Target_3D(Sol(Ti,m).Bound,Sol(Ti,m).Unbound);
          Sol(Ti,m).Skew = skew_probes_Target_3D(Sol(Ti,m).Bound,Sol(Ti,m).Unbound,Sol(Ti,m).MomGen1);
          Sol(Ti,m).Kurt = kurt_probes_Target_3D(Sol(Ti,m).Bound,Sol(Ti,m).Unbound,Sol(Ti,m).MomGen2);
          Sol(Ti,m).Mean0_ON = Mean_G0_3D_MultiProbeConc(ON_IDs,Sol(Ti,m).Mean); 
          Sol(Ti,m).STD0_ON = Var_G0_3D_MultiProbeConc(ON_IDs,Sol(Ti,m).Var);
          Sol(Ti,m).Mean_ON = Mean_G_3D_MultiProbeConc(ON_IDs,Sol(Ti,m).Mean,Sol(Ti,m).Zero); 
          Sol(Ti,m).STD_ON = Var_G_3D_MultiProbeConc(ON_IDs,Sol(Ti,m).Var,Sol(Ti,m).Zero); 
          Sol(Ti,m).sConc_ON = sConc_G_3D_MultiProbeConc(ON_IDs,Sol(Ti,m).Zero);
          Sol(Ti,m).pConc_ON = pConc_G_3D_MultiProbeConc(ON_IDs,Sol(Ti,m).Mean); 
          Sol(Ti,m).sFrac_ON = sFrac_G_3D_MultiProbeConc(ON_IDs,Sol(Ti,m).Zero); 
          Sol(Ti,m).pFrac_ON = pFrac_G_3D_MultiProbeConc(ON_IDs,Sol(Ti,m).Mean); 
          Sol(Ti,m).Mean0_OFF =  Mean_G0_3D_MultiProbeConc(OFF_IDs,Sol(Ti,m).Mean); 
          Sol(Ti,m).STD0_OFF = Var_G0_3D_MultiProbeConc(OFF_IDs,Sol(Ti,m).Var);
          Sol(Ti,m).Mean_OFF = Mean_G_3D_MultiProbeConc(OFF_IDs,Sol(Ti,m).Mean,Sol(Ti,m).Zero); 
          Sol(Ti,m).STD_OFF = Var_G_3D_MultiProbeConc(OFF_IDs,Sol(Ti,m).Var,Sol(Ti,m).Zero); 
          Sol(Ti,m).sConc_OFF = sConc_G_3D_MultiProbeConc(OFF_IDs,Sol(Ti,m).Zero);
          Sol(Ti,m).pConc_OFF = pConc_G_3D_MultiProbeConc(OFF_IDs,Sol(Ti,m).Mean); 
          Sol(Ti,m).sFrac_OFF = sFrac_G_3D_MultiProbeConc(OFF_IDs,Sol(Ti,m).Zero); 
          Sol(Ti,m).pFrac_OFF = pFrac_G_3D_MultiProbeConc(OFF_IDs,Sol(Ti,m).Mean); 
          Sol(Ti,m).PxF = arrayfun(@(c) arrayfun(@(J) F_DiscretePossionBinomial(P_1FLAPs(2)*Sol(Ti,m).Bound(c,Tset(J),find(sum(squeeze(DoesProbeBindSite(Pset,Tset(J),:)),1)>0))),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);       
          Sol(Ti,m).PxF_Matrix = CATnWrapper(arrayfun(@(c) vertcat(Sol(Ti,m).PxF{c}{:}),1:size(EKernel_MultiProbeConc,1),'Un',0),3); 
          Sol(Ti,m).PxF_ExprMatrix = CATnWrapper(arrayfun(@(c) ExpressionMatrix(Tset,c).*squeeze(Sol(Ti,m).PxF_Matrix(:,:,c)),1:size(EKernel_MultiProbeConc,1),'Un',0),3);    
          PxF2D_Weave = arrayfun(@(c) arrayfun(@(J) PxFD_Weave{2}{c}{J}.*PxFD_Weave{1}{c}{J}',1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          PxF2D_Split = arrayfun(@(c) arrayfun(@(J) PxFD_Split{2}{c}{J}.*PxFD_Split{1}{c}{J}',1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          temp_Weave2 = arrayfun(@(c) arrayfun(@(J) PxF2D_Weave{c}{J}./sum(sum(PxF2D_Weave{c}{J},1),2),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          temp_Split2 = arrayfun(@(c) arrayfun(@(J) PxF2D_Split{c}{J}./sum(sum(PxF2D_Split{c}{J},1),2),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          PxF2D_Weave = temp_Weave2;
          PxF2D_Split = temp_Split2;  
          p1_Split = arrayfun(@(c) arrayfun(@(J) Sol(Ti,m).Bound_Split{1}(c,Tset(J),overlap_Split(J)),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          p2_Split = arrayfun(@(c) arrayfun(@(J) Sol(Ti,m).Bound_Split{2}(c,Tset(J),overlap_Split(J)),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          p1_Weave = arrayfun(@(c) arrayfun(@(J) Sol(Ti,m).Bound_Weave{1}(c,Tset(J),overlap_Weave(J)),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          p2_Weave = arrayfun(@(c) arrayfun(@(J) Sol(Ti,m).Bound_Weave{2}(c,Tset(J),overlap_Weave(J)),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          Split_Basis_Matrix = arrayfun(@(c) arrayfun(@(J) [p1_Split{c}{J};p2_Split{c}{J}],1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          Weave_Basis_Matrix = arrayfun(@(c) arrayfun(@(J) [p1_Weave{c}{J};p2_Weave{c}{J}],1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);  
          PxF2D_Split = arrayfun(@(c) arrayfun(@(J) F_DiscretePossionBinomial_Base2(PxF2D_Split{c}{J},Split_Basis_Matrix{c}{J}),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          PxF2D_Weave = arrayfun(@(c) arrayfun(@(J) F_DiscretePossionBinomial_Base2(PxF2D_Weave{c}{J},Weave_Basis_Matrix{c}{J}),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);  
          Dim_Split = arrayfun(@(c) cell2mat(arrayfun(@(J) max(size(PxF2D_Split{c}{J}))+1,1:length(Tset),'Un',0)),1:size(EKernel_MultiProbeConc,1),'Un',0);
          Dim_Weave = arrayfun(@(c) cell2mat(arrayfun(@(J) max(size(PxF2D_Weave{c}{J}))+1,1:length(Tset),'Un',0)),1:size(EKernel_MultiProbeConc,1),'Un',0);
          DimL = max(cell2mat(arrayfun(@(c) max([Dim_Split{c} Dim_Weave{c}]),1:size(EKernel_MultiProbeConc,1),'Un',0)));
          PxF2D_Split = arrayfun(@(c) arrayfun(@(J) padarray(PxF2D_Split{c}{J},[max(size(PxF2D_Split{c}{J}))+1-size(PxF2D_Split{c}{J},1) max(size(PxF2D_Split{c}{J}))+1-size(PxF2D_Split{c}{J},2)],0,'post'),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          PxF2D_Weave = arrayfun(@(c) arrayfun(@(J) padarray(PxF2D_Weave{c}{J},[max(size(PxF2D_Weave{c}{J}))+1-size(PxF2D_Weave{c}{J},1) max(size(PxF2D_Weave{c}{J}))+1-size(PxF2D_Weave{c}{J},2)],0,'post'),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          PxF2D_Expr_Split = arrayfun(@(c) arrayfun(@(J) ExpressionMatrix(Tset(J),c)*PxF2D_Split{c}{J},1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          PxF2D_Expr_Weave = arrayfun(@(c) arrayfun(@(J) ExpressionMatrix(Tset(J),c)*PxF2D_Weave{c}{J},1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          Sol(Ti,m).PxF2D_Split = arrayfun(@(c) arrayfun(@(J) padarray(PxF2D_Split{c}{J},[DimL-size(PxF2D_Split{c}{J},1) DimL-size(PxF2D_Split{c}{J},2)],0,'post'),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          Sol(Ti,m).PxF2D_Expr_Split = arrayfun(@(c) arrayfun(@(J) padarray(PxF2D_Expr_Split{c}{J},[DimL-size(PxF2D_Expr_Split{c}{J},1) DimL-size(PxF2D_Expr_Split{c}{J},2)],0,'post'),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          Sol(Ti,m).PxF2D_Weave = arrayfun(@(c) arrayfun(@(J) padarray(PxF2D_Weave{c}{J},[DimL-size(PxF2D_Weave{c}{J},1) DimL-size(PxF2D_Weave{c}{J},2)],0,'post'),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          Sol(Ti,m).PxF2D_Expr_Weave = arrayfun(@(c) arrayfun(@(J) padarray(PxF2D_Expr_Weave{c}{J},[DimL-size(PxF2D_Expr_Weave{c}{J},1) DimL-size(PxF2D_Expr_Weave{c}{J},2)],0,'post'),1:length(Tset),'Un',0),1:size(EKernel_MultiProbeConc,1),'Un',0);
          %check ability to distiguish isoforms (P(x) similarity)
          if (length(Isoforms)>1)
              for u = 1:length(Isoforms)
                 for v = 1:length(Isoforms)
                     for c = 1:size(EKernel_MultiProbeConc,1)
                     Sol(Ti,m).D_Isoforms{c}(u,v) = 1/2*trapz(0:Smax,abs(Sol(Ti,m).PxF{c}{Corr_ONRNA_IsoAgnostic(u)}-Sol(Ti,m).PxF{c}{Corr_ONRNA_IsoAgnostic(v)}));
                     Sol(Ti,m).D_Isoforms_Split{c}(u,v) = 1/2*trapz(0:DimL-1,trapz(0:DimL-1,abs(Sol(Ti,m).PxF2D_Split{c}{Corr_ONRNA_IsoAgnostic(u)}-Sol(Ti,m).PxF2D_Split{c}{Corr_ONRNA_IsoAgnostic(v)})));
                     Sol(Ti,m).D_Isoforms_Weave{c}(u,v) = 1/2*trapz(0:DimL-1,trapz(0:DimL-1,abs(Sol(Ti,m).PxF2D_Weave{c}{Corr_ONRNA_IsoAgnostic(u)}-Sol(Ti,m).PxF2D_Weave{c}{Corr_ONRNA_IsoAgnostic(v)})));
                     end
                 end
              end
          else
              for c = 1:size(EKernel_MultiProbeConc,1)
              Sol(Ti,m).D_Isoforms{c} = NaN;
              Sol(Ti,m).D_Isoforms_Split{c} = NaN;
              Sol(Ti,m).D_Isoforms_Weave{c} = NaN;
              end
          end
          Sol(Ti,m).MeanIsoIdentifiability = cell2mat(arrayfun(@(c) mean(Sol(Ti,m).D_Isoforms{c}(triu(true(size(Sol(Ti,m).D_Isoforms{c})),1))),1:size(EKernel_MultiProbeConc,1),'Un',0));
          Sol(Ti,m).SplitMeanIsoIdentifiability = cell2mat(arrayfun(@(c) mean(Sol(Ti,m).D_Isoforms_Split{c}(triu(true(size(Sol(Ti,m).D_Isoforms_Split{c})),1))),1:size(EKernel_MultiProbeConc,1),'Un',0));
          Sol(Ti,m).WeaveMeanIsoIdentifiability = cell2mat(arrayfun(@(c) mean(Sol(Ti,m).D_Isoforms_Weave{c}(triu(true(size(Sol(Ti,m).D_Isoforms_Weave{c})),1))),1:size(EKernel_MultiProbeConc,1),'Un',0));  
          %Gibbs Sampling
          W_Norm = squeeze(sum(Sol(Ti,m).PxF_ExprMatrix(:,2:end,:),1));%count dist
          W_ONDNA = squeeze(sum(Sol(Ti,m).PxF_ExprMatrix(Corr_ONDNA,2:end,:),1));%count dist
          W_ONRNA_IsoSpecific = squeeze(sum(Sol(Ti,m).PxF_ExprMatrix(Corr_ONRNA_IsoSpecific,2:end,:),1));%count dist
          W_ONRNA_IsoAgnostic = squeeze(sum(Sol(Ti,m).PxF_ExprMatrix(Corr_ONRNA_IsoAgnostic,2:end,:),1));%count dist
          W_OFFDNA = squeeze(sum(Sol(Ti,m).PxF_ExprMatrix(Corr_OFFDNA,2:end,:),1));%count dist
          W_OFFRNA_withIsos = squeeze(sum(Sol(Ti,m).PxF_ExprMatrix(Corr_OFFRNA_withIsos,2:end,:),1));%count dist
          W_OFFRNA_minusIso = squeeze(sum(Sol(Ti,m).PxF_ExprMatrix(Corr_OFFRNA_minusIso,2:end,:),1));%count dist
          W_OFFRNA_otherIso = squeeze(sum(Sol(Ti,m).PxF_ExprMatrix(Corr_OFFRNA_otherIso,2:end,:),1));%count dist    
          Spot_Detect_Max = sum(W_Norm,1);
          Spot_Detect_ON_IsoSpecific = sum(W_ONDNA,1) + sum(W_ONRNA_IsoSpecific,1);
          Spot_Detect_ON_IsoAgnostic = sum(W_ONDNA,1) + sum(W_ONRNA_IsoAgnostic,1);
          Spot_Detect_OFF_withIsos = sum(W_OFFDNA,1) + sum(W_OFFRNA_withIsos,1);
          Spot_Detect_OFF_minusIso = sum(W_OFFDNA,1) + sum(W_OFFRNA_minusIso,1);
          Spot_Detect_OFF_otherIso = sum(W_OFFDNA,1) + sum(W_OFFRNA_otherIso,1);
          Probe_Detect_Max = cell2mat(arrayfun(@(c) dot(1:size(W_Norm,1),W_Norm(:,c)),1:size(EKernel_MultiProbeConc,1),'Un',0));
          Probe_Detect_ON_IsoSpecific = cell2mat(arrayfun(@(c) dot(1:size(W_ONDNA,1),W_ONDNA(:,c)) + dot(1:size(W_ONRNA_IsoSpecific,1),W_ONRNA_IsoSpecific(:,c)),1:size(EKernel_MultiProbeConc,1),'Un',0));
          Probe_Detect_ON_IsoAgnostic = cell2mat(arrayfun(@(c) dot(1:size(W_ONDNA,1),W_ONDNA(:,c)) + dot(1:size(W_ONRNA_IsoAgnostic,1),W_ONRNA_IsoAgnostic(:,c)),1:size(EKernel_MultiProbeConc,1),'Un',0));
          Probe_Detect_OFF_withIsos = cell2mat(arrayfun(@(c) dot(1:size(W_OFFDNA,1),W_OFFDNA(:,c)) + dot(1:size(W_OFFRNA_withIsos,1),W_OFFRNA_withIsos(:,c)),1:size(EKernel_MultiProbeConc,1),'Un',0));
          Probe_Detect_OFF_minusIso = cell2mat(arrayfun(@(c) dot(1:size(W_OFFDNA,1),W_OFFDNA(:,c)) + dot(1:size(W_OFFRNA_minusIso,1),W_OFFRNA_minusIso(:,c)),1:size(EKernel_MultiProbeConc,1),'Un',0));
          Probe_Detect_OFF_otherIso = cell2mat(arrayfun(@(c) dot(1:size(W_OFFRNA_otherIso,1),W_OFFRNA_otherIso(:,c)),1:size(EKernel_MultiProbeConc,1),'Un',0));  
          Sol(Ti,m).Spot.Detection.Max = Spot_Detect_Max;
          Sol(Ti,m).Spot.Detection.ON_IsoSpecific = Spot_Detect_ON_IsoSpecific;
          Sol(Ti,m).Spot.Detection.ON_IsoAgnostic = Spot_Detect_ON_IsoAgnostic;
          Sol(Ti,m).Spot.Detection.OFF_withIsos = Spot_Detect_OFF_withIsos;
          Sol(Ti,m).Spot.Detection.OFF_minusIso = Spot_Detect_OFF_minusIso;
          Sol(Ti,m).Spot.Detection.OFF_otherIso = Spot_Detect_OFF_otherIso;
          Sol(Ti,m).Spot.Percentage.ON_IsoSpecific = Spot_Detect_ON_IsoSpecific./Spot_Detect_Max;
          Sol(Ti,m).Spot.Percentage.ON_IsoAgnostic = Spot_Detect_ON_IsoAgnostic./Spot_Detect_Max;
          Sol(Ti,m).Spot.Percentage.OFF_withIsos = Spot_Detect_OFF_withIsos./Spot_Detect_Max;
          Sol(Ti,m).Spot.Percentage.OFF_minusIso = Spot_Detect_OFF_minusIso./Spot_Detect_Max;
          Sol(Ti,m).Spot.Percentage.OFF_otherIso = Spot_Detect_OFF_otherIso./Spot_Detect_Max;
          Sol(Ti,m).Probe.Detection.Max = Probe_Detect_Max;
          Sol(Ti,m).Probe.Detection.ON_IsoSpecific = Probe_Detect_ON_IsoSpecific;
          Sol(Ti,m).Probe.Detection.ON_IsoAgnostic = Probe_Detect_ON_IsoAgnostic;
          Sol(Ti,m).Probe.Detection.OFF_withIsos = Probe_Detect_OFF_withIsos;
          Sol(Ti,m).Probe.Detection.OFF_minusIso = Probe_Detect_OFF_minusIso;
          Sol(Ti,m).Probe.Detection.OFF_otherIso = Probe_Detect_OFF_otherIso;
          Sol(Ti,m).Probe.Percentage.ON_IsoSpecific = Probe_Detect_ON_IsoSpecific./Probe_Detect_Max;
          Sol(Ti,m).Probe.Percentage.ON_IsoAgnostic = Probe_Detect_ON_IsoAgnostic./Probe_Detect_Max;
          Sol(Ti,m).Probe.Percentage.OFF_withIsos = Probe_Detect_OFF_withIsos./Probe_Detect_Max;
          Sol(Ti,m).Probe.Percentage.OFF_minusIso = Probe_Detect_OFF_withIsos./Probe_Detect_Max;
          Sol(Ti,m).Probe.Percentage.OFF_otherIso = Probe_Detect_OFF_otherIso./Probe_Detect_Max;
          for groups = 1:5
              switch groups 
                 case 1%RNA_ON_Target
                    Gx = arrayfun(@(x) find(Tset==onRNA(x)),1:length(onRNA),'Un',0);
                    withZero = 1;
                    case 2%RNA_ON_Target_IsoAgnostic
                    Gx = arrayfun(@(x) find(Tset==Isoforms(x)),1:length(Isoforms),'Un',0);
                    withZero = 1;
                    case 3%RNA_OFF_Target_withIso
                    Gx = arrayfun(@(x) find(Tset==offRNA(x)),1:length(offRNA),'Un',0);
                    withZero = 0;
                    case 4%RNA_OFF_Target_withoutIso
                    Gx = arrayfun(@(x) find(Tset==offRNA_minusIso(x)),1:length(offRNA_minusIso),'Un',0);
                    withZero = 0;
                    case 5%RNA_ON_EachIso
                    Gx = arrayfun(@(x) find(Tset==Isoforms(x)),1:length(Isoforms),'Un',0);
                    withZero = 1;                    
                    case 6%DNA_ON_Target
                    Gx = arrayfun(@(x) find(Tset==onDNA(x)),1:length(onDNA),'Un',0);
                    withZero = 1;
                    case 7%DNA_OFF_Target
                    Gx = arrayfun(@(x) find(Tset==offDNA(x)),1:length(offDNA),'Un',0);
                    withZero = 0;
              end
              Gx = [Gx{:}];%ON-DNA, ON-RNA, ON-IsoSpec, ON-IsoAgnostic, OFF-RNAwoIso, OFF-RNAwIso, OFF-DNA    
              if (groups~=5)
                 Sol(Ti,m).CountY(groups,1:Smax+1,1:size(EKernel_MultiProbeConc,1)) = CATnWrapper(arrayfun(@(c)PxF_Count_G(Gx,squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c)))',1:size(EKernel_MultiProbeConc,1),'Un',0),2);
                 Sol(Ti,m).CountCY(groups,1:Smax+1,1:size(EKernel_MultiProbeConc,1)) = CATnWrapper(arrayfun(@(c)cumsum(PxF_Count_G(Gx,squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c))))',1:size(EKernel_MultiProbeConc,1),'Un',0),2);
                 Sol(Ti,m).DistY(groups,1:Smax+1,1:size(EKernel_MultiProbeConc,1)) = CATnWrapper(arrayfun(@(c)PxF_Dist_G(Gx,squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c)))',1:size(EKernel_MultiProbeConc,1),'Un',0),2);
                 Sol(Ti,m).DistCY(groups,1:Smax+1,1:size(EKernel_MultiProbeConc,1)) = CATnWrapper(arrayfun(@(c)cumsum(PxF_Dist_G(Gx,squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c))))',1:size(EKernel_MultiProbeConc,1),'Un',0),2);
                 Sol(Ti,m).DistX(groups,1:Smax+1,1:size(EKernel_MultiProbeConc,1)) = repmat(0:Smax,[size(EKernel_MultiProbeConc,1) 1])';
                 Sol(Ti,m).Mean(groups,1:size(EKernel_MultiProbeConc,1)) = cell2mat(arrayfun(@(c)PxF_Mean_G(Gx,squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c))),1:size(EKernel_MultiProbeConc,1),'Un',0));
                 Sol(Ti,m).Fano(groups,1:size(EKernel_MultiProbeConc,1)) = cell2mat(arrayfun(@(c)PxF_Fano_G(Gx,squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c))),1:size(EKernel_MultiProbeConc,1),'Un',0));
                 Sol(Ti,m).Skew(groups,1:size(EKernel_MultiProbeConc,1)) = cell2mat(arrayfun(@(c)PxF_Moment_G(Gx,squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c)),3),1:size(EKernel_MultiProbeConc,1),'Un',0));
                 Sol(Ti,m).Kurt(groups,1:size(EKernel_MultiProbeConc,1)) = cell2mat(arrayfun(@(c)PxF_Moment_G(Gx,squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c)),4),1:size(EKernel_MultiProbeConc,1),'Un',0));        
             else
                 for vs = 1:length(Gx)
                     Sol(Ti,m).CountY(groups-1+vs,1:Smax+1,1:size(EKernel_MultiProbeConc,1)) = CATnWrapper(arrayfun(@(c) PxF_Count_G(Gx(vs),squeeze(Sol(Ti,m,c).PxF_ExprMatrix(:,:,c)))',1:size(EKernel_MultiProbeConc,1),'Un',0),2);
                     Sol(Ti,m).CountCY(groups-1+vs,1:Smax+1,1:size(EKernel_MultiProbeConc,1)) = CATnWrapper(arrayfun(@(c) cumsum(PxF_Count_G(Gx(vs),squeeze(Sol(Ti,m,c).PxF_ExprMatrix(:,:,c))))',1:size(EKernel_MultiProbeConc,1),'Un',0),2);
                     Sol(Ti,m).DistY(groups-1+vs,1:Smax+1,1:size(EKernel_MultiProbeConc,1)) = CATnWrapper(arrayfun(@(c) PxF_Dist_G(Gx(vs),squeeze(Sol(Ti,m,c).PxF_ExprMatrix))',1:size(EKernel_MultiProbeConc,1),'Un',0),2);
                     Sol(Ti,m).DistCY(groups-1+vs,1:Smax+1,1:size(EKernel_MultiProbeConc,1)) = CATnWrapper(arrayfun(@(c) cumsum(PxF_Dist_G(Gx(vs),squeeze(Sol(Ti,m,c).PxF_ExprMatrix(:,:,c))))',1:size(EKernel_MultiProbeConc,1),'Un',0),2);
                     Sol(Ti,m).DistX(groups-1+vs,1:Smax+1,1:size(EKernel_MultiProbeConc,1)) = repmat(0:Smax,[size(EKernel_MultiProbeConc,1) 1])';
                     Sol(Ti,m).Mean(groups-1+vs,1:size(EKernel_MultiProbeConc,1)) = cell2mat(arrayfun(@(c) PxF_Mean_G(Gx(vs),squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c))),1:size(EKernel_MultiProbeConc,1),'Un',0));   
                     Sol(Ti,m).Fano(groups-1+vs,1:size(EKernel_MultiProbeConc,1)) = cell2mat(arrayfun(@(c) PxF_Fano_G(Gx(vs),squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c))),1:size(EKernel_MultiProbeConc,1),'Un',0));   
                     Sol(Ti,m).Skew(groups-1+vs,1:size(EKernel_MultiProbeConc,1)) = cell2mat(arrayfun(@(c) PxF_Moment_G(Gx(vs),squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c)),3),1:size(EKernel_MultiProbeConc,1),'Un',0));   
                     Sol(Ti,m).Kurt(groups-1+vs,1:size(EKernel_MultiProbeConc,1)) = cell2mat(arrayfun(@(c) PxF_Moment_G(Gx(vs),squeeze(Sol(Ti,m).PxF_ExprMatrix(:,:,c)),4),1:size(EKernel_MultiProbeConc,1),'Un',0));   
                 end
             end                                            
             if (groups~=5)
                 Mat_Split = arrayfun(@(c) sum(CATnWrapper(Sol(Ti,m).PxF2D_Expr_Split{c}{Gx},3),3),1:size(EKernel_MultiProbeConc,1),'Un',0); 
                 Mat_Weave = arrayfun(@(c) sum(CATnWrapper(Sol(Ti,m).PxF2D_Expr_Weave{c}{Gx},3),3),1:size(EKernel_MultiProbeConc,1),'Un',0); 
                 Sol(Ti,m) = SubStats_3D(groups,Mat_Split,Sol(Ti,m),1,withZero,maxProbes);
                 Sol(Ti,m) = SubStats_3D(groups,Mat_Weave,Sol(Ti,m),2,withZero,maxProbes);
             else
                for vs = 1:length(Gx)
                     Mat_Split2 = arrayfun(@(c) Sol(Ti,m).PxF2D_Expr_Split{c}{Gx(vs)},1:size(EKernel_MultiProbeConc,1),'Un',0);
                     Mat_Weave2 = arrayfun(@(c) Sol(Ti,m).PxF2D_Expr_Weave{c}{Gx(vs)},1:size(EKernel_MultiProbeConc,1),'Un',0);
                     Sol(Ti,m) = SubStats_3D(groups-1+vs,Mat_Split2,Sol(Ti,m),1,withZero,maxProbes);
                     Sol(Ti,m) = SubStats_3D(groups-1+vs,Mat_Weave2,Sol(Ti,m),2,withZero,maxProbes);
                end        
             end
          end   
       end       
    end
end

end
function PxF_Group = SubStats_2D(vn,P2D,PxF_Group,type,withZero,maxProbes)
switch type
    case 1%Split
        [XX,YY] = meshgrid(0:size(P2D,1)-1,0:size(P2D,2)-1);
    case 2%Weave
        [XX,YY] = meshgrid(0:size(P2D,1)-1,0:size(P2D,2)-1);
end     
P2Dn = P2D/sum(P2D(:));
P2Dn_MargX = sum(P2Dn,1);
P2Dn_MargY = sum(P2Dn,2);
if (withZero)
    MeanX0 = dot(XX(1,:),P2Dn_MargX);
    VarX0 = dot((XX(1,:)-MeanX0).^2,P2Dn_MargX);
    FanoX0 = VarX0/MeanX0;
    MeanY0 = dot(YY(:,1),P2Dn_MargY);
    VarY0 = dot((YY(:,1)-MeanY0).^2,P2Dn_MargY);
    FanoY0 = VarY0/MeanY0;
    CoVar0 = sum(P2Dn.*(XX-MeanX0).*(YY-MeanY0),'all');
    MeanZ0_ver1 = 0.5*(MeanX0+MeanY0);
    MeanZ0_ver2 = sqrt(MeanX0*MeanY0);
    VarZ0_ver1 = 0.5*(VarX0+VarY0);
    VarZ0_ver2 = sqrt(VarX0*VarY0);
    FanoZ0_ver1 = VarZ0_ver1/MeanZ0_ver1;
    FanoZ0_ver2 = VarZ0_ver2/MeanZ0_ver2;
    Pearson_CC0 = CoVar0/sqrt(VarX0*VarY0);
else
    P2Dn_mod = P2Dn;
    P2Dn_mod(1,1) = 0;
    P2Dn_mod = P2Dn_mod/sum(P2Dn_mod(:));
    MeanX1 = dot(XX(1,2:end),P2Dn_MargX(2:end))/sum(P2Dn_MargX(2:end));
    VarX1 = dot((XX(1,2:end)-MeanX1).^2,P2Dn_MargX(2:end))/sum(P2Dn_MargX(2:end));
    FanoX1 = VarX1/MeanX1;
    MeanY1 = dot(YY(2:end,1),P2Dn_MargY(2:end))/sum(P2Dn_MargY(2:end));
    VarY1 = dot((YY(2:end,1)-MeanY1).^2,P2Dn_MargY(2:end))/sum(P2Dn_MargY(2:end));
    FanoY1 = VarY1/MeanY1;
    CoVar1 = sum(P2Dn_mod.*(XX-MeanX1).*(YY-MeanY1),'all');
    MeanZ1_ver1 = 0.5*(MeanX1+MeanY1);
    MeanZ1_ver2 = sqrt(MeanX1*MeanY1);
    VarZ1_ver1 = 0.5*(VarX1+VarY1);
    VarZ1_ver2 = sqrt(VarX1*VarY1);
    FanoZ1_ver1 = VarZ1_ver1/MeanZ1_ver1;
    FanoZ1_ver2 = VarZ1_ver2/MeanZ1_ver2;
    Pearson_CC1 = CoVar1/sqrt(VarX1*VarY1);
end
r1 = XX./YY;
r2 = YY./XX;
r3 = XX./(YY+XX);
r4 = YY./(YY+XX);
R1 = unique(r1);
R2 = unique(r2);
R3 = unique(r3);
R4 = unique(r4);
P_R1 = zeros(1,length(R1));
P_R2 = zeros(1,length(R2));
P_R3 = zeros(1,length(R3));
P_R4 = zeros(1,length(R4));
%make this step faster
for v = 1:length(R1)
   [xloc,yloc] = find(R1(v)==r1); 
   P_R1(v) = sum(P2Dn(sub2ind(size(P2Dn),xloc,yloc)));
end
for v = 1:length(R2)
   [xloc,yloc] = find(R2(v)==r2); 
   P_R2(v) = sum(P2Dn(sub2ind(size(P2Dn),xloc,yloc)));
end
for v = 1:length(R3)
   [xloc,yloc] = find(R3(v)==r3); 
   P_R3(v) = sum(P2Dn(sub2ind(size(P2Dn),xloc,yloc)));
end
for v = 1:length(R4)
   [xloc,yloc] = find(R4(v)==r4); 
   P_R4(v) = sum(P2Dn(sub2ind(size(P2Dn),xloc,yloc)));
end
R1 = R1(P_R1~=0);
R2 = R2(P_R2~=0);
R3 = R3(P_R3~=0);
R4 = R4(P_R4~=0);
P_R1 = P_R1(P_R1~=0)/sum(P_R1(P_R1~=0));
P_R2 = P_R2(P_R2~=0)/sum(P_R2(P_R2~=0));
P_R3 = P_R3(P_R3~=0)/sum(P_R3(P_R3~=0));
P_R4 = P_R4(P_R4~=0)/sum(P_R4(P_R4~=0));
Mean_PR1 = dot(R1,P_R1);
Mean_PR2 = dot(R2,P_R2);
Mean_PR3 = dot(R3,P_R3);
Mean_PR4 = dot(R4,P_R4);
Var_PR1 = dot((R1-Mean_PR1).^2,P_R1);
Var_PR2 = dot((R2-Mean_PR2).^2,P_R2);
Var_PR3 = dot((R3-Mean_PR3).^2,P_R3);
Var_PR4 = dot((R4-Mean_PR4).^2,P_R4);
Fano_PR1 = Var_PR1/Mean_PR1;
Fano_PR2 = Var_PR2/Mean_PR2;
Fano_PR3 = Var_PR3/Mean_PR3;
Fano_PR4 = Var_PR4/Mean_PR4;
P2D_MargX.X = XX(1,:);
P2D_MargX.XNORM = 2*XX(1,:)/maxProbes;
P2D_MargX.Y = P2Dn_MargX;
P2D_MargY.X = YY(:,1).';
P2D_MargY.XNorm = 2*YY(:,1).'/maxProbes;
P2D_MargY.Y = P2Dn_MargY;
P2D_XY.X = XX(1,:);
P2D_XY.XNORM = 2*XX(1,:)/maxProbes;
P2D_XY.Y = YY(:,1).';
P2D_XY.YNorm = 2*YY(:,1).'/maxProbes;
if (withZero)
    P2D_XY.Z = P2Dn;
    P2D_MargX.Mean = MeanX0;
    P2D_MargX.Var = VarX0;
    P2D_MargX.Fano = FanoX0;
    P2D_MargY.Mean = MeanY0;
    P2D_MargY.Var = VarY0;
    P2D_MargY.Fano = FanoY0;
    P2D_XY.arithMean = MeanZ0_ver1;
    P2D_XY.geoMean = MeanZ0_ver2;
    P2D_XY.arithVar = VarZ0_ver1;
    P2D_XY.geoVar = VarZ0_ver2;
    P2D_XY.arithFano = FanoZ0_ver1;
    P2D_XY.geoFano = FanoZ0_ver2;
    P2D_XY.PearsonCoeff = Pearson_CC0;
else
    P2D_XY.Z = P2Dn_mod;
    P2D_MargX.Mean = MeanX1;
    P2D_MargX.Var = VarX1;
    P2D_MargX.Fano = FanoX1;
    P2D_MargY.Mean = MeanY1;
    P2D_MargY.Var = VarY1;
    P2D_MargY.Fano = FanoY1;
    P2D_XY.arithMean = MeanZ1_ver1;
    P2D_XY.geoMean = MeanZ1_ver2;
    P2D_XY.arithVar = VarZ1_ver1;
    P2D_XY.geoVar = VarZ1_ver2;
    P2D_XY.arithFano = FanoZ1_ver1;
    P2D_XY.geoFano = FanoZ1_ver2;
    P2D_XY.PearsonCoeff = Pearson_CC1;
end
P2D_XdivY.X = R1.';
P2D_XdivY.Y = P_R1;
P2D_XdivY.Mean = Mean_PR1;
P2D_XdivY.Var = Var_PR1;
P2D_XdivY.Fano = Fano_PR1;
P2D_YdivX.X = R2.';
P2D_YdivX.Y = P_R2;
P2D_YdivX.Mean = Mean_PR2;
P2D_YdivX.Var = Var_PR2;
P2D_YdivX.Fano = Fano_PR2;
P2D_XFraction.X = R3.';
P2D_XFraction.Y = P_R3;
P2D_XFraction.Mean = Mean_PR3;
P2D_XFraction.Var = Var_PR3;
P2D_XFraction.Fano = Fano_PR3;
P2D_YFraction.X = R4.';
P2D_YFraction.Y = P_R4;
P2D_YFraction.Mean = Mean_PR4;
P2D_YFraction.Var = Var_PR4;
P2D_YFraction.Fano = Fano_PR4;
switch type
    case 1%Split
        PxF_Group.SplitStats{vn}.XY = P2D_XY;
        PxF_Group.SplitStats{vn}.XMarginal = P2D_MargX;
        PxF_Group.SplitStats{vn}.YMarginal = P2D_MargY;
        PxF_Group.SplitStats{vn}.XFraction = P2D_XdivY;
        PxF_Group.SplitStats{vn}.YdivX = P2D_YdivX;
        PxF_Group.SplitStats{vn}.XdivY = P2D_XFraction;
        PxF_Group.SplitStats{vn}.YFraction = P2D_YFraction;
    case 2%Weave
        PxF_Group.WeaveStats{vn}.XY = P2D_XY;
        PxF_Group.WeaveStats{vn}.XMarginal = P2D_MargX;
        PxF_Group.WeaveStats{vn}.YMarginal = P2D_MargY;
        PxF_Group.WeaveStats{vn}.XFraction = P2D_XdivY;
        PxF_Group.WeaveStats{vn}.YdivX = P2D_YdivX;
        PxF_Group.WeaveStats{vn}.XdivY = P2D_XFraction;
        PxF_Group.WeaveStats{vn}.YFraction = P2D_YFraction;
end 
end
function PxF_Group = SubStats_3D(vn,P2D,PxF_Group,type,withZero,maxProbes)
N_CTypes = length(P2D);
for c = 1:N_CTypes
    switch type
        case 1%Split
            [XX,YY] = meshgrid(0:size(P2D{c},1)-1,0:size(P2D{c},2)-1);
        case 2%Weave
            [XX,YY] = meshgrid(0:size(P2D{c},1)-1,0:size(P2D{c},2)-1);
    end     
P2Dn = P2D{c}/sum(P2D{c}(:));
P2Dn_MargX = sum(P2Dn,1);
P2Dn_MargY = sum(P2Dn,2);
if (withZero)
    MeanX0 = dot(XX(1,:),P2Dn_MargX);
    VarX0 = dot((XX(1,:)-MeanX0).^2,P2Dn_MargX);
    FanoX0 = VarX0/MeanX0;
    MeanY0 = dot(YY(:,1),P2Dn_MargY);
    VarY0 = dot((YY(:,1)-MeanY0).^2,P2Dn_MargY);
    FanoY0 = VarY0/MeanY0;
    CoVar0 = sum(P2Dn.*(XX-MeanX0).*(YY-MeanY0),'all');
    MeanZ0_ver1 = 0.5*(MeanX0+MeanY0);
    MeanZ0_ver2 = sqrt(MeanX0*MeanY0);
    VarZ0_ver1 = 0.5*(VarX0+VarY0);
    VarZ0_ver2 = sqrt(VarX0*VarY0);
    FanoZ0_ver1 = VarZ0_ver1/MeanZ0_ver1;
    FanoZ0_ver2 = VarZ0_ver2/MeanZ0_ver2;
    Pearson_CC0 = CoVar0/sqrt(VarX0*VarY0);
else
    P2Dn_mod = P2Dn;
    P2Dn_mod(1,1) = 0;
    P2Dn_mod = P2Dn_mod/sum(P2Dn_mod(:));
    MeanX1 = dot(XX(1,2:end),P2Dn_MargX(2:end))/sum(P2Dn_MargX(2:end));
    VarX1 = dot((XX(1,2:end)-MeanX1).^2,P2Dn_MargX(2:end))/sum(P2Dn_MargX(2:end));
    FanoX1 = VarX1/MeanX1;
    MeanY1 = dot(YY(2:end,1),P2Dn_MargY(2:end))/sum(P2Dn_MargY(2:end));
    VarY1 = dot((YY(2:end,1)-MeanY1).^2,P2Dn_MargY(2:end))/sum(P2Dn_MargY(2:end));
    FanoY1 = VarY1/MeanY1;
    CoVar1 = sum(P2Dn_mod.*(XX-MeanX1).*(YY-MeanY1),'all');
    MeanZ1_ver1 = 0.5*(MeanX1+MeanY1);
    MeanZ1_ver2 = sqrt(MeanX1*MeanY1);
    VarZ1_ver1 = 0.5*(VarX1+VarY1);
    VarZ1_ver2 = sqrt(VarX1*VarY1);
    FanoZ1_ver1 = VarZ1_ver1/MeanZ1_ver1;
    FanoZ1_ver2 = VarZ1_ver2/MeanZ1_ver2;
    Pearson_CC1 = CoVar1/sqrt(VarX1*VarY1);
end
r1 = XX./YY;
r2 = YY./XX;
r3 = XX./(YY+XX);
r4 = YY./(YY+XX);
R1 = unique(r1);
R2 = unique(r2);
R3 = unique(r3);
R4 = unique(r4);
P_R1 = zeros(1,length(R1));
P_R2 = zeros(1,length(R2));
P_R3 = zeros(1,length(R3));
P_R4 = zeros(1,length(R4));
%make this step faster
for v = 1:length(R1)
   [xloc,yloc] = find(R1(v)==r1); 
   P_R1(v) = sum(P2Dn(sub2ind(size(P2Dn),xloc,yloc)));
end
for v = 1:length(R2)
   [xloc,yloc] = find(R2(v)==r2); 
   P_R2(v) = sum(P2Dn(sub2ind(size(P2Dn),xloc,yloc)));
end
for v = 1:length(R3)
   [xloc,yloc] = find(R3(v)==r3); 
   P_R3(v) = sum(P2Dn(sub2ind(size(P2Dn),xloc,yloc)));
end
for v = 1:length(R4)
   [xloc,yloc] = find(R4(v)==r4); 
   P_R4(v) = sum(P2Dn(sub2ind(size(P2Dn),xloc,yloc)));
end
R1 = R1(P_R1~=0);
R2 = R2(P_R2~=0);
R3 = R3(P_R3~=0);
R4 = R4(P_R4~=0);
P_R1 = P_R1(P_R1~=0)/sum(P_R1(P_R1~=0));
P_R2 = P_R2(P_R2~=0)/sum(P_R2(P_R2~=0));
P_R3 = P_R3(P_R3~=0)/sum(P_R3(P_R3~=0));
P_R4 = P_R4(P_R4~=0)/sum(P_R4(P_R4~=0));
Mean_PR1 = dot(R1,P_R1);
Mean_PR2 = dot(R2,P_R2);
Mean_PR3 = dot(R3,P_R3);
Mean_PR4 = dot(R4,P_R4);
Var_PR1 = dot((R1-Mean_PR1).^2,P_R1);
Var_PR2 = dot((R2-Mean_PR2).^2,P_R2);
Var_PR3 = dot((R3-Mean_PR3).^2,P_R3);
Var_PR4 = dot((R4-Mean_PR4).^2,P_R4);
Fano_PR1 = Var_PR1/Mean_PR1;
Fano_PR2 = Var_PR2/Mean_PR2;
Fano_PR3 = Var_PR3/Mean_PR3;
Fano_PR4 = Var_PR4/Mean_PR4;
P2D_MargX.X(1:size(P2Dn,1),c) = XX(1,:);
P2D_MargX.XNORM(1:size(P2Dn,1),c) = 2*XX(1,:)/maxProbes;
P2D_MargX.Y(1:size(P2Dn,2),c) = P2Dn_MargX;
P2D_MargY.X(1:size(P2Dn,1),c) = YY(:,1).';
P2D_MargY.XNorm(1:size(P2Dn,1),c) = 2*YY(:,1).'/maxProbes;
P2D_MargY.Y(1:size(P2Dn,2),c) = P2Dn_MargY;
P2D_XY.X(1:size(P2Dn,1),c) = XX(1,:);
P2D_XY.XNORM(1:size(P2Dn,1),c) = 2*XX(1,:)/maxProbes;
P2D_XY.Y(1:size(P2Dn,2),c) = YY(:,1).';
P2D_XY.YNorm(1:size(P2Dn,2),c) = 2*YY(:,1).'/maxProbes;
if (withZero)
    P2D_XY.Z(1:size(P2Dn,1),1:size(P2Dn,2),c) = P2Dn;
    P2D_MargX.Mean(c) = MeanX0;
    P2D_MargX.Var(c) = VarX0;
    P2D_MargX.Fano(c) = FanoX0;
    P2D_MargY.Mean(c) = MeanY0;
    P2D_MargY.Var(c) = VarY0;
    P2D_MargY.Fano(c) = FanoY0;
    P2D_XY.arithMean(c) = MeanZ0_ver1;
    P2D_XY.geoMean(c) = MeanZ0_ver2;
    P2D_XY.arithVar(c) = VarZ0_ver1;
    P2D_XY.geoVar(c) = VarZ0_ver2;
    P2D_XY.arithFano(c) = FanoZ0_ver1;
    P2D_XY.geoFano(c) = FanoZ0_ver2;
    P2D_XY.PearsonCoeff(c) = Pearson_CC0;
else
    P2D_XY.Z(1:size(P2Dn,1),1:size(P2Dn,2),c) = P2Dn_mod;
    P2D_MargX.Mean(c) = MeanX1;
    P2D_MargX.Var(c) = VarX1;
    P2D_MargX.Fano(c) = FanoX1;
    P2D_MargY.Mean(c) = MeanY1;
    P2D_MargY.Var(c) = VarY1;
    P2D_MargY.Fano(c) = FanoY1;
    P2D_XY.arithMean(c) = MeanZ1_ver1;
    P2D_XY.geoMean(c) = MeanZ1_ver2;
    P2D_XY.arithVar(c) = VarZ1_ver1;
    P2D_XY.geoVar(c) = VarZ1_ver2;
    P2D_XY.arithFano(c) = FanoZ1_ver1;
    P2D_XY.geoFano(c) = FanoZ1_ver2;
    P2D_XY.PearsonCoeff(c) = Pearson_CC1;
end
P2D_XdivY.X(1:size(P2Dn,1),c) = R1.';
P2D_XdivY.Y(1:size(P2Dn,2),c) = P_R1;
P2D_XdivY.Mean(c) = Mean_PR1;
P2D_XdivY.Var(c) = Var_PR1;
P2D_XdivY.Fano(c) = Fano_PR1;
P2D_YdivX.X(1:size(P2Dn,1),c) = R2.';
P2D_YdivX.Y(1:size(P2Dn,2),c) = P_R2;
P2D_YdivX.Mean(c) = Mean_PR2;
P2D_YdivX.Var(c) = Var_PR2;
P2D_YdivX.Fano(c) = Fano_PR2;
P2D_XFraction.X(1:size(P2Dn,1),c) = R3.';
P2D_XFraction.Y(1:size(P2Dn,2),c) = P_R3;
P2D_XFraction.Mean(c) = Mean_PR3;
P2D_XFraction.Var(c) = Var_PR3;
P2D_XFraction.Fano(c) = Fano_PR3;
P2D_YFraction.X(1:size(P2Dn,1),c) = R4.';
P2D_YFraction.Y(1:size(P2Dn,2),c) = P_R4;
P2D_YFraction.Mean(c) = Mean_PR4;
P2D_YFraction.Var(c) = Var_PR4;
P2D_YFraction.Fano(c) = Fano_PR4;
end
switch type
    case 1%Split
        PxF_Group.SplitStats{vn}.XY = P2D_XY;
        PxF_Group.SplitStats{vn}.XMarginal = P2D_MargX;
        PxF_Group.SplitStats{vn}.YMarginal = P2D_MargY;
        PxF_Group.SplitStats{vn}.XFraction = P2D_XdivY;
        PxF_Group.SplitStats{vn}.YdivX = P2D_YdivX;
        PxF_Group.SplitStats{vn}.XdivY = P2D_XFraction;
        PxF_Group.SplitStats{vn}.YFraction = P2D_YFraction;
    case 2%Weave
        PxF_Group.WeaveStats{vn}.XY = P2D_XY;
        PxF_Group.WeaveStats{vn}.XMarginal = P2D_MargX;
        PxF_Group.WeaveStats{vn}.YMarginal = P2D_MargY;
        PxF_Group.WeaveStats{vn}.XFraction = P2D_XdivY;
        PxF_Group.WeaveStats{vn}.YdivX = P2D_YdivX;
        PxF_Group.WeaveStats{vn}.XdivY = P2D_XFraction;
        PxF_Group.WeaveStats{vn}.YFraction = P2D_YFraction;
end
end
