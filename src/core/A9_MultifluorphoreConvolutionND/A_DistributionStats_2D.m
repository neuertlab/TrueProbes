function [PxF_Group] = A_DistributionStats_2D(I,settings,probes,gene_table,DoesProbeBindSite,Kb,Kb_Complement,ExpressionMatrix,onTypeIndx,offTypeIndx)
ProbeConc =  10^5;P_1FLAPs = [0 1];
LocMax = max(cell2mat(cellfun(@(x) x,{probes{:,3}},'UniformOutput',false)));  
Lpmin = min(cell2mat(cellfun(@length,{probes{:,2}},'UniformOutput',false)));
TargetLength = LocMax + Lpmin - 1;
maxProbes = floor(TargetLength/(Lpmin+settings.ProbeSpacing));
onDNA = onTypeIndx{1};
onRNA = onTypeIndx{2};
offDNA = offTypeIndx{1};
offRNA = offTypeIndx{2};
GeneName = settings.GeneName;
GeneTarget = settings.transcript_IDs;
CellTypeID = settings.CellType_ExprID;
expValType = settings.expressionValType;
isDNA = union(onTypeIndx{1},offTypeIndx{1});
isRNA = union(onTypeIndx{2},offTypeIndx{2});
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
switch expValType
    case 1
        ExprLevels = mean(ExpressionMatrix,2);
    case 2
        ExprLevels = ExpressionMatrix(:,CellTypeID);
    case 3
        ExprLevels = mean(ExpressionMatrix,2);      
    case 4
        ExprLevels = ExpressionMatrix(:,CellTypeID);
end
Isoforms = find(contains(Names,GeneName));
Desired_Isoforms = find(contains(uniNames,extractBefore(GeneTarget{:},'.')));
UnDesired_Isoforms = setdiff(Isoforms,Desired_Isoforms);
offRNA_minusIso = setdiff(offRNA,UnDesired_Isoforms);
clear gene_table2 gene_table3 RNA_ID* contains_RNA MinusStrandedHits RNA_MissedFilteredHits

CTargets_Total_TargetSite = sparse(repmat(ExprLevels,[1 size(Kb,3)]));
CProbes_Total = ProbeConc*ones(size(Kb,1),1);
[K_S,K_CD,~,~] = A_JH_GenerateSecondaryStructureInfo(probes,I,settings);
K_CDeff = K_CD + sparse(diag(diag(full(K_CD)))); %takes care of multiplicity 
GuessConc = 1/ProbeConc;
CProbes_Free = GuessConc*ones(length(I),1);

fun = @(x) dX_ProbeSet(x,I,K_S,K_CDeff,Kb,Kb_Complement,DoesProbeBindSite,isDNA,isRNA,CProbes_Total,CTargets_Total_TargetSite);
try
x1 = CProbes_Free;err = 1;tol = 1e-8;Nruns = 1;
while (err>tol&&Nruns<10)
    xn = fun(x1) + x1;
    err = abs(sum(xn./x1,'omitnan')-length(x1));
    x1 = xn;
    Nruns = Nruns + 1;
end
catch
    x1 = CProbes_Free;err = 1;tol = 1e-8;Nruns = 1;
    while (all(err>tol)&&all(Nruns<10))
        xn = abs(fun(x1) + x1);
        err = abs(sum(xn./x1,'omitnan')-length(x1));
        x1 = abs(xn);
        Nruns = Nruns + 1;
    end   
end
CProbes_Free = x1;
clear x1 xn 

h_S = K_S(I);
h_D = (K_CDeff(I,I)*CProbes_Free);
h_sub_base = sparse(1+squeeze(multiprod(permute(Kb(I,:,:).*DoesProbeBindSite(I,:,:),[2 1 3]),CProbes_Free)));
h_sub_DNA = sparse((sqrt(1+4*CTargets_Total_TargetSite(isDNA,:)./h_sub_base(isDNA,:))-1)./(2*Kb_Complement(isDNA,:)));
h_sub_RNA = CTargets_Total_TargetSite(isRNA,:)./h_sub_base(isRNA,:);
h_sub_DNA_complement = CTargets_Total_TargetSite(isDNA,:)./(1+h_sub_DNA./Kb_Complement(isDNA,:));
h_sub_DNA(isinf(h_sub_DNA))=0;
h_sub_RNA(isinf(h_sub_RNA))=0;
h_sub_DNA_complement(isinf(h_sub_DNA_complement))=0;
h_sub_DNA(isnan(h_sub_DNA))=0;
h_sub_RNA(isnan(h_sub_RNA))=0;
h_DNAx = zeros(length(I),1);
h_RNAx = zeros(length(I),1);
for v = 1:length(I)
    h_DNAx(v) = full(sum(sum(squeeze(Kb(I(v),isDNA,:)).*h_sub_DNA,1,'omitnan'),2,'omitnan'));
    h_RNAx(v) = full(sum(sum(squeeze(Kb(I(v),isRNA,:)).*h_sub_RNA,1,'omitnan'),2,'omitnan'));
end
h_DNA = sparse(h_DNAx);
h_RNA = sparse(h_RNAx);
h = 1 + h_S + h_D + h_DNA + h_RNA;
CProbes_Free = rdivide(CProbes_Total(I),h);%Conservation equilibrium
%Probability of binding sites get p_RNA p_DNA
pRNA_TargetSites_Bound = (squeeze(multiprod(permute(Kb(I,isRNA,:).*DoesProbeBindSite(I,isRNA,:),[2 1 3]),CProbes_Free)).*h_sub_RNA)./...
    (squeeze(multiprod(permute(Kb(I,isRNA,:).*DoesProbeBindSite(I,isRNA,:),[2 1 3]),CProbes_Free)).*h_sub_RNA+h_sub_RNA);
pDNA_TargetSites_Bound = (squeeze(multiprod(permute(Kb(I,isDNA,:).*DoesProbeBindSite(I,isDNA,:),[2 1 3]),CProbes_Free)).*h_sub_DNA)./...
 (squeeze(multiprod(permute(Kb(I,isDNA,:).*DoesProbeBindSite(I,isDNA,:),[2 1 3]),CProbes_Free)).*h_sub_DNA+h_sub_DNA+Kb_Complement(isDNA,:).*h_sub_DNA.*h_sub_DNA_complement);
p_TargetSites_Bound = ndSparse.build([size(Kb,2),size(Kb,3)],0);
p_TargetSites_Bound(isDNA,:) = pDNA_TargetSites_Bound;
p_TargetSites_Bound(isRNA,:) = pRNA_TargetSites_Bound;
p_TargetSites_Bound(isnan(p_TargetSites_Bound)) = 0;
p_TargetSites_Unbound = 1 - p_TargetSites_Bound;

T = find(sum(squeeze(sum(DoesProbeBindSite(I,:,:),1)),2)>0);
Smax = max(full(sum(sum(DoesProbeBindSite(I,T,:),3),1)));
PxF = cell(1,length(T));
for i = 1:length(T)%skip ones fano of NaN 
    Ti = T(i);
    vs = find(sum(squeeze(DoesProbeBindSite(I,Ti,:)),1)>0);
    PxF{i} = 1;
    for j = 1:length(vs)
        v = vs(j);
        PxF{i} = [PxF{i}*full(p_TargetSites_Unbound(Ti,v)+p_TargetSites_Bound(Ti,v)*P_1FLAPs(1)) 0] + [0 PxF{i}*full(p_TargetSites_Bound(Ti,v))*P_1FLAPs(2)];
        PxF{i} = PxF{i}/sum(PxF{i});   
    end
    PxF{i} = [PxF{i} zeros(1,Smax-length(PxF{i})+1)];
end
%two-color tests
PColor2_Split = [1*ones(1,floor(length(I)/2)) 2*ones(1,ceil(length(I)/2))];
if (mod(length(I),2)==0)
    PColor2_Weave = repmat([1 2],[1 floor(length(I)/2)]);
else
    PColor2_Weave = [repmat([1 2],[1 floor(length(I)/2)]) 1];
end
N_Channels = length(unique(PColor2_Split));
Ix_Split = cell(1,N_Channels);
piDNA_TargetSites_Bound_Split = cell(1,N_Channels);
piRNA_TargetSites_Bound_Split = cell(1,N_Channels);
pi_TargetSites_Bound_Split = cell(1,N_Channels);
Ix_Weave = cell(1,N_Channels);
piDNA_TargetSites_Bound_Weave = cell(1,N_Channels);
piRNA_TargetSites_Bound_Weave = cell(1,N_Channels);
pi_TargetSites_Bound_Weave = cell(1,N_Channels);
for k = 1:N_Channels
    Ix_Split{k} = I(PColor2_Split==k);    
    Ix_Weave{k} = I(PColor2_Weave==k);    
    piDNA_TargetSites_Bound_Split{k} = (squeeze(multiprod(permute(Kb(Ix_Split{k},isDNA,:).*DoesProbeBindSite(Ix_Split{k},isDNA,:),[2 1 3]),CProbes_Free(PColor2_Split==k))).*h_sub_DNA)./...
    (squeeze(multiprod(permute(Kb(I,isDNA,:).*DoesProbeBindSite(I,isDNA,:),[2 1 3]),CProbes_Free)).*h_sub_DNA+h_sub_DNA+Kb_Complement(isDNA,:).*h_sub_DNA.*h_sub_DNA_complement);   
    piRNA_TargetSites_Bound_Split{k} = (squeeze(multiprod(permute(Kb(Ix_Split{k},isRNA,:).*DoesProbeBindSite(Ix_Split{k},isRNA,:),[2 1 3]),CProbes_Free(PColor2_Split==k))).*h_sub_RNA)./...
    (squeeze(multiprod(permute(Kb(I,isRNA,:).*DoesProbeBindSite(I,isRNA,:),[2 1 3]),CProbes_Free)).*h_sub_RNA+h_sub_RNA);
    piDNA_TargetSites_Bound_Weave{k} = (squeeze(multiprod(permute(Kb(Ix_Weave{k},isDNA,:).*DoesProbeBindSite(Ix_Weave{k},isDNA,:),[2 1 3]),CProbes_Free(PColor2_Weave==k))).*h_sub_DNA)./...
    (squeeze(multiprod(permute(Kb(I,isDNA,:).*DoesProbeBindSite(I,isDNA,:),[2 1 3]),CProbes_Free)).*h_sub_DNA+h_sub_DNA+Kb_Complement(isDNA,:).*h_sub_DNA.*h_sub_DNA_complement);
    piRNA_TargetSites_Bound_Weave{k} = (squeeze(multiprod(permute(Kb(Ix_Weave{k},isRNA,:).*DoesProbeBindSite(Ix_Weave{k},isRNA,:),[2 1 3]),CProbes_Free(PColor2_Weave==k))).*h_sub_RNA)./...
    (squeeze(multiprod(permute(Kb(I,isRNA,:).*DoesProbeBindSite(I,isRNA,:),[2 1 3]),CProbes_Free)).*h_sub_RNA+h_sub_RNA);
    pi_TargetSites_Bound_Split{k} = ndSparse.build([size(Kb,2),size(Kb,3)],0);
    pi_TargetSites_Bound_Weave{k} = ndSparse.build([size(Kb,2),size(Kb,3)],0);
    pi_TargetSites_Bound_Split{k}(isDNA,:) = piDNA_TargetSites_Bound_Split{k};
    pi_TargetSites_Bound_Split{k}(isRNA,:) = piRNA_TargetSites_Bound_Split{k};
    pi_TargetSites_Bound_Split{k}(isnan(pi_TargetSites_Bound_Split{k})) = 0;
    pi_TargetSites_Bound_Weave{k}(isDNA,:) = piDNA_TargetSites_Bound_Weave{k};
    pi_TargetSites_Bound_Weave{k}(isRNA,:) = piRNA_TargetSites_Bound_Weave{k};
    pi_TargetSites_Bound_Weave{k}(isnan(pi_TargetSites_Bound_Weave{k})) = 0;
end
clear CProbes_Free CProbes_Total CTargets_Total_TargetSite err 
clear h h_D h_DNA h_DNAx h_RNA h_RNAx h_S h_sub_base h_sub_DNA h_sub_RNA h_sub_DNA_complement
PxFD_Split = cell(length(T),N_Channels);PxFD_Weave = cell(length(T),N_Channels);
PxF2D_Split = cell(1,length(T));PxF2D_Weave = cell(1,length(T));
PxF2D_Expr_Split = cell(1,length(T));PxF2D_Expr_Weave = cell(1,length(T));
Dim_Split = zeros(1,length(T));Dim_Weave = zeros(1,length(T));
for i = 1:length(T)
   Ti = T(i);  
   vsi_Split = cell(1,N_Channels);vsi_Weave = cell(1,N_Channels);
   for k = 1:N_Channels
      vsi_Split{k} = find(sum(squeeze(DoesProbeBindSite(Ix_Split{k},Ti,:)),1)>0);
      vsi_Weave{k} = find(sum(squeeze(DoesProbeBindSite(Ix_Weave{k},Ti,:)),1)>0);
   end
   overlap_Split = intersect(vsi_Split{1},vsi_Split{2});
   overlap_Weave = intersect(vsi_Weave{1},vsi_Weave{2});
   for k = 1:N_Channels
       PxFD_Split{i}{k} = 1;PxFD_Weave{i}{k} = 1;
       for j = 1:length(vsi_Split{k})
           v = vsi_Split{k}(j);
           PxFD_Split{i}{k} = [PxFD_Split{i}{k}*full(1-pi_TargetSites_Bound_Split{k}(Ti,v)+pi_TargetSites_Bound_Split{k}(Ti,v)*P_1FLAPs(1)) 0] + [0 PxFD_Split{i}{k}*full(pi_TargetSites_Bound_Split{k}(Ti,v))*P_1FLAPs(2)];
           PxFD_Split{i}{k} = PxFD_Split{i}{k}/sum(PxFD_Split{i}{k});
       end
       for j = 1:length(vsi_Weave{k})
           v = vsi_Weave{k}(j);
           PxFD_Weave{i}{k} = [PxFD_Weave{i}{k}*full(1-pi_TargetSites_Bound_Weave{k}(Ti,v)+pi_TargetSites_Bound_Weave{k}(Ti,v)*P_1FLAPs(1)) 0] + [0 PxFD_Weave{i}{k}*full(pi_TargetSites_Bound_Weave{k}(Ti,v))*P_1FLAPs(2)];
           PxFD_Weave{i}{k} = PxFD_Weave{i}{k}/sum(PxFD_Weave{i}{k});  
       end
   end
   PxF2D_Split{i} = PxFD_Split{i}{2}.*PxFD_Split{i}{1}';PxF2D_Weave{i} = PxFD_Weave{i}{2}.*PxFD_Weave{i}{1}';
   PxF2D_Split{i} = PxF2D_Split{i}./sum(sum(PxF2D_Split{i},1),2);PxF2D_Weave{i} = PxF2D_Weave{i}./sum(sum(PxF2D_Weave{i},1),2);
   for j = 1:length(overlap_Split)
      v = overlap_Split(j);
      p1_Split = full(pi_TargetSites_Bound_Split{1}(Ti,v));%1st dimension k = 1  
      p2_Split = full(pi_TargetSites_Bound_Split{2}(Ti,v));%2nd dimension k = 2
      pub_Split = 1-p1_Split-p2_Split; 
      tmp1_Split = [PxF2D_Split{i}*(pub_Split+p1_Split*P_1FLAPs(1)) ; 0*ones(1,size(PxF2D_Split{i},2))]+[0*ones(1,size(PxF2D_Split{i},2)); PxF2D_Split{i}*(p1_Split*P_1FLAPs(2))];
      tmp2_Split = [PxF2D_Split{i}*(pub_Split+p2_Split*P_1FLAPs(1)) 0*ones(size(PxF2D_Split{i},1),1)]+ [0*ones(size(PxF2D_Split{i},1),1) PxF2D_Split{i}*(p2_Split*P_1FLAPs(2))];
      tmp_Split = [tmp1_Split 0*ones(size(PxF2D_Split{i},1)+1,1)] + [tmp2_Split; 0*ones(1,size(PxF2D_Split{i},2)+1)];
      PxF2D_Split{i} = tmp_Split/sum(sum(tmp_Split,1),2);
   end
   for j = 1:length(overlap_Weave)
      v = overlap_Weave(j);
      p1_Weave = full(pi_TargetSites_Bound_Weave{1}(Ti,v));%1st dimension k = 1  
      p2_Weave = full(pi_TargetSites_Bound_Weave{2}(Ti,v));%2nd dimension k = 2
      pub_Weave = 1-p1_Weave-p2_Weave; 
      tmp1_Weave = [PxF2D_Weave{i}*(pub_Weave+p1_Weave*P_1FLAPs(1)) ; 0*ones(1,size(PxF2D_Weave{i},2))]+[0*ones(1,size(PxF2D_Weave{i},2)); PxF2D_Weave{i}*(p1_Weave*P_1FLAPs(2))];
      tmp2_Weave = [PxF2D_Weave{i}*(pub_Weave+p2_Weave*P_1FLAPs(1)) 0*ones(size(PxF2D_Weave{i},1),1)]+ [0*ones(size(PxF2D_Weave{i},1),1) PxF2D_Weave{i}*(p2_Weave*P_1FLAPs(2))]; 
      tmp_Weave = [tmp1_Weave 0*ones(size(PxF2D_Weave{i},1)+1,1)] + [tmp2_Weave; 0*ones(1,size(PxF2D_Weave{i},2)+1)];
      PxF2D_Weave{i} = tmp_Weave/sum(sum(tmp_Weave,1),2);
   end
   Dim_Split(i) = max(size(PxF2D_Split{i}))+1;
   Dim_Weave(i) = max(size(PxF2D_Weave{i}))+1;
   PxF2D_Split{i} = padarray(PxF2D_Split{i},[max(size(PxF2D_Split{i}))+1-size(PxF2D_Split{i},1) max(size(PxF2D_Split{i}))+1-size(PxF2D_Split{i},2)],0,'post');
   PxF2D_Weave{i} = padarray(PxF2D_Weave{i},[max(size(PxF2D_Weave{i}))+1-size(PxF2D_Weave{i},1) max(size(PxF2D_Weave{i}))+1-size(PxF2D_Weave{i},2)],0,'post');
   PxF2D_Expr_Split{i} = ExprLevels(Ti)*PxF2D_Split{i};
   PxF2D_Expr_Weave{i} = ExprLevels(Ti)*PxF2D_Weave{i};
end
clear pDNA_TargetSites_Bound pRNA_TargetSites_Bound
clear p_TargetSites_Bound p_TargetSites_Unbound
clear piDNA_TargetSites_Bound_Split piDNA_TargetSites_Bound_Weave
clear piRNA_TargetSites_Bound_Split piRNA_TargetSites_Bound_Weave
clear pi_TargetSites_Bound_Split pi_TargetSites_Bound_Weave
clear PxFD_Split PxFD_Weave
DimL = max([Dim_Split Dim_Weave]);
for i = 1:length(T)
   PxF2D_Split{i} = padarray(PxF2D_Split{i},[DimL-size(PxF2D_Split{i},1) DimL-size(PxF2D_Split{i},2)],0,'post');
   PxF2D_Expr_Split{i} = padarray(PxF2D_Expr_Split{i},[DimL-size(PxF2D_Expr_Split{i},1) DimL-size(PxF2D_Expr_Split{i},2)],0,'post');
   PxF2D_Weave{i} = padarray(PxF2D_Weave{i},[DimL-size(PxF2D_Weave{i},1) DimL-size(PxF2D_Weave{i},2)],0,'post');
   PxF2D_Expr_Weave{i} = padarray(PxF2D_Expr_Weave{i},[DimL-size(PxF2D_Expr_Weave{i},1) DimL-size(PxF2D_Expr_Weave{i},2)],0,'post'); 
end
PxF_Matrix = vertcat(PxF{:});
PxF_ExprMatrix = ExprLevels(T).*PxF_Matrix; 
PxF_Group.Fano = zeros(1,6+length(Isoforms));
PxF_Group.Fano = zeros(1,6+length(Isoforms));
PxF_Group.Skew = zeros(1,6+length(Isoforms));
PxF_Group.Kurt = zeros(1,6+length(Isoforms));
PxF_Count_G = @(G) sum(PxF_ExprMatrix(G,2:end),1);
PxF_Dist_G = @(G) sum(PxF_ExprMatrix(G,2:end),1)/sum(sum(PxF_ExprMatrix(G,2:end),1),2);
PxF_Mean_G = @(G) dot(1:Smax,sum(PxF_ExprMatrix(G,2:end),1))/sum(sum(PxF_ExprMatrix(G,2:end),1),2);
PxF_Var_G = @(G) dot((linspace(1,Smax,Smax)-PxF_Mean_G(G)).^2,sum(PxF_ExprMatrix(G,2:end),1))/sum(sum(PxF_ExprMatrix(G,2:end),1),2);
PxF_Fano_G = @(G) PxF_Var_G(G)/PxF_Mean_G(G);
PxF_Moment_G = @(G,n) dot((linspace(1,Smax,Smax)-PxF_Mean_G(G)).^n,sum(PxF_ExprMatrix(G,2:end),1))/(sum(sum(PxF_ExprMatrix(G,2:end),1),2)*PxF_Var_G(G)^(n/2));   

%Gibbs Sampling
Corr_ONDNA = arrayfun(@(x) find(T==offRNA(x)),1:length(onDNA),'Un',0);
Corr_ONRNA_IsoSpecific = arrayfun(@(x) find(T==onRNA(x)),1:length(onRNA),'Un',0);
Corr_ONRNA_IsoAgnostic = arrayfun(@(x) find(T==Isoforms(x)),1:length(Isoforms),'Un',0);
Corr_OFFDNA = arrayfun(@(x) find(T==offDNA(x)),1:length(offDNA),'Un',0);
Corr_OFFRNA_withIsos = arrayfun(@(x) find(T==offRNA(x)),1:length(offRNA),'Un',0);
Corr_OFFRNA_minusIso = arrayfun(@(x) find(T==offRNA_minusIso(x)),1:length(offRNA_minusIso),'Un',0);
Corr_ONDNA = [Corr_ONDNA{:}];
Corr_ONRNA_IsoSpecific = [Corr_ONRNA_IsoSpecific{:}];
Corr_ONRNA_IsoAgnostic = [Corr_ONRNA_IsoAgnostic{:}];
Corr_OFFDNA = [Corr_OFFDNA{:}];
Corr_OFFRNA_withIsos = [Corr_OFFRNA_withIsos{:}];
Corr_OFFRNA_minusIso = [Corr_OFFRNA_minusIso{:}];
Corr_OFFRNA_otherIso = setdiff(Corr_ONRNA_IsoAgnostic,Corr_ONRNA_IsoSpecific);

W_Norm = sum(PxF_ExprMatrix(:,2:end),1);%count dist
W_ONDNA = sum(PxF_ExprMatrix(Corr_ONDNA,2:end),1);%count dist
W_ONRNA_IsoSpecific = sum(PxF_ExprMatrix(Corr_ONRNA_IsoSpecific,2:end),1);%count dist
W_ONRNA_IsoAgnostic = sum(PxF_ExprMatrix(Corr_ONRNA_IsoAgnostic,2:end),1);%count dist
W_OFFDNA = sum(PxF_ExprMatrix(Corr_OFFDNA,2:end),1);%count dist
W_OFFRNA_withIsos = sum(PxF_ExprMatrix(Corr_OFFRNA_withIsos,2:end),1);%count dist
W_OFFRNA_minusIso = sum(PxF_ExprMatrix(Corr_OFFRNA_minusIso,2:end),1);%count dist
W_OFFRNA_otherIso = sum(PxF_ExprMatrix(Corr_OFFRNA_otherIso,2:end),1);%count dist 

%check ability to distiguish isoforms (P(x) similarity)
if (length(Isoforms)>1)
    for u = 1:length(Isoforms)
        for v = 1:length(Isoforms)
            PxF_Group.D_Isoforms(u,v) = 1/2*trapz(0:Smax,abs(PxF{Corr_ONRNA_IsoAgnostic(u)}-PxF{Corr_ONRNA_IsoAgnostic(v)}));
            PxF_Group.D_Isoforms_Split(u,v) = 1/2*trapz(0:DimL-1,trapz(0:DimL-1,abs(PxF2D_Split{Corr_ONRNA_IsoAgnostic(u)}-PxF2D_Split{Corr_ONRNA_IsoAgnostic(v)})));
            PxF_Group.D_Isoforms_Weave(u,v) = 1/2*trapz(0:DimL-1,trapz(0:DimL-1,abs(PxF2D_Weave{Corr_ONRNA_IsoAgnostic(u)}-PxF2D_Weave{Corr_ONRNA_IsoAgnostic(v)})));
        end
    end
else
   PxF_Group.D_Isoforms = NaN;
   PxF_Group.D_Isoforms_Split = NaN;
   PxF_Group.D_Isoforms_Weave = NaN;
end

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
    
    PxF_Group.MeanIsoIdentifiability = mean(PxF_Group.D_Isoforms(triu(true(size(PxF_Group.D_Isoforms)),1)));%Optimize Isoform identifiability (mean,var)
    PxF_Group.SplitMeanIsoIdentifiability = mean(PxF_Group.D_Isoforms_Split(triu(true(size(PxF_Group.D_Isoforms_Split)),1)));%Optimize Isoform identifiability (mean,var)    
    PxF_Group.WeaveMeanIsoIdentifiability = mean(PxF_Group.D_Isoforms_Weave(triu(true(size(PxF_Group.D_Isoforms_Weave)),1)));%Optimize Isoform identifiability (mean,var)    

    PxF_Group.Spot.Detection.Max = Spot_Detect_Max;
    PxF_Group.Spot.Detection.ON_IsoSpecific = Spot_Detect_ON_IsoSpecific;
    PxF_Group.Spot.Detection.ON_IsoAgnostic = Spot_Detect_ON_IsoAgnostic;
    PxF_Group.Spot.Detection.OFF_withIsos = Spot_Detect_OFF_withIsos;
    PxF_Group.Spot.Detection.OFF_minusIso = Spot_Detect_OFF_minusIso;
    PxF_Group.Spot.Detection.OFF_otherIso = Spot_Detect_OFF_otherIso;
    PxF_Group.Spot.Percentage.ON_IsoSpecific = Spot_Detect_ON_IsoSpecific/Spot_Detect_Max;
    PxF_Group.Spot.Percentage.ON_IsoAgnostic = Spot_Detect_ON_IsoAgnostic/Spot_Detect_Max;
    PxF_Group.Spot.Percentage.OFF_withIsos = Spot_Detect_OFF_withIsos/Spot_Detect_Max;
    PxF_Group.Spot.Percentage.OFF_minusIso = Spot_Detect_OFF_minusIso/Spot_Detect_Max;
    PxF_Group.Spot.Percentage.OFF_otherIso = Spot_Detect_OFF_otherIso/Spot_Detect_Max;
    PxF_Group.Probe.Detection.Max = Probe_Detect_Max;
    PxF_Group.Probe.Detection.ON_IsoSpecific = Probe_Detect_ON_IsoSpecific;
    PxF_Group.Probe.Detection.ON_IsoAgnostic = Probe_Detect_ON_IsoAgnostic;
    PxF_Group.Probe.Detection.OFF_withIsos = Probe_Detect_OFF_withIsos;
    PxF_Group.Probe.Detection.OFF_minusIso = Probe_Detect_OFF_minusIso;
    PxF_Group.Probe.Detection.OFF_otherIso = Probe_Detect_OFF_otherIso;
    PxF_Group.Probe.Percentage.ON_IsoSpecific = Probe_Detect_ON_IsoSpecific/Probe_Detect_Max;
    PxF_Group.Probe.Percentage.ON_IsoAgnostic = Probe_Detect_ON_IsoAgnostic/Probe_Detect_Max;
    PxF_Group.Probe.Percentage.OFF_withIsos = Probe_Detect_OFF_withIsos/Probe_Detect_Max;
    PxF_Group.Probe.Percentage.OFF_minusIso = Probe_Detect_OFF_withIsos/Probe_Detect_Max;
    PxF_Group.Probe.Percentage.OFF_otherIso = Probe_Detect_OFF_otherIso/Probe_Detect_Max;

for groups = 1:5
    switch groups 
        %case 1%DNA_ON_Target
        %    Gx = arrayfun(@(x) find(T==onDNA(x)),1:length(onDNA),'Un',0);
        %    withZero = 1;
        case 1%RNA_ON_Target
            Gx = arrayfun(@(x) find(T==onRNA(x)),1:length(onRNA),'Un',0);
            withZero = 1;
        case 2%RNA_ON_Target_IsoAgnostic
            Gx = arrayfun(@(x) find(T==Isoforms(x)),1:length(Isoforms),'Un',0);
            withZero = 1;
        %case 4%DNA_OFF_Target
        %    Gx = arrayfun(@(x) find(T==offDNA(x)),1:length(offDNA),'Un',0);
        %    withZero = 0;
        case 3%RNA_OFF_Target_withIso
            Gx = arrayfun(@(x) find(T==offRNA(x)),1:length(offRNA),'Un',0);
            withZero = 0;
        case 4%RNA_OFF_Target_withoutIso
            Gx = arrayfun(@(x) find(T==offRNA_minusIso(x)),1:length(offRNA_minusIso),'Un',0);
            withZero = 0;
        case 5%RNA_ON_EachIso
            Gx = arrayfun(@(x) find(T==Isoforms(x)),1:length(Isoforms),'Un',0);
            withZero = 1;
    end
    Gx = [Gx{:}];
        %ON-DNA, ON-RNA, ON-IsoSpec, ON-IsoAgnostic, OFF-RNAwoIso, OFF-RNAwIso, OFF-DNA
    if (groups<5)
        PxF_Group.CountY(groups,:) = PxF_Count_G(Gx);
        PxF_Group.CountCY(groups,:) = cumsum(PxF_Count_G(Gx));
        PxF_Group.DistY(groups,:) = PxF_Dist_G(Gx);
        PxF_Group.DistCY(groups,:) = cumsum(PxF_Dist_G(Gx));
        PxF_Group.DistX(groups,:) = 0:Smax;
        PxF_Group.Mean(groups) = PxF_Mean_G(Gx);
        PxF_Group.Fano(groups) = PxF_Fano_G(Gx);
        PxF_Group.Skew(groups) = PxF_Moment_G(Gx,3);
        PxF_Group.Kurt(groups) = PxF_Moment_G(Gx,4);           
    else
        for vs = 1:length(Gx)
            PxF_Group.CountY(groups-1+vs,:) = PxF_Count_G(Gx(vs));
            PxF_Group.CountCY(groups-1+vs,:) = cumsum(PxF_Count_G(Gx(vs)));
            PxF_Group.DistY(groups-1+vs,:) = PxF_Dist_G(Gx(vs)); 
            PxF_Group.DistCY(groups-1+vs,:) = cumsum(PxF_Dist_G(Gx(vs)));
            PxF_Group.DistX(groups-1+vs,:) = 0:Smax;
            PxF_Group.Mean(groups-1+vs) = PxF_Mean_G(Gx(vs));
            PxF_Group.Fano(groups-1+vs) = PxF_Fano_G(Gx(vs));
            PxF_Group.Skew(groups-1+vs) = PxF_Moment_G(Gx(vs),3);
            PxF_Group.Kurt(groups-1+vs) = PxF_Moment_G(Gx(vs),4);
        end
    end      
    if (groups<5)
        Mat_Split = zeros(DimL,DimL);Mat_Weave = zeros(DimL,DimL); 
        for vz = 1:length(Gx)
            Mat_Split = Mat_Split + PxF2D_Expr_Split{Gx(vz)}; 
            Mat_Weave = Mat_Weave + PxF2D_Expr_Weave{Gx(vz)}; 
        end
        PxF_Group = SubStats_2D(groups,Mat_Split,PxF_Group,1,withZero,maxProbes);
        PxF_Group = SubStats_2D(groups,Mat_Weave,PxF_Group,2,withZero,maxProbes);
    else
        for vs = 1:length(Gx)
            PxF_Group = SubStats_2D(groups-1+vs,PxF2D_Expr_Split{Gx(vs)},PxF_Group,1,withZero,maxProbes);
            PxF_Group = SubStats_2D(groups-1+vs,PxF2D_Expr_Weave{Gx(vs)},PxF_Group,2,withZero,maxProbes);
        end        
    end
end

end
function dC = dX_ProbeSet(CProbes_Free,I,K_S,K_CDeff,Kb,Kb_Complement,DoesProbeBindSite,isDNA,isRNA,CProbes_Total,CTargets_Total_TargetSite)
    h_S = K_S(I);
    h_D = (K_CDeff(I,I)*CProbes_Free);
    h_sub_base = sparse(1+squeeze(multiprod(permute(Kb(I,:,:).*DoesProbeBindSite(I,:,:),[2 1 3]),CProbes_Free)));
    h_sub_DNA = sparse((sqrt(1+4*CTargets_Total_TargetSite(isDNA,:)./h_sub_base(isDNA,:))-1)./(2*Kb_Complement(isDNA,:)));
    h_sub_RNA = CTargets_Total_TargetSite(isRNA,:)./h_sub_base(isRNA,:);
    h_sub_DNA(isinf(h_sub_DNA))=0;
    h_sub_RNA(isinf(h_sub_RNA))=0;
    h_sub_DNA(isnan(h_sub_DNA))=0;
    h_sub_RNA(isnan(h_sub_RNA))=0;
    h_DNAx = zeros(length(I),1);
    h_RNAx = zeros(length(I),1);
    for v = 1:length(I)
        h_DNAx(v) = full(sum(sum(squeeze(Kb(I(v),isDNA,:)).*h_sub_DNA,1,'omitnan'),2,'omitnan'));
        h_RNAx(v) = full(sum(sum(squeeze(Kb(I(v),isRNA,:)).*h_sub_RNA,1,'omitnan'),2,'omitnan'));
    end
    h_DNA = sparse(h_DNAx);
    h_RNA = sparse(h_RNAx);
    h = 1 + h_S + h_D + h_DNA + h_RNA;
    dC = rdivide(CProbes_Total(I),h)-CProbes_Free;
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