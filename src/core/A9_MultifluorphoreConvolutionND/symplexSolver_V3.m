function [probes_Target_stats] = symplexSolver_V3(I,Probe_Properties,settings,probes,gene_table,DoesProbeBindSite,Kb,Kb_Complement,ExpressionMatrix,onTypeIndx,offTypeIndx)
N_Channels = 1;ProbeConc =  10^5;
%how do I pick probes.
%matlab function optimizes x = f(x,params),  x is a vector
%plot distribution of different groups
    %ON-DNA, ON-RNA, ON-IsoSpec, ON-IsoAgnostic, OFF-RNAwoIso, OFF-RNAwIso,%OFF-DNA
onDNA = onTypeIndx{1};
onRNA = onTypeIndx{2};
offDNA = offTypeIndx{1};
offRNA = offTypeIndx{2};
GeneName = settings.GeneName;
spacing = settings.ProbeSpacing;
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
        probe_specificities = Probe_Properties.exp_weighted_KON_KOFF_ratio_AllTissueMean;
    case 2
        ExprLevels = ExpressionMatrix(:,CellTypeID);
        probe_specificities = Probe_Properties.exp_weighted_KON_KOFF_ratio_TissueSpecific_KONNorm_Mean(:,1).';
    case 3
        ExprLevels = mean(ExpressionMatrix,2);
        probe_specificities = Probe_Properties.exp_weighted_KON_KOFF_ratio_TissueSpecific_Mean;
    case 4
        ExprLevels = ExpressionMatrix(:,CellTypeID);
        probe_specificities = Probe_Properties.exp_weighted_KON_KOFF_ratio_TissueSpecific(:,CellTypeID);
end

Iz = find(contains(gene_table3.Name,'ribosomal'));
ProbesWithRibosomalHits = unique(gene_table3.ProbeNum(Iz));
if (settings.RemoveProbesWithRibosomalHits)
    ForbiddenProbes1 = ProbesWithRibosomalHits;
else
    ForbiddenProbes1 = [];
end
ForbiddenProbes2 = find(probe_specificities<settings.SpecificityThreshold);
ForbiddenProbes = union(ForbiddenProbes1,ForbiddenProbes2);

Isoforms = find(contains(Names,GeneName));
Desired_Isoforms = find(contains(uniNames,extractBefore(GeneTarget{:},'.')));
UnDesired_Isoforms = setdiff(Isoforms,Desired_Isoforms);
offRNA_minusIso = setdiff(offRNA,UnDesired_Isoforms);
clear gene_table2 gene_table3 RNA_ID* contains_RNA MinusStrandedHits RNA_MissedFilteredHits

AllowableProbes = setdiff(1:size(probes,1),ForbiddenProbes);
AllowMatrix = zeros(length(AllowableProbes),length(AllowableProbes));
for u1 = 1:length(AllowableProbes)
    u = AllowableProbes(u1);
    for v1 = 1:length(AllowableProbes)
       v = AllowableProbes(v1);
       if (u~=v)
       AllowMatrix(u1,v1) = double(length(intersect([probes{u,3}-spacing:probes{u,3}+length(probes{u,2})-1+spacing],[probes{v,3}-spacing:probes{v,3}+length(probes{v,2})-1+spacing]))>spacing);   
       end
    end
end
acceptFun = @(x) double(sum(AllowMatrix(round(x)==1,round(x)==1),'all')==0);


%Mu ON, Mu OFF, Mu_Iso Var ON, Var OFF, Var_Iso
%I = ProbeDesignResults.SoftwareResults(1).Designed_Probes;

CTargets_Total_TargetSite = sparse(repmat(ExprLevels,[1 size(Kb,3)]));
CProbes_Total = ProbeConc*ones(size(Kb,1),1);
[K_S,K_CD,~,~] = A_JH_GenerateSecondaryStructureInfo(probes,I,settings);
K_CDeff = K_CD + sparse(diag(diag(full(K_CD)))); %takes care of multiplicity 
GuessConc = 1/ProbeConc;
CProbes_Free = GuessConc*ones(length(I),1);

fun = @(x) dX_ProbeSet(x,I,K_S,K_CDeff,Kb,Kb_Complement,DoesProbeBindSite,isDNA,isRNA,CProbes_Total,CTargets_Total_TargetSite);
x1 = CProbes_Free;err = 1;tol = 1e-8;Nruns = 1;
while (err>tol&&Nruns<10)
    xn = fun(x1) + x1;
    err = abs(nansum(xn./x1)-length(x1));
    x1 = xn;
    Nruns = Nruns + 1;
end
CProbes_Free = x1;
clear x1 xn 

h_S = K_S(I);
h_D = (K_CDeff(I,I)*CProbes_Free);
h_sub_base = sparse(1+squeeze(multiprod(permute(Kb(I,:,:).*DoesProbeBindSite(I,:,:),[2 1 3]),CProbes_Free)));
h_sub_DNA = sparse((sqrt(1+4*CTargets_Total_TargetSite(isDNA,:)./h_sub_base(isDNA,:))-1)./(2*Kb_Complement(isDNA,:)));
h_sub_RNA = CTargets_Total_TargetSite(isRNA,:)./h_sub_base(isRNA,:);
h_sub_DNA_complement = CTargets_Total_TargetSite(isDNA,:)./(1+h_sub_DNA./Kb_Complement(isDNA,:));
h_DNAx = zeros(length(I),1);
h_RNAx = zeros(length(I),1);
for v = 1:length(I)
    h_DNAx(v) = full(nansum(nansum(squeeze(Kb(I(v),isDNA,:)).*h_sub_DNA,1),2));
    h_RNAx(v) = full(nansum(nansum(squeeze(Kb(I(v),isRNA,:)).*h_sub_RNA,1),2));
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
p_TargetSites_MomGen1 = 1 - 2*p_TargetSites_Bound;
p_TargetSites_MomGen2 = 1 - 6*p_TargetSites_Unbound.*p_TargetSites_Bound;
mean_probes_Target = nansum(p_TargetSites_Bound,2);
var_probes_Target = nansum(p_TargetSites_Unbound.*p_TargetSites_Bound,2);
fano_probes_Target = var_probes_Target./mean_probes_Target;
skew_probes_Target = nansum(p_TargetSites_MomGen1.*p_TargetSites_Unbound.*p_TargetSites_Bound,2);
kurt_probes_Target = nansum(p_TargetSites_MomGen2.*p_TargetSites_Unbound.*p_TargetSites_Bound,2);

probes_Target_stats = full([mean_probes_Target fano_probes_Target skew_probes_Target kurt_probes_Target]);
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
   %two-color test
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
    Ix_Split{k} = I(find(PColor2_Split==k));    
    Ix_Weave{k} = I(find(PColor2_Weave==k));    
    piDNA_TargetSites_Bound_Split{k} = (squeeze(multiprod(permute(Kb(Ix_Split{k},isDNA,:).*DoesProbeBindSite(Ix_Split{k},isDNA,:),[2 1 3]),CProbes_Free(find(PColor2_Split==k)))).*h_sub_DNA)./...
    (squeeze(multiprod(permute(Kb(I,isDNA,:).*DoesProbeBindSite(I,isDNA,:),[2 1 3]),CProbes_Free)).*h_sub_DNA+h_sub_DNA+Kb_Complement(isDNA,:).*h_sub_DNA.*h_sub_DNA_complement);   
    piRNA_TargetSites_Bound_Split{k} = (squeeze(multiprod(permute(Kb(Ix_Split{k},isRNA,:).*DoesProbeBindSite(Ix_Split{k},isRNA,:),[2 1 3]),CProbes_Free(find(PColor2_Split==k)))).*h_sub_RNA)./...
    (squeeze(multiprod(permute(Kb(I,isRNA,:).*DoesProbeBindSite(I,isRNA,:),[2 1 3]),CProbes_Free)).*h_sub_RNA+h_sub_RNA);
    piDNA_TargetSites_Bound_Weave{k} = (squeeze(multiprod(permute(Kb(Ix_Weave{k},isDNA,:).*DoesProbeBindSite(Ix_Weave{k},isDNA,:),[2 1 3]),CProbes_Free(find(PColor2_Weave==k)))).*h_sub_DNA)./...
    (squeeze(multiprod(permute(Kb(I,isDNA,:).*DoesProbeBindSite(I,isDNA,:),[2 1 3]),CProbes_Free)).*h_sub_DNA+h_sub_DNA+Kb_Complement(isDNA,:).*h_sub_DNA.*h_sub_DNA_complement);
    piRNA_TargetSites_Bound_Weave{k} = (squeeze(multiprod(permute(Kb(Ix_Weave{k},isRNA,:).*DoesProbeBindSite(Ix_Weave{k},isRNA,:),[2 1 3]),CProbes_Free(find(PColor2_Weave==k)))).*h_sub_RNA)./...
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
PxF_Var_G = @(G) dot(([1:Smax]-PxF_Mean_G(G)).^2,sum(PxF_ExprMatrix(G,2:end),1))/sum(sum(PxF_ExprMatrix(G,2:end),1),2);
PxF_Fano_G = @(G) PxF_Var_G(G)/PxF_Mean_G(G);
PxF_Moment_G = @(G,n) dot(([1:Smax]-PxF_Mean_G(G)).^n,sum(PxF_ExprMatrix(G,2:end),1))/(sum(sum(PxF_ExprMatrix(G,2:end),1),2)*PxF_Var_G(G)^(n/2));   
for groups = 1:7
    switch groups 
        case 1%DNA_ON_Target
            Gx = arrayfun(@(x) find(T==onDNA(x)),1:length(onDNA),'Un',0);
        case 2%RNA_ON_Target
            Gx = arrayfun(@(x) find(T==onRNA(x)),1:length(onRNA),'Un',0);
        case 3%RNA_ON_Target_IsoAgnostic
            Gx = arrayfun(@(x) find(T==Isoforms(x)),1:length(Isoforms),'Un',0);
        case 4%DNA_OFF_Target
            Gx = arrayfun(@(x) find(T==offDNA(x)),1:length(offDNA),'Un',0);
        case 5%RNA_OFF_Target_withIso
            Gx = arrayfun(@(x) find(T==offRNA(x)),1:length(offRNA),'Un',0);
        case 6%RNA_OFF_Target_withoutIso
            Gx = arrayfun(@(x) find(T==offRNA_minusIso(x)),1:length(offRNA_minusIso),'Un',0);
        case 7%RNA_ON_EachIso
            Gx = arrayfun(@(x) find(T==Isoforms(x)),1:length(Isoforms),'Un',0);
    end
    Gx = [Gx{:}];
        %ON-DNA, ON-RNA, ON-IsoSpec, ON-IsoAgnostic, OFF-RNAwoIso, OFF-RNAwIso, OFF-DNA
    if (groups<7)
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
    if (groups<7)
        Mat_Split = zeros(DimL,DimL);Mat_Weave = zeros(DimL,DimL); 
        for vz = 1:length(Gx)
            Mat_Split = Mat_Split + PxF2D_Expr_Split{Gx(vz)}; 
            Mat_Weave = Mat_Weave + PxF2D_Expr_Weave{Gx(vz)}; 
        end
        PxF_Group.Split{groups} = Mat_Split;
        PxF_Group.Weave{groups} = Mat_Weave;
    else
        for vs = 1:length(Gx)
            PxF_Group.Split{groups-1+vs} = PxF2D_Expr_Split{Gx(vs)};
            PxF_Group.Weave{groups-1+vs} = PxF2D_Expr_Weave{Gx(vs)}; 
        end
    end
end
%Gibbs Sampling
%Obj  - (m-mu)^2/2

%obj(m-m_desired) + obj(s-s_desired),  v = s^2


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
%look at conc accumulating vs with set number    
Conc_ONDNA = sum(W_ONDNA);
Mean_ONDNA = PxF_Mean_G(Corr_ONDNA);
Fano_ONDNA = PxF_Fano_G(Corr_ONDNA);
Conc_ONRNA_IsoSpecific = sum(W_ONRNA_IsoSpecific);
Mean_ONRNA_IsoSpecific = PxF_Mean_G(Corr_ONRNA_IsoSpecific);
Fano_ONRNA_IsoSpecific = PxF_Fano_G(Corr_ONRNA_IsoSpecific);
Conc_ONRNA_IsoAgnostic = sum(W_ONRNA_IsoAgnostic);
Mean_ONRNA_IsoAgnostic = PxF_Mean_G(Corr_ONRNA_IsoAgnostic);
Fano_ONRNA_IsoAgnostic = PxF_Fano_G(Corr_ONRNA_IsoAgnostic);
Conc_OFFDNA = sum(W_OFFDNA);
Mean_OFFDNA = PxF_Mean_G(Corr_OFFDNA);
Fano_OFFDNA = PxF_Fano_G(Corr_OFFDNA);
Conc_OFFRNA_withIsos = sum(W_OFFRNA_withIsos);
Mean_OFFRNA_withIsos = PxF_Mean_G(Corr_OFFRNA_withIsos);
Fano_OFFRNA_withIsos = PxF_Fano_G(Corr_OFFRNA_withIsos);
Conc_OFFRNA_minusIso = sum(W_OFFRNA_minusIso);
Mean_OFFRNA_minusIso = PxF_Mean_G(Corr_OFFRNA_minusIso);
Fano_OFFRNA_minusIso = PxF_Fano_G(Corr_OFFRNA_minusIso);
Conc_OFFRNA_otherIso = sum(W_OFFRNA_otherIso);
Mean_OFFRNA_otherIso = PxF_Mean_G(Corr_OFFRNA_otherIso);
Fano_OFFRNA_otherIso = PxF_Fano_G(Corr_OFFRNA_otherIso);
for v = 1:length(Corr_OFFRNA_otherIso)
    Conc_RNA_otherIso(v) = sum(sum(PxF_ExprMatrix(Corr_OFFRNA_otherIso(v),2:end),1),2);
    Mean_RNA_otherIso(v) = PxF_Mean_G(Corr_OFFRNA_otherIso(v));
    Fano_RNA_otherIso(v) = PxF_Fano_G(Corr_OFFRNA_otherIso(v));
end
%check ability to distiguish isoforms (P(x) similarity)
%check distinguish-ability between isos
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
%Minimize 
    y1 = Fano_ONDNA/(Mean_ONDNA*Conc_ONDNA);%Optimize ON_DNA binding
    y2 = Fano_ONRNA_IsoSpecific/(Mean_ONRNA_IsoSpecific*Conc_ONRNA_IsoSpecific);%Optimize ON_RNA Iso-Specific Binding
    y3 = Fano_ONRNA_IsoAgnostic/(Mean_ONRNA_IsoAgnostic*Conc_ONRNA_IsoAgnostic);%Optimize ON_RNA Iso-Agnostic Binding
    y4 = Conc_OFFDNA*Mean_OFFDNA*Fano_OFFDNA;%Optimize OFF_DNA binding
    y5 = Conc_OFFRNA_withIsos*Mean_OFFRNA_withIsos*Fano_OFFRNA_withIsos;%Optimize OFF_RNA Binding
    y6 = Conc_OFFRNA_minusIso*Mean_OFFRNA_minusIso*Fano_OFFRNA_minusIso;%Optimize OFF_Non-Iso RNA binding
    y7 = Conc_OFFRNA_otherIso*Mean_OFFRNA_otherIso*Fano_OFFRNA_otherIso;%Optimize OFF_Iso     RNA binding
    y8 = Conc_RNA_otherIso./(Conc_RNA_otherIso.*Mean_RNA_otherIso);
    y9 = mean(PxF_Group.D_Isoforms(triu(true(size(PxF_Group.D_Isoforms)),1)));%Optimize Isoform identifiability (mean,var)
    y10 = mean(PxF_Group.D_Isoforms_Split(triu(true(size(PxF_Group.D_Isoforms_Split)),1)));%Optimize Isoform identifiability (mean,var)    
    y11 = mean(PxF_Group.D_Isoforms_Weave(triu(true(size(PxF_Group.D_Isoforms_Weave)),1)));%Optimize Isoform identifiability (mean,var)    
     
    y = (y1+y2+y3+y4+y5+y6+y7+sum(y8)+y9+y10+y11)^2;
    nB = 5;
    %I 
    
    
    %m1, s1,     m1, m2,  s1,  s2 (chain).
    
    
    %return y's and 
    PxF_Group
    probes_Target_stats
    Spot_Backgd_OFF = sum(W_OFFDNA(1:nB)) + sum(W_OFFRNA_withIsos(1:nB));
    Spot_Detect_Max = sum(W_Norm);
    Probe_Backgd_OFF = sum([1:nB].*W_OFFDNA(1:nB))+ sum([1:nB].*W_OFFRNA_withIsos(1:nB));
    Probe_Detect_Max = sum([1:length(W_Norm)].*W_Norm);
    
end
function [dC,J] = dX_ProbeSet_withJacobian(CProbes_Free,I,K_S,K_CDeff,Kb,Kb_Complement,DoesProbeBindSite,isDNA,isRNA,CProbes_Total,CTargets_Total_TargetSite)
    h_S = K_S(I);
    h_D = (K_CDeff(I,I)*CProbes_Free);
    h_sub_base = sparse(1+squeeze(multiprod(permute(Kb(I,:,:).*DoesProbeBindSite(I,:,:),[2 1 3]),CProbes_Free)));
    h_sub_DNA = sparse((sqrt(1+4*CTargets_Total_TargetSite(isDNA,:)./h_sub_base(isDNA,:))-1)./(2*Kb_Complement(isDNA,:)));
    h_sub_RNA = CTargets_Total_TargetSite(isRNA,:)./h_sub_base(isRNA,:);
    h_DNAx = zeros(length(I),1);
    h_RNAx = zeros(length(I),1);
    for v = 1:length(I)
        h_DNAx(v) = full(nansum(nansum(squeeze(Kb(I(v),isDNA,:)).*h_sub_DNA,1),2));
        h_RNAx(v) = full(nansum(nansum(squeeze(Kb(I(v),isRNA,:)).*h_sub_RNA,1),2));
    end
    h_DNA = sparse(h_DNAx);
    h_RNA = sparse(h_RNAx);
    h = 1 + h_S + h_D + h_DNA + h_RNA;
    dC = rdivide(CProbes_Total(I),h)-CProbes_Free;
    hs = h.^2;
    h_top_sub_DNA = sparse(sqrt(1+4*CTargets_Total_TargetSite(isDNA,:)./h_sub_base(isDNA,:)));
    h_top_DNA = 2*(CTargets_Total_TargetSite(isDNA,:)./(h_sub_base(isDNA,:).^2))./h_top_sub_DNA;
    h_top_RNA = CTargets_Total_TargetSite(isRNA,:)./(h_sub_base(isRNA,:).^2);
    h_top_DNAx = zeros(length(I),length(I));
    h_top_RNAx = zeros(length(I),length(I));
    J = zeros(length(I),length(I));
    for a = 1:length(I)
        for b = 1:length(I)
            h_top_DNAx(a,b) = full(nansum(nansum(squeeze(Kb(I(a),isDNA,:)).*squeeze(Kb(I(b),isDNA,:)).*h_top_DNA,1),2));
            h_top_RNAx(a,b) = full(nansum(nansum(squeeze(Kb(I(a),isRNA,:)).*squeeze(Kb(I(b),isRNA,:)).*h_top_RNA,1),2));
            J(a,b) = CProbes_Total(I(a))*(h_top_DNAx(a,b)+h_top_RNAx(a,b)-K_CDeff(I(a),I(b)))/hs(a)-double(a==b);
        end
    end
    
end
function dC = dX_ProbeSet(CProbes_Free,I,K_S,K_CDeff,Kb,Kb_Complement,DoesProbeBindSite,isDNA,isRNA,CProbes_Total,CTargets_Total_TargetSite)
    h_S = K_S(I);
    h_D = (K_CDeff(I,I)*CProbes_Free);
    h_sub_base = sparse(1+squeeze(multiprod(permute(Kb(I,:,:).*DoesProbeBindSite(I,:,:),[2 1 3]),CProbes_Free)));
    h_sub_DNA = sparse((sqrt(1+4*CTargets_Total_TargetSite(isDNA,:)./h_sub_base(isDNA,:))-1)./(2*Kb_Complement(isDNA,:)));
    h_sub_RNA = CTargets_Total_TargetSite(isRNA,:)./h_sub_base(isRNA,:);
    h_DNAx = zeros(length(I),1);
    h_RNAx = zeros(length(I),1);
    for v = 1:length(I)
        h_DNAx(v) = full(nansum(nansum(squeeze(Kb(I(v),isDNA,:)).*h_sub_DNA,1),2));
        h_RNAx(v) = full(nansum(nansum(squeeze(Kb(I(v),isRNA,:)).*h_sub_RNA,1),2));
    end
    h_DNA = sparse(h_DNAx);
    h_RNA = sparse(h_RNAx);
    h = 1 + h_S + h_D + h_DNA + h_RNA;
    dC = rdivide(CProbes_Total(I),h)-CProbes_Free;
end
function J = Jacobian_ProbeSet(CProbes_Free,I,K_S,K_CDeff,Kb,Kb_Complement,DoesProbeBindSite,isDNA,isRNA,CProbes_Total,CTargets_Total_TargetSite)
    h_S = K_S(I);
    h_D = (K_CDeff(I,I)*CProbes_Free);
    h_sub_base = sparse(1+squeeze(multiprod(permute(Kb(I,:,:).*DoesProbeBindSite(I,:,:),[2 1 3]),CProbes_Free)));
    h_sub_DNA = sparse((sqrt(1+4*CTargets_Total_TargetSite(isDNA,:)./h_sub_base(isDNA,:))-1)./(2*Kb_Complement(isDNA,:)));
    h_sub_RNA = CTargets_Total_TargetSite(isRNA,:)./h_sub_base(isRNA,:);
    h_DNAx = zeros(length(I),1);
    h_RNAx = zeros(length(I),1);
    for v = 1:length(I)
        h_DNAx(v) = full(nansum(nansum(squeeze(Kb(I(v),isDNA,:)).*h_sub_DNA,1),2));
        h_RNAx(v) = full(nansum(nansum(squeeze(Kb(I(v),isRNA,:)).*h_sub_RNA,1),2));    
    end
    h_DNA = sparse(h_DNAx);
    h_RNA = sparse(h_RNAx);
    h = 1 + h_S + h_D + h_DNA + h_RNA;
    hs = h.^2;
    h_top_sub_DNA = sparse(sqrt(1+4*CTargets_Total_TargetSite(isDNA,:)./h_sub_base(isDNA,:)));
    h_top_DNA = 2*(CTargets_Total_TargetSite(isDNA,:)./(h_sub_base(isDNA,:).^2))./h_top_sub_DNA;
    h_top_RNA = CTargets_Total_TargetSite(isRNA,:)./(h_sub_base(isRNA,:).^2);
    h_top_DNAx = zeros(length(I),length(I));
    h_top_RNAx = zeros(length(I),length(I));
    J = zeros(length(I),length(I));
    for a = 1:length(I)
        for b = 1:length(I)
            h_top_DNAx(a,b) = full(nansum(nansum(squeeze(Kb(I(a),isDNA,:)).*squeeze(Kb(I(b),isDNA,:)).*h_top_DNA,1),2));
            h_top_RNAx(a,b) = full(nansum(nansum(squeeze(Kb(I(a),isRNA,:)).*squeeze(Kb(I(b),isRNA,:)).*h_top_RNA,1),2));
            J(a,b) = CProbes_Total(I(a))*(h_top_DNAx(a,b)+h_top_RNAx(a,b)-K_CDeff(I(a),I(b)))/hs(a)-double(a==b);
        end
    end
end