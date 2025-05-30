function ModelMetrics = RNAsolver_JH(Pset,settings,probes,gene_table,ExpressionMatrix,DoesProbeBindSite,dHeq_mod,dSeq_mod,dCp_mod,dHeq_Complement,dSeq_Complement,dCp_Complement)   
% This function computes the metrics for RNA-FISH performance.
% The function simulates probe binding dynamics to model probe
% equilibrium binding under a variety of conditions. Varying by
% temperature, expression, probe concentration, etc.
% The distribution of probes bound targets is computed and is 
% used to generate confusion matrix metrics, of probe performance.

% An option is also included to translate this into multiple fluorphores,
% And generate statistics when using any number of fluorphores.
% By Default is option is set to off (AS=0).

%Conditions [Cell-Line (c),NN Model Used (m), probe dilution (d), temperature (t)]. 

R = 0.001987204259;%gas constant [Energy Kcal/mol K
kb = 1.380649*10^-23;%bolzman constant J/K
h = 6.62607015*10^-34;%planks constant J/Hz
Tref = 37+273.15;    

%both isospecific and agnostic output w/wout other isos as off-targets?
GeneName = settings.GeneName;
GeneTarget = settings.transcript_IDs;
gene_table = sortrows(gene_table,[7 6],'ascend');
gene_table2 = gene_table(gene_table.Match>=settings.MinHomologySearchTargetSize,:);
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
DNA_IDs_1 = find(contains(uniNames,'NC_'));%IDs
DNA_IDs_2 = find(contains(uniNames,'NT_'));%IDs
DNA_IDs_3 = find(contains(uniNames,'NW_'));%IDs
DNA_IDs =union(union(DNA_IDs_1,DNA_IDs_2),DNA_IDs_3).';
NonDNA_IDs_1 = find(~contains(uniNames,'NC_'));%IDs
NonDNA_IDs_2 = find(~contains(uniNames,'NT_'));%IDs
NonDNA_IDs_3 = find(~contains(uniNames,'NW_'));%
NonDNA_IDs = intersect(intersect(NonDNA_IDs_1,NonDNA_IDs_2),NonDNA_IDs_3).';
%probe,target,site,model,temp, dilution, cellline  
    %dG = dG - 0.114*N/2*ln(Na+)
    %dS = dG + 0.368*N/2*ln(Na+)
%Add in count statistics
%probe,target,model,temp      
   %nTargets x sites
   %1D

V_Cell = @(R) 4/3*pi*(R^3)/10^15;%um to L
Js = @(x) find(sum(squeeze(sum(DoesProbeBindSite(x,:,:),1)),2)>0);
Js_Sites = @(x) find(sum(sum(DoesProbeBindSite(x,Js(x),:),1),2)>0);
Js_DNA = @(x)DNA_IDs(ismember(DNA_IDs,Js(x)));
Js_RNA = @(x)NonDNA_IDs(ismember(NonDNA_IDs,Js(x)));

%free temp    
    NumNonUniformConditions = size(ExpressionMatrix,2);
    ExpressionMatrix(DNA_IDs,NumNonUniformConditions+1) = 2;
    ExpressionMatrix(NonDNA_IDs,NumNonUniformConditions+1) = 100;

%Setting conditions    
    Rcell = 10;
    PC0 = 5/10^6;%5uM
    Tvec = [37+273.15];
    Mvec = 1:7;
    Dvec = [1 1/2 1/10 1/100 1/1000 1/10^4 1/10^5 1/10^6];%1/10^6 is the point where it drops 
    Cvec = [size(ExpressionMatrix,2) size(ExpressionMatrix,2)];%main expression cell-lines?, and main expr tissue?
    
    Mvec = 1:4;
    Dvec = [1 1/2 1/10];%1/10^6 is the point where it drops 
    GuessConc = 10^-10;
    CProbes_Free = GuessConc*ones(length(Pset),length(Mvec),length(Tvec),length(Dvec),length(Cvec));
    CProbes_Free0 = GuessConc*ones(length(Pset),length(Mvec),length(Tvec),length(Dvec),length(Cvec));
    ProbeConc = PC0*squeeze(permute(repmat(Dvec,[1 1 length(Pset) length(Mvec) length(Tvec) length(Cvec)]),[1 3 4 5 2 6]));
       
    %P M 1 1 D
    nExpressionMatrix = ExpressionMatrix/(V_Cell(Rcell)*6.022*10^23);
	
% ON/OFF-Targets
Isoforms = find(contains(Names,GeneName));
Desired_Isoforms = find(contains(uniNames,extractBefore(GeneTarget{:},'.')));
UnDesired_Isoforms = setdiff(Isoforms,Desired_Isoforms);
ON_RNAIDs = find(strcmp(uniNames,extractBefore(GeneTarget,'.')));
OFF_RNAIDs = setdiff(NonDNA_IDs,ON_RNAIDs);
ON_RNAIDs_Isos = find(contains(Names,GeneName));
Desired_Isoforms = find(contains(uniNames,extractBefore(GeneTarget,'.')));
UnDesired_Isoforms = setdiff(ON_RNAIDs_Isos,Desired_Isoforms);
OFF_RNAIDs_minusIsos = setdiff(OFF_RNAIDs,UnDesired_Isoforms);
ON_IDs = ON_RNAIDs_Isos;
ON_IDs = Desired_Isoforms;
OFF_IDs = OFF_RNAIDs_minusIsos;

%Generate Probe Secondary Structure
    DoesProbeBindTarget = ndSparse(double(sum(DoesProbeBindSite,3)>0));
    DoesTargetHaveSite = ndSparse(double(sum(DoesProbeBindSite,1)>0));
    [SelfSeqParsed,CrossDimerSeqParsed] = A_JH_GenerateSecondaryStructureSeq(probes,Pset,settings);
    crossNotExist = cellfun(@isempty,CrossDimerSeqParsed);
    selfNotExist = cellfun(@isempty,SelfSeqParsed);
    [crossId1,crossId2] = find(crossNotExist==1);
    selfId1 = find(selfNotExist==1);
    for i = 1:length(selfId1)
       SelfSeqParsed{selfId1(i)}=cell(0,2); 
    end
    for i = 1:length(crossId1)
        CrossDimerSeqParsed{crossId1(i),crossId2(i)}=cell(0,2);
    end
    CrossDimerSeqParsed = cellfun(@(x) x(sum(cellfun(@isempty,x),2)'==0,1:2), CrossDimerSeqParsed,'Un',0);
    
    Ns = cellfun(@(x) length(x), SelfSeqParsed);
    Nc = cellfun(@(x) size(x,1), CrossDimerSeqParsed);

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
    Ns_Config = size(Ns_val,2);
    Nc_Config = size(Nc_val,3);
    c_alpha = ones(size(probes,1))+diag(diag(ones(size(probes,1))));
    
    
    [Ks,dHs_eq,dSs_eq,~,~,~,~,...
    dCps_eq,Kd,dHd_eq,dSd_eq,~,~,~,~,dCpd_eq,...
    ~,~] = A_JH_GenerateSecondaryStructureInfo_V2(probes,Pset,settings); 
	
%Format all binding equations to be parallelized for different experimental conditions.	
    %probe,config,model,temp
    DHs_eq = @(T,Tref,p,n,m) repmat(dHs_eq(p,n,m),[1 1 1 length(T)]) + ...
                             repmat(dCps_eq(p,n,m),[1 1 1 length(T)]).*...
                             permute(repmat(T'-Tref,[1 length(p) length(n) length(m)]),[2 3 4 1]);              
    DSs_eq = @(T,Tref,p,n,m) repmat(dSs_eq(p,n,m),[1 1 1 length(T)]) + ...
                             repmat(dCps_eq(p,n,m),[1 1 1 length(T)]).*...
                             permute(repmat(log(T'/Tref),[1 length(p) length(n) length(m)]),[2 3 4 1]);
    Ks_eq = @(T,Tref,p,n,m) repmat(Ns_val(p,n),[1 1 length(m) length(T)]).*...
                            spfun(@exp,-DHs_eq(T,Tref,p,n,m)./...
                            permute(repmat(R*T',[1 length(p) length(n) length(m)]),[2 3 4 1])+...
                            DSs_eq(T,Tref,p,n,m)/R);                          
    %probe,probe,config,model,temp
    DHd_eq = @(T,Tref,p,q,n,m) repmat(dHd_eq(p,q,n,m),[1 1 1 1 length(T)]) + ...
                               repmat(dCpd_eq(p,q,n,m),[1 1 1 1 length(T)]).*...
                               permute(repmat(T'-Tref,[1 length(p) length(q) length(n) length(m)]),[2 3 4 5 1]);
    DSd_eq = @(T,Tref,p,q,n,m) repmat(dSd_eq(p,q,n,m),[1 1 1 1 length(T)]) + ...
                               repmat(dCpd_eq (p,q,n,m),[1 1 1 1 length(T)]).*...
                               permute(repmat(log(T'/Tref),[1 length(p) length(q) length(n) length(m)]),[2 3 4 5 1]);  
    Kd_eq = @(T,Tref,p,q,n,m) repmat(Nc_val(p,q,n),[1 1 1 length(m) length(T)]).*...
                              repmat(c_alpha(p,q),[1 1 length(n) length(m) length(T)]).*...
                              spfun(@exp,-DHd_eq(T,Tref,p,q,n,m)./...
                              permute(repmat(R*T',[1 length(p) length(q) length(n) length(m)]),[2 3 4 5 1])+...
                              DSd_eq(T,Tref,p,q,n,m)/R);                      
    %probe,target,site,model,temp
    DHeq_mod = @(T,Tref,p,q,r,m) repmat(dHeq_mod(p,q,r,m),[1 1 1 1 length(T)]) + ...
                                 repmat(dCp_mod(p,q,r,m),[1 1 1 1 length(T)]).*...
                                 permute(repmat(T'-Tref,[1 length(p) length(q) length(r) length(m)]),[2 3 4 5 1]);
    DSeq_mod = @(T,Tref,p,q,r,m) repmat(dSeq_mod(p,q,r,m),[1 1 1 1 length(T)]) + ...
                                 repmat(dCp_mod(p,q,r,m),[1 1 1 1 length(T)]).*...
                                 permute(repmat(log(T'/Tref),[1 length(p) length(q) length(r) length(m)]),[2 3 4 5 1]);                           
    Keq_mod = @(T,Tref,p,q,r,m) repmat(DoesProbeBindSite(p,q,r),[1 1 1 length(m) length(T)]).*...
                                spfun(@exp,-DHeq_mod(T,Tref,p,q,r,m)./...
                                permute(repmat(R*T',[1 length(p) length(q) length(r) length(m)]),[2 3 4 5 1])...
                                + DSeq_mod(T,Tref,p,q,r,m)/R); 
    %probe,target,model,temp        
    DHeq_Complement = @(T,Tref,q,r,m) repmat(dHeq_Complement(q,r,m),[1 1 1 length(T)]) + ...
                                 repmat(dCp_Complement(q,r,m),[1 1 1 length(T)]).*...
                                 permute(repmat(T'-Tref,[1 length(q) length(r) length(m)]),[2 3 4 1]);
    DSeq_Complement = @(T,Tref,q,r,m) repmat(dSeq_Complement(q,r,m),[1 1 1 length(T)]) + ...
                                 repmat(dCp_Complement(q,r,m),[1 1 1 length(T)]).*...
                                 permute(repmat(log(T'/Tref),[1 length(q) length(r) length(m)]),[2 3 4 1]);                             
    Keq_Complement = @(T,Tref,q,r,m) repmat(DoesTargetHaveSite(q,r),[1 1 length(m) length(T)]).*...
                                spfun(@exp,-DHeq_Complement(T,Tref,q,r,m)./...
                                permute(repmat(R*T',[1 length(q) length(r) length(m)]),[2 3 4 1])...
                                + DSeq_Complement(T,Tref,q,r,m)/R);                                               
    %sub_base 1+KbDPS(P,T,S)*pf    Js_RNA
    sub_base = @(pf,T,Tref,p,q,r,m,d,c) dot(repmat(Keq_mod(T,Tref,p,q,r,m),[1 1 1 1 1 length(d) length(c)]),...
                                         permute(repmat(pf,[1 1 1 1 1 length(q) length(r)]),[1 6 7 2 3 4 5]),1);% [1 x N_Targets x N_sites x N_Models x N_Temperatures x N_Dilutions x N_CellLines]
    sub_base = @(pf,T,Tref,p,q,r,m,d,c) sum(repmat(Keq_mod(T,Tref,p,q,r,m),[1 1 1 1 1 length(d) length(c)]).*...
                                         permute(repmat(pf,[1 1 1 1 1 length(q) length(r)]),[1 6 7 2 3 4 5]),1);% [1 x N_Targets x N_sites x N_Models x N_Temperatures x N_Dilutions x N_CellLines]
  
                                         
    %sub_RNA =   E(:,T,S)/base(P,T,S) Js_RNA                             
    sub_RNA = @(pf,T,Tref,p,q,r,m,d,c,cExpressionMatrix) permute(repmat(cExpressionMatrix(q,c),[1 1 length(p) length(r) length(m) length(T) length(d)]),...
                                        [3 1 4 5 6 7 2])./(1+sub_base(pf,T,Tref,p,q,r,m,d,c));                                          
    sub_RNA = @(pf,T,Tref,p,q,r,m,d,c,cExpressionMatrix) permute(repmat(cExpressionMatrix(q,c),[1 1 length(p) length(r) length(m) length(T) length(d)]),...
                                        [3 1 4 5 6 7 2])./permute(repmat(permute(1+sub_base(pf,T,Tref,p,q,r,m,d,c),[2 3 4 5 6 7 1]),[1 1 1 1 1 1 length(p)]),[7 1 2 3 4 5 6]);                               
    %sub_DNA =    (sqrt(1 + 4*E(:,T,S)/sub_DNA_base(P,T,S)-1)/(2*Kb_Complement(T,S));   
    spfunction = @(x) sqrt(x+1)-1;            
    sub_DNA = @(pf,T,Tref,p,q,r,m,d,c,cExpressionMatrix) 2*spfun(spfunction,...
           4*permute(repmat(cExpressionMatrix(q,c),[1 1 length(p) length(r) length(m) length(T) length(d)]),[3 1 4 5 6 7 2]).*...
           permute(repmat(Keq_Complement(T,Tref,q,r,m),[1 1 1 1 length(p) length(d) length(c)]),[5 1 2 3 4 6 7])./...
              (1+sub_base(pf,T,Tref,p,q,r,m,d,c)))./permute(repmat(Keq_Complement(T,Tref,q,r,m),[1 1 1 1 length(p) length(d) length(c)]),[5 1 2 3 4 6 7]);      
    sub_DNA_Complement = @(pf,T,Tref,p,q,r,m,d,c,cExpressionMatrix) permute(repmat(cExpressionMatrix(q,c),...
                             [1 1 length(p) length(r) length(m) length(T) length(d)]),...
                             [3 1 4 5 6 7 2])./ndSparse(1+sub_DNA(pf,T,Tref,p,q,r,m,d,c,cExpressionMatrix).*...
                             permute(repmat(Keq_Complement(T,Tref,q,r,m),...
                             [1 1 1 1 length(p) length(d) length(c)]),[5 1 2 3 4 6 7]));                                        
    %output pf  probe,model,temp,dilution,cellline      
    %probe,target,site,model,temp, dilution, cellline  
    h_RNA = @(pf,T,Tref,p,q,r,m,d,c)  sum(sum(repmat(Keq_mod(T,Tref,p,q,r,m),[1 1 1 1 1 length(d) length(c)]).*...
                                       sub_RNA(pf,T,Tref,p,q,r,m,d,c),2,'omitnan'),3,'omitnan');                    
    h_DNA = @(pf,T,Tref,p,q,r,m,d,c)  sum(sum(repmat(Keq_mod(T,Tref,p,q,r,m),[1 1 1 1 1 length(d) length(c)]).*...
                                       sub_DNA(pf,T,Tref,p,q,r,m,d,c),2,'omitnan'),3,'omitnan'); 
    %probe,model,temp,dilution,cell-line
    K_S = @(T,Tref,p,n,m,d,c) repmat(squeeze(sum(Ks_eq(T,Tref,p,n,m),2)),[1 1 1 length(d) length(c)]);
    %probe,probe,model,temp
    K_CDeff = @(pf,T,Tref,p,q,n,m,d,c) multiprod(repmat(squeeze(sum(Kd_eq(T,Tref,p,q,n,m),3)),[1 1 1 1 length(d) length(c)]),sum(permute(repmat(pf,[ones(1,length(size(pf))) 2]),[1 1+length(size(pf)) 2:length(size(pf))]),2)/2);

%run solver for finding equilibrium
    vc = 1;
    A1 = squeeze(K_S(Tvec,Tref,Pset,1:Ns_Config,Mvec,Dvec,Cvec));
    A1(isnan(A1)) = 0;
    for i = 1:20
    x1 = squeeze(CProbes_Free);
    A2 = squeeze(K_CDeff(CProbes_Free,Tvec,Tref,Pset,Pset,1:Nc_Config,Mvec,Dvec,Cvec));
    A2(isnan(A2)) = 0;
    A3 = squeeze(h_RNA(CProbes_Free,Tvec,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),Mvec,Dvec,Cvec));
    A3(isnan(A3)) = 0;
    A4 = squeeze(h_DNA(CProbes_Free,Tvec,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),Mvec,Dvec,Cvec));%7-D
    A4(isnan(A4)) = 0;
    CProbes_Free = squeeze(ProbeConc)./(1+A1+A2+A3+A4); 
    xn = CProbes_Free;
    if (ndims(CProbes_Free0)-ndims(xn)>0)
        try
            Ai = arrayfun(@(x)find(size(CProbes_Free0)==x),size(xn));%problem was unique dims
        catch
            Ai = arrayfun(@(x)find(size(CProbes_Free0)==x),size(xn),'Un',0);
            Unique_Ai{1} = Ai{1};
            for ci = 1:length(Ai)-1
                K = length(Unique_Ai);
                if (sum(cell2mat(arrayfun(@(x) isequal(Ai{ci+1},Unique_Ai{x}),1:K,'Un',0)))>0)
                else
                   Unique_Ai{K+1} = Ai{ci+1};    
                end
            end
            lAi = cellfun(@length,Unique_Ai);
            ic = find(lAi>1);
            for kl = 1:length(ic)
                ind = find(cell2mat(arrayfun(@(x) isequal(Unique_Ai{ic(kl)},Ai{x}),1:length(Ai),'Un',0))==1);
                for icc = 1:length(ind)
                    Ai{ind(icc)} = Unique_Ai{ic(kl)}(icc);
                end
            end
            Ai = cell2mat(Ai);    
        end
    Bi_I = setdiff(1:ndims(CProbes_Free0),Ai);
    ReOrder = ones(1,ndims(CProbes_Free0));
    ReOrder(Ai) = 1:ndims(xn);
    ReOrder(Bi_I) = ndims(xn)+1:ndims(CProbes_Free0);
    CProbes_Free = sum(permute(repmat(CProbes_Free,[ones(1,ndims(x1)) 10*ones(1,ndims(CProbes_Free0)-ndims(x1))]),ReOrder),Bi_I)/10^(ndims(CProbes_Free0)-ndims(x1));
    end
    err(vc) = abs(sum(xn./x1,'all','omitnan')-numel(x1));
    vc = vc + 1;
    end    
%Initialize and store target binding probabilities for probe sets under all specified conditions 
p_TargetSites_Bound = ndSparse.build([length(Pset) size(DoesProbeBindSite,2) size(DoesProbeBindSite,3) length(Mvec) length(Tvec) length(Dvec) length(Cvec)],0);    
pRNA_TargetSites_Bound = ...
    (repmat(Keq_mod(Tvec,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),Mvec),[1 1 1 1 1 length(Dvec) length(Cvec)]).*...
    sub_RNA(CProbes_Free,Tvec,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),Mvec,Dvec,Cvec,cExpressionMatrix).*...
    permute(repmat(CProbes_Free,[1 1 1 1 1 length(Js_RNA(Pset)) length(Js_Sites(Pset))]),[1 6 7 2 3 4 5]))./...    
    (sum(repmat(Keq_mod(Tvec,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),Mvec),[1 1 1 1 1 length(Dvec) length(Cvec)]).*...
    sub_RNA(CProbes_Free,Tvec,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),Mvec,Dvec,Cvec,cExpressionMatrix).*...
    permute(repmat(CProbes_Free,[1 1 1 1 1 length(Js_RNA(Pset)) length(Js_Sites(Pset))]),[1 6 7 2 3 4 5]),1)+...
    sub_RNA(CProbes_Free,Tvec,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),Mvec,Dvec,Cvec,cExpressionMatrix));
pRNA_TargetSites_Bound(isnan(pRNA_TargetSites_Bound)) = 0;

p_TargetSites_Bound(:,Js_RNA(Pset),Js_Sites(Pset),:,:,:,:,:) = pRNA_TargetSites_Bound;
try
pDNA_TargetSites_Bound = ...
    repmat(Keq_mod(Tvec,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),Mvec),[1 1 1 1 1 length(Dvec) length(Cvec)]).*...
    sub_DNA(CProbes_Free,Tvec,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),Mvec,Dvec,Cvec,cExpressionMatrix).*...
    permute(repmat(CProbes_Free,[1 1 1 1 1 length(Js_DNA(Pset)) length(Js_Sites(Pset))]),[1 6 7 2 3 4 5])./...
    (sum(repmat(Keq_mod(Tvec,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),Mvec),[1 1 1 1 1 length(Dvec) length(Cvec)]).*...
    sub_DNA(CProbes_Free,Tvec,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),Mvec,Dvec,Cvec,cExpressionMatrix).*...
    permute(repmat(CProbes_Free,[1 1 1 1 1 length(Js_DNA(Pset)) length(Js_Sites(Pset))]),[1 6 7 2 3 4 5]),1)+...
    permute(repmat(Keq_Complement(Tvec,Tref,Js_DNA(Pset),Js_Sites(Pset),Mvec),[1 1 1 1 length(Pset) length(Dvec) length(Cvec)]),[5 1 2 3 4 6 7]).*...
    sub_DNA(CProbes_Free,Tvec,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),Mvec,Dvec,Cvec,cExpressionMatrix).*...
    sub_DNA_Complement(CProbes_Free,Tvec,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),Mvec,Dvec,Cvec,cExpressionMatrix)+...
    sub_DNA(CProbes_Free,Tvec,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),Mvec,Dvec,Cvec,cExpressionMatrix));
pDNA_TargetSites_Bound(isnan(pDNA_TargetSites_Bound)) = 0;
p_TargetSites_Bound(:,Js_DNA(Pset),Js_Sites(Pset),:,:,:,:,:) = pDNA_TargetSites_Bound;
catch
end
p_TargetSites_Bound(isnan(p_TargetSites_Bound)) = 0; 
c_TargetSites_Bound = ndSparse.build([length(Pset) size(DoesProbeBindSite,2) size(DoesProbeBindSite,3) length(Mvec) length(Tvec) length(Dvec) length(Cvec)],0);    
cRNA_TargetSites_Bound = ...
    repmat(Keq_mod(Tvec,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),Mvec),[1 1 1 1 1 length(Dvec) length(Cvec)]).*...
    sub_RNA(CProbes_Free,Tvec,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),Mvec,Dvec,Cvec,cExpressionMatrix).*...
    permute(repmat(CProbes_Free,[1 1 1 1 1 length(Js_RNA(Pset)) length(Js_Sites(Pset))]),[1 6 7 2 3 4 5]);
cRNA_TargetSites_Bound(isnan(cRNA_TargetSites_Bound)) = 0;
c_TargetSites_Bound(:,Js_RNA(Pset),Js_Sites(Pset),:,:,:,:,:) = cRNA_TargetSites_Bound;
try
cDNA_TargetSites_Bound = ...
    repmat(Keq_mod(Tvec,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),Mvec),[1 1 1 1 1 length(Dvec) length(Cvec)]).*...
    sub_DNA(CProbes_Free,Tvec,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),Mvec,Dvec,Cvec,cExpressionMatrix).*...
    permute(repmat(CProbes_Free,[1 1 1 1 1 length(Js_DNA(Pset)) length(Js_Sites(Pset))]),[1 6 7 2 3 4 5]); 
cDNA_TargetSites_Bound(isnan(cDNA_TargetSites_Bound)) = 0;
c_TargetSites_Bound(:,Js_DNA(Pset),Js_Sites(Pset),:,:,:,:,:) = cDNA_TargetSites_Bound;
catch
end
c_TargetSites_Bound(isnan(c_TargetSites_Bound)) = 0;

ModelMetrics.pTargetSites_Bound = p_TargetSites_Bound;
ModelMetrics.cTargetSites_Bound = c_TargetSites_Bound;
pAnyBindSite = sum(p_TargetSites_Bound,1);
Sx = unique(cell2mat(arrayfun(@(x) find(sum(DoesProbeBindSite(Pset,x,:),1)>0)',Js(Pset),'Un',0)));
Smax = length(Sx);%more smax for added to get to zero state min
Index2 = 1:Smax+1;
%ProbeSetMetrics is an CxM struct with 2 field

%Compute Confusion Matrix metrics.
for c=1:length(Cvec)
    for m = 1:length(Mvec)  
    %Index1 = m+(c-1)*length(Mvec); 
        for d = 1:length(Dvec) 
            for t = 1:length(Tvec)              
    pAnyBindSiteCM0 = squeeze(pAnyBindSite(1,:,:,m,t,d,c));   
    pAnyBindSiteCM0(pAnyBindSiteCM0>1) = 1;
    tHit = find(sum(pAnyBindSiteCM0,2)>0);
    if (~isempty(tHit))
    pON_IDs = find(ismember(tHit,ON_IDs));
    pOFF_IDs = find(ismember(tHit,OFF_IDs));
    pAnyBindSiteCM = pAnyBindSiteCM0(tHit,Sx);
    Pt = F_DiscretePossionBinomialMulti(pAnyBindSiteCM);
    Ct = repmat(nExpressionMatrix(tHit,Cvec(c)),[1 length(Sx)+1]).*Pt; 
    OnCounts = sum(Ct(pON_IDs,1:end),1);
    OnDistribution = OnCounts/sum(OnCounts);
    OffCounts = sum(Ct(pOFF_IDs,1:end),1);
    OffDistribution = OffCounts/sum(OffCounts);
    avgCon(c,m,d,t) = sum(OnCounts)/length(Pset);
    avgCoff(c,m,d,t) = sum(OffCounts)/length(Pset);
    Con(c,m,d,t) = sum(OnCounts);
    Coff(c,m,d,t) = sum(OffCounts);   
    FalseNegative_SpotCounts = arrayfun(@(d) sum(Ct(pON_IDs,1:d),'all'),0:Smax);
    TrueNegative_SpotCounts = arrayfun(@(d) sum(Ct(pOFF_IDs,1:d),'all'),0:Smax);
    FalsePositive_SpotCounts = arrayfun(@(d) sum(Ct(pOFF_IDs,d+1:Smax+1),'all'),0:Smax);
    TruePositive_SpotCounts = arrayfun(@(d) sum(Ct(pON_IDs,d+1:Smax+1),'all'),0:Smax);
    FalseNegative_ProbeCounts = arrayfun(@(d) sum(repmat(0:d-1,[length(pON_IDs) 1]).*Ct(pON_IDs,1:d),'all'),0:Smax);
    TrueNegative_ProbeCounts = arrayfun(@(d) sum(repmat(0:d-1,[length(pOFF_IDs) 1]).*Ct(pOFF_IDs,1:d),'all'),0:Smax);
    TruePositive_ProbeCounts = arrayfun(@(d) sum(repmat(d:Smax,[length(pON_IDs) 1]).*Ct(pON_IDs,d+1:Smax+1),'all'),0:Smax);
    FalsePositive_ProbeCounts = arrayfun(@(d) sum(repmat(d:Smax,[length(pOFF_IDs) 1]).*Ct(pOFF_IDs,d+1:Smax+1),'all'),0:Smax);
    True_Spots = TruePositive_SpotCounts + FalseNegative_SpotCounts;
    Removed_Spots = TrueNegative_SpotCounts + FalsePositive_SpotCounts;
    Measured_Spots = TruePositive_SpotCounts + FalsePositive_SpotCounts;
    Negative_Spots = TrueNegative_SpotCounts + FalseNegative_SpotCounts;  
    True_Probe = TruePositive_ProbeCounts + FalseNegative_ProbeCounts;
    Removed_Probe = TrueNegative_ProbeCounts + FalsePositive_ProbeCounts;
    Measured_Probe = TruePositive_ProbeCounts + FalsePositive_ProbeCounts;
    Negative_Probe = TrueNegative_ProbeCounts + FalseNegative_ProbeCounts; 
    Prevalence_Spots = Measured_Spots./(Measured_Spots+Negative_Spots);
    Accuracy_Spots = (TruePositive_SpotCounts+TrueNegative_SpotCounts)./(TruePositive_SpotCounts+TrueNegative_SpotCounts+FalsePositive_SpotCounts+FalseNegative_SpotCounts);
    BalancedAccuracy_Spots = 50*TruePositive_SpotCounts./True_Spots+50*TrueNegative_SpotCounts./Removed_Spots;
    F1_Spots = 2*TruePositive_SpotCounts./(2*TruePositive_SpotCounts+FalsePositive_SpotCounts+FalseNegative_SpotCounts);
    MCC_Spots = (TruePositive_SpotCounts.*TrueNegative_SpotCounts-FalsePositive_SpotCounts.*FalseNegative_SpotCounts)./sqrt(Measured_Spots.*True_Spots.*Removed_Spots.*Negative_Spots);
    Prevalence_Probe = Measured_Probe/(Measured_Probe+Negative_Probe);
    Accuracy_Probe = (TruePositive_ProbeCounts+TrueNegative_ProbeCounts)/(TruePositive_ProbeCounts+TrueNegative_ProbeCounts+FalsePositive_ProbeCounts+FalseNegative_ProbeCounts);
    BalancedAccuracy_Probe = 50*TruePositive_ProbeCounts./True_Probe+50*TrueNegative_ProbeCounts./Removed_Probe;
    F1_Probe = 2*TruePositive_ProbeCounts./(2*TruePositive_ProbeCounts+FalsePositive_ProbeCounts+FalseNegative_ProbeCounts);
    MCC_Probe = (TruePositive_ProbeCounts.*TrueNegative_ProbeCounts-FalsePositive_ProbeCounts.*FalseNegative_ProbeCounts)./sqrt(Measured_Probe.*True_Probe.*Removed_Probe.*Negative_Probe); 
    Spots.TPR(Index2) = 100*TruePositive_SpotCounts./True_Spots;
    Spots.FNR(Index2) = 100*FalseNegative_SpotCounts./True_Spots;   
    Spots.TNR(Index2)= 100*TrueNegative_SpotCounts./Removed_Spots;
    Spots.FPR(Index2) = 100*FalsePositive_SpotCounts./Removed_Spots;   
    Spots.PPV(Index2) = 100*TruePositive_SpotCounts./Measured_Spots;
    Spots.FDR(Index2) = 100*FalsePositive_SpotCounts./Measured_Spots;  
    Spots.NPV(Index2) = 100*TrueNegative_SpotCounts./Negative_Spots;
    Spots.FOR(Index2) = 100*FalseNegative_SpotCounts./Negative_Spots;
    Spots.TP(Index2) = TruePositive_SpotCounts;
    Spots.FN(Index2) = FalseNegative_SpotCounts;  
    Spots.TN(Index2) = TrueNegative_SpotCounts;
    Spots.FP(Index2) = FalsePositive_SpotCounts;   
    Spots.T(Index2) = True_Spots;
    Spots.F(Index2) = Removed_Spots;
    Spots.P(Index2) = Measured_Spots;
    Spots.N(Index2) = Negative_Spots;  
    Spots.PREV(Index2) = Prevalence_Spots;   
    Spots.ACC(Index2) = Accuracy_Spots;   
    Spots.BA(Index2) = BalancedAccuracy_Spots;   
    Spots.F1(Index2) = F1_Spots;   
    Spots.MCC(Index2) = MCC_Spots;  
    Probe.TPR(Index2) = 100*TruePositive_ProbeCounts./True_Probe;
    Probe.FNR(Index2) = 100*FalseNegative_ProbeCounts./True_Probe;   
    Probe.TNR(Index2) = 100*TrueNegative_ProbeCounts./Removed_Probe;
    Probe.FPR(Index2) = 100*FalsePositive_ProbeCounts./Removed_Probe;   
    Probe.PPV(Index2) = 100*TruePositive_ProbeCounts./Measured_Probe;
    Probe.FDR(Index2) = 100*FalsePositive_ProbeCounts./Measured_Probe;  
    Probe.NPV(Index2) = 100*TrueNegative_ProbeCounts./Negative_Probe;
    Probe.FOR(Index2) = 100*FalseNegative_ProbeCounts./Negative_Probe;
    Probe.TP(Index2) = TruePositive_ProbeCounts;
    Probe.FN(Index2) = FalseNegative_ProbeCounts;  
    Probe.TN(Index2) = TrueNegative_ProbeCounts;
    Probe.FP(Index2) = FalsePositive_ProbeCounts;   
    Probe.T(Index2) = True_Probe;
    Probe.F(Index2) = Removed_Probe;
    Probe.P(Index2) = Measured_Probe;
    Probe.N(Index2) = Negative_Probe;  
    Probe.PREV(Index2) = Prevalence_Probe;   
    Probe.ACC(Index2) = Accuracy_Probe;   
    Probe.BA(Index2) = BalancedAccuracy_Probe;   
    Probe.F1(Index2) = F1_Probe;   
    Probe.MCC(Index2) = MCC_Probe; 
    ProbeSetMetrics(c,m,d,t).fields = {'P','N','F','T','FP','FN','TP','TN','FOR','NPV','FDR','PPV','FPR','TPR','TNR','FNR','PREV','ACC','BA','F1','MCC'};   
    ProbeSetMetrics(c,m,d,t).Spots = Spots;
    ProbeSetMetrics(c,m,d,t).Probe = Probe; 
    ProbeSetMetrics(c,m,d,t).Pt = Pt; 
    ProbeSetMetrics(c,m,d,t).Ct = Ct;     
    ProbsTargets(c,m,d,t).x = 0:Smax;
    ProbsTargets(c,m,d,t).xnorm = (0:Smax)/length(Pset);
    ProbsTargets(c,m,d,t).Names = Names(tHit);
    ProbsTargets(c,m,d,t).nExpression = nExpressionMatrix(tHit,Cvec(c));
    ProbsTargets(c,m,d,t).Expresison = ExpressionMatrix(tHit,Cvec(c));
    CountsTargets(c,m,d,t).x = 0:Smax;
    CountsTargets(c,m,d,t).xnorm = (0:Smax)/length(Pset);
    CountsTargets(c,m,d,t).Names = Names(tHit);
    CountsTargets(c,m,d,t).nExpression = nExpressionMatrix(tHit,Cvec(c));
    CountsTargets(c,m,d,t).Expresison = ExpressionMatrix(tHit,Cvec(c));
    ProbsTargets(c,m,d,t).y = F_DiscretePossionBinomialMulti(pAnyBindSiteCM); 
    ProbsTargets(c,m,d,t).mean = sum((0:Smax).*Pt,2);
    ProbsTargets(c,m,d,t).var = sum((ProbsTargets(c,m,d,t).mean-(0:Smax)).^2.*Pt,2);
    ProbsTargets(c,m,d,t).fano = ProbsTargets(c,m,d,t).var./ProbsTargets(c,m,d,t).mean;
    ProbsTargets(c,m,d,t).normmean = sum(ProbsTargets(c,m,d,t).xnorm.*Pt,2);
    ProbsTargets(c,m,d,t).normvar = sum((ProbsTargets(c,m,d,t).normmean-ProbsTargets(c,m,d,t).xnorm).^2.*Pt,2);
    ProbsTargets(c,m,d,t).normfano = ProbsTargets(c,m,d,t).normvar./ProbsTargets(c,m,d,t).normmean;
    CountsTargets(c,m,d,t).y = Ct; 
    CountsTargets(c,m,d,t).mean = sum((0:Smax).*Ct,2)./sum(Ct,2);
    CountsTargets(c,m,d,t).var = sum((CountsTargets(c,m,d,t).mean-(0:Smax)).^2.*Ct,2)./sum(Ct,2);
    CountsTargets(c,m,d,t).fano = CountsTargets(c,m,d,t).var./CountsTargets(c,m,d,t).mean;
    CountsTargets(c,m,d,t).normmean = sum(CountsTargets(c,m,d,t).xnorm.*Ct,2)./sum(Ct,2);
    CountsTargets(c,m,d,t).normvar = sum((CountsTargets(c,m,d,t).normmean-CountsTargets(c,m,d,t).xnorm).^2.*Ct,2)./sum(Ct,2);
    CountsTargets(c,m,d,t).normfano = CountsTargets(c,m,d,t).normvar./CountsTargets(c,m,d,t).normmean;
    ProbsOn(c,m,d,t).x = 0:Smax;
    ProbsOn(c,m,d,t).xnorm = (0:Smax)/length(Pset);
    CountsOn(c,m,d,t).x = 0:Smax;
    CountsOn(c,m,d,t).xnorm = (0:Smax)/length(Pset);
    ProbsOff(c,m,d,t).x = 0:Smax;
    ProbsOff(c,m,d,t).xnorm = (0:Smax)/length(Pset);
    CountsOff(c,m,d,t).x = 0:Smax;
    CountsOff(c,m,d,t).xnorm = (0:Smax)/length(Pset);  
    ProbsOn(c,m,d,t).y = OnDistribution;
    ProbsOn(c,m,d,t).mean = sum((0:Smax).*OnDistribution);
    ProbsOn(c,m,d,t).var = sum((ProbsOn(c,m,d,t).mean-(0:Smax)).^2.*OnDistribution);
    ProbsOn(c,m,d,t).fano = ProbsOn(c,m,d,t).var/ProbsOn(c,m,d,t).mean;
    ProbsOn(c,m,d,t).normmean = sum(ProbsOn(c,m,d,t).xnorm.*OnDistribution);
    ProbsOn(c,m,d,t).normvar = sum((ProbsOn(c,m,d,t).normmean-ProbsOn(c,m,d,t).xnorm).^2.*OnDistribution);
    ProbsOn(c,m,d,t).normfano = ProbsOn(c,m,d,t).normvar/ProbsOn(c,m,d,t).normmean;
    CountsOn(c,m,d,t).y = OnCounts;
    CountsOn(c,m,d,t).mean = sum((0:Smax).*OnCounts)/sum(OnCounts);
    CountsOn(c,m,d,t).var = sum((CountsOn(c,m,d,t).mean-(0:Smax)).^2.*OnCounts)/sum(OnCounts);
    CountsOn(c,m,d,t).fano = CountsOn(c,m,d,t).var/CountsOn(c,m,d,t).mean;
    CountsOn(c,m,d,t).normmean = sum(CountsOn(c,m,d,t).xnorm.*OnCounts)/sum(OnCounts);
    CountsOn(c,m,d,t).normvar = sum((CountsOn(c,m,d,t).normmean-CountsOn(c,m,d,t).xnorm).^2.*OnCounts)/sum(OnCounts);
    CountsOn(c,m,d,t).normfano = CountsOn(c,m,d,t).normvar/CountsOn(c,m,d,t).normmean;
    ProbsOff(c,m,d,t).y = OffDistribution;
    ProbsOff(c,m,d,t).mean = sum((0:Smax).*OffDistribution);
    ProbsOff(c,m,d,t).var = sum((ProbsOff(c,m,d,t).mean-(0:Smax)).^2.*OffDistribution);
    ProbsOff(c,m,d,t).fano = ProbsOff(c,m,d,t).var/ProbsOff(c,m,d,t).mean;
    ProbsOff(c,m,d,t).normmean = sum(ProbsOff(c,m,d,t).xnorm.*OffDistribution);
    ProbsOff(c,m,d,t).normvar = sum((ProbsOff(c,m,d,t).normmean-ProbsOff(c,m,d,t).xnorm).^2.*OffDistribution);
    ProbsOff(c,m,d,t).normfano = ProbsOff(c,m,d,t).normvar/ProbsOff(c,m,d,t).normmean;
    CountsOff(c,m,d,t).y = OffCounts;
    CountsOff(c,m,d,t).mean = sum((0:Smax).*OffCounts)/sum(OffCounts);
    CountsOff(c,m,d,t).var = sum((CountsOff(c,m,d,t).mean-(0:Smax)).^2.*OffCounts)/sum(OffCounts);
    CountsOff(c,m,d,t).fano = CountsOff(c,m,d,t).var/CountsOff(c,m,d,t).mean;  
    CountsOff(c,m,d,t).normmean = sum(CountsOff(c,m,d,t).xnorm.*OffCounts);
    CountsOff(c,m,d,t).normvar = sum((CountsOff(c,m,d,t).normmean-CountsOff(c,m,d,t).xnorm).^2.*OffCounts)/sum(OffCounts);
    CountsOff(c,m,d,t).normfano = CountsOff(c,m,d,t).normvar/CountsOff(c,m,d,t).normmean;
    end
            end
        end
    end
end
LocMax = max(cell2mat(cellfun(@(x) x,{probes{:,3}},'UniformOutput',false)));  
Lpmin = min(cell2mat(cellfun(@length,{probes{:,2}},'UniformOutput',false)));
TargetLength = LocMax + Lpmin - 1;
theoryMaxProbes = floor(TargetLength/(Lpmin+settings.ProbeSpacing));
if (theoryMaxProbes>settings.maxProbes)
    theoryMaxProbes = settings.maxProbes;
end
PackEf = length(Pset)/theoryMaxProbes;

ModelMetrics.ProbeSetMetrics = ProbeSetMetrics;
ModelMetrics.PackEf = PackEf;
ModelMetrics.avgCon = avgCon;
ModelMetrics.avgCoff = avgCoff;
ModelMetrics.Con = Con;
ModelMetrics.Coff = Coff;
ModelMetrics.ProbsOn = ProbsOn;
ModelMetrics.CountsOn = CountsOn;
ModelMetrics.ProbsOff = ProbsOff;
ModelMetrics.CountsOff = CountsOff; 
ModelMetrics.ProbsTargets = ProbsTargets; 
ModelMetrics.CountsTargets = CountsTargets; 
AS = 0;
% save('ModelMetricsTest.mat','ModelMetrics','-v7.3');
if (AS)

N_Colors = 3;
for n = 1:N_Colors
ProbesInColor{n} = Pset(mod(1:length(Pset),N_Colors)==n-1);
end


pAnyBindSite_Color = cell(1,N_Colors);
Tx = find(sum(squeeze(sum(DoesProbeBindSite(Pset,:,:),1)),2)>0);
Tset = Tx;
Sx_Colors = cell(1,N_Colors);
Tx_Colors = cell(1,N_Colors);
Sx_OnlyColor = cell(1,N_Colors);
Tx_OnlyColor = cell(1,N_Colors);
Sx_ShareColors = cell(1,N_Colors);
Tx_ShareColors = cell(1,N_Colors);
for n = 1:N_Colors
Tx_Colors{n} = find(sum(squeeze(sum(DoesProbeBindSite(ProbesInColor{n},:,:),1)),2)>0); 
Sx_Colors{n} = unique(cell2mat(arrayfun(@(x) find(sum(DoesProbeBindSite(ProbesInColor{n},x,:),1)>0)',Js(ProbesInColor{n}),'Un',0)));
end
for n = 1:N_Colors
Tx_OnlyColor{n} = setdiff(Tx_Colors{n},cell2mat(arrayfun(@(x) Tx_Colors{x}',setdiff(1:N_Colors,n),'Un',0)));
Sx_OnlyColor{n} = setdiff(Sx_Colors{n},cell2mat(arrayfun(@(x) Sx_Colors{x}',setdiff(1:N_Colors,n),'Un',0)));
Tx_ShareColors{n} = setdiff(Tx_Colors{n},setdiff(Tx_Colors{n},cell2mat(arrayfun(@(x) Tx_Colors{x}',setdiff(1:N_Colors,n),'Un',0))));
Sx_ShareColors{n} = setdiff(Sx_Colors{n},setdiff(Sx_Colors{n},cell2mat(arrayfun(@(x) Sx_Colors{x}',setdiff(1:N_Colors,n),'Un',0))));
end
SxMC = max(cellfun(@length,Sx_Colors))
Tx_MultiColor = setdiff(Tx,cell2mat(arrayfun(@(x) Tx_OnlyColor{x}',1:N_Colors,'Un',0)));

% find(DoesProbeBindSite(ProbesInColor{n},Tx_MultiColor(x),:))
% intersect(find(DoesProbeBindSite(ProbesInColor{1},Tx_MultiColor(x),:)),find(DoesProbeBindSite(ProbesInColor{2},Tx_MultiColor(x),:)))
% find(DoesProbeBindSite(ProbesInColor{n},Tx_ShareColors{n}(x),:))
% distinct T, distinct S, distinct T/S
% F_DiscretePossionBinomial_Base2(Px2D,p);
% %N_Channels = length(unique(PColor2_Split));
% 
% find(ismember(Tset,Tx_MultiColor))

Ix_Split = cell(1,N_Colors);
vsi_Split = cell(1,N_Colors);
for k = 1:N_Colors
    Ix_Split{k} = ProbesInColor{k};      
    vsi_Split{k} = @(Tx) find(sum(squeeze(DoesProbeBindSite(Ix_Split{k},Tset(Tx),:)),1)>0);
end
vsi_Split_Mul = @(Tx) find(sum(CATnWrapper(arrayfun(@(k) double(sum(squeeze(DoesProbeBindSite(Ix_Split{k},Tset(Tx),:)),1)>0),1:N_Colors,'Un',0),1),1)>1);
vsi_Split_Uni = @(n,Tx) setdiff(vsi_Split{n}(Tx),cell2mat(arrayfun(@(k) vsi_Split{k}(Tx),setdiff(1:N_Colors,n),'Un',0)));

ns_Multi = zeros(length(Tset),1);
for a = 1:length(Tset)
    if (~isempty(vsi_Split_Mul(a)))
        ns_Multi(a) = length(vsi_Split_Mul(a));
    end
end
PxFD_Split = cell(1,N_Channels);
P_1FLAPs = [0 1];

for c=1:length(Cvec)
    for m = 1:length(Mvec)  
        for t = 1:length(Tvec)
            for d = 1:length(Dvec) 
                %tHit = find(sum(pAnyBindSiteCM0,2)>0);%replace Tset with tHit
                pAnyBindSiteCM0_Color = ...
                CATnWrapper(arrayfun(@(n) squeeze(P_1FLAPs(2)*sum(p_TargetSites_Bound(ismember(Pset,ProbesInColor{n}),:,:,Mvec(m),t,d,Cvec(c)),1)),1:N_Colors,'Un',0),3);  
                pAnyBindSiteCM_Color = ...
                CATnWrapper(arrayfun(@(n)cell2mat(arrayfun(@(J) [full(pAnyBindSiteCM0_Color(Tset(J),vsi_Split_Uni(n,J),n)) zeros(1,SxMC-length(vsi_Split_Uni(n,J)))]',1:length(Tset),'Un',0))'...
                ,1:N_Colors,'Un',0),3); 
                PxFD_Split = arrayfun(@(n) F_DiscretePossionBinomialMulti(squeeze(pAnyBindSiteCM_Color(:,:,n))),...
                1:N_Colors,'Un',0);
                PxFD_Split_Permute = arrayfun(@(n) permute(repmat(PxFD_Split{n},[1 1 (SxMC+1)*ones(1,N_Colors-1)]),[1 circshift(2:N_Colors+1,n-1)]),1:N_Colors,'Un',0);% permute T 1 n 1 ... 1
                E = 1;
                for k = 1:length(PxFD_Split_Permute)
                E = E.*PxFD_Split_Permute{k}; 
                end
                temp_Split = E;
                temp_Split2 = temp_Split./repmat(squeeze(sum(temp_Split,2:N_Colors+1)),[1 (SxMC+1)*ones(1,N_Colors)]);
                PxFnd_Split = temp_Split2;
                %need to update adding when same site has different color
                %bound
                if (sum(ns_Multi)>0)
                Split_Basis_Matrix = arrayfun(@(Tx) pAnyBindSiteCM_Color(Tx,vsi_Split_Mul(Tx),:),1:length(Tset),'Un',0);% third dimension basis
                PxFnd_Split = arrayfun(@(J) F_DiscretePossionBinomial_Base2(PxFnd_Split{J},Split_Basis_Matrix{J}),1:length(Tset),'Un',0);
                end
                PxFnd_Expr_Split = repmat(ExpressionMatrix(Tset,c),[1 (SxMC+1)*ones(1,N_Colors)]).*PxFnd_Split;
                %agnostic or not
                pnON_IDs = find(ismember(Tset,ON_IDs));
                pnOFF_IDs = find(ismember(Tset,OFF_IDs));                
                %every combination 1 to Sx_MC up to n times
                %https://www.mathworks.com/matlabcentral/answers/564101-accessing-the-nth-dimension-in-a-variable-sized-multidimensional-array
                subON = repmat({':'},1,ndims(PxFnd_Split));%can make variable number of : for nd
                subOFF = repmat({':'},1,ndims(PxFnd_Split));
                subON{1} = pnON_IDs;
                subOFF{1} = pnOFF_IDs;
                PxFnd_SplitOn = squeeze(sum(PxFnd_Split(subON{:}),1))/sum(PxFnd_Split(subON{:}),'all');
                PxFnd_SplitOff = squeeze(sum(PxFnd_Split(subOFF{:}),1))/sum(PxFnd_Split(subOFF{:}),'all');
                PxFnd_Expr_SplitOn = squeeze(sum(PxFnd_Split(subON{:}),1));
                PxFnd_Expr_SplitOff = squeeze(sum(PxFnd_Split(subOFF{:}),1));
                
                %joint probability any two dimensions
                combos = nchoosek(1:N_Colors,2);%dimension is the combo it is 
                PxFnd_SplitOnJoint.Combos = combos;
                PxFnd_SplitOnJoint.PDF =  CATnWrapper(arrayfun(@(k) squeeze(sum(PxFnd_SplitOn,setdiff(1:N_Colors,combos(k,:)))),1:size(combos,1),'Un',0),3);
                PxFnd_SplitOffJoint.Combos = combos;
                PxFnd_SplitOffJoint.PDF =  CATnWrapper(arrayfun(@(k) squeeze(sum(PxFnd_SplitOff,setdiff(1:N_Colors,combos(k,:)))),1:size(combos,1),'Un',0),3);
                PxFnd_Expr_SplitOnJoint.Combos = combos;
                PxFnd_Expr_SplitOnJoint.PDF =  CATnWrapper(arrayfun(@(k) squeeze(sum(PxFnd_Expr_SplitOn,setdiff(1:N_Colors,combos(k,:)))),1:size(combos,1),'Un',0),3);
                PxFnd_Expr_SplitOffJoint.Combos = combos;
                PxFnd_Expr_SplitOffJoint.PDF =  CATnWrapper(arrayfun(@(k) squeeze(sum(PxFnd_Expr_SplitOff,setdiff(1:N_Colors,combos(k,:)))),1:size(combos,1),'Un',0),3); 
                MultiColorTarget(c,m,d,t).PnD = PxFnd_Split;%Pt_Multi
                MultiColorTarget(c,m,d,t).ExprPnD = PxFnd_Expr_Split;%Ct_Multi
                    
                %check ability to distiguish isoforms (P(x) similarity)
%              if (length(Isoforms)>1)
%                  for u = 1:length(Isoforms)
%                     for v = 1:length(Isoforms)
%                     Sol(Ti,m,c).D_Isoforms(u,v) = 1/2*trapz(0:Smax,abs(Sol(Ti,m,c).PxF{Corr_ONRNA_IsoAgnostic(u)}-Sol(Ti,m,c).PxF{Corr_ONRNA_IsoAgnostic(v)}));
%                     Sol(Ti,m,c).D_Isoforms_Split(u,v) = 1/2*trapz(0:DimL-1,trapz(0:DimL-1,abs(Sol(Ti,m,c).PxF2D_Split{Corr_ONRNA_IsoAgnostic(u)}-Sol(Ti,m,c).PxF2D_Split{Corr_ONRNA_IsoAgnostic(v)})));
%                     Sol(Ti,m,c).D_Isoforms_Weave(u,v) = 1/2*trapz(0:DimL-1,trapz(0:DimL-1,abs(Sol(Ti,m,c).PxF2D_Weave{Corr_ONRNA_IsoAgnostic(u)}-Sol(Ti,m,c).PxF2D_Weave{Corr_ONRNA_IsoAgnostic(v)})));
%                     end
%                  end
%              else

xI = [0:500];
pd1 = makedist('Normal','mu',100,'sigma',10);
pd2 = makedist('Normal','mu',20,'sigma',10);
pd3 = makedist('Normal','mu',30,'sigma',10);
ypd1 = pdf(pd1,xI);
ypd2 = pdf(pd2,xI);
ypd3 = pdf(pd3,xI);
FiPD(1,1,:) = ypd1;
FiPD(2,1,:) = ypd2;
FiPD(3,1,:) = ypd3;
FiPD(1,2,:) = ypd2;
FiPD(2,2,:) = ypd1;
FiPD(3,2,:) = ypd3;
FiPD(1,3,:) = ypd2;
FiPD(2,3,:) = ypd3;
FiPD(3,3,:) = ypd1;
                %intensity
                %FiPD(k,C,:) Intensity of a fluorphores k in channel C
               %v vector of number of fluorophores of each type, given flurophore vector P(intensity in channel C)  
               pIntensity_Pnd = @(C) ifft(prod(CATnWrapper(arrayfun(@(k)permute(repmat(CATnWrapper(arrayfun(@(n) fft(squeeze(FiPD(k,C,:))).^n,0:SxMC,'Un',0),2),[1 1 (SxMC+1)*ones(1,size(FiPD,1)-1)]) ,[1 circshift(2:size(FiPD,1)+1,k-1)]),1:size(FiPD,1),'Un',0),size(FiPD,1)+2),size(FiPD,1)+2));% 
               Rep_Pnd = permute(repmat(Pnd,[ones(1,ndims(Pnd)) size(FiPD,3)]),[1 size(FiPD,1)+2 2:size(FiPD,1)+1]);
               pIntensity_In_C = @(C) sum(Rep_Pnd.*permute(repmat(pIntensity_Pnd(C),[ones(1,ndims(Pnd)) size(Pnd,1)]),[size(FiPD,1)+2 1:size(FiPD,1)+1]),[size(FiPD,1):size(FiPD,1)+2]);
               pndIntensity = prod(CATnWrapper(arrayfun(@(k)permute(repmat(pIntensity_In_C(k),[1 1 size(FiPD,3)*ones(1,size(FiPD,1)-1)]),[1 circshift(2:size(FiPD,1)+1,k-1)]),1:size(FiPD,2),'Un',0),size(FiPD,1)+2),size(FiPD,1)+2);        
               pndIntensity_Expr = repmat(ExpressionMatrix(Tset,c),[1 size(FiPD,3)*ones(1,size(FiPD,1))]).*pndIntensity;
               %probability specific intensity values in n channels as vector. 
               subON = repmat({':'},1,ndims(pndIntensity));%can make variable number of : for nd
               subOFF = repmat({':'},1,ndims(pndIntensity));
               pnON_IDs = find(ismember(Tset,ON_IDs));
               pnOFF_IDs = find(ismember(Tset,OFF_IDs));    
               subON{1} = pnON_IDs;  
               subOFF{1} = pnOFF_IDs; 
               pndIntensityOn = squeeze(sum(pndIntensity(subON{:}),1))/sum(pndIntensity(subON{:}),'all');
               pndIntensityOff = squeeze(sum(pndIntensity(subOFF{:}),1))/sum(pndIntensity(subOFF{:}),'all');
               pndIntensity_ExprOn = squeeze(sum(pndIntensity_Expr(subON{:}),1));
               pndIntensity_ExprOff = squeeze(sum(pndIntensity_Expr(subOFF{:}),1));
               %Nd thresholding positive and negative
               combos = nchoosek(1:size(FiPD,2),2);%dimension is the combo it is 
               pndIntensityOnJoint.Combos = combos;
               pndIntensityOnJoint.PDF =  CATnWrapper(arrayfun(@(k) squeeze(sum(pndIntensityOn,setdiff(1:size(FiPD,2),combos(k,:)))),1:size(combos,1),'Un',0),3);
               pndIntensityOffJoint.Combos = combos;
               pndIntensityOffJoint.PDF =  CATnWrapper(arrayfun(@(k) squeeze(sum(pndIntensityOff,setdiff(1:size(FiPD,2),combos(k,:)))),1:size(combos,1),'Un',0),3);
               pndIntensity_ExprOnJoint.Combos = combos;
               pndIntensity_ExprOnJoint.PDF =  CATnWrapper(arrayfun(@(k) squeeze(sum(pndIntensity_ExprOn,setdiff(1:size(FiPD,2),combos(k,:)))),1:size(combos,1),'Un',0),3);
               pndIntensity_ExprOffJoint.Combos = combos;
               pndIntensity_ExprOffJoint.PDF =  CATnWrapper(arrayfun(@(k) squeeze(sum(pndIntensity_ExprOff,setdiff(1:size(FiPD,2),combos(k,:)))),1:size(combos,1),'Un',0),3); 
               %threshold each channel independently.
               
            end
        end
    end
end

end


end