function [Nvec_RNAmulti,RNAOFF_Score,RNASpecificity_Score,NumRNAOffTargetOptions,Probes_WithNRNAOFF,DNAOFF_Score,DNASpecificity_Score,NumDNAOffTargetOptions,Probes_WithNDNAOFF,Cout] = A0_BasicDesignerStats(targetTypes,removeUndesiredIsos,gene_table,settings,FolderRootName,DoesProbeBindSite,Kon,Kb,Kb_Complement,EKernel)
% This function determine the metrics and statistics used for selecting RNA-FISH probes.
RNASpecificity_Score = zeros(1,size(DoesProbeBindSite,1));RNAOFF_Score = zeros(1,size(DoesProbeBindSite,1));
DNASpecificity_Score = zeros(1,size(DoesProbeBindSite,1));DNAOFF_Score = zeros(1,size(DoesProbeBindSite,1));
NumRNAOffTargetOptions = [];Probes_WithNRNAOFF = [];
NumDNAOffTargetOptions = [];Probes_WithNDNAOFF = [];
isRNA = targetTypes(1);isDNA = targetTypes(2);
AllowableProbes = 1:size(DoesProbeBindSite,1);
Cout = cell(1,2);
Cout{1} = cell(1,7);
Cout{2} = cell(1,9);
gene_table = sortrows(gene_table,[7 6],'ascend');
gene_table = gene_table(gene_table.Match>=settings.MinHomologySearchTargetSize,:);
MinusStrandedHits = find(contains(gene_table.Strand,'Minus'));
gene_table_NamesZ = convertCharsToStrings(gene_table.Name);
contains_RNA = find(ismember(gene_table_NamesZ,settings.RNAdbParser));
RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);
gene_table = gene_table(setdiff(1:size(gene_table,1),RNA_MissedFilteredHits),:);
gene_table.Ax = min(gene_table.SubjectIndices,[],2);
gene_table.Bx = max(gene_table.SubjectIndices,[],2);
gene_table = sortrows(gene_table,[7 13],'ascend');
Names = unique(gene_table.Name);
Names = convertCharsToStrings(Names);
if (and(strcmp(settings.referenceType,"ENSEMBL"),max(double(contains(extractBefore(Names,' '),'ENS')))==0))
    uniNames = extractBefore(Names,' ');
else
    uniNames = extractBefore(Names,'.');
    if (sum(ismissing(uniNames))>0)
        uniNames(ismissing(uniNames)) = extractBefore(Names(ismissing(uniNames)),' ');
    end
end
if (settings.BLASTdna)
DNA_IDs = find(~ismember(Names,settings.DNAdbParser));%IDs
else
DNA_IDs = [];
end
if (settings.BLASTrna)
NonDNA_IDs = find(ismember(Names,settings.RNAdbParser));%IDs
else
NonDNA_IDs =[];
end
ON_RNAIDs = find(ismember(uniNames,extractBefore(settings.transcript_IDs_desired{:},'.')));
OFF_RNAIDs = setdiff(NonDNA_IDs,ON_RNAIDs);
ON_RNAIDs_Isos =  find(ismember(uniNames,extractBefore(settings.transcript_IDs_desired{:},'.')));
Desired_Isoforms =  find(ismember(uniNames,extractBefore(settings.transcript_IDs{:},'.')));
UnDesired_Isoforms = setdiff(ON_RNAIDs_Isos,Desired_Isoforms);
OFF_RNAIDs_minusIsos = setdiff(OFF_RNAIDs,UnDesired_Isoforms);
if (removeUndesiredIsos)
    OFF_RNAIDs = OFF_RNAIDs_minusIsos;
end
%Finds Off-targets and off-target binding sites
Js = @(x) find(sum(squeeze(sum(DoesProbeBindSite(x,:,:),1)),2)>0);
Tp = @(x) find(sum(squeeze(DoesProbeBindSite(:,x,:)),2)>0);
Tx2 =@(y,Z) arrayfun(@(x) find(squeeze(DoesProbeBindSite(x,y,:))==1)',Z,'Un',0);
Sx =@(x,Z) arrayfun(@(y) find(squeeze(DoesProbeBindSite(x,y,:))==1),Z,'Un',0);
Tx =@(y,Z) arrayfun(@(x) find(squeeze(DoesProbeBindSite(x,y,:))==1),Z,'Un',0);
Js_OFFRNA = @(x)OFF_RNAIDs(ismember(OFF_RNAIDs,Js(x)));
Js_OFFRNAi = @(x,y)OFF_RNAIDs(find(cumsum(ismember(OFF_RNAIDs,Js(x)))==y,1));
Js_OFFDNA = @(x)DNA_IDs(ismember(DNA_IDs,Js(x)));
Js_OFFDNAi = @(x,y)DNA_IDs(find(cumsum(ismember(DNA_IDs,Js(x)))==y,1));
% finds list of each off-target both number, site and location, as well as Koff and Kon equilibrium constants
if (isRNA)
    TPvec_RNA0 = cell(1,size(DoesProbeBindSite,1));
    Tvec_RNA = cell(1,size(DoesProbeBindSite,1));
    TSvec_RNA = cell(1,size(DoesProbeBindSite,1));
    TPvec_RNA = cell(1,size(DoesProbeBindSite,1));
    Svec_RNA = cell(1,size(DoesProbeBindSite,1));
    Nvec_RNAsingle = zeros(1,size(DoesProbeBindSite,1));
    Nvec_RNAmulti = zeros(1,size(DoesProbeBindSite,1));
    NTPvec_RNAsingle = zeros(1,size(DoesProbeBindSite,1));
    NTPvec_RNAmulti = zeros(1,size(DoesProbeBindSite,1));
    Tvec_logKOFF_RNA = cell(1,size(DoesProbeBindSite,1));
    Tvec_logKOFFdivON_RNA = cell(1,size(DoesProbeBindSite,1));
    Tvec_logKONdivOFF_RNA = cell(1,size(DoesProbeBindSite,1));
    TPvec_logKOFF_RNA = cell(1,length(OFF_RNAIDs));
    TPvec_logKOFFdivON_RNA = cell(1,length(OFF_RNAIDs));
    TPvec_logKONdivOFF_RNA = cell(1,length(OFF_RNAIDs));
    fprintf("Generating RNA target statistics by probe")
    fprintf('\n')
    fprintf('\n')
    wb = parwaitbar(size(DoesProbeBindSite,1),'WaitMessage','Computing');
    for p = 1:size(DoesProbeBindSite,1)
        Nvec_RNAsingle(p) = length(Js_OFFRNA(p));
        Nvec_RNAmulti(p) = sum(cellfun(@length,Sx(p,Js_OFFRNA(p))));
        Svec_RNA{p} = cell2mat(Sx(p,Js_OFFRNA(p)));%size(N by 1)
        Tvec_RNA{p} = cell2mat(arrayfun(@(x) repmat(Js_OFFRNAi(p,x),[1 cellfun(@length,Sx(p,Js_OFFRNAi(p,x)))]),1:length(Js_OFFRNA(p)),'Un',0));% 1 by N
        Tvec_logKOFF_RNA{p} = log10(diag(full(squeeze(Kb(p,Tvec_RNA{p},Svec_RNA{p}))))');
        Tvec_logKOFFdivON_RNA{p} = log10(diag(full(squeeze(Kb(p,Tvec_RNA{p},Svec_RNA{p}))))'/Kon(p));
        Tvec_logKONdivOFF_RNA{p} = log10(Kon(p)./diag(full(squeeze(Kb(p,Tvec_RNA{p},Svec_RNA{p}))))');
        progress(wb);
    end
    wb.delete();
    fprintf('\n')
    fprintf('\n')
    fprintf("Generating probe statistics by RNA off-targets")
    fprintf('\n')
    fprintf('\n')
    wb = parwaitbar(length(OFF_RNAIDs),'WaitMessage','Computing');
    for t = 1:length(OFF_RNAIDs)
        TPvec_RNA0{t} = Tp(OFF_RNAIDs(t));%does not have multiplicity
        temp_RNAV1b = Tx(OFF_RNAIDs(t),TPvec_RNA0{t});
        temp_RNAV2b = cellfun(@length,temp_RNAV1b);
        NTPvec_RNAsingle(t) = length(TPvec_RNA0{t});
        NTPvec_RNAmulti(t) = sum(temp_RNAV2b);
        TPvec_RNA{t} = cell2mat(arrayfun(@(x) repmat(TPvec_RNA0{t}(x),[1 temp_RNAV2b(x)]),1:length(TPvec_RNA0{t}),'Un',0));
        TSvec_RNA{t} = cell2mat(temp_RNAV1b);%site locations
        TPvec_logKOFF_RNA{t} = log10(diag(full(squeeze(Kb(TPvec_RNA{t},OFF_RNAIDs(t),TSvec_RNA{t}))))');
        TPvec_logKOFFdivON_RNA{t} = log10(diag(full(squeeze(Kb(TPvec_RNA{t},OFF_RNAIDs(t),TSvec_RNA{t}))))'./Kon(TPvec_RNA{t}));
        TPvec_logKONdivOFF_RNA{t} = log10(Kon(TPvec_RNA{t})./diag(full(squeeze(Kb(TPvec_RNA{t},OFF_RNAIDs(t),TSvec_RNA{t}))))');
        progress(wb);
    end
    wb.delete();
    Cout{1}{1} = Tvec_RNA;
    Cout{1}{2} = Svec_RNA;
    Cout{1}{3} = TPvec_RNA;
    Cout{1}{4} = TSvec_RNA;
    Cout{1}{5} = TPvec_logKOFF_RNA;
    Cout{1}{6} = TPvec_logKOFFdivON_RNA;
    Cout{1}{7} = TPvec_logKONdivOFF_RNA;
    Probes_With_RNAOFF = AllowableProbes(Nvec_RNAmulti(AllowableProbes)>0);
    NumRNAOffTargetOptions = unique(Nvec_RNAmulti(AllowableProbes));
    Probes_WithNRNAOFF = arrayfun(@(i) AllowableProbes(Nvec_RNAmulti(AllowableProbes)==NumRNAOffTargetOptions(i)),1:length(NumRNAOffTargetOptions),'Un',0);
    fprintf("Converting RNA Statistics into Probe Specificity and OFF-target Scores")
    fprintf('\n')
    fprintf('\n')
    wb = parwaitbar(length(Probes_With_RNAOFF),'WaitMessage','Converting');
    for v = 1:length(Probes_With_RNAOFF)
        temp_T = Tvec_RNA{Probes_With_RNAOFF(v)};
        temp_KOFF = Tvec_logKOFF_RNA{Probes_With_RNAOFF(v)};
        temp_KOFFdivON = Tvec_logKOFFdivON_RNA{Probes_With_RNAOFF(v)};
        RNASpecificity_Score(Probes_With_RNAOFF(v)) = dot(EKernel(temp_T)',temp_KOFFdivON);
        RNAOFF_Score(Probes_With_RNAOFF(v)) = dot(EKernel(temp_T)',temp_KOFF);
        progress(wb);
    end
    wb.delete();
    fprintf('\n')
    fprintf('\n')
end
if (isDNA)
    TPvec_DNA0 = cell(1,size(DoesProbeBindSite,1));
    Tvec_DNA = cell(1,size(DoesProbeBindSite,1));
    TSvec_DNA = cell(1,size(DoesProbeBindSite,1));
    TPvec_DNA = cell(1,size(DoesProbeBindSite,1));
    Svec_DNA = cell(1,size(DoesProbeBindSite,1));
    Nvec_DNAsingle = zeros(1,size(DoesProbeBindSite,1));
    Nvec_DNAmulti = zeros(1,size(DoesProbeBindSite,1));
    NTPvec_DNAsingle = zeros(1,size(DoesProbeBindSite,1));
    NTPvec_DNAmulti = zeros(1,size(DoesProbeBindSite,1));
    Tvec_logKOFF_DNA = cell(1,size(DoesProbeBindSite,1));
    Tvec_logKOFFdivON_DNA = cell(1,size(DoesProbeBindSite,1));
    Tvec_logKONdivOFF_DNA = cell(1,size(DoesProbeBindSite,1));
    Tvec_logKOFFdivCOMP_DNA = cell(1,size(DoesProbeBindSite,1));
    Tvec_logKCOMPdivOFF_DNA = cell(1,size(DoesProbeBindSite,1));
    TPvec_logKOFF_DNA = cell(1,length(DNA_IDs));
    TPvec_logKOFFdivON_DNA = cell(1,length(DNA_IDs));
    TPvec_logKONdivOFF_DNA = cell(1,length(DNA_IDs));
    TPvec_logKOFFdivCOMP_DNA = cell(1,length(DNA_IDs));
    TPvec_logKCOMPdivOFF_DNA = cell(1,length(DNA_IDs));
    fprintf("Generating DNA target statistics by probe")
    fprintf('\n')
    fprintf('\n')
    wb = parwaitbar(size(DoesProbeBindSite,1),'WaitMessage','Computing');
    for p = 1:size(DoesProbeBindSite,1)
        Nvec_DNAsingle(p) = length(Js_OFFDNA(p));
        Nvec_DNAmulti(p) = sum(cellfun(@length,Sx(p,Js_OFFDNA(p))));
        Svec_DNA{p} =  cell2mat(Sx(p,Js_OFFDNA(p)));
        Tvec_DNA{p} = cell2mat(arrayfun(@(x) repmat(Js_OFFDNAi(p,x),[1 cellfun(@length,Sx(p,Js_OFFDNAi(p,x)))]),1:length(Js_OFFDNA(p)),'Un',0));
        Tvec_logKOFF_DNA{p} = log10(diag(full(squeeze(Kb(p,Tvec_DNA{p},Svec_DNA{p}))))');
        Tvec_logKOFFdivON_DNA{p} = log10(diag(full(squeeze(Kb(p,Tvec_DNA{p},Svec_DNA{p}))))'/Kon(p));
        Tvec_logKONdivOFF_DNA{p} = log10(Kon(p)./diag(full(squeeze(Kb(p,Tvec_DNA{p},Svec_DNA{p}))))');
        Tvec_logKOFFdivCOMP_DNA{p} = log10(diag(full(squeeze(Kb(p,Tvec_DNA{p},Svec_DNA{p}))))'./diag(full(squeeze(Kb_Complement(Tvec_DNA{p},Svec_DNA{p}))))');
        Tvec_logKCOMPdivOFF_DNA{p} = log10(diag(full(squeeze(Kb_Complement(Tvec_DNA{p},Svec_DNA{p}))))'./diag(full(squeeze(Kb(p,Tvec_DNA{p},Svec_DNA{p}))))');
        progress(wb);
    end
    wb.delete();
    fprintf('\n')
    fprintf('\n')
    fprintf("Generating probe statistics by DNA off-targets")
    fprintf('\n')
    fprintf('\n')
    wb = parwaitbar(length(DNA_IDs),'WaitMessage','Computing');
    for t = 1:length(DNA_IDs)
        TPvec_DNA0{t} = Tp(DNA_IDs(t))';
        temp_DNAV1b = Tx2(DNA_IDs(t),TPvec_DNA0{t});
        temp_DNAV2b = cellfun(@length,temp_DNAV1b);
        NTPvec_DNAsingle(t) = length(TPvec_DNA0{t});
        NTPvec_DNAmulti(t) = sum(temp_DNAV2b);
        TPvec_DNA{t} = cell2mat(arrayfun(@(x) repmat(TPvec_DNA0{t}(x),[1 temp_DNAV2b(x)]),1:length(TPvec_DNA0{t}),'Un',0));
        TSvec_DNA{t} = cell2mat(temp_DNAV1b);%site locations
        TPvec_logKOFF_DNA{t} = log10(diag(full(squeeze(Kb(TPvec_DNA{t},DNA_IDs(t),TSvec_DNA{t}))))');
        TPvec_logKOFFdivON_DNA{t} = log10(diag(full(squeeze(Kb(TPvec_DNA{t},DNA_IDs(t),TSvec_DNA{t}))))'./Kon(TPvec_DNA{t}));
        TPvec_logKONdivOFF_DNA{t} = log10(Kon(TPvec_DNA{t})./diag(full(squeeze(Kb(TPvec_DNA{t},DNA_IDs(t),TSvec_DNA{t}))))');
        TPvec_logKOFFdivCOMP_DNA{t} = log10(diag(full(squeeze(Kb(TPvec_DNA{t},DNA_IDs(t),TSvec_DNA{t}))))'./full(Kb_Complement(DNA_IDs(t),TSvec_DNA{t})));
        TPvec_logKCOMPdivOFF_DNA{t} = log10(full(Kb_Complement(DNA_IDs(t),TSvec_DNA{t}))./diag(full(squeeze(Kb(TPvec_DNA{t},DNA_IDs(t),TSvec_DNA{t}))))');
        progress(wb);
    end
    wb.delete();
    fprintf('\n')
    fprintf('\n')
    Cout{2}{1} = Tvec_DNA;
    Cout{2}{2} = Svec_DNA;
    Cout{2}{3} = TPvec_DNA;
    Cout{2}{4} = TSvec_DNA;
    Cout{2}{5} = TPvec_logKOFF_DNA;
    Cout{2}{6} = TPvec_logKOFFdivON_DNA;
    Cout{2}{7} = TPvec_logKONdivOFF_DNA;
    Cout{2}{8} = TPvec_logKOFFdivCOMP_DNA;
    Cout{2}{9} = TPvec_logKCOMPdivOFF_DNA;
    Probes_With_DNAOFF = AllowableProbes(Nvec_DNAmulti(AllowableProbes)>0);
    NumDNAOffTargetOptions = unique(Nvec_DNAmulti(AllowableProbes));
    Probes_WithNDNAOFF = arrayfun(@(i) AllowableProbes(Nvec_DNAmulti(AllowableProbes)==NumDNAOffTargetOptions(i)),1:length(NumDNAOffTargetOptions),'Un',0);
    fprintf("Converting DNA Statistics into Probe Specificity and OFF-target Scores")
    fprintf('\n')
    fprintf('\n')
    wb = parwaitbar(length(Probes_With_DNAOFF),'WaitMessage','Converting');
    for v = 1:length(Probes_With_DNAOFF)
        temp_T = Tvec_DNA{Probes_With_DNAOFF(v)};
        temp_KOFF = Tvec_logKOFF_DNA{Probes_With_DNAOFF(v)};
        temp_KOFFdivON = Tvec_logKOFFdivON_DNA{Probes_With_DNAOFF(v)};
        DNASpecificity_Score(Probes_With_DNAOFF(v)) = dot(EKernel(temp_T)',temp_KOFFdivON);
        DNAOFF_Score(Probes_With_DNAOFF(v)) = dot(EKernel(temp_T)',temp_KOFF);
        progress(wb);
    end
    wb.delete();
    fprintf('\n')
    fprintf('\n')
end
end