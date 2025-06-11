function ModelMetrics = RNAsolver_JH2(Pset,settings,probes,gene_table,ExpressionMatrix,DoesProbeBindSite,dHeq_mod,dSeq_mod,dCp_mod,dHeq_Complement,dSeq_Complement,dCp_Complement)
% This function computes the metrics for RNA-FISH performance.
% The function simulates probe binding dynamics to model probe
% equilibrium binding under a variety of conditions. Varying by
% temperature, expression, probe concentration, etc.
% The distribution of probes bound targets is computed and is
% used to generate confusion matrix metrics, of probe performance.

%% Handle Functions
H2_func = @(z,x) heaviside(z-x');
plus_Random = @(x,y) conv(x,y)/sum(conv(x,y));
minus_Random = @(x,y) conv(x,flip(y))/sum(conv(x,flip(y)));
pIntensity = @(Istep,Mu,Std) pdf('Normal',Istep,Mu,Std)/sum(pdf('Normal',Istep,Mu,Std));
nProbeIntensity = @(Istep,Mu,Std,Nr,x) pdf('Normal',Istep,Mu/Nr*x,Std/Nr*x)/sum(pdf('Normal',Istep,Mu/Nr*x,Std/Nr*x));
V_Cell = @(R) 4/3*pi*(R^3)/10^15;%um to L
z_domain_function = @(z,Pz) z(Pz>0);
Pz_domain_function = @(Pz) Pz(Pz>0)/sum(Pz(Pz>0));

Tvec = settings.SimulationConfiguration.Temperature_Celsius_Model_Vector+273.15;
Mvec = settings.SimulationConfiguration.Gibbs_Model_Vector;
Dvec = settings.SimulationConfiguration.Dilution_Vector;
AutoBackground_Mean = settings.SimulationConfiguration.AutoBackground_Mean;
AutoBackground_STD = settings.SimulationConfiguration.AutoBackground_STD;
NumReferenceProbes = settings.SimulationConfiguration.NumReferenceProbes;
SpotIntensity_Mean = settings.SimulationConfiguration.SpotIntensity_Mean;
SpotIntensity_STD = settings.SimulationConfiguration.SpotIntensity_STD;
Tref = settings.SimulationConfiguration.Tref+273.15;
Mean_Diameter = settings.SimulationConfiguration.Mean_Diameter;
Rcell = settings.SimulationConfiguration.CellRadius;
Rspot = settings.SimulationConfiguration.SpotRadius;
Nstacks = settings.SimulationConfiguration.NumOfZStacks;
GuessConc = settings.SimulationConfiguration.InitialGuessConcc;
PC0 = settings.SimulationConfiguration.ProbeConcentration;
errThreshold = settings.SimulationConfiguration.errThreshold;
MaxIter = settings.SimulationConfiguration.MaxIter;
Diameter_vals = Mean_Diameter;
load('data/TS_DefaultParams.mat','TrueSpotDefaultThParameters')
SI = 0:0.1:3000;%only works if SI
SI_SignalMinusBackgd_I = -max(SI):0.1:1*max(SI);
SI_Signal_wAuto_I = 0:0.1:2*max(SI);
IntensityPerProbe = SpotIntensity_Mean/NumReferenceProbes;
Pr_Auto_I = pIntensity(SI,AutoBackground_Mean,AutoBackground_STD);%I over SI
NumNonUniformConditions = size(ExpressionMatrix,2);
[~, ExpressionMatrix_nTPM] = tmm(ExpressionMatrix,0.3,0.05,-1e10,1);
transcript_Expression_Equal_acrossSample = mean(ExpressionMatrix_nTPM(:));
transcript_Expression_Average_acrossSample = mean(ExpressionMatrix_nTPM,2)';
ExpressionMatrix_nTPM(:,NumNonUniformConditions+1) = transcript_Expression_Equal_acrossSample;
ExpressionMatrix_nTPM(:,NumNonUniformConditions+2) = transcript_Expression_Average_acrossSample;
nExpressionMatrix = ndSparse(ExpressionMatrix_nTPM/(V_Cell(Rcell)*6.022*10^23));
if (settings.ExpressionReferenceForDesigningProbes==0)
Cvec = [size(nExpressionMatrix,2)-1 size(nExpressionMatrix,2)];
else
Cvec = [settings.ExpressionReferenceForDesigningProbes size(nExpressionMatrix,2)-1 size(nExpressionMatrix,2)];
end
GeneName = settings.GeneName;
GeneTarget = settings.transcript_IDs;
gene_table = sortrows(gene_table,[7 6],'ascend');
gene_table = gene_table(gene_table.Match>=settings.MinHomologySearchTargetSize,:);
MinusStrandedHits = find(contains(gene_table.Strand,'Minus'));
if (strcmp(settings.referenceType,'RefSeq'))
    RNA_IDs_1 = find(contains(gene_table.Name,'NM_'));
    RNA_IDs_2 = find(contains(gene_table.Name,'NR_'));
    RNA_IDs_3 = find(contains(gene_table.Name,'XM_'));
    RNA_IDs_4 = find(contains(gene_table.Name,'XR_'));
    contains_RNA = union(union(union(RNA_IDs_1,RNA_IDs_2),RNA_IDs_3),RNA_IDs_4);
elseif (strcmp(settings.referenceType,'ENSEMBL'))
    contains_RNA = find(contains(gene_table.Name,settings.EMBL_RNAparser(Organism)));
end
RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);
gene_table = gene_table(setdiff(1:size(gene_table,1),RNA_MissedFilteredHits),:);
gene_table.Ax = min(gene_table.SubjectIndices,[],2);
gene_table.Bx = max(gene_table.SubjectIndices,[],2);
gene_table = sortrows(gene_table,[7 13],'ascend');
Names = unique(gene_table.Name);
Names = convertCharsToStrings(Names);
uniNames = extractBefore(Names,'.');
if (sum(ismissing(uniNames))>0)
    uniNames(ismissing(uniNames)) = extractBefore(Names(ismissing(uniNames)),' ');
end
ON_IDs_specific = find(contains(uniNames,extractBefore(GeneTarget{:},'.')));
ON_IDs_agnostic = find(contains(Names,GeneName));
OFF_IDs = find(~contains(Names,GeneName));
if (strcmp(settings.referenceType,'RefSeq'))
DNA_IDs_1 = find(contains(uniNames,'NC_'));%IDs
DNA_IDs_2 = find(contains(uniNames,'NT_'));%IDs
DNA_IDs_3 = find(contains(uniNames,'NW_'));%IDs
NonDNA_IDs_1 = find(~contains(uniNames,'NC_'));%IDs
NonDNA_IDs_2 = find(~contains(uniNames,'NT_'));%IDs
NonDNA_IDs_3 = find(~contains(uniNames,'NW_'));%IDs
DNA_IDs =union(union(DNA_IDs_1,DNA_IDs_2),DNA_IDs_3).';
NonDNA_IDs = intersect(intersect(NonDNA_IDs_1,NonDNA_IDs_2),NonDNA_IDs_3).';
elseif (strcmp(settings.referenceType,'ENSEMBL'))
DNA_IDs = find(~contains(uniNames,settings.EMBL_RNAparser(Organism)));%IDs
NonDNA_IDs = find(contains(uniNames,settings.EMBL_RNAparser(Organism)));%IDs
end
if (ndims(dHeq_Complement)~=3)%error quick fix
    dHeq_Complement =  permute(repmat(dHeq_Complement, [1 1 size(dHeq_mod,3)]),[1 3 2]);
    dSeq_Complement = permute(repmat(dSeq_Complement, [1 1 size(dHeq_mod,3)]),[1 3 2]);
    dCp_Complement = permute(repmat(dCp_Complement, [1 1 size(dHeq_mod,3)]),[1 3 2]);
end
LocMax = max(cell2mat(cellfun(@(x) x,{probes{:,3}},'UniformOutput',false)));
Lpmin = min(cell2mat(cellfun(@length,{probes{:,2}},'UniformOutput',false)));
TargetLength = LocMax + Lpmin - 1;
theoryMaxProbes = floor(TargetLength/(Lpmin+settings.ProbeSpacing));
if (theoryMaxProbes>settings.maxProbes)
    theoryMaxProbes = settings.maxProbes;
end
PackEf = length(Pset)/theoryMaxProbes;
ModelMetrics.PackingEfficiency = PackEf;
[Ns_Config,Nc_Config,Js_RNA,Js_DNA,Js_Sites,linearIndexed,MultiDim_PJSMC,MultiDim_PJSMTDC]  = ...
    A_ModelSolverWrapper_V4(probes,Pset,settings,DoesProbeBindSite,DNA_IDs,NonDNA_IDs,dCp_mod,dHeq_mod,dSeq_mod,dCp_Complement,dSeq_Complement,dHeq_Complement);
ProbeSetMetrics.Ns_Config = Ns_Config;
ProbeSetMetrics.Nc_Config = Nc_Config;
ProbeSetMetrics.Js_RNA = Js_RNA;
ProbeSetMetrics.Js_DNA = Js_DNA;
ProbeSetMetrics.Js_Sites = Js_Sites;
ProbeSetMetrics.ModelSolverFunctions_7D = MultiDim_PJSMTDC;
ProbeSetMetrics.ModelSolverFunctions_5D = MultiDim_PJSMC;
ProbeSetMetrics.ModelSolverFunctions_linIndex = linearIndexed;
clear dCp_mod dSeq_mod dHeq_mod dCp_Complement dSeq_Complement dHeq_Complement
[~,m_unique_loc,~] = unique(Mvec);
[~,t_unique_loc,~] = unique(Tvec);
[~,d_unique_loc,~] = unique(Dvec);
ProbeSetMetrics.SolutionModels = Mvec(m_unique_loc);
ProbeSetMetrics.SolutionTemperatures = Tvec(t_unique_loc);
ProbeSetMetrics.SolutionDilutions = Dvec(d_unique_loc);
ProbeSetMetrics.SolutionCellLines = Cvec;
ProbeSetMetrics.iter = zeros(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.err = zeros(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc),length(Cvec));
ProbeSetMetrics.varSSE = zeros(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc),length(Cvec));
ProbeSetMetrics.eqSSE = zeros(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc),length(Cvec));
ProbeSetMetrics.BindingPredictions.CProbes_Free = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.p_TargetSites_Bound_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.c_TargetSites_Bound_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.c_OnOtherOff_nBound_ModelTemperatureDilutionVector =  cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.Con_Distribution_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.Coff_Distribution_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.Cother_Distribution_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.Pon_Distribution_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.Poff_Distribution_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.Pother_Distribution_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.Non_Count_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.Noff_Count_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.Nother_Count_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.Non_history_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.Nother_history_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.Noff_history_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.Noff_tot_history_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.IsoIgnorantConfusion_Probe_P_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.IsoSpecificConfusion_Probe_P_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.IsoAgnosticConfusion_Probe_P_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.IntensityPredictions.Ioff_history_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.IntensityPredictions.Pr_Non_I_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.IntensityPredictions.Pr_Nother_I_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.IntensityPredictions.Pr_Noff_I_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.IntensityPredictions.Non_I_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.IntensityPredictions.Nother_I_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.IntensityPredictions.Noff_I_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.IntensityPredictions.QzSignal_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.IntensityPredictions.PzSignal_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.IntensityPredictions.QzCellBkg_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.IntensityPredictions.PzCellBkg_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.IntensityPredictions.QzSignalMinusBackgd_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.IntensityPredictions.PzSignalMinusBackgd_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.CountPredictions.IsoAgnostic_SpotCountMetrics_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.CountPredictions.IsoSpecific_SpotCountMetrics_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.CountPredictions.IsoAgnostic_SpotCountMetrics_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ModelSolverFunctions = linearIndexed;
for m_unique_loci = 1:length(m_unique_loc)
    for t_unique_loci = 1:length(t_unique_loc)
        for d_unique_loci = 1:length(d_unique_loc)
            CProbes_Free = GuessConc*ones(length(Pset),length(Mvec),length(Cvec));
            CProbes_Free0 = GuessConc*ones(length(Pset),length(Mvec),length(Cvec));
            ProbeConc = PC0*squeeze(permute(repmat(Dvec(d_unique_loc(d_unique_loci)),[1 1 length(Pset) 1 1 length(Cvec)]),[1 3 4 5 2 6]));
            [CProbes_Free,varSSE,err,iter,eqSSE] = A_ModelEquilibriumSolverWrapper_V4(ModelSolverFunctions,MaxIter,errThreshold,Pset,nExpressionMatrix,Tvec,Mvec,Dvec,Cvec,Tref,Ns_Config,Nc_Config,CProbes_Free0,CProbes_Free,ProbeConc,Js_RNA,Js_DNA,Js_Sites,0);
            ProbeSetMetrics.iter(m_unique_loci,t_unique_loci,d_unique_loci)  = iter;
            ProbeSetMetrics.err(m_unique_loci,t_unique_loci,d_unique_loci,1:length(Cvec)) = err;
            ProbeSetMetrics.varSSE(m_unique_loci,t_unique_loci,d_unique_loci,1:length(Cvec))  = varSSE;
            ProbeSetMetrics.eqSSE(m_unique_loci,t_unique_loci,d_unique_loci,1:length(Cvec))  = eqSSE;
            ProbeSetMetrics.CProbes_Free{m_unique_loci,t_unique_loci,d_unique_loci} = CProbes_Free;
            [c_Target_nBound,p_TargetSites_Bound,c_TargetSites_Bound] = A_DetectionSolverWrapper_V4(ModelSolverFunctions,[m_unique_loc(m_unique_loci) t_unique_loc(t_unique_loci) d_unique_loc(d_unique_loci)],Pset,settings,nExpressionMatrix,Tvec,Mvec,Dvec,Cvec,Tref,CProbes_Free,DoesProbeBindSite,Js_RNA,Js_DNA,Js_Sites,Names,ON_IDs_specific,ON_IDs_agnostic,OFF_IDs);
            ProbeSetMetrics.BindingPredictions.p_TargetSites_Bound_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = p_TargetSites_Bound;   %same regardless of t has all t, currently just adds to memory, by duplication
            ProbeSetMetrics.BindingPredictions.c_TargetSites_Bound_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = c_TargetSites_Bound;    %same regardless of t has all t, currently just adds to memory, by duplication
            ProbeSetMetrics.BindingPredictions.c_OnOtherOff_nBound_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = c_Target_nBound;
            %Probe Binding Calculations
            OnSpecificCounts = squeeze(c_Target_nBound(1,:,:));
            OnOtherCounts = squeeze(c_Target_nBound(2,:,:));
            OffCounts = squeeze(c_Target_nBound(3,:,:));
            OnSpecificDistribution = OnSpecificCounts./sum(OnSpecificCounts,1,'omitnan');
            OnOtherDistribution = OnOtherCounts./sum(OnOtherCounts,1,'omitnan');
            OffDistribution = OffCounts./sum(OffCounts,1,'omitnan');
            Non_Counts = sum(OnSpecificCounts,'omitnan')*(V_Cell(Rcell)*6.022*10^23);
            Nother_Counts= sum(OnOtherCounts,'omitnan')*(V_Cell(Rcell)*6.022*10^23);
            Noff_Counts = sum(OffCounts,'omitnan')*(V_Cell(Rcell)*6.022*10^23);
            Non_P = Non_Counts.*OnSpecificDistribution;
            Nother_P = Nother_Counts.*OnOtherDistribution;
            Noff_P = Noff_Counts.*OffDistribution;
            Non_P = Non_P(1:find(squeeze(sum(Non_P,2))>0, 1, 'last'),:);
            Nother_P = Nother_P(1:find(squeeze(sum(Nother_P,2))>0, 1, 'last'),:);
            Noff_P = Noff_P(1:find(squeeze(sum(Noff_P,2))>0, 1, 'last'),:);
            Basic_Noff = sum(Noff_P.*[0:size(Noff_P,1)-1]',1);
            IsoIgnorantConfusion_P = confusionMatrixWrapper_MultiCell(...
                CATnWrapper({permute(Non_P(2:end,:),[2 1]), ndSparse.build([size(Non_P,2) max([size(Non_P,1) size(Noff_P,1) size(Nother_P,1)])-size(Non_P,1)])},2),...
                CATnWrapper({permute(Noff_P(2:end,:),[2 1]), ndSparse.build([size(Noff_P,2) max([size(Non_P,1) size(Noff_P,1) size(Nother_P,1)])-size(Noff_P,1)])},2));
            if (~isempty(Nother_P))
                IsoSpecificConfusion_P = confusionMatrixWrapper_MultiCell(...
                    CATnWrapper({permute(Non_P(2:end,:),[2 1]), ndSparse.build([size(Non_P,2) max([size(Non_P,1) size(Noff_P,1) size(Nother_P,1)])-size(Non_P,1)])},2),...
                    CATnWrapper({permute(Noff_P(2:end,:),[2 1]), ndSparse.build([size(Noff_P,2) max([size(Non_P,1) size(Noff_P,1) size(Nother_P,1)])-size(Noff_P,1)])},2)+...
                    CATnWrapper({permute(Nother_P(2:end,:),[2 1]), ndSparse.build([size(Nother_P,2) max([size(Non_P,1) size(Noff_P,1) size(Nother_P,1)])-size(Nother_P,1)])},2));
                IsoAgnosticConfusion_P = confusionMatrixWrapper_MultiCell(...
                    CATnWrapper({permute(Non_P(2:end,:),[2 1]), ndSparse.build([size(Non_P,2) max([size(Non_P,1) size(Noff_P,1) size(Nother_P,1)])-size(Non_P,1)])},2)+...
                    CATnWrapper({permute(Nother_P(2:end,:),[2 1]), ndSparse.build([size(Nother_P,2) max([size(Non_P,1) size(Noff_P,1) size(Nother_P,1)])-size(Nother_P,1)])},2),...
                    CATnWrapper({permute(Noff_P(2:end,:),[2 1]), ndSparse.build([size(Noff_P,2) max([size(Non_P,1) size(Noff_P,1) size(Nother_P,1)])-size(Noff_P,1)])},2));
            else
                IsoSpecificConfusion_P = [];
                IsoAgnosticConfusion_P = [];
            end
            %store first set of results
            ProbeSetMetrics.BindingPredictions.Con_Distribution_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = OnSpecificCounts;
            ProbeSetMetrics.BindingPredictions.Cother_Distribution_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = OnOtherCounts;
            ProbeSetMetrics.BindingPredictions.Coff_Distribution_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci}= OffCounts;
            ProbeSetMetrics.BindingPredictions.Pon_Distribution_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = OnSpecificDistribution;
            ProbeSetMetrics.BindingPredictions.Pother_Distribution_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = OnOtherDistribution;
            ProbeSetMetrics.BindingPredictions.Poff_Distribution_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = OffDistribution;
            ProbeSetMetrics.BindingPredictions.Non_Count_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Non_Counts;
            ProbeSetMetrics.BindingPredictions.Nother_Count_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Nother_Counts;
            ProbeSetMetrics.BindingPredictions.Noff_Count_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Noff_Counts;
            ProbeSetMetrics.BindingPredictions.Non_history_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Non_P;
            ProbeSetMetrics.BindingPredictions.Nother_history_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci}= Nother_P;
            ProbeSetMetrics.BindingPredictions.Noff_history_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Noff_P;
            ProbeSetMetrics.BindingPredictions.Noff_tot_history_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Basic_Noff;
            ProbeSetMetrics.BindingPredictions.IsoIgnorantConfusion_Probe_P_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = IsoIgnorantConfusion_P;
            ProbeSetMetrics.BindingPredictions.IsoSpecificConfusion_Probe_P_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = IsoSpecificConfusion_P;
            ProbeSetMetrics.BindingPredictions.IsoAgnosticConfusion_Probe_P_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = IsoAgnosticConfusion_P;
            %Intensity Calculations
            Ioff_Pixel = repmat(Basic_Noff,[length(Diameter_vals) 1])*IntensityPerProbe/Nstacks*Rspot^2./(pi/4*Diameter_vals'.^2);%for all diameter entries
            Pn_XtoI_ON =  [CATnWrapper(arrayfun(@(x) nProbeIntensity(SI,SpotIntensity_Mean,SpotIntensity_STD,NumReferenceProbes,x),1:size(Non_P,1)-1,'Un',0),1)];
            Pn_XtoI_OTHER =  [CATnWrapper(arrayfun(@(x) nProbeIntensity(SI,SpotIntensity_Mean,SpotIntensity_STD,NumReferenceProbes,x),1:size(Nother_P,1)-1,'Un',0),1)];
            Pn_XtoI_OFF =  [CATnWrapper(arrayfun(@(x) nProbeIntensity(SI,SpotIntensity_Mean,SpotIntensity_STD,NumReferenceProbes,x),1:size(Noff_P,1)-1,'Un',0),1)];
            Pon_I = permute(squeeze(sum(permute(repmat(Non_P(2:end,:)',[1 1 length(SI)]),[2 3 1]).*Pn_XtoI_ON,1))./sum(Non_P(2:end,:),1),[2 1]);%I over SI
            Non_I = sum(Non_P(2:end,:),1)'.*Pon_I;
            if (~isempty(Nother_P))
                Pother_I = permute(squeeze(sum(permute(repmat(Nother_P(2:end,:)',[1 1 length(SI)]),[2 3 1]).*Pn_XtoI_OTHER,1))./sum(Nother_P(2:end,:),1),[2 1]);%I over SI
                Nother_I = sum(Nother_P(2:end,:),1)'.*Pother_I;
            else
                Pother_I = [];
                Nother_I = [];
            end
            Poff_I = permute(squeeze(sum(permute(repmat(Noff_P(2:end,:)',[1 1 length(SI)]),[2 3 1]).*Pn_XtoI_OFF,1))./sum(Noff_P(2:end,:),1),[2 1]);%I over SI
            Noff_I = sum(Noff_P(2:end,:),1)'.*Poff_I;
            Pon_wAuto_I = CATnWrapper(arrayfun(@(x) plus_Random(full(squeeze(Pon_I(x,:))),Pr_Auto_I),1:size(Pon_I,1),'Un',0),1);
            if (~isempty(Nother_P))
                Pother_wAuto_I = CATnWrapper(arrayfun(@(x) plus_Random(full(squeeze(Pother_I(x,:))),Pr_Auto_I),1:size(Pother_I,1),'Un',0),1);
            else
                Pother_wAuto_I = [];
            end
            Poff_NonAverage_wAuto_I = CATnWrapper(arrayfun(@(x) plus_Random(full(squeeze(Poff_I(x,:))),Pr_Auto_I),1:size(Poff_I,1),'Un',0),1);
            IsoIgnorantConfusion_I = confusionMatrixWrapper_MultiCell(Non_I,Noff_I);
            if (~isempty(Nother_P))
                IsoSpecificConfusion_I = confusionMatrixWrapper_MultiCell(Non_I,Noff_I+Nother_I);
                IsoAgnosticConfusion_I = confusionMatrixWrapper_MultiCell(Non_I+Nother_I,Noff_I);
            else
                IsoSpecificConfusion_I = [];
                IsoAgnosticConfusion_I = [];
            end
            %% Itensity with Cell Size Calculations
            %Signal Prediction Calculation
            x =  SI_Signal_wAuto_I;
            Px = Pon_wAuto_I;
            Py = Pz_domain_function(Pr_Auto_I);
            Q_Signal_Func = @(c) sum(H2_func(SI,z_domain_function(x,squeeze(Px(c,:)))).*Pz_domain_function(squeeze(Px(c,:)))',1);
            Q_Backgd_Func = @(c,ci) sum(H2_func(SI,z_domain_function(SI+squeeze(Ioff_Pixel(ci,c)),Pr_Auto_I)).*Py',1);
            Qsignal = CATnWrapper(arrayfun(@(nth_cell) Q_Signal_Func(nth_cell),1:size(Px,1),'Un',0),1);
            Psignal = CATnWrapper(arrayfun(@(nth_cell) gradient(squeeze(Qsignal(nth_cell,:)),SI)/(max(squeeze(Qsignal(nth_cell,:)))-min(squeeze(Qsignal(nth_cell,:)))),1:length(Cvec),'Un',0),1);
            Qbackgd = CATnWrapper(arrayfun(@(nth_cell) Q_Backgd_Func(nth_cell,1),1:length(Cvec),'Un',0),1);
            Pbackgd = CATnWrapper(arrayfun(@(nth_cell) gradient(squeeze(Qbackgd(nth_cell,:)),SI)/(max(squeeze(Qbackgd(nth_cell,:)))-min(squeeze(Qbackgd(nth_cell,:)))),1:length(Cvec),'Un',0),1);
            Psignal_minus_backgd = CATnWrapper(arrayfun(@(nth_cell) minus_Random(squeeze(Psignal(nth_cell,:)),squeeze(Pbackgd(nth_cell,:))),1:length(Cvec),'Un',0),1);
            Qsignal_minus_backgd = CATnWrapper(arrayfun(@(nth_cell) cumsum(squeeze(Psignal_minus_backgd(nth_cell,:))),1:length(Cvec),'Un',0),1);
            Psignal_minus_backgd = CATnWrapper(arrayfun(@(nth_cell) gradient(squeeze(Qsignal_minus_backgd(nth_cell,:)),SI_SignalMinusBackgd_I)/(max(squeeze(Qsignal_minus_backgd(nth_cell,:)))-min(squeeze(Qsignal_minus_backgd(nth_cell,:)))),1:length(Cvec),'Un',0),1);
            ProbeSetMetrics.IntensityPredictions.Ioff_history_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Ioff_Pixel;
            ProbeSetMetrics.IntensityPredictions.Non_I_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Non_I;%I over SI
            ProbeSetMetrics.IntensityPredictions.Nother_I_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Nother_I;%I over SI
            ProbeSetMetrics.IntensityPredictions.Noff_I_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Noff_I;%I over SI
            ProbeSetMetrics.IntensityPredictions.P_Non_I_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Pon_I;%I over SI
            ProbeSetMetrics.IntensityPredictions.P_Nother_I_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Pother_I;%I over SI
            ProbeSetMetrics.IntensityPredictions.P_Noff_I_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Poff_I;%I over SI
            ProbeSetMetrics.IntensityPredictions.P_Signal_wAuto_I_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Pon_wAuto_I;
            ProbeSetMetrics.IntensityPredictions.P_SignalOther_wAuto_I_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Pother_wAuto_I;
            ProbeSetMetrics.IntensityPredictions.P_SignalOffNonAverage_wAuto_I_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Poff_NonAverage_wAuto_I;
            ProbeSetMetrics.IntensityPredictions.IsoIgnorantConfusion_Intensity_I_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = IsoIgnorantConfusion_I;
            ProbeSetMetrics.IntensityPredictions.IsoSpecificConfusion_Intensity_I_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = IsoSpecificConfusion_I;
            ProbeSetMetrics.IntensityPredictions.IsoAgnosticConfusion_Intensity_I_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = IsoAgnosticConfusion_I;
            ProbeSetMetrics.IntensityPredictions.QzCellBkg_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Qbackgd;
            ProbeSetMetrics.IntensityPredictions.PzCellBkg_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Pbackgd;
            ProbeSetMetrics.IntensityPredictions.QzSignalMinusBackgd_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Qsignal_minus_backgd;
            ProbeSetMetrics.IntensityPredictions.PzSignalMinusBackgd_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Psignal_minus_backgd;
            ProbeSetMetrics.IntensityPredictions.QzSignal_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Qsignal;
            ProbeSetMetrics.IntensityPredictions.PzSignal_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Psignal;
            %Computational Threshold Calculations
            vars = {'AP','PP','TP','FP','FN','FPR','FNR','TPR','TNR','FDR','FOR','PPV','NPV','Accuracy','BalancedAccuracy','F1','P4','MCC','CKC','FowlkesMallowsIndex','Accuracy'};
            Positive_SpotCounts_Matrix = IsoIgnorantConfusion_I.PP;
            F1_ScoreCurve_Matrix = IsoIgnorantConfusion_I.F1;
            F2_ScoreCurve_Matrix =  IsoIgnorantConfusion_I.Fbeta(2);
            Fhalf_ScoreCurve_Matrix =  IsoIgnorantConfusion_I.Fbeta(0.5);
            Filtered_SpotCountCurves = arrayfun(@(nth_cell)  unique(round(Positive_SpotCounts_Matrix(nth_cell,:)),'stable'),1:length(Cvec),'Un',0);
            Filtered_SpotCountIndexes = arrayfun(@(nth_cell) arrayfun(@(z) find(round(Positive_SpotCounts_Matrix(nth_cell,:))==z,1),Filtered_SpotCountCurves{nth_cell}),1:length(Cvec),'Un',0);
            param_struct_vector = arrayfun(@(x) RNAThreshold.genEmptyThresholdParamStruct(),1:length(Cvec),'Un',0);
            for nn = 1:length(Cvec)
                param_struct_vector{nn} = TrueSpotDefaultThParameters;
                param_struct_vector{nn}.sample_spot_table = [[1:length(Filtered_SpotCountCurves{nn})]' Filtered_SpotCountCurves{nn}'];
            end
            scThresholdSuggestions = arrayfun(@(nn) RNAThreshold.scoreThresholdSuggestions(RNAThreshold.estimateThreshold(param_struct_vector{nn})),1:length(Cvec),'Un',0);
            scThresholdSuggestions =  [scThresholdSuggestions{:}];
            subfield_groups = {'pool','thstats'};
            for v = 1:length(subfield_groups)
                subfields = fieldnames(scThresholdSuggestions(1).(subfield_groups{v}));
                for subf = 1:length(subfields)
                    subf_vals = arrayfun(@(nn) scThresholdSuggestions(nn).(subfield_groups{v}).(subfields{subf}),1:size(scThresholdSuggestions,2),'UniformOutput',0);
                    [scThresholdSuggestions.(subfields{subf})] = subf_vals{:};
                end
            end
            scThresholdSuggestions = rmfield(scThresholdSuggestions,subfield_groups);
            TrueSpot_FilteredThreshold = [scThresholdSuggestions.threshold];
            TrueSpot_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) Filtered_SpotCountIndexes{nth_cell}(TrueSpot_FilteredThreshold(nth_cell)),1:length(Cvec),'Un',0));
            F1_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(F1_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(F1_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
            F1_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) Filtered_SpotCountIndexes{nth_cell}(find(F1_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(F1_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
            F2_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(F2_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(F2_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
            F2_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) Filtered_SpotCountIndexes{nth_cell}(find(F2_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(F2_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
            Fhalf_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(Fhalf_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(Fhalf_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
            Fhalf_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) Filtered_SpotCountIndexes{nth_cell}(find(Fhalf_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(Fhalf_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
            TS_CountMetricFunction = @(k) IsoIgnorantConfusion_I.(vars{k})(sub2ind(size(IsoIgnorantConfusion_I.(vars{k})),1:length(Cvec),TrueSpot_ThresholdLocations));
            F1_CountMetricFunction = @(k) IsoIgnorantConfusion_I.(vars{k})(sub2ind(size(IsoIgnorantConfusion_I.(vars{k})),1:length(Cvec),F1_ThresholdLocations));
            F2_CountMetricFunction = @(k) IsoIgnorantConfusion_I.(vars{k})(sub2ind(size(IsoIgnorantConfusion_I.(vars{k})),1:length(Cvec),F2_ThresholdLocations));
            Fhalf_CountMetricFunction = @(k) IsoIgnorantConfusion_I.(vars{k})(sub2ind(size(IsoIgnorantConfusion_I.(vars{k})),1:length(Cvec),Fhalf_ThresholdLocations));
            SpotCountMetrics = [];
            SpotCountMetrics.TrueSpot.Thresholds = TrueSpot_FilteredThreshold;
            SpotCountMetrics.F1.Thresholds = F1_FilteredThreshold;
            SpotCountMetrics.F2.Thresholds = F2_FilteredThreshold;
            SpotCountMetrics.Fhalf.Thresholds = Fhalf_FilteredThreshold;
            for k = 1:length(vars)
                SpotCountMetrics.TrueSpot.(vars{k}) = TS_CountMetricFunction(k);
                SpotCountMetrics.F1.(vars{k}) = F1_CountMetricFunction(k);
                SpotCountMetrics.F2.(vars{k}) = F2_CountMetricFunction(k);
                SpotCountMetrics.Fhalf.(vars{k}) = Fhalf_CountMetricFunction(k);
            end
            ProbeSetMetrics.CountPredictions.IsoIgnorant_SpotCountMetrics_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = SpotCountMetrics;
            if (~isempty(Nother_P))
                IsoSpecificConfusion_I = confusionMatrixWrapper_MultiCell(Non_I,Noff_I+Nother_I);
                IsoAgnosticConfusion_I = confusionMatrixWrapper_MultiCell(Non_I+Nother_I,Noff_I);
                Positive_SpotCounts_Matrix = IsoSpecificConfusion_I.PP;
                F1_ScoreCurve_Matrix = IsoSpecificConfusion_I.F1;
                F2_ScoreCurve_Matrix =  IsoSpecificConfusion_I.Fbeta(2);
                Fhalf_ScoreCurve_Matrix =  IsoSpecificConfusion_I.Fbeta(0.5);
                Filtered_SpotCountCurves = arrayfun(@(nth_cell)  unique(round(Positive_SpotCounts_Matrix(nth_cell,:)),'stable'),1:length(Cvec),'Un',0);
                Filtered_SpotCountIndexes = arrayfun(@(nth_cell) arrayfun(@(z) find(round(Positive_SpotCounts_Matrix(nth_cell,:))==z,1),Filtered_SpotCountCurves{nth_cell}),1:length(Cvec),'Un',0);
                param_struct_vector = arrayfun(@(x) RNAThreshold.genEmptyThresholdParamStruct(),1:length(Cvec),'Un',0);
                for nn = 1:length(Cvec)
                    param_struct_vector{nn} = TrueSpotDefaultThParameters;
                    param_struct_vector{nn}.sample_spot_table = [[1:length(Filtered_SpotCountCurves{nn})]' Filtered_SpotCountCurves{nn}'];
                end
                scThresholdSuggestions = arrayfun(@(nn) RNAThreshold.scoreThresholdSuggestions(RNAThreshold.estimateThreshold(param_struct_vector{nn})),1:length(Cvec),'Un',0);
                scThresholdSuggestions =  [scThresholdSuggestions{:}];
                subfield_groups = {'pool','thstats'};
                for v = 1:length(subfield_groups)
                    subfields = fieldnames(scThresholdSuggestions(1).(subfield_groups{v}));
                    for subf = 1:length(subfields)
                        subf_vals = arrayfun(@(nn) scThresholdSuggestions(nn).(subfield_groups{v}).(subfields{subf}),1:size(scThresholdSuggestions,2),'UniformOutput',0);
                        [scThresholdSuggestions.(subfields{subf})] = subf_vals{:};
                    end
                end
                scThresholdSuggestions = rmfield(scThresholdSuggestions,subfield_groups);
                TrueSpot_FilteredThreshold = [scThresholdSuggestions.threshold];
                TrueSpot_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) Filtered_SpotCountIndexes{nth_cell}(TrueSpot_FilteredThreshold(nth_cell)),1:length(Cvec),'Un',0));
                F1_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(F1_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(F1_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                F1_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) Filtered_SpotCountIndexes{nth_cell}(find(F1_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(F1_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                F2_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(F2_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(F2_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                F2_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) Filtered_SpotCountIndexes{nth_cell}(find(F2_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(F2_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                Fhalf_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(Fhalf_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(Fhalf_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                Fhalf_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) Filtered_SpotCountIndexes{nth_cell}(find(Fhalf_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(Fhalf_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                TS_CountMetricFunction = @(k) IsoSpecificConfusion_I.(vars{k})(sub2ind(size(IsoSpecificConfusion_I.(vars{k})),1:length(Cvec),TrueSpot_ThresholdLocations));
                F1_CountMetricFunction = @(k) IsoSpecificConfusion_I.(vars{k})(sub2ind(size(IsoSpecificConfusion_I.(vars{k})),1:length(Cvec),F1_ThresholdLocations));
                F2_CountMetricFunction = @(k) IsoSpecificConfusion_I.(vars{k})(sub2ind(size(IsoSpecificConfusion_I.(vars{k})),1:length(Cvec),F2_ThresholdLocations));
                Fhalf_CountMetricFunction = @(k) IsoSpecificConfusion_I.(vars{k})(sub2ind(size(IsoSpecificConfusion_I.(vars{k})),1:length(Cvec),Fhalf_ThresholdLocations));
                IsoSpecific_SpotCountMetrics = [];
                IsoSpecific_SpotCountMetrics.TrueSpot.Thresholds = TrueSpot_FilteredThreshold;
                IsoSpecific_SpotCountMetrics.F1.Thresholds = F1_FilteredThreshold;
                IsoSpecific_SpotCountMetrics.F2.Thresholds = F2_FilteredThreshold;
                IsoSpecific_SpotCountMetrics.Fhalf.Thresholds = Fhalf_FilteredThreshold;
                for k = 1:length(vars)
                    IsoSpecific_SpotCountMetrics.TrueSpot.(vars{k}) = TS_CountMetricFunction(k);
                    IsoSpecific_SpotCountMetrics.F1.(vars{k}) = F1_CountMetricFunction(k);
                    IsoSpecific_SpotCountMetrics.F2.(vars{k}) = F2_CountMetricFunction(k);
                    IsoSpecific_SpotCountMetrics.Fhalf.(vars{k}) = Fhalf_CountMetricFunction(k);
                end
                ProbeSetMetrics.CountPredictions.IsoSpecific_SpotCountMetrics{m_unique_loci,t_unique_loci,d_unique_loci} = IsoSpecific_SpotCountMetrics;
                Positive_SpotCounts_Matrix = IsoAgnosticConfusion_I.PP;
                F1_ScoreCurve_Matrix = IsoAgnosticConfusion_I.F1;
                F2_ScoreCurve_Matrix =  IsoAgnosticConfusion_I.Fbeta(2);
                Fhalf_ScoreCurve_Matrix =  IsoAgnosticConfusion_I.Fbeta(0.5);
                Filtered_SpotCountCurves = arrayfun(@(nth_cell)  unique(round(Positive_SpotCounts_Matrix(nth_cell,:)),'stable'),1:length(Cvec),'Un',0);
                Filtered_SpotCountIndexes = arrayfun(@(nth_cell) arrayfun(@(z) find(round(Positive_SpotCounts_Matrix(nth_cell,:))==z,1),Filtered_SpotCountCurves{nth_cell}),1:length(Cvec),'Un',0);
                param_struct_vector = arrayfun(@(x) RNAThreshold.genEmptyThresholdParamStruct(),1:length(Cvec),'Un',0);
                for nn = 1:length(Cvec)
                    param_struct_vector{nn} = TrueSpotDefaultThParameters;
                    param_struct_vector{nn}.sample_spot_table = [[1:length(Filtered_SpotCountCurves{nn})]' Filtered_SpotCountCurves{nn}'];
                end
                scThresholdSuggestions = arrayfun(@(nn) RNAThreshold.scoreThresholdSuggestions(RNAThreshold.estimateThreshold(param_struct_vector{nn})),1:length(Cvec),'Un',0);
                scThresholdSuggestions =  [scThresholdSuggestions{:}];
                subfield_groups = {'pool','thstats'};
                for v = 1:length(subfield_groups)
                    subfields = fieldnames(scThresholdSuggestions(1).(subfield_groups{v}));
                    for subf = 1:length(subfields)
                        subf_vals = arrayfun(@(nn) scThresholdSuggestions(nn).(subfield_groups{v}).(subfields{subf}),1:size(scThresholdSuggestions,2),'UniformOutput',0);
                        [scThresholdSuggestions.(subfields{subf})] = subf_vals{:};
                    end
                end
                scThresholdSuggestions = rmfield(scThresholdSuggestions,subfield_groups);
                TrueSpot_FilteredThreshold = [scThresholdSuggestions.threshold];
                TrueSpot_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) Filtered_SpotCountIndexes{nth_cell}(TrueSpot_FilteredThreshold(nth_cell)),1:length(Cvec),'Un',0));
                F1_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(F1_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(F1_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                F1_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) Filtered_SpotCountIndexes{nth_cell}(find(F1_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(F1_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                F2_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(F2_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(F2_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                F2_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) Filtered_SpotCountIndexes{nth_cell}(find(F2_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(F2_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                Fhalf_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(Fhalf_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(Fhalf_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                Fhalf_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) Filtered_SpotCountIndexes{nth_cell}(find(Fhalf_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})==max(Fhalf_ScoreCurve_Matrix(nth_cell,Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                TS_CountMetricFunction = @(k) IsoAgnosticConfusion_I.(vars{k})(sub2ind(size(IsoAgnosticConfusion_I.(vars{k})),1:length(Cvec),TrueSpot_ThresholdLocations));
                F1_CountMetricFunction = @(k) IsoAgnosticConfusion_I.(vars{k})(sub2ind(size(IsoAgnosticConfusion_I.(vars{k})),1:length(Cvec),F1_ThresholdLocations));
                F2_CountMetricFunction = @(k) IsoAgnosticConfusion_I.(vars{k})(sub2ind(size(IsoAgnosticConfusion_I.(vars{k})),1:length(Cvec),F2_ThresholdLocations));
                Fhalf_CountMetricFunction = @(k) IsoAgnosticConfusion_I.(vars{k})(sub2ind(size(IsoAgnosticConfusion_I.(vars{k})),1:length(Cvec),Fhalf_ThresholdLocations));
                IsoAgnostic_SpotCountMetrics = [];
                IsoAgnostic_SpotCountMetrics.TrueSpot.Thresholds = TrueSpot_FilteredThreshold;
                IsoAgnostic_SpotCountMetrics.F1.Thresholds = F1_FilteredThreshold;
                IsoAgnostic_SpotCountMetrics.F2.Thresholds = F2_FilteredThreshold;
                IsoAgnostic_SpotCountMetrics.Fhalf.Thresholds = Fhalf_FilteredThreshold;
                for k = 1:length(vars)
                    IsoAgnostic_SpotCountMetrics.TrueSpot.(vars{k}) = TS_CountMetricFunction(k);
                    IsoAgnostic_SpotCountMetrics.F1.(vars{k}) = F1_CountMetricFunction(k);
                    IsoAgnostic_SpotCountMetrics.F2.(vars{k}) = F2_CountMetricFunction(k);
                    IsoAgnostic_SpotCountMetrics.Fhalf.(vars{k}) = Fhalf_CountMetricFunction(k);
                end
                ProbeSetMetrics.CountPredictions.IsoAgnostic_SpotCountMetrics_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = IsoAgnostic_SpotCountMetrics;
            else
                ProbeSetMetrics.CountPredictions.IsoSpecific_SpotCountMetrics_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                ProbeSetMetrics.CountPredictions.IsoAgnostic_SpotCountMetrics_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
            end
        end
    end
end
ModelMetrics.ProbeSetMetrics = ProbeSetMetrics;
end