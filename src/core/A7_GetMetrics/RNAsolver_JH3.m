function ModelMetrics = RNAsolver_JH3(Pset,settings,probes,gene_table,ExpressionMatrix,DoesProbeBindSite,dHeq_mod,dSeq_mod,dCp_mod,dHeq_Complement,dSeq_Complement,dCp_Complement)
%% This function computes the metrics for RNA-FISH performance.
% The function simulates probe binding dynamics to model probe
% equilibrium binding under a variety of conditions. Varying by
% temperature, expression, probe concentration, etc.
% The distribution of probes bound targets is computed and is
% used to generate confusion matrix metrics, of probe performance.
N_Fluorophores = 2;
loaded_TrueSpot_default_params = 0;
if (isfile('data/TS_DefaultParams.mat'))
    load('data/TS_DefaultParams.mat','TrueSpotDefaultThParameters');
    loaded_TrueSpot_default_params = 1;
end
%% Handle Functions
if ~(ismcc || isdeployed)
    %#exclude heaviside
    H2_func = @(z,x) heaviside(z-x');
else
    H2_func = @(z,x) deployable_heaviside(z-x');
end
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
GuessConc = settings.SimulationConfiguration.InitialGuessConc;
PC0 = settings.SimulationConfiguration.ProbeConcentration;
errThreshold = settings.SimulationConfiguration.errThreshold;
MaxIter = settings.SimulationConfiguration.MaxIter;
Diameter_vals = Mean_Diameter;
SI_StepSize = settings.SimulationConfiguration.Signal_StepSize;
SI_MaxSignal = settings.SimulationConfiguration.Signal_MaxValue;
logRatioTrim = settings.TMM.logRatioTrim;
sumTrim = settings.TMM.sumTrim;
Acutoff = settings.TMM.Acutoff;
doWeighting = settings.TMM.doWeighting;
SI = 0:SI_StepSize:SI_MaxSignal;%only works if SI
SI_SignalMinusBackgd_I = -max(SI):SI_StepSize:1*max(SI);
SI_Signal_wAuto_I = 0:SI_StepSize:2*max(SI);
IntensityPerProbe_Mean = SpotIntensity_Mean/NumReferenceProbes;
IntensityPerProbe_STD = SpotIntensity_STD/sqrt(NumReferenceProbes);

Pr_Auto_I = pIntensity(SI,AutoBackground_Mean,AutoBackground_STD);%I over SI
ExpressionMatrix(isnan(ExpressionMatrix)) = 0;
NumNonUniformConditions = size(ExpressionMatrix,2);


[~, ExpressionMatrix_nTPM] = tmm(double(ExpressionMatrix),logRatioTrim,sumTrim,Acutoff,doWeighting);
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


gene_table = sortrows(gene_table,[7 6],'ascend');
gene_table = gene_table(gene_table.Match>=settings.MinHomologySearchTargetSize,:);
MinusStrandedHits = find(contains(gene_table.Strand,'Minus'));
gene_table_NamesZ = convertCharsToStrings(gene_table.Name);
contains_RNA = find(ismember(gene_table_NamesZ,settings.RNAdbParser));clear gene_table_NamesZ
RNA_MissedFilteredHits = intersect(MinusStrandedHits,contains_RNA);clear contains_RNA
gene_table = gene_table(setdiff(1:size(gene_table,1),RNA_MissedFilteredHits),:);clear RNA_MissedFilteredHits
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
    DNA_IDs = find(ismember(Names,settings.DNAdbParser));%IDs
else
    DNA_IDs = [];
end
if (settings.BLASTrna)
    NonDNA_IDs = find(ismember(Names,settings.RNAdbParser));%IDs
else
    NonDNA_IDs =[];
end
ON_IDs_specific = find(ismember(uniNames,extractBefore(settings.transcript_IDs{:},'.')));
ON_IDs_agnostic = find(ismember(uniNames,extractBefore(settings.transcript_IDs_desired{:},'.')));
OFF_IDs = find(~ismember(uniNames,extractBefore(settings.transcript_IDs_desired{:},'.')));
if (ndims(dHeq_Complement)~=3)%error quick fix
    dHeq_Complement =  permute(repmat(dHeq_Complement, [1 1 size(dHeq_mod,3)]),[1 3 2]);
    dSeq_Complement = permute(repmat(dSeq_Complement, [1 1 size(dHeq_mod,3)]),[1 3 2]);
    dCp_Complement = permute(repmat(dCp_Complement, [1 1 size(dHeq_mod,3)]),[1 3 2]);
end
fprintf("Computing Probe Design Target Packing Efficiency")
fprintf('\n')
fprintf('\n')
LocMax = max(cell2mat(cellfun(@(x) x,probes(:,3),'UniformOutput',false)));
Lpmin = min(cell2mat(cellfun(@length,probes(:,2),'UniformOutput',false)));
TargetLength = LocMax + Lpmin - 1;
theoryMaxProbes = floor(TargetLength/(Lpmin+settings.ProbeSpacing));
if (theoryMaxProbes>settings.maxProbes)
    theoryMaxProbes = settings.maxProbes;
end
PackEf = length(Pset)/theoryMaxProbes;
ModelMetrics.PackingEfficiency = PackEf;
fprintf("Generating thermodynamic-kinetic model structure")
fprintf('\n')
fprintf('\n')
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
ProbeSetMetrics.BindingPredictions.p_IndividualTargets_nBound_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
ProbeSetMetrics.BindingPredictions.c_IndividualTargets_nBound_ModelTemperatureDilutionVector = cell(length(m_unique_loc),length(t_unique_loc),length(d_unique_loc));
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
fprintf("Solving model for all specified solution configurations")
fprintf('\n')
fprintf('\n')
for m_unique_loci = 1:length(m_unique_loc)
    for t_unique_loci = 1:length(t_unique_loc)
        for d_unique_loci = 1:length(d_unique_loc)
            CProbes_Free = GuessConc*ones(length(Pset),length(Mvec(m_unique_loc(m_unique_loci))),length(Cvec));
            CProbes_Free0 = GuessConc*ones(length(Pset),length(Mvec(m_unique_loc(m_unique_loci))),length(Cvec));
            ProbeConc = PC0*squeeze(permute(repmat(Dvec(d_unique_loc(d_unique_loci)),[1 1 length(Pset) 1 1 length(Cvec)]),[1 3 4 5 2 6]));
            fprintf("Solving model for free-steady state probe concentrations")
            fprintf('\n')
            fprintf(strcat("Initial Probe Concentration = ",string(PC0*Dvec(d_unique_loc(d_unique_loci))),"μM, Gibbs Model = ",string(Mvec(m_unique_loc(m_unique_loci))),",Temperature = ",string(Tvec(t_unique_loc(t_unique_loci))-273.15),"°C"))
            fprintf('\n')
            fprintf('\n')
            [CProbes_Free,varSSE,err,iter,eqSSE] = A_ModelEquilibriumSolverWrapper_V4(ModelSolverFunctions,MaxIter,errThreshold,Pset,nExpressionMatrix,Tvec(t_unique_loc(t_unique_loci)),Mvec(m_unique_loc(m_unique_loci)),Dvec(d_unique_loc(d_unique_loci)),Cvec,Tref,Ns_Config,Nc_Config,CProbes_Free0,CProbes_Free,ProbeConc,Js_RNA,Js_DNA,Js_Sites,0);
            ProbeSetMetrics.iter(m_unique_loci,t_unique_loci,d_unique_loci)  = iter;
            ProbeSetMetrics.err(m_unique_loci,t_unique_loci,d_unique_loci,1:length(Cvec)) = err;
            ProbeSetMetrics.varSSE(m_unique_loci,t_unique_loci,d_unique_loci,1:length(Cvec))  = varSSE;
            ProbeSetMetrics.eqSSE(m_unique_loci,t_unique_loci,d_unique_loci,1:length(Cvec))  = eqSSE;
            ProbeSetMetrics.CProbes_Free{m_unique_loci,t_unique_loci,d_unique_loci} = CProbes_Free;
            fprintf("Computing model steady-state equilibrium probe-target duplex concentrations")
            fprintf('\n')
            fprintf(strcat("Initial Probe Concentration = ",string(PC0*Dvec(d_unique_loc(d_unique_loci))),"μM, Gibbs Model = ",string(Mvec(m_unique_loc(m_unique_loci))),",Temperature = ",string(Tvec(t_unique_loc(t_unique_loci))-273.15),"°C"))
            fprintf('\n')
            fprintf('\n')
            if (nnz(dHeq_mod(:,:,:,m_unique_loc(m_unique_loci)))+nnz(dSeq_mod(:,:,:,m_unique_loc(m_unique_loci)))>0)
                [c_Target_nBound,p_TargetSites_Bound,c_TargetSites_Bound,c_IndividualTargets_nBound,p_IndividualTargets_nBound,tHit] = ...
                    A_DetectionSolverWrapper_V4(ModelSolverFunctions,...
                    [1 1 1],Pset,settings,nExpressionMatrix,Tvec(t_unique_loc(t_unique_loci)),Mvec(m_unique_loc(m_unique_loci)),Dvec(d_unique_loc(d_unique_loci)),Cvec,Tref,CProbes_Free,DoesProbeBindSite,Js_RNA,Js_DNA,Js_Sites,Names,ON_IDs_specific,ON_IDs_agnostic,OFF_IDs);
                ProbeSetMetrics.BindingPredictions.p_IndividualTargets_nBound_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = p_IndividualTargets_nBound;   %same regardless of t has all t, currently just adds to memory, by duplication
                ProbeSetMetrics.BindingPredictions.c_IndividualTargets_nBound_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = c_IndividualTargets_nBound;    %same regardless of t has all t, currently just adds to memory, by duplication
                ProbeSetMetrics.BindingPredictions.p_TargetSites_Bound_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = p_TargetSites_Bound;   %same regardless of t has all t, currently just adds to memory, by duplication
                ProbeSetMetrics.BindingPredictions.c_TargetSites_Bound_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = c_TargetSites_Bound;    %same regardless of t has all t, currently just adds to memory, by duplication
                ProbeSetMetrics.BindingPredictions.c_OnOtherOff_nBound_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = c_Target_nBound;
                pON_IDs = ismember(tHit,ON_IDs_specific);
                pOFF_IDs = ismember(tHit,OFF_IDs);
                pON_IDs_other = ismember(tHit,setdiff(ON_IDs_agnostic,ON_IDs_specific));
                if (~isempty(pON_IDs_other))
                    p_Isoforms_1D = p_IndividualTargets_nBound([find(pON_IDs); find(pON_IDs_other)],:,:);
                    p_Isoform_Identifiability_P = CATnWrapper(arrayfun(@(c) arrayfun(@(x,y) 1/2*sum(abs(p_Isoforms_1D(x,:,c)-p_Isoforms_1D(y,:,c)),'all'), ...
                        meshgrid(1: length(ON_IDs_agnostic),1:length(ON_IDs_agnostic)),...
                        meshgrid(1:length(ON_IDs_agnostic),1: length(ON_IDs_agnostic))'),1:length(Cvec),'UniformOutput',0),3);
                    ProbeSetMetrics.BindingPredictions.p_Isoforms_P_1D_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = p_Isoforms_1D;
                    ProbeSetMetrics.BindingPredictions.p_Isoform_Identifiability_P_1D_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = p_Isoform_Identifiability_P;
                else
                    ProbeSetMetrics.BindingPredictions.p_Isoforms_P_1D_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                    ProbeSetMetrics.BindingPredictions.p_Isoform_Identifiability_P_1D_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                end
                if (N_Fluorophores>1)
                    tilePosition = cell2mat(probes(Pset,3));
                    [~,probe_individual_binding_tile_order] = sort(tilePosition,'ascend');
                    PsetTilePositionOrder = 1:length(Pset);
                    % Get even-valued entries or odd
                    Fluorophore_labeled_probes = cell(1,N_Fluorophores);
                    pAnyBindSiteCM0_Fluorophores = cell(1,N_Fluorophores);
                    for nflo = 1:N_Fluorophores
                        Fluorophore_labeled_probes{nflo} = probe_individual_binding_tile_order(PsetTilePositionOrder(mod(PsetTilePositionOrder, N_Fluorophores) == nflo-1));
                    end
                    %subset of targets to actually binded
                    %P,T,S,M, .....
                    Js = @(x) find(sum(squeeze(sum(DoesProbeBindSite(x,:,:),1)),2)>0);
                    Sx = unique(cell2mat(arrayfun(@(x) find(sum(DoesProbeBindSite(Pset,x,:),1)>0)',Js(Pset),'Un',0)));
                    %different for multi-weave and split or probes having any number of labels
                    switch ModelSolverFunctions.solverType
                        case 0
                            if (length(Pset)>1)
                                for nflo = 1:N_Fluorophores
                                    pAnyBindSiteCM0_Fluorophores{nflo} = squeeze(sum(squeeze(p_TargetSites_Bound(Fluorophore_labeled_probes{nflo},:,Sx,1,:))));
                                end
                            else
                                for nflo = 1:N_Fluorophores
                                    pAnyBindSiteCM0_Fluorophores{nflo} = squeeze(p_TargetSites_Bound(Fluorophore_labeled_probes{nflo},:,Sx,1,:));
                                end
                            end
                        case 1
                            if (length(Pset)>1)
                                for nflo = 1:N_Fluorophores
                                    pAnyBindSiteCM0_Fluorophores{nflo} = squeeze(sum(squeeze(p_TargetSites_Bound(Fluorophore_labeled_probes{nflo},:,Sx,1,:)),1));
                                end
                            else
                                for nflo = 1:N_Fluorophores
                                    pAnyBindSiteCM0_Fluorophores{nflo} = squeeze(p_TargetSites_Bound(Fluorophore_labeled_probes{nflo},:,Sx,1,:));
                                end
                            end
                        case 2
                            if (length(Pset)>1)
                                for nflo = 1:N_Fluorophores
                                    pAnyBindSiteCM0_Fluorophores{nflo} = squeeze(sum(squeeze(p_TargetSites_Bound(Fluorophore_labeled_probes{nflo},:,Sx,1,1,1,:)),1));
                                end
                            else
                                for nflo = 1:N_Fluorophores
                                    pAnyBindSiteCM0_Fluorophores{nflo} = squeeze(p_TargetSites_Bound(Fluorophore_labeled_probes{nflo},:,Sx,1,1,1,:));
                                end
                            end
                    end
                    tHit = unique(cell2mat(arrayfun(@(nflo) reshape(find(squeeze(sum(sum(pAnyBindSiteCM0_Fluorophores{nflo},2),3))>0),1,[]),1:N_Fluorophores,'Un',0)));
                    pON_IDs = ismember(tHit,ON_IDs_specific);
                    pOFF_IDs = ismember(tHit,OFF_IDs);
                    pON_IDs_other = ismember(tHit,setdiff(ON_IDs_agnostic,ON_IDs_specific));
                    p_multi_single_labeled_Fluorophores = permute(CATnWrapper(pAnyBindSiteCM0_Fluorophores,4),[2 4 1 3]);
                    p_multi_nBound = F_DiscretePoissonMultinomialMultiAny_V2(p_multi_single_labeled_Fluorophores(:,:,tHit,:),'LS');
                    if (~isempty(pON_IDs_other))
                        front_indices = repmat({':'}, 1, N_Fluorophores);
                        handle_ab = @(a,b) [{a},front_indices, {b}];
                        handle_slice_indexing = @(A) A{:};
                        dimAgnosticON_IDs_indices = [find(pON_IDs) find(pON_IDs_other)];
                        end_indices = repmat({':'}, 1, 1);
                        all_indices_AgnosticON_IDs = [front_indices, {dimAgnosticON_IDs_indices},end_indices];
                        p_Isoforms_ND = permute(p_multi_nBound(all_indices_AgnosticON_IDs{:}),[N_Fluorophores+1 1:N_Fluorophores N_Fluorophores+2]);
                        p_Isoform_Identifiability_ND_P = CATnWrapper(arrayfun(@(c) arrayfun(@(x,y) 1/2*sum(squeeze(abs(p_Isoforms_ND(handle_slice_indexing(handle_ab(x,c)))-p_Isoforms_ND(handle_slice_indexing(handle_ab(y,c))))),'all'), ...
                            meshgrid(1: length(ON_IDs_agnostic),1:length(ON_IDs_agnostic)),...
                            meshgrid(1:length(ON_IDs_agnostic),1: length(ON_IDs_agnostic))'),1:length(Cvec),'UniformOutput',0),N_Fluorophores+1);
                        ProbeSetMetrics.BindingPredictions.p_Isoforms_P_ND_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = p_Isoforms_ND;
                        ProbeSetMetrics.BindingPredictions.p_Isoform_Identifiability_P_ND_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = p_Isoform_Identifiability_ND_P;
                    else
                        ProbeSetMetrics.BindingPredictions.p_Isoforms_P_ND_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                        ProbeSetMetrics.BindingPredictions.p_Isoform_Identifiability_P_ND_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                    end
                    multi_dims = size(p_multi_nBound);
                    Nt_multi_nBound = repmat(reshape(nExpressionMatrix(tHit,Cvec),[1 1 length(tHit) length(Cvec)]),[multi_dims(1:N_Fluorophores) 1 1]).*p_multi_nBound;
                    front_indices = repmat({':'}, 1, N_Fluorophores);
                    dimON_IDs_indices = find(pON_IDs);
                    dimON_IDs_other_indices = find(pON_IDs_other);
                    dimOFF_IDs_indices = find(pOFF_IDs);
                    end_indices = repmat({':'}, 1, 1);
                    all_indices_pON_IDs = [front_indices, {dimON_IDs_indices},end_indices];
                    all_indices_pON_IDs_other = [front_indices, {dimON_IDs_other_indices},end_indices];
                    all_indices_pOFF_IDs = [front_indices, {dimOFF_IDs_indices},end_indices];
                    NtON_multi_specific = squeeze(sum(Nt_multi_nBound(all_indices_pON_IDs{:}),N_Fluorophores+1,'omitnan'));
                    NtON_multi_other = squeeze(sum(Nt_multi_nBound(all_indices_pON_IDs_other{:}),N_Fluorophores+1,'omitnan'));
                    NtOFF_multi = squeeze(sum(Nt_multi_nBound(all_indices_pOFF_IDs{:}),N_Fluorophores+1,'omitnan'));
                    NtN_multi = permute(CATnWrapper({NtON_multi_specific,NtON_multi_other,NtOFF_multi},N_Fluorophores+2),[N_Fluorophores+2 1:N_Fluorophores+1]);
                    NtTotal_multi = squeeze(sum(NtN_multi,1));
                    ProbeSetMetrics.BindingPredictions.p_IndividualTargets_nBoundMulti_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = p_multi_nBound;
                    ProbeSetMetrics.BindingPredictions.c_IndividualTargets_nBoundMulti_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = Nt_multi_nBound;    %same regardless of t has all t, currently just adds to memory, by duplication
                    ProbeSetMetrics.BindingPredictions.c_OnOtherOff_nBoundMulti_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = NtN_multi;
                    ProbeSetMetrics.BindingPredictions.c_Total_nBoundMulti_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = NtTotal_multi;
                    ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_Ind_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = F_PMF_generalMetrics(p_multi_nBound, [1:N_Fluorophores], [], [], [], N_Fluorophores + [1 2], [], [],[]);
                    ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_Total_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = F_PMF_generalMetrics(NtTotal_multi , [1:N_Fluorophores], N_Fluorophores + [1], [], [], [], [], [],[]);
                    ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_ON_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = F_PMF_generalMetrics(NtON_multi_specific, [1:N_Fluorophores], N_Fluorophores + [1], [], [], [], [], [],[]);
                    ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_OFF_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} =  F_PMF_generalMetrics(NtOFF_multi, [1:N_Fluorophores], N_Fluorophores + [1], [], [], [], [], [],[]);
                    if ~isempty(find(pON_IDs_other, 1))
                        ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_OTHER_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = F_PMF_generalMetrics(NtON_multi_other, [1:N_Fluorophores], N_Fluorophores + [1], [], [], [], [], [],[]);
                        ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_AgnosticON_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = F_PMF_generalMetrics(NtON_multi_other+NtOFF_multi, [1:N_Fluorophores], N_Fluorophores + [1], [], [], [], [], [],[]);
                        ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_SpecificOFF_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = F_PMF_generalMetrics(NtON_multi_other+NtON_multi_specific, [1:N_Fluorophores], N_Fluorophores + [1], [], [], [], [], [],[]);
                    else
                        ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_OTHER_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                        ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_AgnosticON_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                        ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_SpecificOFF_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                    end
                else
                    ProbeSetMetrics.BindingPredictions.p_Isoforms_P_ND_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                    ProbeSetMetrics.BindingPredictions.p_Isoform_Identifiability_P_ND_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                    ProbeSetMetrics.BindingPredictions.p_IndividualTargets_nBoundMulti_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                    ProbeSetMetrics.BindingPredictions.c_IndividualTargets_nBoundMulti_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                    ProbeSetMetrics.BindingPredictions.c_OnOtherOff_nBoundMulti_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                    ProbeSetMetrics.BindingPredictions.c_Total_nBoundMulti_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                    ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_Ind_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                    ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_Total_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci}  = [];
                    ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_ON_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci}  = [];
                    ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_OFF_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                    ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_OTHER_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                    ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_AgnosticON_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                    ProbeSetMetrics.BindingPredictions.MultiFluoroMetrics_P_SpecificOFF_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                end
                %Probe Binding Calculations
                fprintf("Computing probe set on/off-target binding model metrics")
                fprintf('\n')
                fprintf(strcat("Initial Probe Concentration = ",string(PC0*Dvec(d_unique_loc(d_unique_loci))),"μM, Gibbs Model = ",string(Mvec(m_unique_loc(m_unique_loci))),",Temperature = ",string(Tvec(t_unique_loc(t_unique_loci))-273.15),"°C"))
                fprintf('\n')
                fprintf(strcat("Cell Radius = ",string(Rcell),"μm"))
                fprintf('\n')
                fprintf('\n')
                % if (settings.BLASTdna)
                %     DNA_IDs = find(ismember(Names,settings.DNAdbParser));%IDs
                % else
                %     DNA_IDs = [];
                % end
                % if (settings.BLASTrna)
                %     NonDNA_IDs = find(ismember(Names,settings.RNAdbParser));%IDs
                % else
                %     NonDNA_IDs =[];
                % end
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
                Non_P(isnan(Non_P)) = 0;
                Noff_P(isnan(Noff_P)) = 0;
                Nother_P(isnan(Nother_P)) = 0;
                Non_P = Non_P(1:find(squeeze(sum(Non_P,2))>0, 1, 'last'),:);
                Nother_P = Nother_P(1:find(squeeze(sum(Nother_P,2))>0, 1, 'last'),:);
                Noff_P = Noff_P(1:find(squeeze(sum(Noff_P,2))>0, 1, 'last'),:);
                output_exists = 0;
                if (isempty(Non_P))
                    Non_P = ndSparse.build([max([size(Non_P,1) size(Noff_P,1) size(Nother_P,1)]) length(Cvec)],0);
                    output_exists = 1;
                end
                if (isempty(Nother_P))
                    Nother_P = ndSparse.build([max([size(Non_P,1) size(Noff_P,1) size(Nother_P,1)]) length(Cvec)],0);
                    output_exists = 1;
                end
                if (isempty(Noff_P))
                    Noff_P = ndSparse.build([max([size(Non_P,1) size(Noff_P,1) size(Nother_P,1)])  length(Cvec)],0);
                    output_exists = 1;
                end
                if (output_exists)
                    Basic_Noff = sum(Noff_P.*(0:size(Noff_P,1)-1)',1);
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
                    fprintf("Computing predicted probe set spot intensity and spot detection model metrics")
                    fprintf('\n')
                    fprintf(strcat("Initial Probe Concentration = ",string(PC0*Dvec(d_unique_loc(d_unique_loci))),"μM, Gibbs Model = ",string(Mvec(m_unique_loc(m_unique_loci))),",Temperature = ",string(Tvec(t_unique_loc(t_unique_loci))-273.15),"°C"))
                    fprintf('\n')
                    fprintf(strcat("Reference Mean Autofluorescence Background Intensity = ",string(AutoBackground_Mean)," ","a.u."))
                    fprintf('\n')
                    fprintf(strcat("Reference Autofluorescence Background Intensity Standard Deviation = ",string(AutoBackground_STD)))
                    fprintf('\n')
                    fprintf(strcat("Reference Mean Single Probe Intensity = ",string(IntensityPerProbe_Mean)," ","a.u."))
                    fprintf('\n')
                    fprintf(strcat("Reference Single Probe Intensity Standard Deviation = ",string(IntensityPerProbe_STD)))
                    fprintf('\n')
                    fprintf(strcat("Cell Pixel Diameter = ",string(Mean_Diameter),"px"))
                    fprintf('\n')
                    fprintf(strcat("Spot Radius = ",string(Rspot),"px"))
                    fprintf('\n')
                    fprintf('\n')
                    %Intensity Calculations
                    Ioff_Pixel = repmat(Basic_Noff,[length(Diameter_vals) 1])*IntensityPerProbe_Mean/Nstacks*Rspot^2./(pi/4*Diameter_vals'.^2);%for all diameter entries
                    Pn_XtoI_ON =  [CATnWrapper(arrayfun(@(x) nProbeIntensity(SI,SpotIntensity_Mean,SpotIntensity_STD,NumReferenceProbes,x),1:size(Non_P,1)-1,'Un',0),1)];
                    Pn_XtoI_OTHER =  [CATnWrapper(arrayfun(@(x) nProbeIntensity(SI,SpotIntensity_Mean,SpotIntensity_STD,NumReferenceProbes,x),1:size(Nother_P,1)-1,'Un',0),1)];
                    Pn_XtoI_OFF =  [CATnWrapper(arrayfun(@(x) nProbeIntensity(SI,SpotIntensity_Mean,SpotIntensity_STD,NumReferenceProbes,x),1:size(Noff_P,1)-1,'Un',0),1)];
                    Pon_I = permute(squeeze(sum(permute(repmat(Non_P(2:end,:)',[1 1 length(SI)]),[2 3 1]).*Pn_XtoI_ON,1))./sum(Non_P(2:end,:),1),[2 1]);%I over SI
                    Pon_I(isnan(Pon_I)) = 0;
                    if (sum(ismember(sum(Non_P,1),0))>0)
                        Pon_I(sum(Non_P,1)==0,1) = 1;
                    end
                    Non_I = sum(Non_P(2:end,:),1)'.*Pon_I;
                    if (~isempty(Nother_P))
                        Pother_I = permute(squeeze(sum(permute(repmat(Nother_P(2:end,:)',[1 1 length(SI)]),[2 3 1]).*Pn_XtoI_OTHER,1))./sum(Nother_P(2:end,:),1),[2 1]);%I over SI
                        Pother_I(isnan(Pother_I)) = 0;
                        if (sum(ismember(sum(Nother_P,1),0))>0)
                            Pother_I(sum(Nother_P,1)==0,1) = 1;
                        end
                        Nother_I = sum(Nother_P(2:end,:),1)'.*Pother_I;
                    else
                        Pother_I = [];
                        Nother_I = [];
                    end
                    Poff_I = permute(squeeze(sum(permute(repmat(Noff_P(2:end,:)',[1 1 length(SI)]),[2 3 1]).*Pn_XtoI_OFF,1))./sum(Noff_P(2:end,:),1),[2 1]);%I over SI
                    Poff_I(isnan(Poff_I)) = 0;
                    if (sum(ismember(sum(Noff_P,1),0))>0)
                        Poff_I(sum(Noff_P,1)==0,1) = 1;
                    end
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
                    Positive_SpotCounts_IsoIgnorant_I_Matrix = IsoIgnorantConfusion_I.PP;
                    F1_ScoreCurve_IsoIgnorant_I_Matrix = IsoIgnorantConfusion_I.F1;
                    F2_ScoreCurve_IsoIgnorant_I_Matrix =  IsoIgnorantConfusion_I.Fbeta(2);
                    Fhalf_ScoreCurve_IsoIgnorant_I_Matrix =  IsoIgnorantConfusion_I.Fbeta(0.5);
                    MCC_ScoreCurve_IsoIgnorant_I_Matrix = IsoIgnorantConfusion_I.MCC;
                    P4_ScoreCurve_IsoIgnorant_I_Matrix = IsoIgnorantConfusion_I.P4;
                    %filtered curves
                    IsoIgnorant_I_Filtered_SpotCountCurves = arrayfun(@(nth_cell)  unique(round(Positive_SpotCounts_IsoIgnorant_I_Matrix(nth_cell,:)),'stable'),1:length(Cvec),'Un',0);
                    IsoIgnorant_I_Filtered_SpotCountIndexes = arrayfun(@(nth_cell) arrayfun(@(z) find(round(Positive_SpotCounts_IsoIgnorant_I_Matrix(nth_cell,:))==z,1),IsoIgnorant_I_Filtered_SpotCountCurves{nth_cell}),1:length(Cvec),'Un',0);
                    % %unfiltered curves
                    % IsoIgnorant_I_Filtered_SpotCountCurves = arrayfun(@(nth_cell)  Positive_SpotCounts_IsoIgnorant_I_Matrix(nth_cell,:),1:length(Cvec),'Un',0);
                    % IsoIgnorant_I_Filtered_SpotCountIndexes = arrayfun(@(nth_cell) 1:size(Positive_SpotCounts_IsoIgnorant_I_Matrix,2),1:length(Cvec),'Un',0);
                    IsoIgnorant_I_param_struct_vector = arrayfun(@(x) RNAThreshold.genEmptyThresholdParamStruct(),1:length(Cvec),'Un',0);
                    for nn = 1:length(Cvec)
                        if (loaded_TrueSpot_default_params)
                            IsoIgnorant_I_param_struct_vector{nn} = TrueSpotDefaultThParameters;
                        end
                        IsoIgnorant_I_param_struct_vector{nn}.sample_spot_table = [(1:length(IsoIgnorant_I_Filtered_SpotCountCurves{nn}))' IsoIgnorant_I_Filtered_SpotCountCurves{nn}'];
                    end
                    IsoIgnorant_I_scThresholdSuggestions = arrayfun(@(nn) RNAThreshold.scoreThresholdSuggestions(RNAThreshold.estimateThreshold(IsoIgnorant_I_param_struct_vector{nn})),1:length(Cvec),'Un',0);
                    IsoIgnorant_I_scThresholdSuggestions =  [IsoIgnorant_I_scThresholdSuggestions{:}];
                    subfield_groups = {'pool','thstats'};
                    for v = 1:length(subfield_groups)
                        subfields = fieldnames(IsoIgnorant_I_scThresholdSuggestions(1).(subfield_groups{v}));
                        for subf = 1:length(subfields)
                            subf_vals = arrayfun(@(nn) IsoIgnorant_I_scThresholdSuggestions(nn).(subfield_groups{v}).(subfields{subf}),1:size(IsoIgnorant_I_scThresholdSuggestions,2),'UniformOutput',0);
                            [IsoIgnorant_I_scThresholdSuggestions.(subfields{subf})] = subf_vals{:};
                        end
                    end
                    IsoIgnorant_I_scThresholdSuggestions = rmfield(IsoIgnorant_I_scThresholdSuggestions,subfield_groups);
                    IsoIgnorant_I_TrueSpot_FilteredThreshold = [IsoIgnorant_I_scThresholdSuggestions.threshold];
                    IsoIgnorant_I_TrueSpot_FilteredThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell}(IsoIgnorant_I_TrueSpot_FilteredThreshold(nth_cell)),1:length(Cvec),'Un',0));
                    IsoIgnorant_I_TrueSpot_SuggFilteredThreshold = vertcat(IsoIgnorant_I_scThresholdSuggestions.sugg_m);
                    IsoIgnorant_I_TrueSpot_SuggFilteredThresholdLocations = CATnWrapper(arrayfun(@(nth_cell) IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell}(IsoIgnorant_I_TrueSpot_SuggFilteredThreshold(nth_cell,:)),1:length(Cvec),'Un',0),1);
                    IsoIgnorant_I_F1_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(F1_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})==max(F1_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                    IsoIgnorant_I_F1_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell}(find(F1_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})==max(F1_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                    IsoIgnorant_I_F2_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(F2_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})==max(F2_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                    IsoIgnorant_I_F2_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell}(find(F2_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})==max(F2_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                    IsoIgnorant_I_Fhalf_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(Fhalf_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})==max(Fhalf_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                    IsoIgnorant_I_Fhalf_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell}(find(Fhalf_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})==max(Fhalf_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                    IsoIgnorant_I_P4_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(P4_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})==max(P4_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                    IsoIgnorant_I_P4_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell}(find(P4_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})==max(P4_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                    IsoIgnorant_I_MCC_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(MCC_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})==max(MCC_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                    IsoIgnorant_I_MCC_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell}(find(MCC_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})==max(MCC_ScoreCurve_IsoIgnorant_I_Matrix(nth_cell,IsoIgnorant_I_Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                    IsoIgnorant_I_TS_CountMetricFunction = @(k) IsoIgnorantConfusion_I.(vars{k})(sub2ind(size(IsoIgnorantConfusion_I.(vars{k})),1:length(Cvec),IsoIgnorant_I_TrueSpot_FilteredThresholdLocations));
                    IsoIgnorant_I_TSS_CountMetricFunction = @(k) reshape(IsoIgnorantConfusion_I.(vars{k})(sub2ind(size(IsoIgnorantConfusion_I.(vars{k})),...
                        reshape(ones(size(IsoIgnorant_I_TrueSpot_SuggFilteredThresholdLocations)).*[1:length(Cvec)]',1,[]),...
                        reshape(IsoIgnorant_I_TrueSpot_SuggFilteredThresholdLocations,1,[]))),[length(Cvec) size(IsoIgnorant_I_TrueSpot_SuggFilteredThreshold,2)]);
                    IsoIgnorant_I_F1_CountMetricFunction = @(k) IsoIgnorantConfusion_I.(vars{k})(sub2ind(size(IsoIgnorantConfusion_I.(vars{k})),1:length(Cvec),IsoIgnorant_I_F1_ThresholdLocations));
                    IsoIgnorant_I_F2_CountMetricFunction = @(k) IsoIgnorantConfusion_I.(vars{k})(sub2ind(size(IsoIgnorantConfusion_I.(vars{k})),1:length(Cvec),IsoIgnorant_I_F2_ThresholdLocations));
                    IsoIgnorant_I_Fhalf_CountMetricFunction = @(k) IsoIgnorantConfusion_I.(vars{k})(sub2ind(size(IsoIgnorantConfusion_I.(vars{k})),1:length(Cvec),IsoIgnorant_I_Fhalf_ThresholdLocations));
                    IsoIgnorant_I_P4_CountMetricFunction = @(k) IsoIgnorantConfusion_I.(vars{k})(sub2ind(size(IsoIgnorantConfusion_I.(vars{k})),1:length(Cvec),IsoIgnorant_I_P4_ThresholdLocations));
                    IsoIgnorant_I_MCC_CountMetricFunction = @(k) IsoIgnorantConfusion_I.(vars{k})(sub2ind(size(IsoIgnorantConfusion_I.(vars{k})),1:length(Cvec),IsoIgnorant_I_MCC_ThresholdLocations));
                    IsoIgnorant_SpotCountMetrics = [];
                    IsoIgnorant_SpotCountMetrics.TrueSpot.Thresholds = IsoIgnorant_I_TrueSpot_FilteredThreshold;
                    IsoIgnorant_SpotCountMetrics.TrueSpot_Sugg.Thresholds = IsoIgnorant_I_TrueSpot_SuggFilteredThreshold;
                    IsoIgnorant_SpotCountMetrics.F1.Thresholds = IsoIgnorant_I_F1_FilteredThreshold;
                    IsoIgnorant_SpotCountMetrics.F2.Thresholds = IsoIgnorant_I_F2_FilteredThreshold;
                    IsoIgnorant_SpotCountMetrics.Fhalf.Thresholds = IsoIgnorant_I_Fhalf_FilteredThreshold;
                    IsoIgnorant_SpotCountMetrics.MCC.Thresholds = IsoIgnorant_I_MCC_FilteredThreshold;
                    IsoIgnorant_SpotCountMetrics.P4.Thresholds = IsoIgnorant_I_P4_FilteredThreshold;
                    for k = 1:length(vars)
                        IsoIgnorant_SpotCountMetrics.TrueSpot.(vars{k}) = IsoIgnorant_I_TS_CountMetricFunction(k);
                        IsoIgnorant_SpotCountMetrics.TrueSpot_Sugg.(vars{k}) = IsoIgnorant_I_TSS_CountMetricFunction(k);
                        IsoIgnorant_SpotCountMetrics.F1.(vars{k}) = IsoIgnorant_I_F1_CountMetricFunction(k);
                        IsoIgnorant_SpotCountMetrics.F2.(vars{k}) = IsoIgnorant_I_F2_CountMetricFunction(k);
                        IsoIgnorant_SpotCountMetrics.Fhalf.(vars{k}) = IsoIgnorant_I_Fhalf_CountMetricFunction(k);
                        IsoIgnorant_SpotCountMetrics.MCC.(vars{k}) = IsoIgnorant_I_MCC_CountMetricFunction(k);
                        IsoIgnorant_SpotCountMetrics.P4.(vars{k}) = IsoIgnorant_I_P4_CountMetricFunction(k);
                    end
                    ProbeSetMetrics.CountPredictions.IsoIgnorant_SpotCountMetrics_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = IsoIgnorant_SpotCountMetrics;
                    if (~isempty(Nother_P))
                        IsoSpecificConfusion_I = confusionMatrixWrapper_MultiCell(Non_I,Noff_I+Nother_I);
                        IsoAgnosticConfusion_I = confusionMatrixWrapper_MultiCell(Non_I+Nother_I,Noff_I);
                        Positive_SpotCounts_IsoSpecific_I_Matrix = IsoSpecificConfusion_I.PP;
                        F1_ScoreCurve_IsoSpecific_I_Matrix = IsoSpecificConfusion_I.F1;
                        F2_ScoreCurve_IsoSpecific_I_Matrix =  IsoSpecificConfusion_I.Fbeta(2);
                        Fhalf_ScoreCurve_IsoSpecific_I_Matrix =  IsoSpecificConfusion_I.Fbeta(0.5);
                        MCC_ScoreCurve_IsoSpecific_I_Matrix = IsoSpecificConfusion_I.MCC;
                        P4_ScoreCurve_IsoSpecific_I_Matrix = IsoSpecificConfusion_I.P4;
                        IsoSpecific_I_Filtered_SpotCountCurves = arrayfun(@(nth_cell)  unique(round(Positive_SpotCounts_IsoSpecific_I_Matrix(nth_cell,:)),'stable'),1:length(Cvec),'Un',0);
                        IsoSpecific_I_Filtered_SpotCountIndexes = arrayfun(@(nth_cell) arrayfun(@(z) find(round(Positive_SpotCounts_IsoSpecific_I_Matrix(nth_cell,:))==z,1),IsoSpecific_I_Filtered_SpotCountCurves{nth_cell}),1:length(Cvec),'Un',0);
                        IsoSpecific_I_param_struct_vector = arrayfun(@(x) RNAThreshold.genEmptyThresholdParamStruct(),1:length(Cvec),'Un',0);
                        for nn = 1:length(Cvec)
                            if (loaded_TrueSpot_default_params)
                                IsoSpecific_I_param_struct_vector{nn} = TrueSpotDefaultThParameters;
                            end
                            IsoSpecific_I_param_struct_vector{nn}.sample_spot_table = [(1:length(IsoSpecific_I_Filtered_SpotCountCurves{nn}))' IsoSpecific_I_Filtered_SpotCountCurves{nn}'];
                        end
                        IsoSpecific_I_scThresholdSuggestions = arrayfun(@(nn) RNAThreshold.scoreThresholdSuggestions(RNAThreshold.estimateThreshold(IsoSpecific_I_param_struct_vector{nn})),1:length(Cvec),'Un',0);
                        IsoSpecific_I_scThresholdSuggestions =  [IsoSpecific_I_scThresholdSuggestions{:}];
                        subfield_groups = {'pool','thstats'};
                        for v = 1:length(subfield_groups)
                            subfields = fieldnames(IsoSpecific_I_scThresholdSuggestions(1).(subfield_groups{v}));
                            for subf = 1:length(subfields)
                                subf_vals = arrayfun(@(nn) IsoSpecific_I_scThresholdSuggestions(nn).(subfield_groups{v}).(subfields{subf}),1:size(IsoSpecific_I_scThresholdSuggestions,2),'UniformOutput',0);
                                [IsoSpecific_I_scThresholdSuggestions.(subfields{subf})] = subf_vals{:};
                            end
                        end
                        IsoSpecific_I_scThresholdSuggestions = rmfield(IsoSpecific_I_scThresholdSuggestions,subfield_groups);
                        IsoSpecific_I_TrueSpot_FilteredThreshold = [IsoSpecific_I_scThresholdSuggestions.threshold];
                        IsoSpecific_I_TrueSpot_FilteredThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell}(IsoSpecific_I_TrueSpot_FilteredThreshold(nth_cell)),1:length(Cvec),'Un',0));
                        IsoSpecific_I_TrueSpot_SuggFilteredThreshold = vertcat(IsoSpecific_I_scThresholdSuggestions.sugg_m);
                        IsoSpecific_I_TrueSpot_SuggFilteredThresholdLocations = CATnWrapper(arrayfun(@(nth_cell) IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell}(IsoSpecific_I_TrueSpot_SuggFilteredThreshold(nth_cell,:)),1:length(Cvec),'Un',0),1);
                        IsoSpecific_I_F1_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(F1_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})==max(F1_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                        IsoSpecific_I_F1_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell}(find(F1_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})==max(F1_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                        IsoSpecific_I_F2_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(F2_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})==max(F2_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                        IsoSpecific_I_F2_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell}(find(F2_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})==max(F2_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                        IsoSpecific_I_Fhalf_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(Fhalf_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})==max(Fhalf_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                        IsoSpecific_I_Fhalf_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell}(find(Fhalf_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})==max(Fhalf_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                        IsoSpecific_I_P4_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(P4_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})==max(P4_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                        IsoSpecific_I_P4_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell}(find(P4_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})==max(P4_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                        IsoSpecific_I_MCC_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(MCC_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})==max(MCC_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                        IsoSpecific_I_MCC_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell}(find(MCC_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})==max(MCC_ScoreCurve_IsoSpecific_I_Matrix(nth_cell,IsoSpecific_I_Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                        IsoSpecific_I_TS_CountMetricFunction = @(k) IsoSpecificConfusion_I.(vars{k})(sub2ind(size(IsoSpecificConfusion_I.(vars{k})),1:length(Cvec),IsoSpecific_I_TrueSpot_FilteredThresholdLocations));
                        IsoSpecific_I_TSS_CountMetricFunction = @(k) reshape(IsoSpecificConfusion_I.(vars{k})(sub2ind(size(IsoSpecificConfusion_I.(vars{k})),...
                            reshape(ones(size(IsoSpecific_I_TrueSpot_SuggFilteredThresholdLocations)).*[1:length(Cvec)]',1,[]),...
                            reshape(IsoSpecific_I_TrueSpot_SuggFilteredThresholdLocations,1,[]))),[length(Cvec) size(IsoSpecific_I_TrueSpot_SuggFilteredThreshold,2)]);
                        IsoSpecific_I_F1_CountMetricFunction = @(k) IsoSpecificConfusion_I.(vars{k})(sub2ind(size(IsoSpecificConfusion_I.(vars{k})),1:length(Cvec),IsoSpecific_I_F1_ThresholdLocations));
                        IsoSpecific_I_F2_CountMetricFunction = @(k) IsoSpecificConfusion_I.(vars{k})(sub2ind(size(IsoSpecificConfusion_I.(vars{k})),1:length(Cvec),IsoSpecific_I_F2_ThresholdLocations));
                        IsoSpecific_I_Fhalf_CountMetricFunction = @(k) IsoSpecificConfusion_I.(vars{k})(sub2ind(size(IsoSpecificConfusion_I.(vars{k})),1:length(Cvec),IsoSpecific_I_Fhalf_ThresholdLocations));
                        IsoSpecific_I_P4_CountMetricFunction = @(k) IsoSpecificConfusion_I.(vars{k})(sub2ind(size(IsoSpecificConfusion_I.(vars{k})),1:length(Cvec),IsoSpecific_I_P4_ThresholdLocations));
                        IsoSpecific_I_MCC_CountMetricFunction = @(k) IsoSpecificConfusion_I.(vars{k})(sub2ind(size(IsoSpecificConfusion_I.(vars{k})),1:length(Cvec),IsoSpecific_I_MCC_ThresholdLocations));
                        IsoSpecific_SpotCountMetrics = [];
                        IsoSpecific_SpotCountMetrics.TrueSpot.Thresholds = IsoSpecific_I_TrueSpot_FilteredThreshold;
                        IsoSpecific_SpotCountMetrics.TrueSpot_Sugg.Thresholds = IsoSpecific_I_TrueSpot_SuggFilteredThreshold;
                        IsoSpecific_SpotCountMetrics.F1.Thresholds = IsoSpecific_I_F1_FilteredThreshold;
                        IsoSpecific_SpotCountMetrics.F2.Thresholds = IsoSpecific_I_F2_FilteredThreshold;
                        IsoSpecific_SpotCountMetrics.Fhalf.Thresholds = IsoSpecific_I_Fhalf_FilteredThreshold;
                        IsoSpecific_SpotCountMetrics.MCC.Thresholds = IsoSpecific_I_MCC_FilteredThreshold;
                        IsoSpecific_SpotCountMetrics.P4.Thresholds = IsoSpecific_I_P4_FilteredThreshold;
                        for k = 1:length(vars)
                            IsoSpecific_SpotCountMetrics.TrueSpot.(vars{k}) = IsoSpecific_I_TS_CountMetricFunction(k);
                            IsoSpecific_SpotCountMetrics.TrueSpot_Sugg.(vars{k}) = IsoSpecific_I_TSS_CountMetricFunction(k);
                            IsoSpecific_SpotCountMetrics.F1.(vars{k}) = IsoSpecific_I_F1_CountMetricFunction(k);
                            IsoSpecific_SpotCountMetrics.F2.(vars{k}) = IsoSpecific_I_F2_CountMetricFunction(k);
                            IsoSpecific_SpotCountMetrics.Fhalf.(vars{k}) = IsoSpecific_I_Fhalf_CountMetricFunction(k);
                            IsoSpecific_SpotCountMetrics.MCC.(vars{k}) = IsoSpecific_I_MCC_CountMetricFunction(k);
                            IsoSpecific_SpotCountMetrics.P4.(vars{k}) = IsoSpecific_I_P4_CountMetricFunction(k);
                        end
                        ProbeSetMetrics.CountPredictions.IsoSpecific_SpotCountMetrics_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = IsoSpecific_SpotCountMetrics;
                        Positive_SpotCounts_IsoAgnostic_I_Matrix = IsoAgnosticConfusion_I.PP;
                        F1_ScoreCurve_IsoAgnostic_I_Matrix = IsoAgnosticConfusion_I.F1;
                        F2_ScoreCurve_IsoAgnostic_I_Matrix =  IsoAgnosticConfusion_I.Fbeta(2);
                        Fhalf_ScoreCurve_IsoAgnostic_I_Matrix =  IsoAgnosticConfusion_I.Fbeta(0.5);
                        MCC_ScoreCurve_IsoAgnostic_I_Matrix = IsoAgnosticConfusion_I.MCC;
                        P4_ScoreCurve_IsoAgnostic_I_Matrix = IsoAgnosticConfusion_I.P4;
                        IsoAgnostic_I_Filtered_SpotCountCurves = arrayfun(@(nth_cell)  unique(round(Positive_SpotCounts_IsoAgnostic_I_Matrix(nth_cell,:)),'stable'),1:length(Cvec),'Un',0);
                        IsoAgnostic_I_Filtered_SpotCountIndexes = arrayfun(@(nth_cell) arrayfun(@(z) find(round(Positive_SpotCounts_IsoAgnostic_I_Matrix(nth_cell,:))==z,1),IsoAgnostic_I_Filtered_SpotCountCurves{nth_cell}),1:length(Cvec),'Un',0);
                        IsoAgnostic_I_param_struct_vector = arrayfun(@(x) RNAThreshold.genEmptyThresholdParamStruct(),1:length(Cvec),'Un',0);
                        for nn = 1:length(Cvec)
                            if (loaded_TrueSpot_default_params)
                                IsoAgnostic_I_param_struct_vector{nn} = TrueSpotDefaultThParameters;
                            end
                            IsoAgnostic_I_param_struct_vector{nn}.sample_spot_table = [(1:length(IsoAgnostic_I_Filtered_SpotCountCurves{nn}))' IsoAgnostic_I_Filtered_SpotCountCurves{nn}'];
                        end
                        IsoAgnostic_I_scThresholdSuggestions = arrayfun(@(nn) RNAThreshold.scoreThresholdSuggestions(RNAThreshold.estimateThreshold(IsoAgnostic_I_param_struct_vector{nn})),1:length(Cvec),'Un',0);
                        IsoAgnostic_I_scThresholdSuggestions =  [IsoAgnostic_I_scThresholdSuggestions{:}];
                        subfield_groups = {'pool','thstats'};
                        for v = 1:length(subfield_groups)
                            subfields = fieldnames(IsoAgnostic_I_scThresholdSuggestions(1).(subfield_groups{v}));
                            for subf = 1:length(subfields)
                                subf_vals = arrayfun(@(nn) IsoAgnostic_I_scThresholdSuggestions(nn).(subfield_groups{v}).(subfields{subf}),1:size(IsoAgnostic_I_scThresholdSuggestions,2),'UniformOutput',0);
                                [IsoAgnostic_I_scThresholdSuggestions.(subfields{subf})] = subf_vals{:};
                            end
                        end
                        IsoAgnostic_I_scThresholdSuggestions = rmfield(IsoAgnostic_I_scThresholdSuggestions,subfield_groups);
                        IsoAgnostic_I_TrueSpot_FilteredThreshold = [IsoAgnostic_I_scThresholdSuggestions.threshold];
                        IsoAgnostic_I_TrueSpot_FilteredThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell}(IsoAgnostic_I_TrueSpot_FilteredThreshold(nth_cell)),1:length(Cvec),'Un',0));
                        IsoAgnostic_I_TrueSpot_SuggFilteredThreshold = vertcat(IsoAgnostic_I_scThresholdSuggestions.sugg_m);
                        IsoAgnostic_I_TrueSpot_SuggFilteredThresholdLocations = CATnWrapper(arrayfun(@(nth_cell) IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell}(IsoAgnostic_I_TrueSpot_SuggFilteredThreshold(nth_cell,:)),1:length(Cvec),'Un',0),1);
                        IsoAgnostic_I_F1_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(F1_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})==max(F1_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                        IsoAgnostic_I_F1_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell}(find(F1_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})==max(F1_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                        IsoAgnostic_I_F2_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(F2_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})==max(F2_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                        IsoAgnostic_I_F2_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell}(find(F2_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})==max(F2_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                        IsoAgnostic_I_Fhalf_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(Fhalf_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})==max(Fhalf_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                        IsoAgnostic_I_Fhalf_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell}(find(Fhalf_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})==max(Fhalf_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                        IsoAgnostic_I_P4_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(P4_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})==max(P4_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                        IsoAgnostic_I_P4_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell}(find(P4_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})==max(P4_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                        IsoAgnostic_I_MCC_FilteredThreshold = cell2mat(arrayfun(@(nth_cell) find(MCC_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})==max(MCC_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})),1),1:length(Cvec),'Un',0));
                        IsoAgnostic_I_MCC_ThresholdLocations = cell2mat(arrayfun(@(nth_cell) IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell}(find(MCC_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})==max(MCC_ScoreCurve_IsoAgnostic_I_Matrix(nth_cell,IsoAgnostic_I_Filtered_SpotCountIndexes{nth_cell})),1)),1:length(Cvec),'Un',0));
                        IsoAgnostic_I_TS_CountMetricFunction = @(k) IsoAgnosticConfusion_I.(vars{k})(sub2ind(size(IsoAgnosticConfusion_I.(vars{k})),1:length(Cvec),IsoAgnostic_I_TrueSpot_FilteredThresholdLocations));
                        IsoAgnostic_I_TSS_CountMetricFunction = @(k) reshape(IsoAgnosticConfusion_I.(vars{k})(sub2ind(size(IsoAgnosticConfusion_I.(vars{k})),...
                            reshape(ones(size(IsoAgnostic_I_TrueSpot_SuggFilteredThresholdLocations)).*[1:length(Cvec)]',1,[]),...
                            reshape(IsoAgnostic_I_TrueSpot_SuggFilteredThresholdLocations,1,[]))),[length(Cvec) size(IsoAgnostic_I_TrueSpot_SuggFilteredThreshold,2)]);
                        IsoAgnostic_I_F1_CountMetricFunction = @(k) IsoAgnosticConfusion_I.(vars{k})(sub2ind(size(IsoAgnosticConfusion_I.(vars{k})),1:length(Cvec),IsoAgnostic_I_F1_ThresholdLocations));
                        IsoAgnostic_I_F2_CountMetricFunction = @(k) IsoAgnosticConfusion_I.(vars{k})(sub2ind(size(IsoAgnosticConfusion_I.(vars{k})),1:length(Cvec),IsoAgnostic_I_F2_ThresholdLocations));
                        IsoAgnostic_I_Fhalf_CountMetricFunction = @(k) IsoAgnosticConfusion_I.(vars{k})(sub2ind(size(IsoAgnosticConfusion_I.(vars{k})),1:length(Cvec),IsoAgnostic_I_Fhalf_ThresholdLocations));
                        IsoAgnostic_I_P4_CountMetricFunction = @(k) IsoAgnosticConfusion_I.(vars{k})(sub2ind(size(IsoAgnosticConfusion_I.(vars{k})),1:length(Cvec),IsoAgnostic_I_P4_ThresholdLocations));
                        IsoAgnostic_I_MCC_CountMetricFunction = @(k) IsoAgnosticConfusion_I.(vars{k})(sub2ind(size(IsoAgnosticConfusion_I.(vars{k})),1:length(Cvec),IsoAgnostic_I_MCC_ThresholdLocations));
                        IsoAgnostic_SpotCountMetrics = [];
                        IsoAgnostic_SpotCountMetrics.TrueSpot.Thresholds = IsoAgnostic_I_TrueSpot_FilteredThreshold;
                        IsoAgnostic_SpotCountMetrics.TrueSpot_Sugg.Thresholds = IsoAgnostic_I_TrueSpot_SuggFilteredThreshold;
                        IsoAgnostic_SpotCountMetrics.F1.Thresholds = IsoAgnostic_I_F1_FilteredThreshold;
                        IsoAgnostic_SpotCountMetrics.F2.Thresholds = IsoAgnostic_I_F2_FilteredThreshold;
                        IsoAgnostic_SpotCountMetrics.Fhalf.Thresholds = IsoAgnostic_I_Fhalf_FilteredThreshold;
                        IsoAgnostic_SpotCountMetrics.MCC.Thresholds = IsoAgnostic_I_MCC_FilteredThreshold;
                        IsoAgnostic_SpotCountMetrics.P4.Thresholds = IsoAgnostic_I_P4_FilteredThreshold;
                        for k = 1:length(vars)
                            IsoAgnostic_SpotCountMetrics.TrueSpot.(vars{k}) = IsoAgnostic_I_TS_CountMetricFunction(k);
                            IsoAgnostic_SpotCountMetrics.TrueSpot_Sugg.(vars{k}) = IsoAgnostic_I_TSS_CountMetricFunction(k);
                            IsoAgnostic_SpotCountMetrics.F1.(vars{k}) = IsoAgnostic_I_F1_CountMetricFunction(k);
                            IsoAgnostic_SpotCountMetrics.F2.(vars{k}) = IsoAgnostic_I_F2_CountMetricFunction(k);
                            IsoAgnostic_SpotCountMetrics.Fhalf.(vars{k}) = IsoAgnostic_I_Fhalf_CountMetricFunction(k);
                            IsoAgnostic_SpotCountMetrics.MCC.(vars{k}) = IsoAgnostic_I_MCC_CountMetricFunction(k);
                            IsoAgnostic_SpotCountMetrics.P4.(vars{k}) = IsoAgnostic_I_P4_CountMetricFunction(k);
                        end
                        ProbeSetMetrics.CountPredictions.IsoAgnostic_SpotCountMetrics_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = IsoAgnostic_SpotCountMetrics;
                    else
                        ProbeSetMetrics.CountPredictions.IsoSpecific_SpotCountMetrics_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                        ProbeSetMetrics.CountPredictions.IsoAgnostic_SpotCountMetrics_ModelTemperatureDilutionVector{m_unique_loci,t_unique_loci,d_unique_loci} = [];
                    end

                end
            end
        end
    end
end
ModelMetrics.ProbeSetMetrics = ProbeSetMetrics;
end
function x = deployable_heaviside(x)
x(x>=0) = 1;
x(x<0) = 0;
end