function [CProbes_Free,varSSE,err,i,eqSSE] = A_ModelEquilibriumSolverWrapper_V4(ModelSolverFunctions,MaxIter,errThreshold,Pset,nExpressionMatrix,Tvec,Mvec,Dvec,Cvec,Tref,Ns_Config,Nc_Config,CProbes_Free0,CProbes_Free,ProbeConc,Js_RNA,Js_DNA,Js_Sites,i)
%% This function iteratively updates probe binding equilibrium values using fixed point iteration
K_S = ModelSolverFunctions.K_S;
K_CDeff = ModelSolverFunctions.K_CDeff;
h_RNA = ModelSolverFunctions.h_RNA;
h_DNA = ModelSolverFunctions.h_DNA;
if (ModelSolverFunctions.solverType==0)
    tsList_RNA = ModelSolverFunctions.Paired_RNATargetSites;
    tsList_DNA = ModelSolverFunctions.Paired_DNATargetSites;    
end
x1 = squeeze(CProbes_Free);
switch ModelSolverFunctions.solverType
    case 0
        A1 = squeeze(K_S(Tvec,Tref,Pset,1:Ns_Config,Mvec,Cvec));
        A1(isnan(A1)) = 0;
        A2 = squeeze(K_CDeff(CProbes_Free,Tvec,Tref,Pset,Pset,1:Nc_Config,Mvec,Cvec));
        A2(isnan(A2)) = 0;
        A3 = squeeze(h_RNA(CProbes_Free,Tvec,Tref,Pset,tsList_RNA,Mvec,Cvec,nExpressionMatrix));
        A3(isnan(A3)) = 0;
        A4 = squeeze(h_DNA(CProbes_Free,Tvec,Tref,Pset,tsList_DNA,Mvec,Cvec,nExpressionMatrix));%7-D
        A4(isnan(A4)) = 0;
    case 1
        A1 = squeeze(K_S(Tvec,Tref,Pset,1:Ns_Config,Mvec,Cvec));
        A1(isnan(A1)) = 0;
        x1 = squeeze(CProbes_Free);
        A2 = squeeze(K_CDeff(CProbes_Free,Tvec,Tref,Pset,Pset,1:Nc_Config,Mvec,Cvec));
        A2(isnan(A2)) = 0;
        A3 = squeeze(h_RNA(CProbes_Free,Tvec,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),Mvec,Cvec,nExpressionMatrix));
        A3(isnan(A3)) = 0;
        A4 = squeeze(h_DNA(CProbes_Free,Tvec,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),Mvec,Cvec,nExpressionMatrix));%7-D
        A4(isnan(A4)) = 0;
    case 2
        A1 = squeeze(K_S(Tvec,Tref,Pset,1:Ns_Config,Mvec,Dvec,Cvec));
        A1(isnan(A1)) = 0;
        x1 = squeeze(CProbes_Free);
        A2 = squeeze(K_CDeff(CProbes_Free,Tvec,Tref,Pset,Pset,1:Nc_Config,Mvec,Dvec,Cvec));
        A2(isnan(A2)) = 0;
        A3 = squeeze(h_RNA(CProbes_Free,Tvec,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),Mvec,Dvec,Cvec,nExpressionMatrix));
        A3(isnan(A3)) = 0;
        A4 = squeeze(h_DNA(CProbes_Free,Tvec,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),Mvec,Dvec,Cvec,nExpressionMatrix));%7-D
        A4(isnan(A4)) = 0;
end
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
err = abs(sum(xn./x1,1,'omitnan')-size(x1,1));
z0 = x1;
z1 = xn;
z0(isnan(z0)) = 0;
z1(isnan(z1)) = 0;
varSSE = sqrt(sum((z1-z0).^2,1));
switch ModelSolverFunctions.solverType
    case 0
        A1 = squeeze(K_S(Tvec,Tref,Pset,1:Ns_Config,Mvec,Cvec));
        A1(isnan(A1)) = 0;
        A2 = squeeze(K_CDeff(CProbes_Free,Tvec,Tref,Pset,Pset,1:Nc_Config,Mvec,Cvec));
        A2(isnan(A2)) = 0;
        A3 = squeeze(h_RNA(CProbes_Free,Tvec,Tref,Pset,tsList_RNA,Mvec,Cvec,nExpressionMatrix));
        A3(isnan(A3)) = 0;
        A4 = squeeze(h_DNA(CProbes_Free,Tvec,Tref,Pset,tsList_DNA,Mvec,Cvec,nExpressionMatrix));%7-D
        A4(isnan(A4)) = 0;
    case 1
        A1 = squeeze(K_S(Tvec,Tref,Pset,1:Ns_Config,Mvec,Cvec));
        A1(isnan(A1)) = 0;
        A2 = squeeze(K_CDeff(CProbes_Free,Tvec,Tref,Pset,Pset,1:Nc_Config,Mvec,Cvec));
        A2(isnan(A2)) = 0;
        A3 = squeeze(h_RNA(CProbes_Free,Tvec,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),Mvec,Cvec,nExpressionMatrix));
        A3(isnan(A3)) = 0;
        A4 = squeeze(h_DNA(CProbes_Free,Tvec,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),Mvec,Cvec,nExpressionMatrix));%7-D
        A4(isnan(A4)) = 0;
    case 2
        A1 = squeeze(K_S(Tvec,Tref,Pset,1:Ns_Config,Mvec,Dvec,Cvec));
        A1(isnan(A1)) = 0;
        A2 = squeeze(K_CDeff(CProbes_Free,Tvec,Tref,Pset,Pset,1:Nc_Config,Mvec,Dvec,Cvec));
        A2(isnan(A2)) = 0;
        A3 = squeeze(h_RNA(CProbes_Free,Tvec,Tref,Pset,Js_RNA(Pset),Js_Sites(Pset),Mvec,Dvec,Cvec,nExpressionMatrix));
        A3(isnan(A3)) = 0;
        A4 = squeeze(h_DNA(CProbes_Free,Tvec,Tref,Pset,Js_DNA(Pset),Js_Sites(Pset),Mvec,Dvec,Cvec,nExpressionMatrix));%7-D
        A4(isnan(A4)) = 0;
end
eqSSE = sum((squeeze(ProbeConc)-squeeze(CProbes_Free).*(1+A1+A2+A3+A4)).^2,1,'omitnan');
i = i + 1;
if (double(i>MaxIter)+double(max(err)<errThreshold)==0)
    [CProbes_Free,varSSE,err,i,eqSSE] = A_ModelEquilibriumSolverWrapper_V4(ModelSolverFunctions,MaxIter,errThreshold,Pset,nExpressionMatrix,Tvec,Mvec,Dvec,Cvec,Tref,Ns_Config,Nc_Config,CProbes_Free0,CProbes_Free,ProbeConc,Js_RNA,Js_DNA,Js_Sites,i);
end
end