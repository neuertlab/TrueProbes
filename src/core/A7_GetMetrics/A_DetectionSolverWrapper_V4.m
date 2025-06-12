function [CtN,p_TargetSites_Bound,c_TargetSites_Bound] = A_DetectionSolverWrapper_V4(ModelSolverFunctions,v,Pset,settings,nExpressionMatrix,Tvec,Mvec,Dvec,Cvec,Tref,CProbes_Free,DoesProbeBindSite,Js_RNA,Js_DNA,Js_Sites,Names,ON_IDs_specific,ON_IDs_agnostic,OFF_IDs)
m = v(1);
t = v(2);
d = v(3);
Js = @(x) find(sum(squeeze(sum(DoesProbeBindSite(x,:,:),1)),2)>0);
Sx = unique(cell2mat(arrayfun(@(x) find(sum(DoesProbeBindSite(Pset,x,:),1)>0)',Js(Pset),'Un',0)));
pRNA_TargetSites_Bound_Func = ModelSolverFunctions.pRNA_TargetSites_Bound_Func;
pDNA_TargetSites_Bound_Func = ModelSolverFunctions.pDNA_TargetSites_Bound_Func;
cRNA_TargetSites_Bound_Func= ModelSolverFunctions.cRNA_TargetSites_Bound_Func;
cDNA_TargetSites_Bound_Func = ModelSolverFunctions.cDNA_TargetSites_Bound_Func;
if (ModelSolverFunctions.solverType==0)
    tsList_RNA = ModelSolverFunctions.Paired_RNATargetSites;
    tsList_DNA = ModelSolverFunctions.Paired_DNATargetSites;    
end
switch ModelSolverFunctions.solverType
    case 0
        pRNA_TargetSites_Bound = pRNA_TargetSites_Bound_Func(CProbes_Free,Tvec,Tref,Pset,tsList_RNA,Mvec,Cvec,nExpressionMatrix);
        pDNA_TargetSites_Bound = pDNA_TargetSites_Bound_Func(CProbes_Free,Tvec,Tref,Pset,tsList_DNA,Mvec,Cvec,nExpressionMatrix);
        cRNA_TargetSites_Bound = cRNA_TargetSites_Bound_Func(CProbes_Free,Tvec,Tref,Pset,tsList_RNA,Mvec,Cvec,nExpressionMatrix);
        cDNA_TargetSites_Bound = cDNA_TargetSites_Bound_Func(CProbes_Free,Tvec,Tref,Pset,tsList_DNA,Mvec,Cvec,nExpressionMatrix);
    case 1
        pRNA_TargetSites_Bound = pRNA_TargetSites_Bound_Func(CProbes_Free,Tvec,Tref,Pset,Mvec,Cvec,nExpressionMatrix);
        pDNA_TargetSites_Bound = pDNA_TargetSites_Bound_Func(CProbes_Free,Tvec,Tref,Pset,Mvec,Cvec,nExpressionMatrix);
        cRNA_TargetSites_Bound = cRNA_TargetSites_Bound_Func(CProbes_Free,Tvec,Tref,Pset,Mvec,Cvec,nExpressionMatrix);
        cDNA_TargetSites_Bound = cDNA_TargetSites_Bound_Func(CProbes_Free,Tvec,Tref,Pset,Mvec,Cvec,nExpressionMatrix);
    case 2
        pRNA_TargetSites_Bound = pRNA_TargetSites_Bound_Func(CProbes_Free,Tvec,Tref,Pset,Mvec,Dvec,Cvec,nExpressionMatrix);
        pDNA_TargetSites_Bound = pDNA_TargetSites_Bound_Func(CProbes_Free,Tvec,Tref,Pset,Mvec,Dvec,Cvec,nExpressionMatrix);
        cRNA_TargetSites_Bound = cRNA_TargetSites_Bound_Func(CProbes_Free,Tvec,Tref,Pset,Mvec,Dvec,Cvec,nExpressionMatrix);
        cDNA_TargetSites_Bound = cDNA_TargetSites_Bound_Func(CProbes_Free,Tvec,Tref,Pset,Mvec,Dvec,Cvec,nExpressionMatrix);
end
pRNA_TargetSites_Bound(isnan(pRNA_TargetSites_Bound)) = 0;
pDNA_TargetSites_Bound(isnan(pDNA_TargetSites_Bound)) = 0;
cRNA_TargetSites_Bound(isnan(cRNA_TargetSites_Bound)) = 0;
cDNA_TargetSites_Bound(isnan(cDNA_TargetSites_Bound)) = 0;
switch ModelSolverFunctions.solverType
    case 0
        p_TargetSites_Bound = ndSparse.build([length(Pset) size(DoesProbeBindSite,2) size(DoesProbeBindSite,3) length(Mvec) length(Cvec)],0);
        c_TargetSites_Bound = ndSparse.build([length(Pset) size(DoesProbeBindSite,2) size(DoesProbeBindSite,3) length(Mvec) length(Cvec)],0);
        if (~isempty(Js_RNA(Pset)))
            [pList_pRNA,subtsList_pRNA,mList_pRNA,cList_pRNA] = ind2sub(size(pRNA_TargetSites_Bound),find(pRNA_TargetSites_Bound>0));
            [pList_cRNA,subtsList_cRNA,mList_cRNA,cList_cRNA] = ind2sub(size(cRNA_TargetSites_Bound),find(cRNA_TargetSites_Bound>0));
            p_TargetSites_Bound(sub2ind(size(p_TargetSites_Bound),...
                pList_pRNA,tsList_RNA(1,subtsList_pRNA)',tsList_RNA(2,subtsList_pRNA)',...
                mList_pRNA,cList_pRNA)) = pRNA_TargetSites_Bound(sub2ind(size(pRNA_TargetSites_Bound),pList_pRNA,subtsList_pRNA,mList_pRNA,cList_pRNA));
            c_TargetSites_Bound(sub2ind(size(c_TargetSites_Bound),...
                pList_cRNA,tsList_RNA(1,subtsList_cRNA)',tsList_RNA(2,subtsList_cRNA)',...
                mList_cRNA,cList_cRNA)) = cRNA_TargetSites_Bound(sub2ind(size(cRNA_TargetSites_Bound),pList_cRNA,subtsList_cRNA,mList_cRNA,cList_cRNA));
        end
        if (~isempty(Js_DNA(Pset)))
            [pList_pDNA,subtsList_pDNA,mList_pDNA,cList_pDNA] = ind2sub(size(pDNA_TargetSites_Bound),find(pDNA_TargetSites_Bound>0));
            [pList_cDNA,subtsList_cDNA,mList_cDNA,cList_cDNA] = ind2sub(size(cDNA_TargetSites_Bound),find(cDNA_TargetSites_Bound>0));
            p_TargetSites_Bound(sub2ind(size(p_TargetSites_Bound),...
                pList_pDNA,tsList_DNA(1,subtsList_pDNA)',tsList_DNA(2,subtsList_pDNA)',...
                mList_pDNA,cList_pDNA)) = pDNA_TargetSites_Bound(sub2ind(size(pDNA_TargetSites_Bound),pList_pDNA,subtsList_pDNA,mList_pDNA,cList_pDNA));
            c_TargetSites_Bound(sub2ind(size(c_TargetSites_Bound),...
                pList_cDNA,tsList_DNA(1,subtsList_cDNA)',tsList_DNA(2,subtsList_cDNA)',...
                mList_cDNA,cList_cDNA)) = cDNA_TargetSites_Bound(sub2ind(size(cDNA_TargetSites_Bound),pList_cDNA,subtsList_cDNA,mList_cDNA,cList_cDNA));
        end
    case 1
        p_TargetSites_Bound = ndSparse.build([length(Pset) size(DoesProbeBindSite,2) size(DoesProbeBindSite,3) length(Mvec) length(Cvec)],0);
        c_TargetSites_Bound = ndSparse.build([length(Pset) size(DoesProbeBindSite,2) size(DoesProbeBindSite,3) length(Mvec) length(Cvec)],0);
        if (~isempty(Js_RNA(Pset)))
            p_TargetSites_Bound(:,Js_RNA(Pset),Js_Sites(Pset),:,:) = pRNA_TargetSites_Bound;
            c_TargetSites_Bound(:,Js_RNA(Pset),Js_Sites(Pset),:,:) = cRNA_TargetSites_Bound;
        end
        if (~isempty(Js_DNA(Pset)))
            p_TargetSites_Bound(:,Js_DNA(Pset),Js_Sites(Pset),:,:) = pDNA_TargetSites_Bound;
            c_TargetSites_Bound(:,Js_DNA(Pset),Js_Sites(Pset),:,:) = cDNA_TargetSites_Bound;
        end
    case 2
        p_TargetSites_Bound = ndSparse.build([length(Pset) size(DoesProbeBindSite,2) size(DoesProbeBindSite,3) length(Mvec) length(Tvec) length(Dvec) length(Cvec)],0);
        c_TargetSites_Bound = ndSparse.build([length(Pset) size(DoesProbeBindSite,2) size(DoesProbeBindSite,3) length(Mvec) length(Tvec) length(Dvec) length(Cvec)],0);
        if (~isempty(Js_RNA(Pset)))
            p_TargetSites_Bound(:,Js_RNA(Pset),Js_Sites(Pset),:,:,:,:,:) = pRNA_TargetSites_Bound;
            c_TargetSites_Bound(:,Js_RNA(Pset),Js_Sites(Pset),:,:,:,:,:) = cRNA_TargetSites_Bound;
        end
        if (~isempty(Js_DNA(Pset)))
            p_TargetSites_Bound(:,Js_DNA(Pset),Js_Sites(Pset),:,:,:,:,:) = pDNA_TargetSites_Bound;
            c_TargetSites_Bound(:,Js_DNA(Pset),Js_Sites(Pset),:,:,:,:,:) = cDNA_TargetSites_Bound;
        end
end
p_TargetSites_Bound(isnan(p_TargetSites_Bound)) = 0;
c_TargetSites_Bound(isnan(c_TargetSites_Bound)) = 0;
switch ModelSolverFunctions.solverType
    case 0
        if (length(Pset)>1)
            pAnyBindSiteCM0 = squeeze(sum(squeeze(p_TargetSites_Bound(:,:,Sx,m,:)),1));
        else
            pAnyBindSiteCM0 = squeeze(p_TargetSites_Bound(:,:,Sx,m,:));
        end
    case 1
        if (length(Pset)>1)
            pAnyBindSiteCM0 = squeeze(sum(squeeze(p_TargetSites_Bound(:,:,Sx,m,:)),1));
        else
            pAnyBindSiteCM0 = squeeze(p_TargetSites_Bound(:,:,Sx,m,:));
        end
    case 2
        if (length(Pset)>1)
            pAnyBindSiteCM0 = squeeze(sum(squeeze(p_TargetSites_Bound(:,:,Sx,m,t,d,:)),1));
        else
            pAnyBindSiteCM0 = squeeze(p_TargetSites_Bound(:,:,Sx,m,t,d,:));
        end
end
pAnyBindSiteCM0(pAnyBindSiteCM0>1) = 1;
tHit = find(squeeze(sum(sum(pAnyBindSiteCM0,2),3))>0);
pON_IDs = find(ismember(tHit,ON_IDs_specific));
pOFF_IDs = find(ismember(tHit,OFF_IDs));
pON_IDs_other = find(ismember(tHit,setdiff(ON_IDs_agnostic,ON_IDs_specific)));
if (~isempty(tHit))  
    Pt_ModelCellVector = F_DiscretePossionBinomialMultiCell(pAnyBindSiteCM0(tHit,:,:));
    Ct = permute(repmat(nExpressionMatrix(tHit,Cvec),[1 1 size(Pt_ModelCellVector,2)]),[1 3 2]).*Pt_ModelCellVector;
    CtON_specific = squeeze(sum(Ct(pON_IDs,:,:),1,'omitnan'));
    CtON_other = squeeze(sum(Ct(pON_IDs_other,:,:),1,'omitnan'));
    CtOFF = squeeze(sum(Ct(pOFF_IDs,:,:),1,'omitnan'));
    CtN = permute(CATnWrapper({CtON_specific,CtON_other,CtOFF},3),[3 1 2]);
end
end