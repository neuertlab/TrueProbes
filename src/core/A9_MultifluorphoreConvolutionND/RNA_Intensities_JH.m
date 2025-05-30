function IntensityMetrics = RNA_Intensities_JH(Pset,settings,probes,gene_table,ExpressionMatrix,DoesProbeBindSite,dHeq_mod,dSeq_mod,dCp_mod,dHeq_Complement,dSeq_Complement,dCp_Complement)   
R = 0.001987204259;%gas constant [Energy Kcal/mol K
kb = 1.380649*10^-23;%bolzman constant J/K
h = 6.62607015*10^-34;%planks constant J/Hz
Tref = 37+273.15;    
V_Cell = @(R) 4/3*pi*(R^3)/10^15;%um to L
Js = @(x) find(sum(squeeze(sum(DoesProbeBindSite(x,:,:),1)),2)>0);
Js_Sites = @(x) find(sum(sum(DoesProbeBindSite(x,Js(x),:),1),2)>0);

    




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

Sx = unique(cell2mat(arrayfun(@(x) find(sum(DoesProbeBindSite(Pset,x,:),1)>0)',Js(Pset),'Un',0)));
Smax = length(Sx);
Index2 = 1:Smax+1;
%ProbeSetMetrics is an CxM struct with 2 field

pnON_IDs = find(ismember(Tset,ON_IDs));
pnOFF_IDs = find(ismember(Tset,OFF_IDs)); 
ypd = squeeze(FiPD(1,1,:));
subIntDist = cell2mat(arrayfun(@(n) ifft(fft(ypd).^n),0:Smax,'Un',0));%P(Intensity I in channel C given n probes with A Fluorphore) 
Index2 = 1:length(xI);
for c=1:length(Cvec)
    for m = 1:length(Mvec)  
        for d = 1:length(Dvec) 
            for t = 1:length(Tvec)              
                Pd = ModelMetrics.ProbsTargets(c,m,d,t).y;
                pAnyBindSiteCM0 = squeeze(pAnyBindSite(1,:,:,Mvec(m),t,d,Cvec(c)));   
                tHit = find(sum(pAnyBindSiteCM0,2)>0);
                if (~isempty(tHit))
                pON_IDs = find(ismember(tHit,ON_IDs));
                pOFF_IDs = find(ismember(tHit,OFF_IDs));
                pdIntensity = cell2mat(arrayfun(@(a)sum(repmat(Pd(a,:),[size(subIntDist,1) 1]).*subIntDist,2),1:size(Pd,1),'Un',0))';%probability specific intensity values in n channels as vector. 
                subON = repmat({':'},1,ndims(pdIntensity));%can make variable number of : for nd
                subOFF = repmat({':'},1,ndims(pdIntensity));   
                subON{1} = pnON_IDs;  
                subOFF{1} = pnOFF_IDs; 
                pdIntensityOn = squeeze(sum(pdIntensity(subON{:}),1))/sum(pdIntensity(subON{:}),'all');
                pdIntensityOff = squeeze(sum(pdIntensity(subOFF{:}),1))/sum(pdIntensity(subOFF{:}),'all');
                pdIntensity_ExprOn = squeeze(sum(pdIntensity_Expr(subON{:}),1));
                pdIntensity_ExprOff = squeeze(sum(pdIntensity_Expr(subOFF{:}),1));
                It = pdIntensity;               
                It = It./repmat(sum(It,2),[1 size(It,2)]);
                Qi = cumsum(It,2);%cumultative probability of intensity per target.
                I_Detected = repmat(nExpressionMatrix(tHit,Cvec(c)),[1 size(It,2)]).*(1-Qi);%expression 
                I_Missed = repmat(nExpressionMatrix(tHit,Cvec(c)),[1 size(It,2)]).*Qi;%expression 
                Num_Detected = sum(I_Detected,1);
                Num_Missed = sum(I_Missed,1);
                FalseNegative_SpotCounts = sum(I_Missed(pON_IDs,:),1);
                TrueNegative_SpotCounts = sum(I_Missed(pOFF_IDs,:),1);
                FalsePositive_SpotCounts = sum(I_Detected(pOFF_IDs,:),1);
                TruePositive_SpotCounts = sum(I_Detected(pON_IDs,:),1);
                True_SpotIntensities = TruePositive_SpotCounts + FalseNegative_SpotCounts;
                Removed_SpotIntensities = TrueNegative_SpotCounts + FalsePositive_SpotCounts;
                Measured_SpotIntensities = TruePositive_SpotCounts + FalsePositive_SpotCounts;
                Negative_SpotIntensities = TrueNegative_SpotCounts + FalseNegative_SpotCounts;  
                Prevalence_SpotIntensities = Measured_SpotIntensities./(Measured_SpotIntensities+Negative_SpotIntensities);
                Accuracy_SpotIntensities = (TruePositive_SpotCounts+TrueNegative_SpotCounts)./(TruePositive_SpotCounts+TrueNegative_SpotCounts+FalsePositive_SpotCounts+FalseNegative_SpotCounts);
                BalancedAccuracy_SpotIntensities = 50*TruePositive_SpotCounts./True_SpotIntensities+50*TrueNegative_SpotCounts./Removed_SpotIntensities;
                F1_SpotIntensities = 2*TruePositive_SpotCounts./(2*TruePositive_SpotCounts+FalsePositive_SpotCounts+FalseNegative_SpotCounts);
                MCC_SpotIntensities = (TruePositive_SpotCounts.*TrueNegative_SpotCounts-FalsePositive_SpotCounts.*FalseNegative_SpotCounts)./sqrt(Measured_SpotIntensities.*True_SpotIntensities.*Removed_SpotIntensities.*Negative_SpotIntensities);
                SpotIntensities.TPR(Index2) = 100*TruePositive_SpotCounts./True_SpotIntensities;
                SpotIntensities.FNR(Index2) = 100*FalseNegative_SpotCounts./True_SpotIntensities;   
                SpotIntensities.TNR(Index2)= 100*TrueNegative_SpotCounts./Removed_SpotIntensities;
                SpotIntensities.FPR(Index2) = 100*FalsePositive_SpotCounts./Removed_SpotIntensities;   
                SpotIntensities.PPV(Index2) = 100*TruePositive_SpotCounts./Measured_SpotIntensities;
                SpotIntensities.FDR(Index2) = 100*FalsePositive_SpotCounts./Measured_SpotIntensities;  
                SpotIntensities.NPV(Index2) = 100*TrueNegative_SpotCounts./Negative_SpotIntensities;
                SpotIntensities.FOR(Index2) = 100*FalseNegative_SpotCounts./Negative_SpotIntensities;
                SpotIntensities.TP(Index2) = TruePositive_SpotCounts;
                SpotIntensities.FN(Index2) = FalseNegative_SpotCounts;  
                SpotIntensities.TN(Index2) = TrueNegative_SpotCounts;
                SpotIntensities.FP(Index2) = FalsePositive_SpotCounts;   
                SpotIntensities.T(Index2) = True_SpotIntensities;
                SpotIntensities.F(Index2) = Removed_SpotIntensities;
                SpotIntensities.P(Index2) = Measured_SpotIntensities;
                SpotIntensities.N(Index2) = Negative_SpotIntensities;  
                SpotIntensities.PREV(Index2) = Prevalence_SpotIntensities;   
                SpotIntensities.ACC(Index2) = Accuracy_SpotIntensities;   
                SpotIntensities.BA(Index2) = BalancedAccuracy_SpotIntensities;   
                SpotIntensities.F1(Index2) = F1_SpotIntensities;   
                SpotIntensities.MCC(Index2) = MCC_SpotIntensities;  
                ProbeSetMetrics(c,m,d,t).fields = {'P','N','F','T','FP','FN','TP','TN','FOR','NPV','FDR','PPV','FPR','TPR','TNR','FNR','PREV','ACC','BA','F1','MCC'};   
                ProbeSetMetrics(c,m,d,t).SpotIntensities = SpotIntensities;
                IntensityTarget(c,m,d,t).pdIntensity = pdIntensity;  
                IntensityTarget(c,m,d,t).pdIntensity_Expr = pdIntensity_Expr;  
                IntensityTarget(c,m,d,t).pdIntensityOn = pdIntensityOn;  
                IntensityTarget(c,m,d,t).pdIntensityOff = pdIntensityOff;  
                IntensityTarget(c,m,d,t).pdIntensity_ExprOn = pdIntensity_ExprOn;  
                IntensityTarget(c,m,d,t).pdIntensity_ExprOff = pdIntensity_ExprOff;  
                end
            end
        end
    end
end



    
for c=1:length(Cvec)
    for m = 1:length(Mvec)  
        for t = 1:length(Tvec)
            for d = 1:length(Dvec)            
                Pnd = MultiColorTarget(c,m,d,t).PnD;   
                pIntensity_Pnd = @(C) ifft(prod(CATnWrapper(arrayfun(@(k)permute(repmat(CATnWrapper(arrayfun(@(n) fft(squeeze(FiPD(k,C,:))).^n,0:SxMC,'Un',0),2),[1 1 (SxMC+1)*ones(1,size(FiPD,1)-1)]) ,[1 circshift(2:size(FiPD,1)+1,k-1)]),1:size(FiPD,1),'Un',0),size(FiPD,1)+2),size(FiPD,1)+2));% 
                Rep_Pnd = permute(repmat(Pnd,[ones(1,ndims(Pnd)) size(FiPD,3)]),[1 size(FiPD,1)+2 2:size(FiPD,1)+1]);
                pIntensity_In_C = @(C) sum(Rep_Pnd.*permute(repmat(pIntensity_Pnd(C),[ones(1,ndims(Pnd)) size(Pnd,1)]),[size(FiPD,1)+2 1:size(FiPD,1)+1]),[size(FiPD,1):size(FiPD,1)+2]);
                pndIntensity = prod(CATnWrapper(arrayfun(@(k)permute(repmat(pIntensity_In_C(k),[1 1 size(FiPD,3)*ones(1,size(FiPD,1)-1)]),[1 circshift(2:size(FiPD,1)+1,k-1)]),1:size(FiPD,2),'Un',0),size(FiPD,1)+2),size(FiPD,1)+2);        
                pndIntensity_Expr = repmat(ExpressionMatrix(Tset,c),[1 size(FiPD,3)*ones(1,size(FiPD,1))]).*pndIntensity;
                %probability specific intensity values in n channels as vector. 
                subON = repmat({':'},1,ndims(pndIntensity));%can make variable number of : for nd
                subOFF = repmat({':'},1,ndims(pndIntensity));
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
                pndIntensityOnJoint.Combos = combos;
                pndIntensityOffJoint.PDF =  CATnWrapper(arrayfun(@(k) squeeze(sum(pndIntensityOff,setdiff(1:size(FiPD,2),combos(k,:)))),1:size(combos,1),'Un',0),3);
                pndIntensity_ExprOnJoint.Combos = combos;
                pndIntensity_ExprOnJoint.PDF =  CATnWrapper(arrayfun(@(k) squeeze(sum(pndIntensity_ExprOn,setdiff(1:size(FiPD,2),combos(k,:)))),1:size(combos,1),'Un',0),3);
                pndIntensity_ExprOffJoint.Combos = combos;
                pndIntensity_ExprOffJoint.PDF =  CATnWrapper(arrayfun(@(k) squeeze(sum(pndIntensity_ExprOff,setdiff(1:size(FiPD,2),combos(k,:)))),1:size(combos,1),'Un',0),3); 
                %threshold each channel independently.   
               
                
                
                combos = nchoosek(1:size(FiPD,2),1);%dimension is the combo it is 
                pndIntensityOnJoint.Combos = combos;
                pndIntensityOnJoint.PDF =  CATnWrapper(arrayfun(@(k) squeeze(sum(pndIntensityOn,setdiff(1:size(FiPD,2),combos(k,:)))),1:size(combos,1),'Un',0),3);
                pndIntensityOnJoint.Combos = combos;
                pndIntensityOffJoint.PDF =  CATnWrapper(arrayfun(@(k) squeeze(sum(pndIntensityOff,setdiff(1:size(FiPD,2),combos(k,:)))),1:size(combos,1),'Un',0),3);
                pndIntensity_ExprOnJoint.Combos = combos;
                pndIntensity_ExprOnJoint.PDF =  CATnWrapper(arrayfun(@(k) squeeze(sum(pndIntensity_ExprOn,setdiff(1:size(FiPD,2),combos(k,:)))),1:size(combos,1),'Un',0),3);
                pndIntensity_ExprOffJoint.Combos = combos;
                pndIntensity_ExprOffJoint.PDF =  CATnWrapper(arrayfun(@(k) squeeze(sum(pndIntensity_ExprOff,setdiff(1:size(FiPD,2),combos(k,:)))),1:size(combos,1),'Un',0),3); 
             
                
                MultiIntensityTarget(c,m,d,t).pndIntensity = pndIntensity;  
                MultiIntensityTarget(c,m,d,t).pndIntensity_Expr = pndIntensity_Expr;  
                MultiIntensityTarget(c,m,d,t).pndIntensityOn = pndIntensityOn;  
                MultiIntensityTarget(c,m,d,t).pndIntensityOff = pndIntensityOff;  
                MultiIntensityTarget(c,m,d,t).pndIntensity_ExprOn = pndIntensity_ExprOn;  
                MultiIntensityTarget(c,m,d,t).pndIntensity_ExprOff = pndIntensity_ExprOff;  
                MultiIntensityTarget(c,m,d,t).pndIntensityOnJoint = pndIntensityOnJoint;  
                MultiIntensityTarget(c,m,d,t).pndIntensityOffJoint = pndIntensityOffJoint;  
                MultiIntensityTarget(c,m,d,t).pndIntensity_ExprOnJoint = pndIntensity_ExprOnJoint;  
                MultiIntensityTarget(c,m,d,t).pndIntensity_ExprOffJoint = pndIntensity_ExprOffJoint;         
            end
        end
    end
end

end
IntensityMetrics.MultiIntensityTarget = MultiIntensityTarget;

end