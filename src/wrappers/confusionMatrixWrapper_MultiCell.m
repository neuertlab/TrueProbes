<<<<<<< HEAD
<<<<<<< HEAD
function Spots = confusionMatrixWrapper_MultiCell(Non_MultiCellVector,Noff_MultiCellVector)
Ncells = size(Non_MultiCellVector,1);
TruePositive_SpotCounts = CATnWrapper({sum(full(Non_MultiCellVector(:,1:end)),2),cumsum(full(Non_MultiCellVector(:,1:end)),2,'reverse'),zeros([size(Non_MultiCellVector,1) 1])},2);
FalsePositive_SpotCounts = CATnWrapper({sum(full(Noff_MultiCellVector(:,1:end)),2),cumsum(full(Noff_MultiCellVector(:,1:end)),2,'reverse'),zeros([size(Noff_MultiCellVector,1) 1])},2);
FalseNegative_SpotCounts = CATnWrapper({zeros([size(Non_MultiCellVector,1) 1]),cumsum(full(Non_MultiCellVector(:,1:end)),2),sum(full(Non_MultiCellVector(:,1:end)),2)},2);
TrueNegative_SpotCounts = CATnWrapper({zeros([size(Noff_MultiCellVector,1) 1]),cumsum(full(Noff_MultiCellVector(:,1:end)),2),sum(full(Noff_MultiCellVector(:,1:end)),2)},2);
True_Classified_Spots = TruePositive_SpotCounts + TrueNegative_SpotCounts;%True Classified
False_Classified_Spots = FalsePositive_SpotCounts+FalseNegative_SpotCounts;%Negative Classified
Actual_Positive_Spots = TruePositive_SpotCounts + FalseNegative_SpotCounts;%Actual Positive
Actual_Negative_Spots = TrueNegative_SpotCounts + FalsePositive_SpotCounts;% Actual Negative
Predicted_Positive_Spots = TruePositive_SpotCounts + FalsePositive_SpotCounts;%Predicted Positive
Predicted_Negative_Spots = TrueNegative_SpotCounts + FalseNegative_SpotCounts;% Predicted Negative
Prevalence_Spots = Actual_Positive_Spots./(Predicted_Positive_Spots+Predicted_Negative_Spots);
Accuracy_Spots = (TruePositive_SpotCounts+TrueNegative_SpotCounts)./(TruePositive_SpotCounts+TrueNegative_SpotCounts+FalsePositive_SpotCounts+FalseNegative_SpotCounts);
BalancedAccuracy_Spots = 50*TruePositive_SpotCounts./Actual_Positive_Spots+50*TrueNegative_SpotCounts./Actual_Negative_Spots;
P4_Spots = 4*TruePositive_SpotCounts.*TrueNegative_SpotCounts./(4*TruePositive_SpotCounts.*TrueNegative_SpotCounts+True_Classified_Spots.*False_Classified_Spots);
F1_Spots = 2*TruePositive_SpotCounts./(2*TruePositive_SpotCounts+FalsePositive_SpotCounts+FalseNegative_SpotCounts);
MCC_Spots = (TruePositive_SpotCounts.*TrueNegative_SpotCounts-FalsePositive_SpotCounts.*FalseNegative_SpotCounts)./sqrt(Predicted_Positive_Spots.*Actual_Positive_Spots.*Actual_Negative_Spots.*Predicted_Negative_Spots);
CKC_Spots = 2*(TruePositive_SpotCounts.*TrueNegative_SpotCounts-FalseNegative_SpotCounts.*FalsePositive_SpotCounts)./(Predicted_Positive_Spots.*Actual_Negative_Spots+Actual_Positive_Spots.*Predicted_Negative_Spots);
Spots.T = True_Classified_Spots;
Spots.F = False_Classified_Spots;
Spots.AP = Actual_Positive_Spots;
Spots.AN = Actual_Negative_Spots;
Spots.PP = Predicted_Positive_Spots;
Spots.PN = Predicted_Negative_Spots;
Spots.TP = TruePositive_SpotCounts;
Spots.FN = FalseNegative_SpotCounts;
Spots.TN = TrueNegative_SpotCounts;
Spots.FP = FalsePositive_SpotCounts;
Spots.TPR = 100*TruePositive_SpotCounts./Actual_Positive_Spots;
Spots.FNR = 100*FalseNegative_SpotCounts./Actual_Positive_Spots;
Spots.TNR= 100*TrueNegative_SpotCounts./Actual_Negative_Spots;
Spots.FPR = 100*FalsePositive_SpotCounts./Actual_Negative_Spots;
Spots.PPV = 100*TruePositive_SpotCounts./Predicted_Positive_Spots;
Spots.FDR = 100*FalsePositive_SpotCounts./Predicted_Positive_Spots;
Spots.NPV = 100*TrueNegative_SpotCounts./Predicted_Negative_Spots;
Spots.FOR = 100*FalseNegative_SpotCounts./Predicted_Negative_Spots;
Spots.PLR = Spots.TPR./Spots.FPR;
Spots.NLR = Spots.FNR./Spots.TNR;
Spots.DOR = Spots.PLR./Spots.NLR;
Spots.Markedness = (Spots.PPV+Spots.NPV)/100-1;
Spots.Informedness = (Spots.TPR+Spots.TNR)/100-1;
Spots.PREV_THRESH = (sqrt(Spots.TPR.*Spots.FPR)-Spots.FPR)./(Spots.TPR-Spots.FPR);
Spots.Prevalence = Prevalence_Spots;
Spots.Accuracy = Accuracy_Spots;
Spots.BalancedAccuracy = BalancedAccuracy_Spots;
Spots.FowlkesMallowsIndex = sqrt(Spots.PPV/100.*Spots.TPR/100);
Spots.JaccardIndex = Spots.TP./(Spots.TP+Spots.FN+Spots.FP);
Spots.F1 = F1_Spots;
Spots.Fbeta = @(beta) (1+beta^2)*(TruePositive_SpotCounts)./((1+beta^2)*TruePositive_SpotCounts+beta^2*FalseNegative_SpotCounts+FalsePositive_SpotCounts);
Spots.P4 = P4_Spots;
Spots.MCC = MCC_Spots;
Spots.CKC = CKC_Spots;
Spots.TOC_AUC = cell2mat(arrayfun(@(nth_cell) trapz(fliplr(Spots.TP(nth_cell,:)+Spots.FP(nth_cell,:)),fliplr(Spots.TP(nth_cell,:)))/max(Spots.TP(nth_cell,:).*Spots.PP(nth_cell,:)),1:Ncells,'Un',0));
Spots.ROC_AUC = cell2mat(arrayfun(@(nth_cell) trapz(fliplr(Spots.FPR(nth_cell,~isnan(Spots.FPR(nth_cell,:)).*~isnan(Spots.TPR(nth_cell,:))==1)),...
                                                                                fliplr(Spots.TPR(nth_cell,~isnan(Spots.FPR(nth_cell,:)).*~isnan(Spots.TPR(nth_cell,:))==1))),1:Ncells,'Un',0))/100^2;
Spots.PR_AUC = cell2mat(arrayfun(@(nth_cell) trapz(fliplr(Spots.TPR(nth_cell,~isnan(Spots.TPR(nth_cell,:)).*~isnan(Spots.PPV(nth_cell,:))==1)),...
                                                                                fliplr(Spots.PPV(nth_cell,~isnan(Spots.TPR(nth_cell,:)).*~isnan(Spots.PPV(nth_cell,:))==1))),1:Ncells,'Un',0))/100^2;
%TOC_Curve shows Hits (TP), Misses (FN), false alarm (FP), correct rejection (TN)
%TOC_Curve (TP+FP,TP), also plot parallelogram of (0,0) (N,0) (N+P,P) and (P,P)
%ROC_Curve (FPR, TPR)
%PR Curve (Recall, Precision)/ (TPR, PPV)
%Average Precision pr curve precision(recall). integral recall 0 to 1.
Spots.ConfusionMatrixFunction = @(d) permute(CATnWrapper(arrayfun(@(z) [TruePositive_SpotCounts(z,d(z)); FalsePositive_SpotCounts(z,d(z)); FalseNegative_SpotCounts(z,d(z)); TrueNegative_SpotCounts(z,d(z))],1:Ncells,'Un',0),3),[3 1 2]);
=======
function Spots = confusionMatrixWrapper_MultiCell(Non_MultiCellVector,Noff_MultiCellVector)
Ncells = size(Non_MultiCellVector,1);
TruePositive_SpotCounts = CATnWrapper({sum(full(Non_MultiCellVector(:,1:end)),2),cumsum(full(Non_MultiCellVector(:,1:end)),2,'reverse'),zeros([size(Non_MultiCellVector,1) 1])},2);
FalsePositive_SpotCounts = CATnWrapper({sum(full(Noff_MultiCellVector(:,1:end)),2),cumsum(full(Noff_MultiCellVector(:,1:end)),2,'reverse'),zeros([size(Noff_MultiCellVector,1) 1])},2);
FalseNegative_SpotCounts = CATnWrapper({zeros([size(Non_MultiCellVector,1) 1]),cumsum(full(Non_MultiCellVector(:,1:end)),2),sum(full(Non_MultiCellVector(:,1:end)),2)},2);
TrueNegative_SpotCounts = CATnWrapper({zeros([size(Noff_MultiCellVector,1) 1]),cumsum(full(Noff_MultiCellVector(:,1:end)),2),sum(full(Noff_MultiCellVector(:,1:end)),2)},2);
True_Classified_Spots = TruePositive_SpotCounts + TrueNegative_SpotCounts;%True Classified
False_Classified_Spots = FalsePositive_SpotCounts+FalseNegative_SpotCounts;%Negative Classified
Actual_Positive_Spots = TruePositive_SpotCounts + FalseNegative_SpotCounts;%Actual Positive
Actual_Negative_Spots = TrueNegative_SpotCounts + FalsePositive_SpotCounts;% Actual Negative
Predicted_Positive_Spots = TruePositive_SpotCounts + FalsePositive_SpotCounts;%Predicted Positive
Predicted_Negative_Spots = TrueNegative_SpotCounts + FalseNegative_SpotCounts;% Predicted Negative
Prevalence_Spots = Actual_Positive_Spots./(Predicted_Positive_Spots+Predicted_Negative_Spots);
Accuracy_Spots = (TruePositive_SpotCounts+TrueNegative_SpotCounts)./(TruePositive_SpotCounts+TrueNegative_SpotCounts+FalsePositive_SpotCounts+FalseNegative_SpotCounts);
BalancedAccuracy_Spots = 50*TruePositive_SpotCounts./Actual_Positive_Spots+50*TrueNegative_SpotCounts./Actual_Negative_Spots;
P4_Spots = 4*TruePositive_SpotCounts.*TrueNegative_SpotCounts./(4*TruePositive_SpotCounts.*TrueNegative_SpotCounts+True_Classified_Spots.*False_Classified_Spots);
F1_Spots = 2*TruePositive_SpotCounts./(2*TruePositive_SpotCounts+FalsePositive_SpotCounts+FalseNegative_SpotCounts);
MCC_Spots = (TruePositive_SpotCounts.*TrueNegative_SpotCounts-FalsePositive_SpotCounts.*FalseNegative_SpotCounts)./sqrt(Predicted_Positive_Spots.*Actual_Positive_Spots.*Actual_Negative_Spots.*Predicted_Negative_Spots);
CKC_Spots = 2*(TruePositive_SpotCounts.*TrueNegative_SpotCounts-FalseNegative_SpotCounts.*FalsePositive_SpotCounts)./(Predicted_Positive_Spots.*Actual_Negative_Spots+Actual_Positive_Spots.*Predicted_Negative_Spots);
Spots.T = True_Classified_Spots;
Spots.F = False_Classified_Spots;
Spots.AP = Actual_Positive_Spots;
Spots.AN = Actual_Negative_Spots;
Spots.PP = Predicted_Positive_Spots;
Spots.PN = Predicted_Negative_Spots;
Spots.TP = TruePositive_SpotCounts;
Spots.FN = FalseNegative_SpotCounts;
Spots.TN = TrueNegative_SpotCounts;
Spots.FP = FalsePositive_SpotCounts;
Spots.TPR = 100*TruePositive_SpotCounts./Actual_Positive_Spots;
Spots.FNR = 100*FalseNegative_SpotCounts./Actual_Positive_Spots;
Spots.TNR= 100*TrueNegative_SpotCounts./Actual_Negative_Spots;
Spots.FPR = 100*FalsePositive_SpotCounts./Actual_Negative_Spots;
Spots.PPV = 100*TruePositive_SpotCounts./Predicted_Positive_Spots;
Spots.FDR = 100*FalsePositive_SpotCounts./Predicted_Positive_Spots;
Spots.NPV = 100*TrueNegative_SpotCounts./Predicted_Negative_Spots;
Spots.FOR = 100*FalseNegative_SpotCounts./Predicted_Negative_Spots;
Spots.PLR = Spots.TPR./Spots.FPR;
Spots.NLR = Spots.FNR./Spots.TNR;
Spots.DOR = Spots.PLR./Spots.NLR;
Spots.Markedness = (Spots.PPV+Spots.NPV)/100-1;
Spots.Informedness = (Spots.TPR+Spots.TNR)/100-1;
Spots.PREV_THRESH = (sqrt(Spots.TPR.*Spots.FPR)-Spots.FPR)./(Spots.TPR-Spots.FPR);
Spots.Prevalence = Prevalence_Spots;
Spots.Accuracy = Accuracy_Spots;
Spots.BalancedAccuracy = BalancedAccuracy_Spots;
Spots.FowlkesMallowsIndex = sqrt(Spots.PPV/100.*Spots.TPR/100);
Spots.JaccardIndex = Spots.TP./(Spots.TP+Spots.FN+Spots.FP);
Spots.F1 = F1_Spots;
Spots.Fbeta = @(beta) (1+beta^2)*(TruePositive_SpotCounts)./((1+beta^2)*TruePositive_SpotCounts+beta^2*FalseNegative_SpotCounts+FalsePositive_SpotCounts);
Spots.P4 = P4_Spots;
Spots.MCC = MCC_Spots;
Spots.CKC = CKC_Spots;
Spots.TOC_AUC = cell2mat(arrayfun(@(nth_cell) trapz(fliplr(Spots.TP(nth_cell,:)+Spots.FP(nth_cell,:)),fliplr(Spots.TP(nth_cell,:)))/max(Spots.TP(nth_cell,:).*Spots.PP(nth_cell,:)),1:Ncells,'Un',0));
Spots.ROC_AUC = cell2mat(arrayfun(@(nth_cell) trapz(fliplr(Spots.FPR(nth_cell,~isnan(Spots.FPR(nth_cell,:)).*~isnan(Spots.TPR(nth_cell,:))==1)),...
                                                                                fliplr(Spots.TPR(nth_cell,~isnan(Spots.FPR(nth_cell,:)).*~isnan(Spots.TPR(nth_cell,:))==1))),1:Ncells,'Un',0))/100^2;
Spots.PR_AUC = cell2mat(arrayfun(@(nth_cell) trapz(fliplr(Spots.TPR(nth_cell,~isnan(Spots.TPR(nth_cell,:)).*~isnan(Spots.PPV(nth_cell,:))==1)),...
                                                                                fliplr(Spots.PPV(nth_cell,~isnan(Spots.TPR(nth_cell,:)).*~isnan(Spots.PPV(nth_cell,:))==1))),1:Ncells,'Un',0))/100^2;
%TOC_Curve shows Hits (TP), Misses (FN), false alarm (FP), correct rejection (TN)
%TOC_Curve (TP+FP,TP), also plot parallelogram of (0,0) (N,0) (N+P,P) and (P,P)
%ROC_Curve (FPR, TPR)
%PR Curve (Recall, Precision)/ (TPR, PPV)
%Average Precision pr curve precision(recall). integral recall 0 to 1.
Spots.ConfusionMatrixFunction = @(d) permute(CATnWrapper(arrayfun(@(z) [TruePositive_SpotCounts(z,d(z)); FalsePositive_SpotCounts(z,d(z)); FalseNegative_SpotCounts(z,d(z)); TrueNegative_SpotCounts(z,d(z))],1:Ncells,'Un',0),3),[3 1 2]);
>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
=======
function Spots = confusionMatrixWrapper_MultiCell(Non_MultiCellVector,Noff_MultiCellVector)
Ncells = size(Non_MultiCellVector,1);
TruePositive_SpotCounts = CATnWrapper({sum(full(Non_MultiCellVector(:,1:end)),2),cumsum(full(Non_MultiCellVector(:,1:end)),2,'reverse'),zeros([size(Non_MultiCellVector,1) 1])},2);
FalsePositive_SpotCounts = CATnWrapper({sum(full(Noff_MultiCellVector(:,1:end)),2),cumsum(full(Noff_MultiCellVector(:,1:end)),2,'reverse'),zeros([size(Noff_MultiCellVector,1) 1])},2);
FalseNegative_SpotCounts = CATnWrapper({zeros([size(Non_MultiCellVector,1) 1]),cumsum(full(Non_MultiCellVector(:,1:end)),2),sum(full(Non_MultiCellVector(:,1:end)),2)},2);
TrueNegative_SpotCounts = CATnWrapper({zeros([size(Noff_MultiCellVector,1) 1]),cumsum(full(Noff_MultiCellVector(:,1:end)),2),sum(full(Noff_MultiCellVector(:,1:end)),2)},2);
True_Classified_Spots = TruePositive_SpotCounts + TrueNegative_SpotCounts;%True Classified
False_Classified_Spots = FalsePositive_SpotCounts+FalseNegative_SpotCounts;%Negative Classified
Actual_Positive_Spots = TruePositive_SpotCounts + FalseNegative_SpotCounts;%Actual Positive
Actual_Negative_Spots = TrueNegative_SpotCounts + FalsePositive_SpotCounts;% Actual Negative
Predicted_Positive_Spots = TruePositive_SpotCounts + FalsePositive_SpotCounts;%Predicted Positive
Predicted_Negative_Spots = TrueNegative_SpotCounts + FalseNegative_SpotCounts;% Predicted Negative
Prevalence_Spots = Actual_Positive_Spots./(Predicted_Positive_Spots+Predicted_Negative_Spots);
Accuracy_Spots = (TruePositive_SpotCounts+TrueNegative_SpotCounts)./(TruePositive_SpotCounts+TrueNegative_SpotCounts+FalsePositive_SpotCounts+FalseNegative_SpotCounts);
BalancedAccuracy_Spots = 50*TruePositive_SpotCounts./Actual_Positive_Spots+50*TrueNegative_SpotCounts./Actual_Negative_Spots;
P4_Spots = 4*TruePositive_SpotCounts.*TrueNegative_SpotCounts./(4*TruePositive_SpotCounts.*TrueNegative_SpotCounts+True_Classified_Spots.*False_Classified_Spots);
F1_Spots = 2*TruePositive_SpotCounts./(2*TruePositive_SpotCounts+FalsePositive_SpotCounts+FalseNegative_SpotCounts);
MCC_Spots = (TruePositive_SpotCounts.*TrueNegative_SpotCounts-FalsePositive_SpotCounts.*FalseNegative_SpotCounts)./sqrt(Predicted_Positive_Spots.*Actual_Positive_Spots.*Actual_Negative_Spots.*Predicted_Negative_Spots);
CKC_Spots = 2*(TruePositive_SpotCounts.*TrueNegative_SpotCounts-FalseNegative_SpotCounts.*FalsePositive_SpotCounts)./(Predicted_Positive_Spots.*Actual_Negative_Spots+Actual_Positive_Spots.*Predicted_Negative_Spots);
Spots.T = True_Classified_Spots;
Spots.F = False_Classified_Spots;
Spots.AP = Actual_Positive_Spots;
Spots.AN = Actual_Negative_Spots;
Spots.PP = Predicted_Positive_Spots;
Spots.PN = Predicted_Negative_Spots;
Spots.TP = TruePositive_SpotCounts;
Spots.FN = FalseNegative_SpotCounts;
Spots.TN = TrueNegative_SpotCounts;
Spots.FP = FalsePositive_SpotCounts;
Spots.TPR = 100*TruePositive_SpotCounts./Actual_Positive_Spots;
Spots.FNR = 100*FalseNegative_SpotCounts./Actual_Positive_Spots;
Spots.TNR= 100*TrueNegative_SpotCounts./Actual_Negative_Spots;
Spots.FPR = 100*FalsePositive_SpotCounts./Actual_Negative_Spots;
Spots.PPV = 100*TruePositive_SpotCounts./Predicted_Positive_Spots;
Spots.FDR = 100*FalsePositive_SpotCounts./Predicted_Positive_Spots;
Spots.NPV = 100*TrueNegative_SpotCounts./Predicted_Negative_Spots;
Spots.FOR = 100*FalseNegative_SpotCounts./Predicted_Negative_Spots;
Spots.PLR = Spots.TPR./Spots.FPR;
Spots.NLR = Spots.FNR./Spots.TNR;
Spots.DOR = Spots.PLR./Spots.NLR;
Spots.Markedness = (Spots.PPV+Spots.NPV)/100-1;
Spots.Informedness = (Spots.TPR+Spots.TNR)/100-1;
Spots.PREV_THRESH = (sqrt(Spots.TPR.*Spots.FPR)-Spots.FPR)./(Spots.TPR-Spots.FPR);
Spots.Prevalence = Prevalence_Spots;
Spots.Accuracy = Accuracy_Spots;
Spots.BalancedAccuracy = BalancedAccuracy_Spots;
Spots.FowlkesMallowsIndex = sqrt(Spots.PPV/100.*Spots.TPR/100);
Spots.JaccardIndex = Spots.TP./(Spots.TP+Spots.FN+Spots.FP);
Spots.F1 = F1_Spots;
Spots.Fbeta = @(beta) (1+beta^2)*(TruePositive_SpotCounts)./((1+beta^2)*TruePositive_SpotCounts+beta^2*FalseNegative_SpotCounts+FalsePositive_SpotCounts);
Spots.P4 = P4_Spots;
Spots.MCC = MCC_Spots;
Spots.CKC = CKC_Spots;
Spots.TOC_AUC = cell2mat(arrayfun(@(nth_cell) trapz(fliplr(Spots.TP(nth_cell,:)+Spots.FP(nth_cell,:)),fliplr(Spots.TP(nth_cell,:)))/max(Spots.TP(nth_cell,:).*Spots.PP(nth_cell,:)),1:Ncells,'Un',0));
Spots.ROC_AUC = cell2mat(arrayfun(@(nth_cell) trapz(fliplr(Spots.FPR(nth_cell,~isnan(Spots.FPR(nth_cell,:)).*~isnan(Spots.TPR(nth_cell,:))==1)),...
                                                                                fliplr(Spots.TPR(nth_cell,~isnan(Spots.FPR(nth_cell,:)).*~isnan(Spots.TPR(nth_cell,:))==1))),1:Ncells,'Un',0))/100^2;
Spots.PR_AUC = cell2mat(arrayfun(@(nth_cell) trapz(fliplr(Spots.TPR(nth_cell,~isnan(Spots.TPR(nth_cell,:)).*~isnan(Spots.PPV(nth_cell,:))==1)),...
                                                                                fliplr(Spots.PPV(nth_cell,~isnan(Spots.TPR(nth_cell,:)).*~isnan(Spots.PPV(nth_cell,:))==1))),1:Ncells,'Un',0))/100^2;
%TOC_Curve shows Hits (TP), Misses (FN), false alarm (FP), correct rejection (TN)
%TOC_Curve (TP+FP,TP), also plot parallelogram of (0,0) (N,0) (N+P,P) and (P,P)
%ROC_Curve (FPR, TPR)
%PR Curve (Recall, Precision)/ (TPR, PPV)
%Average Precision pr curve precision(recall). integral recall 0 to 1.
Spots.ConfusionMatrixFunction = @(d) permute(CATnWrapper(arrayfun(@(z) [TruePositive_SpotCounts(z,d(z)); FalsePositive_SpotCounts(z,d(z)); FalseNegative_SpotCounts(z,d(z)); TrueNegative_SpotCounts(z,d(z))],1:Ncells,'Un',0),3),[3 1 2]);
>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
end