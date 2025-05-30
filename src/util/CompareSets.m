<<<<<<< HEAD
<<<<<<< HEAD
function CompareSets(Pset1,Pset2)
Pbasis = union(Pset1,Pset2);
Pshared = intersect(Pset1,Pset2);%Probes Shared Between Set 1 and Set 2
Pset1_exclusive = setdiff(Pset1,Pset2);%Probes that are Unique Set 1
Pset2_exclusive = setdiff(Pset2,Pset1);%Probes that are Unique Set 2
Prob_Shared = 100*length(Pshared)/length(Pbasis);
Prob_Exclusive_Set1 = 100*length(Pset1_exclusive)/length(Pbasis);
Prob_Exclusive_Set2 = 100*length(Pset2_exclusive)/length(Pbasis);

Tset_RNA_shared = Js_OFFRNA(Pshared);%Targets from probes that are shared
Tset1_RNA_exclusive = Js_OFFRNA(Pset1_exclusive);%Targets that are in unique Set 1 probes
Tset2_RNA_exclusive = Js_OFFRNA(Pset2_exclusive);%Targets that are in unique Set 2 probes
    Tset_shared_logKOFF_RNA = cell2mat(arrayfun(@(x) Tvec_logKOFF_RNA{x},Pshared,'Un',0));
    Tset_shared_logKOFFdivON_RNA = cell2mat(arrayfun(@(x) Tvec_logKOFFdivON_RNA{x},Pshared,'Un',0));
    Tset_shared_logKONdivOFF_RNA = cell2mat(arrayfun(@(x) Tvec_logKONdivOFF_RNA{x},Pshared,'Un',0));
Tset1_RNA_unique = setdiff(Tset1_RNA_exclusive,Tset2_RNA_exclusive);
    Tset1_Pset1_unique_logKOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_RNA_unique,OFF_RNAIDs)),'Un',0));
    Tset1_Pset1_unique_logKOFFdivON_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_RNA_unique,OFF_RNAIDs)),'Un',0));
    Tset1_Pset1_unique_logKONdivOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_RNA_unique,OFF_RNAIDs)),'Un',0)); 
Tset2_RNA_unique = setdiff(Tset2_RNA_exclusive,Tset1_RNA_exclusive);
    Tset2_Pset2_unique_logKOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset2_RNA_unique,OFF_RNAIDs)),'Un',0));
    Tset2_Pset2_unique_logKOFFdivON_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset2_RNA_unique,OFF_RNAIDs)),'Un',0));
    Tset2_Pset2_unique_logKONdivOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset2_RNA_unique,OFF_RNAIDs)),'Un',0)); 
Tset1_Tset2_RNA_union = intersect(Tset1_RNA_exclusive,Tset2_RNA_exclusive);%targets they share, are they dif or unique sites
    Tset1_Tset2_RNA_union_Pset1_exclusive_logKOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_RNA_union_Pset1_exclusive_logKOFFdivON_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_RNA_union_Pset1_exclusive_logKONdivOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_RNA_union_Pset2_exclusive_logKOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_RNA_union_Pset2_exclusive_logKOFFdivON_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_RNA_union_Pset2_exclusive_logKONdivOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));

    
    
Tset_DNA_shared = Js_OFFDNA(Pshared);%Targets from probes that are shared
Tset1_DNA_exclusive = Js_OFFDNA(Pset1_exclusive);%Targets that are in unique Set 1 probes
Tset2_DNA_exclusive = Js_OFFDNA(Pset2_exclusive);%Targets that are in unique Set 2 probes
    Tset_shared_logKOFF_DNA = cell2mat(arrayfun(@(x) TPvec_logKOFF_DNA{x},find(ismember(Tset_DNA_shared,DNA_IDs)),'Un',0));
    Tset_shared_logKOFFdivON_DNA = cell2mat(arrayfun(@(x) TPvec_logKOFFdivON_DNA{x},find(ismember(Tset_DNA_shared,DNA_IDs)),'Un',0));
    Tset_shared_logKONdivOFF_DNA = cell2mat(arrayfun(@(x) TPvec_logKONdivOFF_DNA{x},find(ismember(Tset_DNA_shared,DNA_IDs)),'Un',0));
    Tset_shared_logKOFFdivCOMP_DNA = cell2mat(arrayfun(@(x) TPvec_logKOFFdivCOMP_DNA{x},find(ismember(Tset_DNA_shared,DNA_IDs)),'Un',0));
    Tset_shared_logKCOMPdivOFF_DNA = cell2mat(arrayfun(@(x) TPvec_logKCOMPdivOFF_DNA{x},find(ismember(Tset_DNA_shared,DNA_IDs)),'Un',0));
Tset1_DNA_unique = setdiff(Tset1_DNA_exclusive,Tset2_DNA_exclusive);
    Tset1_Pset1_unique_logKOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_DNA_unique,DNA_IDs)),'Un',0));
    Tset1_Pset1_unique_logKOFFdivON_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_DNA_unique,DNA_IDs)),'Un',0));
    Tset1_Pset1_unique_logKONdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_DNA_unique,DNA_IDs)),'Un',0)); 
    Tset1_Pset1_unique_logKOFFdivCOMP_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivCOMP_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_DNA_unique,DNA_IDs)),'Un',0));
    Tset1_Pset1_unique_logKCOMPdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKCOMPdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_DNA_unique,DNA_IDs)),'Un',0)); 
Tset2_DNA_unique = setdiff(Tset2_DNA_exclusive,Tset1_DNA_exclusive);
    Tset2_Pset2_unique_logKOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset2_DNA_unique,DNA_IDs)),'Un',0));
    Tset2_Pset2_unique_logKOFFdivON_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset2_DNA_unique,DNA_IDs)),'Un',0));
    Tset2_Pset2_unique_logKONdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset2_DNA_unique,DNA_IDs)),'Un',0));
    Tset2_Pset2_unique_logKOFFdivCOMP_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivCOMP_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset2_DNA_unique,DNA_IDs)),'Un',0));
    Tset2_Pset2_unique_logKCOMPdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKCOMPdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset2_DNA_unique,DNA_IDs)),'Un',0));
Tset1_Tset2_DNA_union = intersect(Tset1_DNA_exclusive,Tset2_DNA_exclusive);%targets they share, are they dif or unique sites
    Tset1_Tset2_DNA_union_Pset1_exclusive_logKOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset1_exclusive_logKOFFdivON_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset1_exclusive_logKONdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset1_exclusive_logKOFFdivCOMP_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivCOMP_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset1_exclusive_logKCOMPdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKCOMPdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset2_exclusive_logKOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_DNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset2_exclusive_logKOFFdivON_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset2_exclusive_logKONdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset2_exclusive_logKOFFdivCOMP_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivCOMP_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset2_exclusive_logKCOMPdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKCOMPdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));

=======
function CompareSets(Pset1,Pset2)
Pbasis = union(Pset1,Pset2);
Pshared = intersect(Pset1,Pset2);%Probes Shared Between Set 1 and Set 2
Pset1_exclusive = setdiff(Pset1,Pset2);%Probes that are Unique Set 1
Pset2_exclusive = setdiff(Pset2,Pset1);%Probes that are Unique Set 2
Prob_Shared = 100*length(Pshared)/length(Pbasis);
Prob_Exclusive_Set1 = 100*length(Pset1_exclusive)/length(Pbasis);
Prob_Exclusive_Set2 = 100*length(Pset2_exclusive)/length(Pbasis);

Tset_RNA_shared = Js_OFFRNA(Pshared);%Targets from probes that are shared
Tset1_RNA_exclusive = Js_OFFRNA(Pset1_exclusive);%Targets that are in unique Set 1 probes
Tset2_RNA_exclusive = Js_OFFRNA(Pset2_exclusive);%Targets that are in unique Set 2 probes
    Tset_shared_logKOFF_RNA = cell2mat(arrayfun(@(x) Tvec_logKOFF_RNA{x},Pshared,'Un',0));
    Tset_shared_logKOFFdivON_RNA = cell2mat(arrayfun(@(x) Tvec_logKOFFdivON_RNA{x},Pshared,'Un',0));
    Tset_shared_logKONdivOFF_RNA = cell2mat(arrayfun(@(x) Tvec_logKONdivOFF_RNA{x},Pshared,'Un',0));
Tset1_RNA_unique = setdiff(Tset1_RNA_exclusive,Tset2_RNA_exclusive);
    Tset1_Pset1_unique_logKOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_RNA_unique,OFF_RNAIDs)),'Un',0));
    Tset1_Pset1_unique_logKOFFdivON_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_RNA_unique,OFF_RNAIDs)),'Un',0));
    Tset1_Pset1_unique_logKONdivOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_RNA_unique,OFF_RNAIDs)),'Un',0)); 
Tset2_RNA_unique = setdiff(Tset2_RNA_exclusive,Tset1_RNA_exclusive);
    Tset2_Pset2_unique_logKOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset2_RNA_unique,OFF_RNAIDs)),'Un',0));
    Tset2_Pset2_unique_logKOFFdivON_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset2_RNA_unique,OFF_RNAIDs)),'Un',0));
    Tset2_Pset2_unique_logKONdivOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset2_RNA_unique,OFF_RNAIDs)),'Un',0)); 
Tset1_Tset2_RNA_union = intersect(Tset1_RNA_exclusive,Tset2_RNA_exclusive);%targets they share, are they dif or unique sites
    Tset1_Tset2_RNA_union_Pset1_exclusive_logKOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_RNA_union_Pset1_exclusive_logKOFFdivON_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_RNA_union_Pset1_exclusive_logKONdivOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_RNA_union_Pset2_exclusive_logKOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_RNA_union_Pset2_exclusive_logKOFFdivON_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_RNA_union_Pset2_exclusive_logKONdivOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));

    
    
Tset_DNA_shared = Js_OFFDNA(Pshared);%Targets from probes that are shared
Tset1_DNA_exclusive = Js_OFFDNA(Pset1_exclusive);%Targets that are in unique Set 1 probes
Tset2_DNA_exclusive = Js_OFFDNA(Pset2_exclusive);%Targets that are in unique Set 2 probes
    Tset_shared_logKOFF_DNA = cell2mat(arrayfun(@(x) TPvec_logKOFF_DNA{x},find(ismember(Tset_DNA_shared,DNA_IDs)),'Un',0));
    Tset_shared_logKOFFdivON_DNA = cell2mat(arrayfun(@(x) TPvec_logKOFFdivON_DNA{x},find(ismember(Tset_DNA_shared,DNA_IDs)),'Un',0));
    Tset_shared_logKONdivOFF_DNA = cell2mat(arrayfun(@(x) TPvec_logKONdivOFF_DNA{x},find(ismember(Tset_DNA_shared,DNA_IDs)),'Un',0));
    Tset_shared_logKOFFdivCOMP_DNA = cell2mat(arrayfun(@(x) TPvec_logKOFFdivCOMP_DNA{x},find(ismember(Tset_DNA_shared,DNA_IDs)),'Un',0));
    Tset_shared_logKCOMPdivOFF_DNA = cell2mat(arrayfun(@(x) TPvec_logKCOMPdivOFF_DNA{x},find(ismember(Tset_DNA_shared,DNA_IDs)),'Un',0));
Tset1_DNA_unique = setdiff(Tset1_DNA_exclusive,Tset2_DNA_exclusive);
    Tset1_Pset1_unique_logKOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_DNA_unique,DNA_IDs)),'Un',0));
    Tset1_Pset1_unique_logKOFFdivON_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_DNA_unique,DNA_IDs)),'Un',0));
    Tset1_Pset1_unique_logKONdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_DNA_unique,DNA_IDs)),'Un',0)); 
    Tset1_Pset1_unique_logKOFFdivCOMP_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivCOMP_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_DNA_unique,DNA_IDs)),'Un',0));
    Tset1_Pset1_unique_logKCOMPdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKCOMPdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_DNA_unique,DNA_IDs)),'Un',0)); 
Tset2_DNA_unique = setdiff(Tset2_DNA_exclusive,Tset1_DNA_exclusive);
    Tset2_Pset2_unique_logKOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset2_DNA_unique,DNA_IDs)),'Un',0));
    Tset2_Pset2_unique_logKOFFdivON_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset2_DNA_unique,DNA_IDs)),'Un',0));
    Tset2_Pset2_unique_logKONdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset2_DNA_unique,DNA_IDs)),'Un',0));
    Tset2_Pset2_unique_logKOFFdivCOMP_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivCOMP_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset2_DNA_unique,DNA_IDs)),'Un',0));
    Tset2_Pset2_unique_logKCOMPdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKCOMPdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset2_DNA_unique,DNA_IDs)),'Un',0));
Tset1_Tset2_DNA_union = intersect(Tset1_DNA_exclusive,Tset2_DNA_exclusive);%targets they share, are they dif or unique sites
    Tset1_Tset2_DNA_union_Pset1_exclusive_logKOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset1_exclusive_logKOFFdivON_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset1_exclusive_logKONdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset1_exclusive_logKOFFdivCOMP_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivCOMP_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset1_exclusive_logKCOMPdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKCOMPdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset2_exclusive_logKOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_DNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset2_exclusive_logKOFFdivON_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset2_exclusive_logKONdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset2_exclusive_logKOFFdivCOMP_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivCOMP_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset2_exclusive_logKCOMPdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKCOMPdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));

>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
=======
function CompareSets(Pset1,Pset2)
Pbasis = union(Pset1,Pset2);
Pshared = intersect(Pset1,Pset2);%Probes Shared Between Set 1 and Set 2
Pset1_exclusive = setdiff(Pset1,Pset2);%Probes that are Unique Set 1
Pset2_exclusive = setdiff(Pset2,Pset1);%Probes that are Unique Set 2
Prob_Shared = 100*length(Pshared)/length(Pbasis);
Prob_Exclusive_Set1 = 100*length(Pset1_exclusive)/length(Pbasis);
Prob_Exclusive_Set2 = 100*length(Pset2_exclusive)/length(Pbasis);

Tset_RNA_shared = Js_OFFRNA(Pshared);%Targets from probes that are shared
Tset1_RNA_exclusive = Js_OFFRNA(Pset1_exclusive);%Targets that are in unique Set 1 probes
Tset2_RNA_exclusive = Js_OFFRNA(Pset2_exclusive);%Targets that are in unique Set 2 probes
    Tset_shared_logKOFF_RNA = cell2mat(arrayfun(@(x) Tvec_logKOFF_RNA{x},Pshared,'Un',0));
    Tset_shared_logKOFFdivON_RNA = cell2mat(arrayfun(@(x) Tvec_logKOFFdivON_RNA{x},Pshared,'Un',0));
    Tset_shared_logKONdivOFF_RNA = cell2mat(arrayfun(@(x) Tvec_logKONdivOFF_RNA{x},Pshared,'Un',0));
Tset1_RNA_unique = setdiff(Tset1_RNA_exclusive,Tset2_RNA_exclusive);
    Tset1_Pset1_unique_logKOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_RNA_unique,OFF_RNAIDs)),'Un',0));
    Tset1_Pset1_unique_logKOFFdivON_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_RNA_unique,OFF_RNAIDs)),'Un',0));
    Tset1_Pset1_unique_logKONdivOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_RNA_unique,OFF_RNAIDs)),'Un',0)); 
Tset2_RNA_unique = setdiff(Tset2_RNA_exclusive,Tset1_RNA_exclusive);
    Tset2_Pset2_unique_logKOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset2_RNA_unique,OFF_RNAIDs)),'Un',0));
    Tset2_Pset2_unique_logKOFFdivON_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset2_RNA_unique,OFF_RNAIDs)),'Un',0));
    Tset2_Pset2_unique_logKONdivOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset2_RNA_unique,OFF_RNAIDs)),'Un',0)); 
Tset1_Tset2_RNA_union = intersect(Tset1_RNA_exclusive,Tset2_RNA_exclusive);%targets they share, are they dif or unique sites
    Tset1_Tset2_RNA_union_Pset1_exclusive_logKOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_RNA_union_Pset1_exclusive_logKOFFdivON_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_RNA_union_Pset1_exclusive_logKONdivOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_RNA{x}(ismember(TPvec_RNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_RNA_union_Pset2_exclusive_logKOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_RNA_union_Pset2_exclusive_logKOFFdivON_RNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_RNA_union_Pset2_exclusive_logKONdivOFF_RNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_RNA{x}(ismember(TPvec_RNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_RNA_union,OFF_RNAIDs)),'Un',0));

    
    
Tset_DNA_shared = Js_OFFDNA(Pshared);%Targets from probes that are shared
Tset1_DNA_exclusive = Js_OFFDNA(Pset1_exclusive);%Targets that are in unique Set 1 probes
Tset2_DNA_exclusive = Js_OFFDNA(Pset2_exclusive);%Targets that are in unique Set 2 probes
    Tset_shared_logKOFF_DNA = cell2mat(arrayfun(@(x) TPvec_logKOFF_DNA{x},find(ismember(Tset_DNA_shared,DNA_IDs)),'Un',0));
    Tset_shared_logKOFFdivON_DNA = cell2mat(arrayfun(@(x) TPvec_logKOFFdivON_DNA{x},find(ismember(Tset_DNA_shared,DNA_IDs)),'Un',0));
    Tset_shared_logKONdivOFF_DNA = cell2mat(arrayfun(@(x) TPvec_logKONdivOFF_DNA{x},find(ismember(Tset_DNA_shared,DNA_IDs)),'Un',0));
    Tset_shared_logKOFFdivCOMP_DNA = cell2mat(arrayfun(@(x) TPvec_logKOFFdivCOMP_DNA{x},find(ismember(Tset_DNA_shared,DNA_IDs)),'Un',0));
    Tset_shared_logKCOMPdivOFF_DNA = cell2mat(arrayfun(@(x) TPvec_logKCOMPdivOFF_DNA{x},find(ismember(Tset_DNA_shared,DNA_IDs)),'Un',0));
Tset1_DNA_unique = setdiff(Tset1_DNA_exclusive,Tset2_DNA_exclusive);
    Tset1_Pset1_unique_logKOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_DNA_unique,DNA_IDs)),'Un',0));
    Tset1_Pset1_unique_logKOFFdivON_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_DNA_unique,DNA_IDs)),'Un',0));
    Tset1_Pset1_unique_logKONdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_DNA_unique,DNA_IDs)),'Un',0)); 
    Tset1_Pset1_unique_logKOFFdivCOMP_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivCOMP_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_DNA_unique,DNA_IDs)),'Un',0));
    Tset1_Pset1_unique_logKCOMPdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKCOMPdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_DNA_unique,DNA_IDs)),'Un',0)); 
Tset2_DNA_unique = setdiff(Tset2_DNA_exclusive,Tset1_DNA_exclusive);
    Tset2_Pset2_unique_logKOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset2_DNA_unique,DNA_IDs)),'Un',0));
    Tset2_Pset2_unique_logKOFFdivON_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset2_DNA_unique,DNA_IDs)),'Un',0));
    Tset2_Pset2_unique_logKONdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset2_DNA_unique,DNA_IDs)),'Un',0));
    Tset2_Pset2_unique_logKOFFdivCOMP_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivCOMP_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset2_DNA_unique,DNA_IDs)),'Un',0));
    Tset2_Pset2_unique_logKCOMPdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKCOMPdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset2_DNA_unique,DNA_IDs)),'Un',0));
Tset1_Tset2_DNA_union = intersect(Tset1_DNA_exclusive,Tset2_DNA_exclusive);%targets they share, are they dif or unique sites
    Tset1_Tset2_DNA_union_Pset1_exclusive_logKOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset1_exclusive_logKOFFdivON_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset1_exclusive_logKONdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset1_exclusive_logKOFFdivCOMP_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivCOMP_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset1_exclusive_logKCOMPdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKCOMPdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset1_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset2_exclusive_logKOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_DNA_union,OFF_RNAIDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset2_exclusive_logKOFFdivON_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivON_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset2_exclusive_logKONdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKONdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset2_exclusive_logKOFFdivCOMP_DNA = cell2mat(arrayfun(@(x)TPvec_logKOFFdivCOMP_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));
    Tset1_Tset2_DNA_union_Pset2_exclusive_logKCOMPdivOFF_DNA = cell2mat(arrayfun(@(x)TPvec_logKCOMPdivOFF_DNA{x}(ismember(TPvec_DNA{x},Pset2_exclusive)),find(ismember(Tset1_Tset2_DNA_union,DNA_IDs)),'Un',0));

>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
end