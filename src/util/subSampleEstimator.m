<<<<<<< HEAD
<<<<<<< HEAD
function [MeanDifference, StdDifference, isNormal,h] = subSampleEstimator(Var,N_pop,Size)
True_Mean = mean(Var);
True_Std = std(Var);
MeanDifference = zeros(1,length(N_pop));
StdDifference = zeros(1,length(N_pop));
isNormal = zeros(1,length(N_pop));
for i = 1:length(N_pop)
    RandomSample_Perms = cell2mat(arrayfun(@(x) randperm(length(Var)).'*length(x),1:Size,'Un',0)).';
    RandomSample_Indexes = RandomSample_Perms(:,1:N_pop(i));
    if (N_pop(i)~=1)
        Sample_Means = mean(Var(RandomSample_Indexes),2);
    else
        Sample_Means = Var(RandomSample_Indexes);  
    end
    MeanDifference(i) = mean(Sample_Means)-True_Mean;
    StdDifference(i) = std(Sample_Means)-True_Std;
    isNormal(i) = adtest(Sample_Means);
    h(i) = vartest(Sample_Means,var(Sample_Means));
    clear RandomSample_Perms RandomSample_Indexes Sample_Means
end

end

%number of variation to run is be high 1000
%repeat on full distribution or overlap cut-off 250 genes.

=======
function [MeanDifference, StdDifference, isNormal,h] = subSampleEstimator(Var,N_pop,Size)
True_Mean = mean(Var);
True_Std = std(Var);
MeanDifference = zeros(1,length(N_pop));
StdDifference = zeros(1,length(N_pop));
isNormal = zeros(1,length(N_pop));
for i = 1:length(N_pop)
    RandomSample_Perms = cell2mat(arrayfun(@(x) randperm(length(Var)).'*length(x),1:Size,'Un',0)).';
    RandomSample_Indexes = RandomSample_Perms(:,1:N_pop(i));
    if (N_pop(i)~=1)
        Sample_Means = mean(Var(RandomSample_Indexes),2);
    else
        Sample_Means = Var(RandomSample_Indexes);  
    end
    MeanDifference(i) = mean(Sample_Means)-True_Mean;
    StdDifference(i) = std(Sample_Means)-True_Std;
    isNormal(i) = adtest(Sample_Means);
    h(i) = vartest(Sample_Means,var(Sample_Means));
    clear RandomSample_Perms RandomSample_Indexes Sample_Means
end

end

%number of variation to run is be high 1000
%repeat on full distribution or overlap cut-off 250 genes.

>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
=======
function [MeanDifference, StdDifference, isNormal,h] = subSampleEstimator(Var,N_pop,Size)
True_Mean = mean(Var);
True_Std = std(Var);
MeanDifference = zeros(1,length(N_pop));
StdDifference = zeros(1,length(N_pop));
isNormal = zeros(1,length(N_pop));
for i = 1:length(N_pop)
    RandomSample_Perms = cell2mat(arrayfun(@(x) randperm(length(Var)).'*length(x),1:Size,'Un',0)).';
    RandomSample_Indexes = RandomSample_Perms(:,1:N_pop(i));
    if (N_pop(i)~=1)
        Sample_Means = mean(Var(RandomSample_Indexes),2);
    else
        Sample_Means = Var(RandomSample_Indexes);  
    end
    MeanDifference(i) = mean(Sample_Means)-True_Mean;
    StdDifference(i) = std(Sample_Means)-True_Std;
    isNormal(i) = adtest(Sample_Means);
    h(i) = vartest(Sample_Means,var(Sample_Means));
    clear RandomSample_Perms RandomSample_Indexes Sample_Means
end

end

%number of variation to run is be high 1000
%repeat on full distribution or overlap cut-off 250 genes.

>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
