function ProbeSets = randomProbeSets(probes,Length_Range,Spacing_Range,Number_of_sets,N_Keep)
% This function can generate random probe sets against a target gene given the probe length range,
% and the range of spacing between probes with a maximum number of probes to design.

%Note added correction for fake time out genes that do not have any recorded 
Lpmin = Length_Range(1);
Lpmax = Length_Range(2);
spacing_req = Spacing_Range(1);
probe_poses = zeros(size(probes,1),1);
for pos1 = 1:size(probes,1)
   probe_poses(pos1) = probes{pos1,3};
end
ProbeSets = []; 
for v = 1:Number_of_sets
    %Pick Subset of probes
    probe_IDs = [];
    
    spacing_matrix = zeros(1,spacing_req*2+Lpmin*2+ max(probe_poses(:,1)));
   try
        x = str2num(getenv('SLURM_ARRAY_TASK_ID'));        
   catch
        try
            x = getenv('SLURM_ARRAY_TASK_ID');    
        catch
            x = feature('getpid');
        end
   end
    s = RandStream('mlfg6331_64','seed',sum(100*(1+labindex/numlabs)*clock*x));
    random_order = randperm(s,length(probes));
    p = 1;
    for i = 1:length(probes)
       probe_num = random_order(i);
       Lp = length(probes{probe_num,2});
       if (spacing_matrix(probes{probe_num,3}+spacing_req+Lp-1) == 0)
           spacing_matrix(probes{probe_num,3}:probes{probe_num,3} + spacing_req*2+Lp+ Lpmin-1) = 1;  
           probe_IDs(p) = probe_num;
           p = p + 1;
       end
    end
    if (N_Keep<=length(probe_IDs))
        random_order2 = randperm(length(probe_IDs));
        ProbeSets{v} = sort(probe_IDs(random_order2(1:N_Keep)),'ascend');
    else
        ProbeSets{v} = sort(probe_IDs,'ascend');
    end
end
end