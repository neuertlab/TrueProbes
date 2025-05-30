function [FITS_PREDICTIONS_SENS] = integrator(nodeID)
    warning('off','all');

    traing_dataset_id = 2;
    kinetics_name = ['TRAINING', num2str(traing_dataset_id,'%02d'),'_t0stt1t2t5'];

    cluster_vs_local = "cluster"; %"local", "cluster";
    use_parallel = "no"; % "yes", "no"

    pop_or_singlecell = "PopAverage"; % "PopAverage", "SingleCells"
    do_sens_analysis = "yes";

    choose_train_data = [6 16 31 38 62 100];
    number_of_MLEs = 50;
    number_of_singlecells = 100;

%     models_to_fit=[1 5 9 13];
    choose_cond = ceil(nodeID/number_of_MLEs); %choose_cond = models_to_fit(choose_con);  % choose models 1, 5, 9, 13%
    cell_id = mod(nodeID,number_of_singlecells); if cell_id==0; cell_id=100; end
    %randomize the seed for random numbers generator
    rng('shuffle');
    seeds = randi(2^32,3e3,1);
    rng(seeds(nodeID));

    % get a refrence for the speed of the node.
    for getnodespead=1:10
        tic
        sum(sum(ones(10000,10000)));
        nodespeads(getnodespead)=toc;
    end
    nodespeadref = mean(nodespeads);

    MLE_id = mod(nodeID,number_of_MLEs);
    if MLE_id==0; MLE_id=number_of_MLEs;end

    global dirName; dirName =  [kinetics_name, '/', 'MODEL',num2str(choose_cond,'%03d'), '/', 'MLE',num2str(MLE_id,'%04d')]; mkdir(dirName);
    global min_OBJ; min_OBJ = inf;
    global best_pars; best_pars = [];
    global model_index; model_index = choose_cond;

    % load models
    MMM = load('../../HJ_202108_HOGSignalingModeling/HJ_20210817_HOGModels/HOGmodels_version_04/Models.mat');
    Model = MMM.Model;

    % load data
    DDD = load('../../HJ_202108_HOGSignalingModeling/HJ_20210830_HOGDATA/HogData/HogDataSet.mat');
    HogDataSet = DDD.HogDataSet;
    DataSet1 = HogDataSet.DataSet1; % 100 monotonic stimulations of 25min
    DataSet3 = HogDataSet.DataSet3; % staircases

    DataSet1_IDs = [1:length(DataSet1)]; % DataSet1 indexes
    train_data_IDs = [DataSet1_IDs(choose_train_data)]; % [in1 in2 etc] % training data IDs
    n_traindata=length(train_data_IDs);

    test1_data_IDs = DataSet1_IDs; % test data IDs
    for ssss=train_data_IDs
        test1_data_IDs=test1_data_IDs(test1_data_IDs~=ssss);
    end

    train_inputs = cell(1,n_traindata);
%     train_inputs = cell(1,n_traindata+1);
%     DataSet3{2}.mean_sHog1nuc(end:end+10) = DataSet3{2}.mean_sHog1nuc(end);
%     train_inputs{1} = DataSet3{2}; % include staircase(0.2M,25min) in training dataset
    for i=1:length(train_data_IDs)
        train_inputs{i} = DataSet1{train_data_IDs(i)};
    end
    n_traindata = length(train_inputs);

    model = Model{model_index};
    train_models = cell(1,n_traindata);
    for trainmodels=1:length(train_inputs)
        train_models{trainmodels}.model = Get_ODE(model,train_inputs{trainmodels}.stimulus);
    end
    model=train_models{1}.model;

    % optimization algorithm params
    iMax = 110; % Number of GA-fmean search optimization
    TPdrop = 0.25; % Number of data points to drop
    NsPOP = [260 260]; % Population size
    Ngen = 20; % # of generations through G.A.
    nElites = 10; % Number of elites
    N_Hv = 50; % Number of best parameter sets to move forward through optimization.

    ParamsLimits = [-1 1];  % bounds on params

    % fit optimization condition
    FitConditionsVec.nodeID = nodeID;
    FitConditionsVec.Npop = NsPOP;
    FitConditionsVec.Ngen = Ngen;
    FitConditionsVec.iMax = iMax;
    FitConditionsVec.TPdrop = TPdrop;
    FitConditionsVec.ParamsLimits = ParamsLimits;
    FitConditionsVec.nElites = nElites;
    FitConditionsVec.randnumRef = rand;
    FitConditionsVec.nodespeadref = nodespeadref;

    %% G.A. Search - logarithmic then Fminnsearch
    constrains_log10.LB=ParamsLimits(1)*ones(1,model.n_params);  constrains_log10.UB=ParamsLimits(2)*ones(1,model.n_params);
    constrains.LB = 10.^constrains_log10.LB; constrains.UB = 10.^constrains_log10.UB;
    train_models{1}.constrains = constrains;
    train_models{1}.cell_id = cell_id; train_models{1}.pop_or_singlecell = pop_or_singlecell;
    trainDataModel.train_models = train_models; trainDataModel.train_inputs = train_inputs;

    Mut_Fun = @(parents,options,nvars,FitnessFcn,state,thisScore,thisPopulation)GA_mutation_function(parents,options,nvars,FitnessFcn,state,thisScore,thisPopulation,constrains_log10);
    gaopt = gaoptimset('Display','none','MutationFcn',Mut_Fun);
    gaopt = gaoptimset(gaopt,'EliteCount',nElites,'Generations',Ngen);
    opt = optimset('display','none','MaxIter',500);

    if use_parallel=="yes"

        switch cluster_vs_local
        case "local"
            myCluster = parcluster('local');
            n_cores = myCluster.NumWorkers; % call # of cores on the local machine
        case "cluster"
            n_cores = str2num(getenv('SLURM_JOB_CPUS_PER_NODE')); % call # of cores on the cluster (SLURM)
        end
        parpool(n_cores); % start parallel pool

        gaopt = gaoptimset('useparallel',1);
        opt = optimset('useparallel',1);
        spmd
            warning('off')
        end

    end

    OBJ_GA_full = @(X)objective_function(10.^X', train_inputs, train_models, 0, 0);
    OBJ_fmin_full = @(X)objective_function(10.^X, train_inputs, train_models, 0, 0);

    MLEs_OBJs = []; Hv = []; ee = [];
    % get best parameter sets from mean pop optimization to use for single cell opt
    switch pop_or_singlecell
        case "SingleCells"
            [~, MLEs_OBJs, ee, Hv] = best_params_PopAve_Opt(200,kinetics_name,choose_cond); % load MLEs
    end
    try
    nL=exp(-ee)/nansum(exp(-ee));
    [~, IDs] = sort(ee,'ascend'); Hv50 = Hv(flip(IDs(1:N_Hv)),:); OBJs50 = MLEs_OBJs(flip(IDs(1:N_Hv)));
    end
    % get the mean and covariance of the best parameters sets
    mu0 = (nL*Hv);  % weighted mean of best parameters sets
    ss0 = nL.*(Hv-mu0)'*(Hv-mu0); % weighted covariance of best parameters sets

    Hv = []; OOBJs=[]; PARAMS = []; bestPARAMS = []; OBJs = NaN(iMax,5); sim_time = NaN(iMax,1);
    zetas = 10.^(-2 + (-1-(-2))*rand(iMax,1)); zetas = sort(zetas,'descend'); % coeeficients to rescale parameter range for G.A. resampling [1e-1 to 1e-2]
    etas = 10.^(-4 + (-2-(-4))*rand(iMax,1)); etas = sort(etas,'descend'); % permute collected parameter sets [1e-2 to 1e-4]
    randss = rand(iMax,1); randss = sort(randss,'descend');

    for i=1:iMax
        tic
        disp(['i = ', num2str(i), ' started.'])
        zeta = zetas(i); eta = etas(i);

        FITS_PREDICTIONS_SENS.Hv = Hv;
        FITS_PREDICTIONS_SENS.best_pars = best_pars;
        FITS_PREDICTIONS_SENS.bestPARAMS = bestPARAMS;
        OBJ = min_OBJ; FITS_PREDICTIONS_SENS.OBJ = OBJ;
        save([dirName,'/FITS_PREDICTIONS_SENS'], 'FITS_PREDICTIONS_SENS');

        if (i<=25 | mod(i,10)==0)
            Npop = NsPOP(1);
            switch pop_or_singlecell
                case "PopAverage"
                    disp(' G.A. initial population uniform random sampling')
                    gaopt.InitialPopulation = ParamsLimits(1) + range(ParamsLimits)*rand(Npop,model.n_params); % [-L, L]
                case "SingleCells"
                    disp(' G.A. initial population sampled as mvnrnd with mean and cov from best param sets of mean populations')
                    mu = mu0; ss = ss0;
                    gaopt.InitialPopulation = mvnrnd(mu, ss, Npop);%[fig] = plot_MVN(gaopt.InitialPopulation,Hv,mu0,ss0,free_params_IDs,2,dirName);
            end
        else
            Npop = NsPOP(2);
            switch pop_or_singlecell
                case "PopAverage"
                    disp(' G.A. initial population permuted from previous steps (rand)')
                    Resample_GA_InitialPopulation = zeta*ParamsLimits(1) + zeta*range(ParamsLimits)*rand(Npop,length(model.n_params)); % 0.1*[-L, L]
                case "SingleCells"
                    disp(' G.A. initial population permuted from previous steps (mvrnd)')
                    mu = (nL*Hv); ss = nL.*(Hv-mu)'*(Hv-mu);
                    Resample_GA_InitialPopulation  = zeta*mvnrnd(mu, ss, Npop);
            end
            gaopt.InitialPopulation = Resample_GA_InitialPopulation + repmat(H,Npop,1);
            randomization = eta*Resample_GA_InitialPopulation(1:size(Hv,1),:);
            gaopt.InitialPopulation(1:size(Hv,1),:) = Hv+randomization;
        end
        % time points
        gaopt = gaoptimset(gaopt,'PopulationSize',Npop);
        wwt=cell(1,n_traindata);
        for w=1:n_traindata
            wt = rand(size(train_inputs{w}.mHog1nuc));
            wt(wt<=TPdrop)=0;
            wt(wt>TPdrop)=1; wt(1:12)=1;
            wwt{w}=wt;
        end
        disp('runing GA - red')
        OBJ_GA_red = @(X)objective_function(10.^X', train_inputs, train_models, 1, wwt);
        H = ga(OBJ_GA_red,model.n_params,gaopt); % Run the G.A.
%         PARAMS = [PARAMS H'];
        OBJs(i, 1) = OBJ_fmin_full(H');
        disp(['final obj = ',num2str(OBJs(i, 1))])

        disp('runing FMIN - red')
        OBJ_fmin_red = @(X)objective_function(10.^X, train_inputs, train_models, 1, wwt);
        H = fminsearch(OBJ_fmin_red,H',opt)';
%         PARAMS = [PARAMS H'];
        if size(Hv,1)>50; Hv(1,:)=[];OOBJs(1)=[];end; Hv = [Hv;H];
        OBJs(i, 2) = OBJ_fmin_full(H');
        OOBJs = [OOBJs OBJs(i, 2)]; nOOBJs = OOBJs/nanmax(OOBJs); nL=exp(-nOOBJs)/nansum(exp(-nOOBJs));
        disp(['final obj = ',num2str(OBJs(i, 2))])

        zeta = randss(i)*zeta; eta = randss(i)*eta;
        switch pop_or_singlecell
            case "PopAverage"
                disp('runing GA - full (initial population sampled rand)')
                Resample_GA_InitialPopulation = zeta*ParamsLimits(1) + zeta*range(ParamsLimits)*rand(Npop,length(model.n_params)); % 0.1*[-L, L]
            case "SingleCells"
                disp('runing GA - full (initial population sampled mvnrnd)')
                mu = (nL*Hv); ss = nL.*(Hv-mu)'*(Hv-mu);
                Resample_GA_InitialPopulation  = zeta*mvnrnd(mu, ss, Npop);
        end
        gaopt.InitialPopulation = Resample_GA_InitialPopulation + repmat(H,Npop,1);
        randomization = eta*Resample_GA_InitialPopulation(1:size(Hv,1),:);
        gaopt.InitialPopulation(1:size(Hv,1),:) = Hv+randomization;

        H = ga(OBJ_GA_full,model.n_params,gaopt); % Run the G.A.
        PARAMS = [PARAMS H'];
        OBJs(i, 3) = OBJ_fmin_full(H');
        disp(['final obj = ',num2str(OBJs(i, 3))])

        disp('runing FMIN - full')
        H = fminsearch(OBJ_fmin_full,H',opt)';
        PARAMS = [PARAMS H'];
        if size(Hv,1)>50; Hv(1,:)=[];OOBJs(1)=[];end; Hv = [Hv;H];
        OBJs(i, 4) = OBJ_fmin_full(H');
        OOBJs = [OOBJs OBJs(i, 4)]; nOOBJs = OOBJs/nanmax(OOBJs); nL=exp(-nOOBJs)/nansum(exp(-nOOBJs));

        disp(['final obj = ',num2str(OBJs(i, 4))]);
        OBJs(i, 5) = min_OBJ;
        disp(num2str(OBJs(i,:)));

        sim_time(i) = toc;
        disp(['time = ', num2str(sim_time(i)/60), ' | i = ', num2str(i), ' finished.'])
        disp(' ')
        bestPARAMS = [bestPARAMS best_pars];
    end
    sim_time = sim_time/60;
    OBJ = min_OBJ;

    FITS_PREDICTIONS_SENS.dirName = dirName;
    FITS_PREDICTIONS_SENS.sim_time = sim_time;
    FITS_PREDICTIONS_SENS.FitConditionsVec = FitConditionsVec;
    FITS_PREDICTIONS_SENS.sim_data_indx = DataSet1_IDs;
    FITS_PREDICTIONS_SENS.train_data_indx = train_data_IDs;
    FITS_PREDICTIONS_SENS.test_data_indx = test1_data_IDs;
    FITS_PREDICTIONS_SENS.model_index = model_index;
    FITS_PREDICTIONS_SENS.trainDataModel = trainDataModel;
    FITS_PREDICTIONS_SENS.best_pars = best_pars;
    FITS_PREDICTIONS_SENS.bestPARAMS = bestPARAMS;
    FITS_PREDICTIONS_SENS.OBJ = OBJ;
    FITS_PREDICTIONS_SENS.OBJs = OBJs;
    FITS_PREDICTIONS_SENS.Hv = Hv;
    FITS_PREDICTIONS_SENS.PARAMS = PARAMS;
    save([dirName,'/FITS_PREDICTIONS_SENS'], 'FITS_PREDICTIONS_SENS');

    % get predictions upon each stimuli input for all mutant strains
    disp(['get predictions upon each stimuli input for all mutant strains  '])
    FITS_PREDICTIONS_SENS = get_predictions(FITS_PREDICTIONS_SENS,Model,model,best_pars,HogDataSet,train_data_IDs,test1_data_IDs,pop_or_singlecell,cell_id);

    % get model sensitivity to model parameters
    if do_sens_analysis=="yes"
        npoints = 101; perc = 0.1; % number of sampled parameter values, percentage within which to vary each parameter
        [OBJs0,MODEL_SENSITIVITY] = model_sensitivity_to_parameters(model,best_pars,npoints,perc,DataSet1,pop_or_singlecell,cell_id);
        FITS_PREDICTIONS_SENS.MODEL_SENSITIVITY = MODEL_SENSITIVITY;
        FITS_PREDICTIONS_SENS.OBJs0 = OBJs0;
        save([dirName,'/FITS_PREDICTIONS_SENS'], 'FITS_PREDICTIONS_SENS');
    end

end
