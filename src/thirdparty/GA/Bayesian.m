function [FITS_PREDICTIONS_SENS] = Bayesian(local_vs_cluser,training_condition,train_data_IDs,model_version,model_index,chain,prior,resample_dist,re_sampling_outliers,pop_or_singlecell,do_sens_analysis,GIBBS,doMetHast)
    %% Bayesian algorithm with lognormal parameter sampling for training models
    warning('off','all');
    restoredefaultpath;
    
    if GIBBS.do=="yes"
        switch local_vs_cluser
        case 'local'
            myCluster = parcluster('local');
            n_cores = myCluster.NumWorkers; % call # of cores on the local machine
        case 'cluster'
            n_cores = str2num(getenv('SLURM_JOB_CPUS_PER_NODE')); % call # of cores on the cluster (SLURM)
            restoredefaultpath;
        end

        parpool(n_cores); % start parallel pool
        restoredefaultpath;
        spmd
            warning('off','all')
        end
    end

    %% Bayesian algorithm variable and model constrains
    % Gibbs sampling 
    bounds1 = [-3 +3]; % parameters upper and lower bound in log10; uniform prior [-2, +2] initialization then normal sampling
    bounds2 = [-3 +3]; % more constrained bounds in log10 for certian modules;
    delta = 0.075; % percentage of the parameters limits range (U-L) to resample parameters in
    nSamples0 = 1e4; % number of parameter sets to be sampled at 1st iteration (~ .33% SUCCESS rate; 3e4 will result in ~100 SUCCESS)
    nSamples = 1e3; % number of parameter sets to be sampled at each iteration
    flattenOBJ = [.01*nSamples 3 1]; % flatten [1% 3 none] lowest OBJs to their max
    iMax = 5e3; % number of iterations
    eta = 0.5;  % update inertia
    beta = 0; %1e-2; % coefficient of the identity matrix added to the covariance matrix (widen sampling)
    zeta = 0; % coefficient to sample % of parameter sets from prior at each iteration
    gamma = 1.0*iMax; % the percentage of iterations to do mvnrnd sampling before disperse sampling
    lin_dist_slope = 2; % for re_sampling the outliers, the slope of linear distribution, 1 or 2
    save_samples = "yes"; % "yes"  "no"
    % Metropolis–Hastings MCMC sampling
    MetHastPars.N_Chain = 5e5;
    MetHastPars.N_Burn = -1;
    MetHastPars.N_Thin  = 1;
    MetHastPars.perc = .50; % percent of the parameters to alter during Met Hast search
    MetHastPars.delta = 1e-3;% coefficient to alter the parameters during Met Hast search
    MetHastPars.intervals = 400;
    MetHastPars.re_sampling_outliers = re_sampling_outliers; 
    
    % Minimum signal criteria on nodes to make all sub-branches and components contribute to Hog1 activation
    % [[1:length(model.DynamicNodes)]' model.DynamicNodes]
    MinimalSignalNodes = struct(... NodeName,Threshold value
        'Sln1',.04, 'Ssk222p',.04,... SLN1
        'Hkr1a',.04, 'Msb2a',.04, 'Sho1a',.04, 'Ste50pSte11',.04, 'Ste50Ste11p',.04,... SHO1
        'Pbs2Sp',.02, 'Pbs2Tp',.02, 'Pbs2SpTp',.02, 'X',.02,... PBS2
        'Hog1cytoYp',.02, 'Hog1cytoTp',.02,'Hog1cytoYpTp',.02, 'Hog1nuc',.02, 'Hog1nucYp',.02, 'Hog1nucTp',.02, 'Hog1nucYpTp',.02,... HOG1
        'Rgc12p',.02, 'Gpd12trans',.02, 'Glyc__synth',.02... Glyc
        );     

    %% load models and data and set seeds
%     dirName =  ['chain',num2str(chain,'%04d'),'/']; mkdir(dirName);
    dirName =  [training_condition, '/', 'chain',num2str(chain,'%04d'),'/']; mkdir(dirName);

        %load models
    MMM = load(['../../HJ_202203_HOGSignalingModeling/HJ_20220301_HOGModels/',model_version,'/Models.mat']);
    Model = MMM.Model;

        % load data
    DDD = load('../../HJ_202203_HOGSignalingModeling/HJ_20220301_HOGDATA/HogDataSet/HogDataSet.mat');
    HogDataSet = DDD.HogDataSet;
    DataSet1 = HogDataSet.DataSet1; % 100 monotonic stimulations of 25min
    DataSet3 = HogDataSet.DataSet3; % staircases

        % randomize the seed for rng (important to initiate different chains)
    rng('shuffle');
    seeds = randi(2^32,1e4,1);
    rng(seeds(randi(1e4,1)));

        % training data and training model ODEs
    DataSet1_IDs = [1:length(DataSet1)]; % DataSet1 indexes
    n_traindata=length(train_data_IDs);
    test1_data_IDs = DataSet1_IDs; % test data IDs
    for ssss=train_data_IDs
        test1_data_IDs=test1_data_IDs(test1_data_IDs~=ssss);
    end

    train_inputs = cell(1,n_traindata);
%     train_inputs = cell(1,n_traindata+1);
%     DataSet3{2}.mean_sHog1nuc(end:end+10) = DataSet3{2}.mean_sHog1nuc(end);
%     train_inputs{1} = DataSet3{2}; % include staircase(0.2M,25min) in training dataset
    dd100tt=DataSet1{train_data_IDs(end)}.tt; len0=length(dd100tt); lenn=len0+100;  
    for i=1:length(train_data_IDs)
        dd = DataSet1{train_data_IDs(i)};
        %  include post-stimuli response in training data
        if length(dd.tt)<lenn
            dd.tt = [dd.tt; dd.tt(end)+[1:100]'];
        end
        dd.mHog1nuc = [dd.mHog1nuc; dd.mHog1nuc(1:100)];
        dd.mean_sHog1nuc = [dd.mean_sHog1nuc; dd.mean_sHog1nuc(1:100)];
        if i==1
            dd.HOG1P5 = [dd.mHog1_non_phospho(:,1) dd.mHog1Y176_phospho(:,1) dd.mHog1T174_phospho(:,1) dd.mHog1YpTp(:,1) dd.mHog1_phospho(:,1)];
            INDEXp=[101:length(dd.tt)]';
            HOG1P5 = dd.HOG1P5(INDEXp,:);
            wtt0=1e-2*(.1*sum(~isnan(HOG1P5))).*ones(100,5); 
            dd.wtt_p=[wtt0;~isnan(HOG1P5)]/sum(sum(([wtt0;~isnan(HOG1P5)])));
        end
        dd.nV=HogDataSet.nucV(1); 
        dd.bHog1n = HogDataSet.bHog1nuc(1); 
        train_inputs{i} = dd;
    end
    n_traindata = length(train_inputs);

    model = Model{model_index};
    train_models = cell(1,n_traindata);
    for i=1:length(train_inputs)
        train_models{i}.model = Get_ODE(model,train_inputs{i}.stimulus);
    end
    model=train_models{1}.model; n_params=model.n_params;

    number_of_singlecells = 100;
    cell_id = mod(chain,number_of_singlecells); if cell_id==0; cell_id=100; end

    %% set training constarins (parameter bounds, nodes ID and nodes minimal signals conditions)
    constrains.delta=delta;
    constrains.log10LB=ones(1,n_params);  constrains.log10UB=ones(1,n_params); % lower and upper bound initiation
    bounds1_IDs = [model.PBS2paramsID' model.HOG1paramsID'];
    bounds2_IDs = [model.SLN1paramsID' model.SHO1paramsID' model.Glyc_paramsID'];

    constrains.log10LB(bounds1_IDs)=bounds1(1); constrains.log10UB(bounds1_IDs)=bounds1(2);
    constrains.log10LB(bounds2_IDs)=bounds2(1); constrains.log10UB(bounds2_IDs)=bounds2(2);
    constrains.LB = 10.^constrains.log10LB; constrains.UB = 10.^constrains.log10UB;

        % initialize parameters with uniform prior [LB, UB]
    InitialPopulation = NaN(nSamples0,n_params); % NaN
    InitialPopulation(:,bounds1_IDs) = bounds1(1) + range(bounds1)*rand(nSamples0,length(bounds1_IDs)); % uniform prior [LB, UB]
    InitialPopulation(:,bounds2_IDs) = bounds2(1) + range(bounds2)*rand(nSamples0,length(bounds2_IDs)); % uniform prior [LB, UB]
          
    MinimalSignalNodesNames = fieldnames(MinimalSignalNodes);
    for i=1:length(MinimalSignalNodesNames)
        MinimalSignalNodesID(i) = find(model.DynamicNodes==sym(MinimalSignalNodesNames{i}));
%         MinimalSignalNodesValues(i) = getfield(MinimalSignalNodes,MinimalSignalNodesNames{i});
        MinimalSignalNodesValues(i) = model.DynamicNodesIC(MinimalSignalNodesID(i),2)+getfield(MinimalSignalNodes,MinimalSignalNodesNames{i});
    end
    constrains.MinimalSignalNodesID = MinimalSignalNodesID; % model.DynamicNodes(constrains.MinimalSignalNodesID)
    constrains.MinimalSignalNodesValues = MinimalSignalNodesValues; % signal range is [0 1], criteria: a min 10%, 5%, or 3% signal over time for the specified nodes
        
    train_models{1}.constrains = constrains;
    train_models{1}.cell_id = cell_id; train_models{1}.pop_or_singlecell = pop_or_singlecell;
    trainDataModel.train_models = train_models; trainDataModel.train_inputs = train_inputs;

    % parameters mean and cov
    muu0 = mean(InitialPopulation); % mean parameters
    covv0 = cov(InitialPopulation); % covariance parameters

    % load posterior from optimization using mean population data to use for optimization using single cell data
    switch pop_or_singlecell
        case "SingleCells"
            % load muu and covv
    end

    % define the objective function for model and training data
    [OBJ] = @(paramset) objective_function(train_models,train_inputs,10.^paramset);

    %% optimization
    if save_samples == "yes"
        posterior_params = NaN([iMax 100 n_params]); % nSamples
    end
    posterior_Likelihoods = NaN([iMax nSamples]);
    posterior_mu = NaN([iMax n_params]);
    posterior_cov = NaN([iMax n_params n_params]);
    posterior_OBJs = NaN([1 iMax]);

    switch prior
        case "uniform" % 1st interation uniform prior
            paramsets = InitialPopulation;
        case "lognormal" % 1st interation lognormal prior
            [paramsets] = sample_parameters(muu0,covv0,nSamples0,constrains,[],[],"mvnrnd","","");
    end
    
    if GIBBS.do=="yes"
    muu=muu0; covv= covv0; posterior=[];
    for i = 1:iMax
        clear t0; t0=tic;
        disp(['i = ', num2str(i), ' started.'])
        constrains.delta=delta*iMax/(iMax+i);

        % sample parameter sets
        if i<=gamma; sampling = "mvnrnd"; else; sampling = "disperse"; end
        if i>1
            if zeta~=0;[paramsets0] = sample_parameters(muu0,covv0,(zeta)*nSamples,constrains,[],[],"mvnrnd","","");else;paramsets0=[];end
            [paramsets_] = sample_parameters(muu,covv,(1-zeta)*nSamples,constrains,posterior,lin_dist_slope,sampling,resample_dist,re_sampling_outliers);
            paramsets = [paramsets0; paramsets_];
        end
%         [fig] = plot_MVN(paramsets,free_parameters,i,dirName);

        % solve the model and calculate OBJ for each parameter set
        OBJs = inf(1,size(paramsets,1));
        parfor j=1:size(paramsets,1)
            OBJs(j) = OBJ(paramsets(j,:)');
        end
        OBJs(isnan(OBJs)==1)=inf; % replace any NaN with inf
        discards = length(find(isinf(OBJs)==1)); % quantify the # parameter set that ODE solver was unsuccesful or nodes min signal criteria did not met.
        SUCCESSi = size(paramsets,1) - discards;
        SUCCESS(i) = SUCCESSi;

        [sOBJs, INXs] = sort(OBJs, 'ascend'); % sort OBJs
        % flatten the lowest [1% 3 none] of OBJs to the highest OBJ among them
        if i<=.2*iMax; flattenOBJ_ = flattenOBJ(1); elseif i<=.9*iMax; flattenOBJ_ = flattenOBJ(2);else flattenOBJ_ = flattenOBJ(3);end 
        if SUCCESSi>(10*flattenOBJ_)
            OBJs(INXs(1:flattenOBJ_)) = OBJs(INXs(flattenOBJ_));
%         elseif SUCCESSi>20
%             OBJs(INXs(1:3)) = OBJs(INXs(3));
        elseif SUCCESSi>10
            OBJs(INXs(1:2)) = OBJs(INXs(2));
        end
        nL=exp(-OBJs)/nansum(exp(-OBJs));  % likelihood of parameter sets
        posterior_Likelihoods(i,:)=nL(INXs(1:nSamples));
        if save_samples == "yes"
            posterior_params(i,:,:) = paramsets(INXs(1:100),:); %nSamples
        end
        posterior.PARAS=paramsets(INXs(1:nSamples),:); posterior.LIKELIHOOD=nL(INXs(1:nSamples));

        % weighted mean and covariance based on the likelihood values
        this_mu = (nL*paramsets);  % weighted mean of parameters
        this_cov = nL.*(paramsets-this_mu)'*(paramsets-this_mu);  % weighted coveriance of parameters

        % add inertia terms to mu and cov
        muu = eta*this_mu+(1-eta)*muu;
		covv = eta*this_cov+(1-eta)*covv;
		if i<iMax/10
			covv = covv + beta*eye(n_params);
        end
        covv = .5*(covv+covv');
        posterior.mu = muu; posterior.cov=covv; 

        % collect mu and cov through the iterations
        posterior_mu(i,:) = muu;
        posterior_cov(i,:,:) =covv;
        posterior_OBJs(i)=sOBJs(1); % min OBJ this iteration
        sim_times(i) = toc(t0)/60;
        disp(['success: ', num2str(SUCCESSi), ', best OBJ value is ', num2str(sOBJs(1))]);
        disp(['time = ', num2str(toc(t0)/60), 'min | i = ', num2str(i), ' finished.'])
        disp(' ')
        if SUCCESSi==0
            return
        end
    end        
    %% save results if successful  % i=1;
    if i==iMax
        try
            sim_time = nanmean(sim_times);
            OBJv = sOBJs(1);
            best_pars = 10.^paramsets(INXs(1),:)';

            FITS_PREDICTIONS_SENS.dirName = dirName;
            FITS_PREDICTIONS_SENS.sim_time = sim_time;
            FITS_PREDICTIONS_SENS.sim_data_indx = DataSet1_IDs;
            FITS_PREDICTIONS_SENS.train_data_indx = train_data_IDs;
            FITS_PREDICTIONS_SENS.test_data_indx = test1_data_IDs;
            FITS_PREDICTIONS_SENS.model_version = model_version;
            FITS_PREDICTIONS_SENS.model_index = model_index;
            FITS_PREDICTIONS_SENS.trainDataModel = trainDataModel;
            FITS_PREDICTIONS_SENS.SUCCESS = SUCCESS;
            FITS_PREDICTIONS_SENS.posterior_Likelihoods = posterior_Likelihoods;
            FITS_PREDICTIONS_SENS.OBJ = OBJv;
            FITS_PREDICTIONS_SENS.OBJs = posterior_OBJs;
            FITS_PREDICTIONS_SENS.posterior_mu = posterior_mu;
            FITS_PREDICTIONS_SENS.posterior_cov = posterior_cov;
            FITS_PREDICTIONS_SENS.best_pars = best_pars;
            if save_samples == "yes"
                FITS_PREDICTIONS_SENS.posterior_params = posterior_params;
            end

            save([dirName,'FITS_PREDICTIONS_SENS'], 'FITS_PREDICTIONS_SENS');
        end
    end
    end
    
    %% Run Markov chain Monte Carlo (MCMC) Metropolis–Hastings search. Start: best_pars from previous optimization
    if doMetHast=="yes"
            if GIBBS.do=="yes"
                initial_guess = log10(best_pars);
            elseif GIBBS.do~="yes"
                % process Gibss results
                [initial_guess,covv,covv1] = processGIBBSresults(n_params,GIBBS,training_condition);                
                FITS_PREDICTIONS_SENS.GIBBS=GIBBS; 
                FITS_PREDICTIONS_SENS.start_MH = initial_guess;
                FITS_PREDICTIONS_SENS.cov = covv; 
                FITS_PREDICTIONS_SENS.cov1 = covv1;  
%                 FITS_PREDICTIONS_SENS.OBJs = FPS.FITS_PREDICTIONS_SENS.OBJs;                
                i=iMax;             
            end
        if i==iMax
            clear t0; t0=tic;
            disp([' Metropolis–Hastings search started ...'])
            [MetHastParChain,MetHastFunChain] = MetHast(OBJ,initial_guess,covv,MetHastPars,constrains);
            [mMetHastFunChain,INDX] = min(MetHastFunChain);
            best_pars = MetHastParChain(:,INDX);
            FITS_PREDICTIONS_SENS.best_pars_MH = best_pars;
            FITS_PREDICTIONS_SENS.MetHastParChain = MetHastParChain;
            FITS_PREDICTIONS_SENS.MetHastFunChain = MetHastFunChain;
            save([dirName,'FITS_PREDICTIONS_SENS'], 'FITS_PREDICTIONS_SENS');
            disp(['Metropolis–Hastings search finished (time = ', num2str(toc(t0)/60), 'min).'])
            disp(['The best [obj, id] is [',num2str(mMetHastFunChain), ', ',num2str(INDX),'].'])
        end
    end

    %% get predictions if successful
    if i==iMax        
        % get predictions upon each stimuli input for all mutant strains
        disp(['get predictions upon each stimuli input for all mutant strains  '])
        FITS_PREDICTIONS_SENS = get_predictions(dirName,FITS_PREDICTIONS_SENS,Model,model,best_pars,HogDataSet,train_data_IDs,test1_data_IDs,pop_or_singlecell,cell_id);

        % get model sensitivity to model parameters
        if do_sens_analysis=="yes"
            npoints = 101; perc = 0.1; % number of sampled parameter values, percentage within which to vary each parameter
            clear dd; dd.nV=HogDataSet.nucV(1); dd.bHog1n = HogDataSet.bHog1nuc(1); %set yeast cells nuc/total size and pre-stimuli Hog1nuc from data
            [OBJs0,MODEL_SENSITIVITY] = model_sensitivity_to_parameters(model,best_pars,npoints,perc,DataSet1,pop_or_singlecell,cell_id,dd);
            FITS_PREDICTIONS_SENS.MODEL_SENSITIVITY = MODEL_SENSITIVITY;
            FITS_PREDICTIONS_SENS.OBJs0 = OBJs0;
            save([dirName,'FITS_PREDICTIONS_SENS'], 'FITS_PREDICTIONS_SENS');
        end
        disp(['   job succeeded.'])
    elseif i<iMax
        try
        rmdir(dirName,'s'); 
        end
        disp(['   job terminated.'])
    end
    % shutdown parallel pool
    if doGIBBS=="yes"
        poolobj = gcp('nocreate');
        delete(poolobj);
    end
end
