function [OBJ] = objective_function(params,train_inputs,train_models,sparsobj,wwt)

constrains = train_models{1}.constrains; 
cell_id = train_models{1}.cell_id; pop_or_singlecell = train_models{1}.pop_or_singlecell; 

params = max([params'; constrains.LB])';
params = min([params'; constrains.UB])';

mm = train_models{1}; 
DynamicNodesMaxVal = mm.model.DynamicNodesMaxVal(:,2); % Dynamic Nodes Max Values (normalized to Hog1), <moleculaes/cell> 
ConstNodesVal = mm.model.ConstNodesVal(:,2);   % Constitute Nodes Values (normalized to Hog1), <moleculaes/cell> 
BasalNodesVal = mm.model.BasalNodesVal(:,2);   % Basal Nodes Values (normalized to Hog1), <moleculaes/cell> 
TunableParamsVal = mm.model.TunableParamsVal;  % Tunable paramters Values (binding affinity) 
DynamicNodesIC = mm.model.DynamicNodesIC(:,2); % Dynamic Nodes Initial Conditions Values (normalized to Hog1), <moleculaes/cell>
    
OBJ = 0; 
for i=1:length(train_inputs)-1

    dd = train_inputs{i}; 
    mm = train_models{i}; 
    
    switch pop_or_singlecell
        case "PopAverage"
            data_mean = dd.mHog1nuc;
        case "SingleCells"
            data_mean = dd.scHog1nuc(:,cell_id);
    end
    data_sigma = dd.mean_sHog1nuc;
    
    % set model ODEs for the data stimulus input profile, s(t)
    ODE = @(t,x)mm.model.ODE(t,x,params,BasalNodesVal,DynamicNodesMaxVal,ConstNodesVal,TunableParamsVal);
    JAC = @(t,x)mm.model.Jacobian(t,x,params,BasalNodesVal,DynamicNodesMaxVal,ConstNodesVal,TunableParamsVal);
    options = odeset('Jacobian',JAC);

    clear time_obj; time_obj=tic; 
    [~,yout] = ode23s(ODE,dd.tt*60,DynamicNodesIC,options); % ode23s ode15s
    if toc(time_obj)<0.5 % avoid potential ode23s solver stock
        NodesDynamics = real(yout); % size(NodesDynamics)
        model_mean = mm.model.Observables(NodesDynamics);
        
        if sparsobj==0
            wt = ones(size(data_mean));
        else  % Use only a sparse set of the data in each evaluation.
            wt=wwt{i};
        end
        wt = wt/sum(wt);

        try
            OBJ = OBJ + sum(((model_mean-data_mean)./data_sigma).^2.*wt);
        catch % ode23 failor 
            disp('ode23 failor')
            check_vec_sizes = [size(model_mean); size(data_mean); size(data_sigma)] 
            OBJ = OBJ + Inf;
            [sort_params, pars_ids] = sort(abs(params), 'ascend'); 
            disp(['4mins and 4maxs |params|:  ', num2str(sort_params([1:4 end-3:end])')]) 
            disp(['4mins and 4maxs |params| ids:  ', num2str(pars_ids([1:4 end-3:end])')]) 
        end

    else
        OBJ = OBJ + Inf; 
    end
end 

%% add control data to OBJ; all nodes must stay unchanged (at their steady state IC) upon a constant basal input of u = 0.1
i = length(train_inputs); 
dd = train_inputs{i}; 
mm = train_models{i}; 

switch pop_or_singlecell
    case "PopAverage"
        basal_response = dd.mHog1nuc;
    case "SingleCells"
        basal_response = dd.scHog1nuc(:,cell_id);
end
data_sigma = dd.mean_sHog1nuc;

% set model ODEs for the data stimulus input profile, s(t)
ODE = @(t,x)mm.model.ODE(t,x,params,BasalNodesVal,DynamicNodesMaxVal,ConstNodesVal,TunableParamsVal);
JAC = @(t,x)mm.model.Jacobian(t,x,params,BasalNodesVal,DynamicNodesMaxVal,ConstNodesVal,TunableParamsVal);
options = odeset('Jacobian',JAC);

clear time_obj; time_obj=tic; 
[~,yout] = ode23s(ODE,dd.tt*60,DynamicNodesIC,options); % ode23s ode15s
if toc(time_obj)<0.5 % avoid potential ode23s solver stock
    
    if sparsobj==0
        wt = ones(size(data_mean));
    else  % Use only a sparse set of the data in each evaluation.
        wt=wwt{i}; 
    end
    wt = wt/sum(wt);

    try
        for node = 1:length(mm.model.DynamicNodes)
            model_mean  = real(yout(:,node));
            data_mean = .1*basal_response + DynamicNodesIC(node); 
            OBJ = OBJ + (1/length(mm.model.DynamicNodes))*sum(((model_mean-data_mean)./data_sigma).^2.*wt);
        end

    catch % ode23 failor 
        disp('ode23s failure SS')
        check_vec_sizes = [size(model_mean); size(data_mean); size(data_sigma)] 
        OBJ = OBJ + Inf;
        [sort_params, pars_ids] = sort(abs(params), 'ascend'); 
        disp(['5mins and 5maxs params:  ', num2str(sort_params([1:5 end-4:end])')]) 
        disp(['5mins and 5maxs |params| ids:  ', num2str(pars_ids([1:5 end-4:end])')]) 
    end
else
    OBJ = OBJ + Inf; 
end
OBJ = OBJ/length(train_inputs); 

global min_OBJ
if OBJ < min_OBJ 
    global best_pars
    min_OBJ = OBJ; 
    best_pars = params; 
end 
end

