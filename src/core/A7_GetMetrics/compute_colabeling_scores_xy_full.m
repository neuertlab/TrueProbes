function scores = compute_colabeling_scores_xy_full(xvals, yvals, p, opts)
% Computes S_A, S_B, S_C, S_D scores, expectation, variance, Fano factor, and full distribution
% xvals, yvals, p: N x 1 vectors
% opts: alpha,beta,rho,gamma,kappa,w00,kernel,eps, nbins

if ~isfield(opts,'nbins'), opts.nbins = 100; end
% --- normalize ---
if ~isfield(opts,'xmin'), opts.xmin = min(xvals); end
if ~isfield(opts,'xmax'), opts.xmax = max(xvals); end
if ~isfield(opts,'ymin'), opts.ymin = min(yvals); end
if ~isfield(opts,'ymax'), opts.ymax = max(yvals); end
if ~isfield(opts,'eps'), opts.eps = 1e-9; end
if ~isfield(opts,'alpha'), opts.alpha = 1; end
if ~isfield(opts,'beta'), opts.beta = 1; end
if ~isfield(opts,'rho'), opts.rho = 1; end
if ~isfield(opts,'gamma'), opts.gamma = 2; end
if ~isfield(opts,'kappa'), opts.kappa = 4; end
if ~isfield(opts,'w00'), opts.w00 = 0; end
if ~isfield(opts,'kernel'), opts.kernel='exp'; end

u = (xvals - opts.xmin)./(opts.xmax-opts.xmin);
v = (yvals - opts.ymin)./(opts.ymax-opts.ymin);
u = min(1,max(0,u)); v = min(1,max(0,v));

m = 0.5*(u+v); b = 1-abs(u-v); sig = (2*u-1).*(2*v-1);
logr = log( (xvals+opts.eps)./(yvals+opts.eps) );
c_ratio = exp(-opts.gamma*abs(logr));

% --- per-point scores ---
% --- A ---
SA_point = sig .* (2*(m.^opts.alpha).*(b.^opts.beta)-1);
% --- B ---
d11 = sqrt((u-1).^2 + (v-1).^2)/sqrt(2);
d10 = sqrt((u-1).^2 + (v-0).^2)/sqrt(2);
d01 = sqrt((u-0).^2 + (v-1).^2)/sqrt(2);
d00 = sqrt((u-0).^2 + (v-0).^2)/sqrt(2);
z11 = 1-d11; z10=1-d10; z01=1-d01; z00=1-d00;
a11 = apply_kernel(z11,opts.kernel,opts.kappa);
a10 = apply_kernel(z10,opts.kernel,opts.kappa);
a01 = apply_kernel(z01,opts.kernel,opts.kappa);
a00 = apply_kernel(z00,opts.kernel,opts.kappa);
SB_point = a11 - max(a10,a01) - opts.w00*a00;
% --- C ---
SC_point = sig .* (2*(m.^opts.alpha).*(c_ratio.^opts.rho)-1);
% --- D ---
L = a11 ./ (a11 + max(a10,a01) + opts.w00*a00 + 1e-12);
SD_point = (2*(m.^opts.alpha).*(b.^opts.beta).*(c_ratio.^opts.rho)-1) .* L;

% --- Expectation under pmf ---
points = {SA_point, SB_point, SC_point, SD_point};
fields = {'SA','SB','SC','SD'};

p = p/sum(p,'all');
for f = 1:length(fields)
    s_point = points{f};
    
    scores.(fields{f}).P_per_point = p;
    scores.(fields{f}).S_per_point = s_point;
    scores.(fields{f}).S = sum(p .* s_point,'all');
    
    % Binned Full distribution
    [S_unique_binned, P_S_binned] = score_distribution(s_point, p, opts.nbins);
    
    % Moments
    mu_binned = sum(S_unique_binned .* P_S_binned);
    var_s_binned = sum((S_unique_binned - mu_binned).^2 .* P_S_binned);
    fano_binned = var_s_binned / (mu_binned + eps);
    
    scores.(fields{f}).binned.mean = mu_binned;
    scores.(fields{f}).binned.var = var_s_binned;
    scores.(fields{f}).binned.fano = fano_binned;
    scores.(fields{f}).binned.S_unique = S_unique_binned;
    scores.(fields{f}).binned.PPF_S = P_S_binned;
    scores.(fields{f}).binned.CDF_S = cumsum(P_S_binned);
    scores.(fields{f}).binned.Survivor_S = 1-cumsum(P_S_binned);
    scores.(fields{f}).binned.CDF_Percentile_To_Score = @(x) interp1(scores.(fields{f}).binned.CDF_S,scores.(fields{f}).binned.S_unique,x);
    scores.(fields{f}).binned.Score_To_CDF_Percentile = @(x) interp1(scores.(fields{f}).binned.S_unique,scores.(fields{f}).binned.CDF_S,x);
    scores.(fields{f}).binned.Survivor_Percentile_To_Score = @(x) interp1(scores.(fields{f}).binned.Survivor_S,scores.(fields{f}).binned.S_unique,x);
    scores.(fields{f}).binned.Score_To_Survivor_Percentile = @(x) interp1(scores.(fields{f}).binned.S_unique,scores.(fields{f}).binned.Survivor_S,x);

    %Exact Full distribution
    [S_unique_exact, ~, idx] = unique(s_point(:));
    P_S_exact = accumarray(idx, sparse(p(:)));

    % Moments
    mu_exact = sum(S_unique_exact .* P_S_exact);
    var_s_exact = sum((S_unique_exact - mu_exact).^2 .* P_S_exact);
    fano_exact = var_s_exact / (mu_exact + eps);
    
    scores.(fields{f}).exact.mean = mu_exact;
    scores.(fields{f}).exact.var = var_s_exact;
    scores.(fields{f}).exact.fano = fano_exact;
    scores.(fields{f}).exact.S_unique = S_unique_exact;
    scores.(fields{f}).exact.PPF_S = P_S_exact;
    scores.(fields{f}).exact.CDF_S = cumsum(P_S_exact);
    scores.(fields{f}).exact.Survivor_S = 1-cumsum(P_S_exact);
    scores.(fields{f}).exact.CDF_Percentile_To_Score = @(x) interp1(scores.(fields{f}).exact.CDF_S,scores.(fields{f}).exact.S_unique,x);
    scores.(fields{f}).exact.Score_To_CDF_Percentile = @(x) interp1(scores.(fields{f}).exact.S_unique,scores.(fields{f}).exact.CDF_S,x);
    scores.(fields{f}).exact.Survivor_Percentile_To_Score = @(x) interp1(scores.(fields{f}).exact.Survivor_S,scores.(fields{f}).exact.S_unique,x);
    scores.(fields{f}).exact.Score_To_Survivor_Percentile = @(x) interp1(scores.(fields{f}).exact.S_unique,scores.(fields{f}).exact.Survivor_S,x);
end

scores.corners.a11 = a11; scores.corners.a10 = a10;
scores.corners.a01 = a01; scores.corners.a00 = a00;
end

%% Helper functions
function a = apply_kernel(z,type,kappa)
    z = min(1,max(0,z));
    switch lower(type)
        case 'linear', a=z;
        case 'exp'
            if kappa<=0, a=z; else a=(1-exp(-kappa*z))/(1-exp(-kappa)); end
        case 'gauss'
            if kappa<=0, a=z; else a=(exp(-kappa*(1-z).^2)-exp(-kappa))/(1-exp(-kappa)); end
        otherwise, a=z;
    end
end

function [S_unique, P_S] = score_distribution(S_point, p, nbins)
Smin = min(S_point(:)); Smax = max(S_point(:));
edges = linspace(Smin, Smax, nbins+1);
S_unique = zeros(nbins,1); P_S = zeros(nbins,1);
for i=1:nbins
    idx = S_point >= edges(i) & S_point < edges(i+1);
    P_S(i) = sum(p(idx));
    S_unique(i) = 0.5*(edges(i)+edges(i+1));
end
P_S = P_S / sum(P_S); % normalize

end
