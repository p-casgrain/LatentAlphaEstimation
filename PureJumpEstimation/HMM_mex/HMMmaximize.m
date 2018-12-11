function [OutParams,post,HMM_l] = HMMmaximize(X, DX, InitParams, reltol, objtol, maxIt,options)
% HMMmaximize - Implementation of Expectation-Maximization algorithm
%
%  in : X = Price Paths
%       DX = Price Path Increments
%       InitParams = struct with all of the starting parameters
%       reltol = Relative Tolerance of optimization procedure
%       objtol = Absolute Tolerance of optimization procedure
%       maxIt = Maximum number of Iterations
%       options = extra options for intermediate optimization steps
%
% out : OutParams = struct with all of the estimated parameters
%       l = log-likelihood of observed sequence
%
% By Philippe Casgrain, 2016

timestamp = datetime('now','Format','yyyyMMddhhmmss');

% Based on the number of input arguments, set default values for
% error tolerance level and max number of iterations allowed
if nargin<6, maxIt = 100; end
if nargin<5, objtol = 1e-2; end
if nargin<4, reltol = 1e-2; end
if nargin<7, options.optimskip=1;options.maxIter=30; end
if ( any(size(InitParams.mu)~=size(InitParams.kappa)) ) || ...
        ( any(size(InitParams.kappa)~=size(InitParams.ThetaValues)) ) ||  ...
        ( size(InitParams.mu,1)~=1 )
    error('Mismatch in Initial Parameter Dimensions');
end

% Set extra optimization params
optimskip = options.optimskip; maxIter = options.maxIter;

% Load in initial parameters
mu = InitParams.mu; kappa = InitParams.kappa;
theta = InitParams.ThetaValues; Delta = InitParams.Delta;
Q = InitParams.Q; nu = InitParams.nu;

% Get number of indepent paths, path length and number of possible states
M = size(X,1);
Ndt = size(X,2);
K = numel(theta);

% Allocate Memory for 'current' & 'old' Q and nu result
it = 0; oldQ = Q; oldnu=nu;

% Create 'paramarray' vector, with all parameters in sequence
% plus an 'old' version.
paramarray = [mu,kappa,theta];
oldparamarray = paramarray;

% Count number of times numerical optim occured
optimcount = 0;

% Initialize SGD parameters
SGDbetam = 0.95;
SGDbetav = 0.99;
SGDit = 0;
SGDm = zeros(size(paramarray)).';
SGDv = zeros(size(paramarray)).';
alpha = 4e-3; % Base learning rate

% Constraint Matrix
ConMat = -[diag(ones(1,2*K)),zeros(2*K,K)];
ConVec = zeros(1,2*K) - 1e-4;

% Storing Log-Likelihood
HMM_l = 0;
HMM_l_old = 0;

%% Begin E-M Loop

while relTol( reltol,[Q(:).',nu(:).',paramarray], [oldQ(:).',oldnu(:).',oldparamarray] ) ...
        && (it<maxIt) && (norm(HMM_l_old - HMM_l) > objtol) || (it<=1)
    
    it = it + 1;
    
    %% Compute Emission Matrix With Current Parameters
    g = Pois_Transition_Prob(mu,kappa,theta,Delta,X,DX);
    
    %% Re-Estimate Q and nu by using the Baum-Welch algorithm
    
    N = zeros(size(Q));
    nu_m = zeros(size(nu));
    c = zeros(M,Ndt); post = zeros(M,K,Ndt);
    
    for m=1:M
        g_m = squeeze(g(m,:,:));
        
        % compute the posterior distribution for the current parameters
        [phi_m, c_m] = HMMfilter_C(nu, Q, g_m);
        beta_m = HMMsmoother_C(Q, g_m, c_m);
        post_m = phi_m .* beta_m;
        
        % expectation of the number of transitions under the current parameters
        N = N + Q.*(phi_m(:, 1:(end-1))*(beta_m(:, 2:end).*g_m(:, 2:end)./(ones(K, 1)*c_m(2:end)))');
        
        nu_m = nu_m + min(post_m(:,1).',1);
        
        c(m,:) = c_m;
        
        post(m,:,:) = post_m;
    end
    
    % re-estimation
    oldQ = Q; oldnu = nu;
    
    Q = N ./ (sum(N, 2) * ones(1, K));
    nu = nu_m/sum(nu_m);
    
    
    %% Re-Estimate mu, kappa and theta via numerical optimization
    
    if  ( ~mod(it,optimskip) ) || (it == 1)
        
        optimcount = optimcount+1;
        
        oldparamarray = paramarray;
        
        loglik_handle = @(y) Pois_Lik_Objective_Merged(y(1:K),y((K+1):(2*K)),y((2*K+1):end),InitParams.Delta,X,DX,post);
        
        options = optimoptions( 'fmincon',...
            'SpecifyObjectiveGradient',true,...
            'CheckGradients',false,...
            'Display','off',...
            'MaxIterations',maxIter);
        
        paramarray = fmincon(loglik_handle,paramarray,ConMat,ConVec,[],[],[],[],[],options);
        
        % Storing current log-likelihood
        HMM_l_old = HMM_l;
        HMM_l = sum(log(c(:)));
        
        % Reassign the individual parameters to make the code more readable
        mu = paramarray(1:K); kappa = paramarray((K+1):(2*K));
        theta = paramarray((2*K+1):end);
        
        % Display new parameter value & log likelihood
%         disp(strcat(num2str(paramarray),' - Log-Likelihood: ',num2str(sum(log(c(:))))));
        
    end
    
    %%  Store Values to Disk every 25 optim. steps
%     if ~mod(optimcount,25)
%         OutParams.mu = mu; OutParams.kappa = kappa;
%         OutParams.ThetaValues = theta;
%         OutParams.Q = Q; OutParams.nu = nu;
%         OutParams.Delta = Delta;
%         OutParams.loglik = sum(log(c(:))) ;
%         
%         filename = strcat('./TempData/TempOutParams-',datestr(timestamp),'.mat');
%         save(filename,'OutParams');
%     end
   
end

%% Final Parameter Output & Result Saving

if (it+1)>=maxIt
    disp('Reached Maximum Number of Iterations');
end

disp(['Number of iterations taken:',num2str(it)]);

OutParams.mu = mu; OutParams.kappa = kappa;
OutParams.ThetaValues = theta;
OutParams.Q = Q; OutParams.nu = nu;
OutParams.Delta = Delta;

% filename = strcat('./TempData/TempOutParams-',datestr(timestamp),'.mat');
% save(filename,'OutParams');


HMM_l = sum(log(c(:))); % Log likelihood of observed data


%% Projection Helper Function

function out = proj(X)
    NX = [max(X(1:(2*K)),1e-10) , X((2*K+1):end)];
    out = NX/norm(NX)*norm(X);
end


end
