% This script simulates the Pure-Jump Latent alpha process and tests to see
% if the EM procedure can accurately recover the parameters that generated
% it.

%% Add to path
addpath('./Tools/');
addpath('./PoissonFunctions_mex/');
addpath('./HMM_mex/');

%% Number of Independent Paths to Simulate and their length
Ndt = 3601;
dt = 1;
Npaths = 249;
Ntrials = 1;

%% Load in Example Parameters for Simulation
load('../Example_Params.mat');

% r_spec = @(x) round(x,2,'significant');
r_spec = @(x) arrayfun(@(y) round(y,2,'significant') , x);
TrueParams = OutParams;
TrueParams.nu = r_spec(TrueParams.nu);
TrueParams.Q = r_spec(TrueParams.Q);
TrueParams.Q = TrueParams.Q ./ sum(TrueParams.Q,2);
TrueParams.mu = r_spec(TrueParams.mu);
TrueParams.kappa = r_spec(TrueParams.kappa);
TrueParams.ThetaValues = r_spec(TrueParams.ThetaValues);

%% Simulate Processes & Assign to Appropriately sized variables.
[S,Z_ind] = SimulatePJProcess(TrueParams,Ndt,dt,Npaths*Ntrials);
X = reshape( S(:,1:(end-1)) , [Ntrials,Npaths,Ndt-1] ); 
DX = reshape( diff(S,1,2) , [Ntrials,Npaths,Ndt-1] );

%% Run Viterbi Algorithm to Obtain Most Likely Path

X = squeeze(X);
DX = squeeze(DX);

[psi,del] = HMMviterbi(X,DX,TrueParams);

idx = 15;

figure(200);
clf;
subplot(1,3,1)
plot(psi(idx,:));
title('Viterbi');
subplot(1,3,2);
plot(Z_ind(idx,:));
title('Real')
subplot(1,3,3);
plot(psi(idx,:)-Z_ind(idx,1:(end-1)));
title('Difference');


%% Initial Guess for Parameters

Nstates = numel(TrueParams.mu);
InitParams = preCalibrate(X,Nstates,dt); % Generate Initial Guesses

InitParams.mu = InitParams.mu + normrnd(0,0.01,1,Nstates); % Perturb guesses with noise
InitParams.kappa = InitParams.kappa + normrnd(0,0.01,1,Nstates);

%% Create Empty Cell to fill with Estimated Parameters
CalibratedParams = cell(1,Ntrials);

%% Run Repeated Estimation
options.maxIter=30; options.optimskip=50;

timestamp = datetime('now','Format','yyyyMMddhhmmss');
string1 = strcat(num2str(Nstates),'state-');
filename = strcat('./TempData/SimTest-',string1,datestr(timestamp),'.mat');

tic;
for k=1:Ntrials
    CalibratedParams{k} = ...
        HMMmaximize(X,DX,InitParams,1e-5,1e-3,1e3,options);
end
toc;

save(filename,'CalibratedParams');


