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

%% Load in Simulation Parameters
load('./TempData/FinalResults/Direct-OutParams-2state-03-Mar-2017 11:44:47.mat');

% r_spec = @(x) round(x,2,'significant');
r_spec = @(x) arrayfun(@(y) round(y,2,'significant') , x);
TrueParams = OutParams;
% TrueParams.nu = r_spec(TrueParams.nu);
% TrueParams.Q = r_spec(TrueParams.Q);
TrueParams.mu = r_spec(TrueParams.mu);
TrueParams.kappa = r_spec(TrueParams.kappa);
TrueParams.ThetaValues = r_spec(TrueParams.ThetaValues);

%% Simulate Processes & Assign to Appropriately sized variables.
[S,Z_ind] = SimulatePJProcess(TrueParams,Ndt,dt,Npaths*Ntrials);
X = reshape( S(:,1:(end-1)) , [Ntrials,Npaths,Ndt-1] ); 
DX = reshape( diff(S,1,2) , [Ntrials,Npaths,Ndt-1] );

%% Run Viterbi Algorithm

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


% %% Load in Initial Guess for Parameters (calibrated Parameters for one dimension less)
% load('./TempData/FinalResults/Direct-OutParams-1state-10-Mar-2017.mat');
% 
% ParamsLast = OutParams;
% Nstates = numel(TrueParams.nu);
% 
% noise = (rand(1,3) - 0.5)*2e-2;
% InitParams.mu = [ ParamsLast.mu , mean(ParamsLast.mu) + noise(1)]/2;
% InitParams.kappa = [ ParamsLast.kappa , mean(ParamsLast.kappa) + noise(2)]*2;
% InitParams.ThetaValues = [ ParamsLast.ThetaValues , mean(ParamsLast.ThetaValues) + noise(3)];
% InitParams.nu = ones(1,Nstates)/Nstates;
% 
% diagProb = 0.995; offdiagProb = (1-diagProb)/(Nstates-1);
% InitParams.Q = diag(diagProb*ones(1,Nstates)-offdiagProb) + offdiagProb*ones(Nstates);
% 
% InitParams.Delta = dt;
% 
% 
% %% Create Empty Cell to fill with Calibrated Parameters
% CalibratedParams = cell(1,Ntrials);
% 
% %% Run Repeated Calibration - Parallelized
% options.maxIter=50; options.optimskip=1;
% 
% timestamp = datetime('now','Format','yyyyMMddhhmmss');
% string1 = strcat(num2str(Nstates),'state-');
% filename = strcat('./TempData/SimTest-',string1,datestr(timestamp),'.mat');
% 
% tic;
% parfor k=1:Ntrials
%     CalibratedParams{k} = HMMmaximize(squeeze(X(k,:,:)),squeeze(DX(k,:,:)),InitParams,1e-3,0.25*1e-3,1e5,options);
% end
% toc;
% 
% save(filename,'CalibratedParams');
% 

