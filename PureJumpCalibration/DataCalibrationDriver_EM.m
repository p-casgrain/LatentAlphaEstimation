%% Add to path
addpath('./Tools/');
addpath('./PoissonFunctions_mex/');
addpath('./HMM_mex/');
addpath('./MiscCalibration/');


%% Define Global Data Parameters
dt = 1; % Size of time intervals e.g. in minutes
b = 0.01; % Tick size in dollars

%% Import Pre-Parsed Example Data
% Data Arrays X and X_normalized are imported
% X consists of 249 price paths of length 3600
load('../Example_Data.mat') 

DX = diff(X,1,2); % Compute Price Increments

%% Pre-Calibration Step
% This part computes a good guess for the initial parameter values the EM
% procedure should start at.

InitParams = preCalibrate(X,1,dt);

%% No Hidden State MLE Estimation
% This part provides another method for computing good initial parameter
% guesses.

InitParams1D = InitParams;
InitParams1D.theta = mean(InitParams.ThetaValues);

% InitParams.mu = 0.0333;
% InitParams.kappa = 0.0783;

OutParams1D = onestateMLE(X,DX,InitParams1D);

%% Run EM Algorithm
% This section runs the EM algorithm starting at the initial parameter
% guesses contained in the cell array InitParams, and prints the results.

tic;
[OutParams,OutPost,OutLogLikelihood] = HMMmaximize(X,DX,InitParams,1e-5,1e-3,1e3);
toc;

disp('Calibrated Transition Matrix')
disp(OutParams.Q);
disp('Calibrated Parameters')
disp(OutParams);


