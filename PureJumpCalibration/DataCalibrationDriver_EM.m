%% Add to path
addpath('./Tools/');
addpath('./PoissonFunctions_mex/');
addpath('./HMM_mex/');
addpath('./MiscCalibration/');


%% Import LOB Data & Parse
if ~exist('LOB','var')
    load('../../Data/Parsed/LOB_Parsed_Data_10-11AM.mat')
end

% dt = 1/3600; % Data time interval in hours (essentially only needed for mean parameter scale)
dt = 1;
b = 0.01; % Tick size in dollars
startTime = 10; % Window from 10am
endTime = 11;   % to 11am

if ~exist('X','var') || 1
    % Get all of the midprices between the start and end time
    
    startInd = (startTime-9.5)*3600;
    endInd = (endTime-9.5)*3600;
    TimeIntInd = startInd:endInd;

    X = zeros(numel(LOB),numel(TimeIntInd));
    
    for m=1:numel(LOB)
        MidPrice = 0.5*LOB{m}.BuyPrice(TimeIntInd) ...
                 + 0.5*LOB{m}.SellPrice(TimeIntInd) ;
        X(m,:) = MidPrice.'*1e-4;
    end
    
    DX = diff(X,1,2);
    X = X(:,1:end-1);
    
    X = bsxfun(@minus,X,X(:,1));
end


%% Pre-Calibration

InitParams = preCalibrate(X,1,dt);

%% No Hidden State MLE Calibration

InitParams1D = InitParams;
InitParams.mu = 0.0333;
InitParams1D.kappa = 0.0783;
InitParams1D.theta = mean(InitParams.ThetaValues);

OutParams1D = onestateMLE(X,DX,InitParams1D);

%% Run Baum-Welch type of thing
% load('./TempData/TempOutParams-02-Mar-2017 13:40:47.mat');
InitParams=OutParams;

% InitParams.mu = 0.1.*[1+1e-2,1-1e-2,1];
% InitParams.kappa = mean(InitParams1D.kappa).*[1,1,1];
% InitParams.ThetaValues = mean(InitParams1D.ThetaValues).*[1,1,1];
% InitParams.Q = ones(3)/3;
% InitParams.nu = (1/3)*[1,1,1];

% InitParams.mu = [0.07,0.08];
% InitParams.kappa = [0.068,0.088];
% InitParams.ThetaValues = [0.02,0.02];
% InitParams.Q = [0.997,1-0.997;1-0.997,0.997];
% InitParams.nu = [0.2,0.8];

% InitParams = preCalibrate(X,numel(OutParams.mu)+1,dt);
% InitParams.mu = [OutParams.mu , mean(OutParams.mu)];
% InitParams.kappa = [OutParams.kappa , mean(OutParams.kappa)];
% InitParams.ThetaValues = [OutParams.ThetaValues , mean(OutParams.ThetaValues)];

tic;
[OutParams,OutPost,OutLogLikelihood] = HMMmaximize(X,DX,InitParams,1e-5,1e-3,1e3);
toc;

disp('Calibrated Transition Matrix')
disp(OutParams.Q);
disp('Calibrated Parameters')
disp(OutParams);


