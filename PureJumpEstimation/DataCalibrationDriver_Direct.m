%% Add to path
addpath('./Tools/');
addpath('./PoissonFunctions_mex/');
addpath('./HMM_mex/');
addpath('./MiscCalibration/');


%% Import LOB Data & Parse

% dt = 1/3600; % Data time interval in hours (essentially only needed for mean parameter scale)
dt = 1;
b = 0.01; % Tick size in dollars
startTime = 10; % Window from 10am
endTime = 11;   % to 11am

if ~exist('X','var')
    
    if ~exist('LOB','var')
        load('../../Data/Parsed/LOB_Parsed_Data_10-11AM.mat')
    end
 
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

% InitParams = preCalibrate(X,6,dt);

%% No Hidden State MLE Calibration

% InitParams1D = InitParams;
% InitParams.mu = 0.0333;
% InitParams1D.kappa = 0.0783;
% InitParams1D.theta = mean(InitParams.ThetaValues);
% 
% OutParams1D = onestateMLE(X,DX,InitParams1D);

%% Run Constrained Direct Maximization

% Starting Values
InitParams = OutParams;
% P0 =  [InitParams.Q(:).',InitParams.nu,InitParams.mu,InitParams.kappa,InitParams.ThetaValues];
P0 = deObjectify(InitParams);

% Constraint Matrices
K = numel(InitParams.nu);
N = K*(K+4);

IneqConA = [zeros(K,K^2+2*K),-diag(ones(K,1)),zeros(K)];
IneqConA = [IneqConA;zeros(K,K^2+K),-diag(ones(K,1)),zeros(K,2*K)];
IneqConA = [IneqConA;-diag(ones(1,K^2+K)),zeros(K^2+K,N-K^2-K)];
IneqConA = [IneqConA;diag(ones(1,K^2+K)),zeros(K^2+K,N-K^2-K)];

IneqConb = [zeros(1,K + K^2 + K + K),ones(1,K^2+K)];

EqConb = ones(1,K+1);
EqConA = zeros(K+1,N);
for i=1:(K+1)
    EqConA( i,K*(i-1)+(1:K) ) = ones(1,K);
end

% Handle
loglikObj = @(P) -HMMObjective(X,DX,P,dt);

options = optimoptions( 'fmincon',...
                        'Display','iter-detailed',...
                        'MaxIterations',1e0,...
                        'MaxFunctionEvaluations',1e0,...
                        'StepTolerance',1e-10,...
                        'ConstraintTolerance',1e-10,...
                        'FunctionTolerance',1e-3  ,...
                        'algorithm','sqp');

[param,nloglik,~,~,~,~,hess] = fmincon(loglikObj,P0,IneqConA,IneqConb,EqConA,EqConb,[],[],[],options);

% Play sound
load gong.mat;
sound(y,3*Fs);

OutParams = Objectify(param,-nloglik,dt); % Store results in struct

% Save Result to Disk
timestamp = datetime('now','Format','yyyyMMddhhmmss');
string1 = strcat(num2str(K),'state-');
filename = strcat('./TempData/Direct-OutParams-',string1,datestr(timestamp),'.mat');
save(filename,'OutParams');


%% Helper Function

invarMean = @(P) [P.nu*P.Q^(1e3)*P.mu.',P.nu*P.Q^(1e3)*P.kappa.',P.nu*P.Q^(1e3)*P.ThetaValues.'];

BIC = @(P,X) P.loglik - log(numel(X))*( (numel(P.mu)*(numel(P.mu)+4)) - (numel(P.mu)+1) );

