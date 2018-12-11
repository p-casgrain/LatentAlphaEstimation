addpath('./Tools/');
addpath('./PoissonFunctions_mex/');
addpath('./HMM_mex/');
addpath('./MiscCalibration/');


%% Import LOB Data

% dt = 1/3600; % Data time interval in hours (essentially only needed for mean parameter scale)
dt = 1;
b = 0.01; % Tick size in dollars


%% Import Pre-Parsed Example Data
% Data Arrays X and X_normalized are imported
% X consists of 249 price paths of length 3600
load('../Example_Data.mat') 

DX = diff(X,1,2); % Compute Price Increments



%% Compare BIC and ICL

filedir = './TempData/FinalResults';
matfiles = dir(fullfile(filedir, '*.mat'));
nfiles = length(matfiles);
Params = cell(1,nfiles);
for i = 1 : nfiles
   fid = load(fullfile(filedir, matfiles(i).name));
   Params{i} = fid.OutParams;
end

BIC = @(P,X) P.loglik - 0.5*log(numel(X))*( (numel(P.mu)*(numel(P.mu)+4)) - (numel(P.mu)+1) );
BICarray = cellfun( @(y) BIC(y,X),Params);

ICLarray = cellfun( @(y) ICLCompute(X,DX,y), Params );

[~,BICchoice] = max(BICarray);
[~,ICLchoice] = max(ICLarray);


%% Reorder & Present Parameters

[~,muorder] = sort(OutParams.mu);
muorder = muorder(end:-1:1);
OutParams1 = reorderParams(OutParams,muorder);

paramMat1 = [ OutParams1.nu.',OutParams1.mu.',OutParams1.kappa.',OutParams1.ThetaValues.'];
paramMat2 = logm( bsxfun(@rdivide,OutParams1.Q,sum(OutParams1.Q,2)) );

paramMat = [paramMat1 paramMat2];


