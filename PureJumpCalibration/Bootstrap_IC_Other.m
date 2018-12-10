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

if ~exist('X','var') || 1
    
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

% % Get NsimsX(K*3) array of parameters
% [~,paramarray]=CalibAvg(CalibratedParams,TrueParams);
% 
% % Remove Outliers - Do this by using K-Means
% StatArray = zeros(2,size(paramarray,2));
% 
% for pdim = 1:size(paramarray,2)
%     p_idx = kmeans(paramarray(:,pdim),3); % Use K-Means to split the groups
%     idx_max = mode(p_idx); % Choose group with most points in it
%     
%     sub_param = paramarray(p_idx==idx_max,pdim); % Get parameters from largest group
%     
%     StatArray(1,pdim) = mean(sub_param);
%     StatArray(2,pdim) = sqrt(var(sub_param));
%     
% end
% 
% disp(StatArray);
% 

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


