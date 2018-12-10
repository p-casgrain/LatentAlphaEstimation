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

% Display Paths
figure(1);
clf;
for n=1:size(X,1)
    plot(X(n,:));
    ylim([-0.5,0.5]);
    title(num2str(n));
    pause(0.3);
end


