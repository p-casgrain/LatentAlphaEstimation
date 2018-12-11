function [ OutParams ] = preCalibrate( X , numstates , dt )
%PRE_CALIBRATE This function 'pre-calibrates' the pure jump model that will
% be fully calibrated with the EM algorithm. The pre-calibration
% uses least squares linear regression and the k-means algorithm to compute
% the parameters of an approximate model, which will serve as a starting
% point for the EM algorithm
% X is an MxN array. Each X(m,:) is an independent sequence of observed
% prices. numstates is a positive integer representing the number of
% possible states for the HMM.

M = size(X,1);
N = size(X,2);

OutParams.Delta = dt;

b = abs(max(X(1,2:end)-X(1,1:end-1)));

%% Obtain most likely sequence and most likely state values
[Tind,Tvalues] = kmeans(X(:),numstates); % Tind -> Which state it's in.
                                         % Tvalues -> State means.

Tpath = Tvalues(Tind); % Get path of theta series mean.
Tpath = reshape(Tpath,N,M).'; % Match X's size conventions.

Tind = reshape(Tind,N,M).'; 


OutParams.ThetaValues = sort(Tvalues).'; % Assign sorted state values to the output

%% Obtain values for mu and kappa 
warning('off','all')

mu = 0; kappa = 0;
for m=1:M
    [a,~,~,~,stats] = regress(X(m,2:end).'-Tpath(m,1:end-1).',...
                                    X(m,1:end-1).'-Tpath(m,1:end-1).');
                                
    kappa_m = -log(a)/(b*dt); % Compute kappa for sequence m
    mu_m = stats(4)*kappa_m/(1-a^2)/b; % Compute mu for sequence m
    
    mu = ( mu*(m-1) + mu_m )/m ; % Average out mu
    kappa = ( kappa*(m-1) + kappa_m )/m ; % Average out kappa
end

warning('on','all');

OutParams.mu = mu*ones(1,numstates); OutParams.kappa = kappa*ones(1,numstates); % Assign to output

%% Obtain starting values Q
Q = zeros(numstates,numstates);
for m=1:M
    Q_m = zeros(size(Q));
    for t=2:N
        % Count and add detected switch events to Q
        Q_m(Tind(m,t),Tind(m,t-1)) = Q_m(Tind(m,t),Tind(m,t-1)) + 1;
    end
    Q = Q + Q_m; % Add these to Q
end

Q = Q ./ (sum(Q, 2) * ones(1,numstates)); % Normalize Q
OutParams.Q = Q; % Assign to output

%% Obtain starting values for nu.

nu = zeros(1,numstates);
for k=1:numstates
    nu(k) = sum(Tind(:,1)==k);
end
nu = nu/sum(nu);

OutParams.nu = nu;
OutParams.Delta = dt;

