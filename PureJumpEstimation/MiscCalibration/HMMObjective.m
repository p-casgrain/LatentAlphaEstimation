function out = HMMObjective(X, DX, P, Delta)
% HMMmaximize - Implementation of Expectation-Maximization algorithm
%
%  in : InitParams = struct with all of the starting parameters
% out : OutParams = struct with all of the estimated parameters
%       l = log-likelihood of observed sequence
%
% By Philippe Casgrain, 2016
% Based on code from Aurelien Garivier, CNRS & Telecom ParisTech

% Get number of indepent paths, path length and number of possible states
M = size(X,1);
Ndt = size(X,2);
K = sqrt(numel(P)+4)-2;

% Load in initial parameters
% mu = Params.mu; kappa = Params.kappa;
% theta = Params.ThetaValues; Delta = Params.Delta;
% Q = Params.Q; nu = Params.nu;


Q = reshape(P(1:K^2),K,K);
Q = Q ./ repmat(sum(Q,2),[1,K]);
pos = K^2;

nu = P(pos + (1:K));
nu=nu./sum(nu);
pos = pos + K;

mu = P(pos + (1:K));
pos = pos + K;

kappa = P(pos + (1:K));
pos = pos + K;

theta = P(pos + (1:K));


    
%% Compute Emission Matrix With Current Parameters
g = Pois_Transition_Prob(mu,kappa,theta,Delta,X,DX);

%% Re-Estimate Q and nu by using the Baum-Welch algorithm

% N = zeros(size(Q));
% nu_m = zeros(size(nu));
c = zeros(M,Ndt); 
% post = zeros(M,K,Ndt);

for m=1:M
    g_m = squeeze(g(m,:,:));

    [~, c_m] = HMMfilter_C(nu, Q, g_m);

    % compute the posterior distribution for the current parameters
%     [phi_m, c_m] = HMMfilter_C(nu, Q, g_m);
%     beta_m = HMMsmoother_C(Q, g_m, c_m);
%     post_m = phi_m .* beta_m;
% 
%     % expectation of the number of transitions under the current parameters
%     N = N + Q.*(phi_m(:, 1:(end-1))*(beta_m(:, 2:end).*g_m(:, 2:end)./(ones(K, 1)*c_m(2:end)))');
% 
%     nu_m = nu_m + min(post_m(:,1).',1);

    c(m,:) = c_m;

%     post(m,:,:) = post_m;
end

out = sum(log(c(:)));

end