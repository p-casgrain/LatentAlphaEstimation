function [ ICL ] = ICLCompute( X,DX,Params )
%ICLCOMPUTE This fucntion computes the ICL for a given set of parameters

% Load in initial parameters
mu = Params.mu; kappa = Params.kappa;
theta = Params.ThetaValues; Delta = Params.Delta;
Q = Params.Q; nu = Params.nu;

K = numel(mu);
N = size(X,2);
D = size(X,1);

seq = HMMviterbi(X,DX,Params);

if K==1
    placeholder_post = ones(D,N,K);
else
    placeholder_post = zeros(D,N,K);
    idx = sub2ind([D,N,K],repmat([1:D].',[1,N]),repmat(1:N,[D,1]),seq);
    placeholder_post(idx) = 1;
end

% for idx1 = 1:D
%     for idx2 = 1:N
%         placeholder_post(idx1,idx2,seq(idx1,idx2)) = 1;
%     end
% end

lik = -Pois_Lik_Objective_Merged(mu,kappa,theta,Delta,X,DX,placeholder_post);

ICL = lik - 0.5*log(N*D)*( K*(K+4) - (K+1) );


end

