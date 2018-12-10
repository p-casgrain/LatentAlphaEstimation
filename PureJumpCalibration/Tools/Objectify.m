function [ out ] = Objectify( P , loglik, Delta )
%OBJECTIFY Summary of this function goes here
%   Detailed explanation goes here

K = sqrt(numel(P)+4)-2;

% Load in initial parameters
% out.mu = Params.out.mu; out.kappa = Params.out.kappa;
% out.ThetaValues = Params.ThetaValues; Delta = Params.Delta;
% Q = Params.Q; out.nu = Params.out.nu;


out.Q = reshape(P(1:K^2),K,K);
out.Q = out.Q ./ repmat(sum(out.Q,2),[1,K]);
pos = K^2;

out.nu = P(pos + (1:K));
out.nu=out.nu./sum(out.nu);
pos = pos + K;

out.mu = P(pos + (1:K));
pos = pos + K;

out.kappa = P(pos + (1:K));
pos = pos + K;

out.ThetaValues = P(pos + (1:K));

out.Delta = Delta;

out.loglik = loglik;

end

