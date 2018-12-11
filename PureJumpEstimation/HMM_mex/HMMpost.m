function [ post ] = HMMpost( nu, Q, g )
%HMMPOST
% Returns the posterior probability of the hidden portion of the HMM being
% in each state using the Forward-Backward algorithm

[phi, c] = HMMfilter_C(nu, Q, g);
beta = HMMsmoother_C(Q, g, c);
post = phi .* beta;

end

