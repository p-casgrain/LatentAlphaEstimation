function [ OutParams ] = onestateMLE(X, DX,InitParams)
%1STATEMLE Computes the MLE estimates of the mean-reverting pure-jump model
% with the assumption that there are no hidden states

Delta = InitParams.Delta;
paramarray = [InitParams.mu,InitParams.kappa,InitParams.theta]; % mu, kappa, 

loglik_handle = @(y) Pois_Lik_Objective_Merged(y(1),y(2),y(3),Delta,X,DX,ones(size(X)));

options = optimoptions( 'fminunc',...
    'SpecifyObjectiveGradient',true,...
    'CheckGradients',false,...
    'Display','off',...
    'Algorithm','trust-region',...
    'HessianFcn','objective',...
    'StepTolerance',1e-10,...
    'OptimalityTolerance',1e-10);

[paramarray,loglik] = fminunc(loglik_handle,paramarray,options);

OutParams.Delta = InitParams.Delta;
OutParams.mu = paramarray(1);
OutParams.kappa = paramarray(2);
OutParams.theta = paramarray(3);
OutParams.loglik = -loglik; 

end

