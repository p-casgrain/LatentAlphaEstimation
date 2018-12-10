function [ OutParams ] = reorderParams( Params, order )
%REORDERPARAMS This function re-orders the parameter vectors in the struct
%Params according to the vector order

% Load in initial parameters
mu = Params.mu; kappa = Params.kappa;
theta = Params.ThetaValues; Delta = Params.Delta;
Q = Params.Q; nu = Params.nu;

K = numel(nu);

if numel(order)~=K
    error('Order Vector is not of length K');
elseif numel(unique(order))~=K
    error('Order Vector does not contain all possible indices');
elseif any(order<=0) || any(order>K)
    error('Index out of range');
end

mu = mu(order); kappa = kappa(order); theta = theta(order);
nu = nu(order); Q = Q(order,order);

OutParams = Params;
OutParams.mu = mu; OutParams.kappa = kappa; OutParams.ThetaValues = theta;
OutParams.nu = nu; OutParams.Q = Q;

end

