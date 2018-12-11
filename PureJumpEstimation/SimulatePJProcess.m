function [ S , Theta_Chain_Ind ] = SimulatePJProcess( Params, Ndt, dt, Nsims )
%SIMULATEFUN Function that simulates pure jump process given a number of
%parameters.

% rng('shuffle');
% rng(522542);

%% Defining Global Simulation Parameters

S_0 = 0; % Initial asset midprice value

%% Defining Model Parameters
Theta_Values = Params.ThetaValues;
pi_0 = Params.nu;
dP = Params.Q;
mu_0 = Params.mu;
kappa = Params.kappa;

lambda_p = @(S,j) mu_0(j) + kappa(j).*max( 0, Theta_Values(j) - S.' ) ;
lambda_m = @(S,j) mu_0(j) + kappa(j).*max( 0, S.' - Theta_Values(j) ) ;

b = 0.01; % Tick Size

%% Initializing Simulation Arrays
S = S_0*ones(Nsims, Ndt);  % Midprice

%% Generating the Random Bits
Unif_dNP = rand(Nsims,Ndt);
Unif_dNM = rand(Nsims,Ndt);

%% Generate Hidden Markov Chain Path
Theta_Chain_Ind = HMMsample(pi_0,dP,Ndt,Nsims);

for t=2:Ndt    
    % Advancing S
    dNP = log(Unif_dNP(:,t)) > -dt*lambda_p( S(:,t-1), Theta_Chain_Ind(:,t-1) ).' ;
    dNM = log(Unif_dNM(:,t)) > -dt*lambda_m( S(:,t-1), Theta_Chain_Ind(:,t-1) ).' ;
    
    S(:,t) = S(:,t-1) + b*( dNP - dNM );
end


end

