function [psi,del1] = HMMviterbi( X, DX, Params )
%HMMVITERBI Viterbi Algorithm for a sequence X.

% Load in initial parameters
mu = Params.mu; kappa = Params.kappa;
theta = Params.ThetaValues; Delta = Params.Delta;
Q = Params.Q; nu = Params.nu;

K = numel(mu);
N = size(X,2);
D = size(X,1);

if K==1
    psi=ones([D,N]);
    del1=zeros([D,N]);
else

    Q = max(Q,0);
    Q = Q ./ (sum(Q, 2) * ones(1, K));
    nu = nu/sum(nu);

    % Compute Emission Probabilities & other
    f = Pois_Transition_Prob(mu,kappa,theta,Delta,X,DX);
    f = permute(f,[1,3,2]);
    l_f = log(f);

    l_Q = repmat(log(Q.'),[1,1,D]);
    l_Q = permute(l_Q,[3,1,2]);

    % Commence Viterbi Algorithm
    del1 = ones(D,N,K);
    del2 = ones(D,N,K);
    psi = zeros(D,N);

    del1(:,1,:) = log(repmat(nu,[D,1]));
    del2(:,1,:) = 0;

    for n=2:N
        for k=1:K
           [del1(:,n,k),del2(:,n,k)] = max( del1(:,n-1,:) + l_f(:,n-1,:) + l_Q(:,k,:) ,[], 3 );
        end
    end

    [~,psi(:,N)] = max( del1(:,N,:),[],3);

    for n=(N-1):-1:1
        idxs = sub2ind([D,N,K],1:D,(n+1)*ones([1,D]),psi(:,n+1).');
        psi(:,n) = del2(idxs);
    end        
    
end

end

