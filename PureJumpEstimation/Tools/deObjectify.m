function [ out ] = deObjectify( InitParams )

out = [InitParams.Q(:).',InitParams.nu,InitParams.mu,InitParams.kappa,InitParams.ThetaValues];

end

