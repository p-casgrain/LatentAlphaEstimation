function [ deriv ] = FDDeriv( handle , x , dim , eps )
%FDDERIV Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    eps = 1e-5;
end

unit = zeros(size(x));
unit(dim) = 1;

deriv = (handle(x + eps*unit) - handle(x - eps*unit))/(2*eps);

end

