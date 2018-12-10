function [ out ] = relTol( tol, X, Y )
% Binary Funtion. Evaluates to 1 if the relative difference
% between numeric vectors of size N, X & Y, exceed the
% tolerance level (positive scalar tol).

out = norm( (X - Y) ./ ( min(abs(X), abs(Y)) + 1e-8 ) ) > tol;

end

