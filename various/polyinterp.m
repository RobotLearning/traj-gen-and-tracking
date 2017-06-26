% Polynomial interpolation using function values and derivatives
% Using Vandermonde matrix
%
% Returns the coefficients of the polynomial
%
% Compared with hermite function and they match

function p = polyinterp(x,y,yprime)

% build the polynomial matrix
% assuming 2n-1 order polynomial (goes through each point exactly)

% number of data points
n = length(x);
M = zeros(2*n);
for i = 2*n-1:-1:0
    M(1:n,2*n-i) = x.^i;
    M(n+1:2*n,2*n-i) = i*x.^(i-1);    
end

vals = [y(:);yprime(:)];
p = M \ vals;

end