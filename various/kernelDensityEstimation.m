% Kernel multivariate density estimation using
% Silverman's rule of thumb to determine the widths 

function density = kernelDensityEstimation(x)

% dimension of samples
d = size(x,1);
% num of samples
n = size(x,2);

mu = sum(x,2)/n;
res = x - repmat(mu,1,n);
Sigma = res * res' / (n-1);

% Silverman's rule of thumb to determine the width matrix H
h = ((4/(d+2))^(1/(d+4)) * n^(-1/(d+4)))^2 * diag(Sigma);
H = diag(h);

% Kernel function
kern = @(x,mu) (2*pi)^(-d/2) * det(H)^(-0.5) * exp(-0.5*(x-mu)'*(H\(x-mu)));

density = @(z) formDensity(kern,z,x);

end

function out = formDensity(kern,z,xs)

n = size(xs,2);
out = 0;
for i = 1:n
     out = out + kern(z,xs(:,i));
end
out = out/n;


end