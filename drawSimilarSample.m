% Draw similar samples
% Fitting squared exponentials with fixed width on each sample point
% According to a rule of thumb (check Wikipedia kernel density est.)

function sample = drawSimilarSample(x)

% dimension of samples
d = size(x,1);
% num of samples
n = size(x,2);

mu = sum(x,2)/n;
res = x - repmat(mu,1,n);
Sigma = res * res' / (n-1);

% Silverman's rule of thumb to determine the widths
s = ((4/(d+2))^(1/(d+4)) * n^(-1/(d+4))) * sqrt(diag(Sigma));
%s = zeros(6,1);
% H = diag(h);

% Kernel function
%kern = @(x,mu) (2*pi)^(-d/2) * det(H)^(-0.5) * exp(-0.5*(x-mu)'*(H\(x-mu)));

a = randi(n);
mu_samp = x(:,a);
s_samp = s;

sample = mu_samp + diag(s_samp) * randn(d,1);
