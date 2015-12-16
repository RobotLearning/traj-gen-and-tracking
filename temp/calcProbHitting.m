%% Calculate probability of hitting 
% for a particular ball and racket trajectory and orientation
clc; clear; close all;

n = 10;
mu = [1, -1];
SIGMA = [0.9, 0.4; 0.4, 0.3];
X1 = mvnrnd(mu, SIGMA, n);
p = mvnpdf(X1,mu,SIGMA);
xl = [-2,-2];
xu = [2,2];
mvncdf(xl,xu,mu,SIGMA)

mu = [1;-1];
d = length(mu);
X2 = repmat(mu,1,n) + chol(SIGMA)*randn(2,n);
fun = @(x,y) exp(-(1/2).*([x;y]-mu)'*(SIGMA\([x;y]-mu))) ./ ((2*pi*abs(det(SIGMA)))^(d/2));
xmin = [-2;-2];
xmax = [2;2];
prob = quad2d(fun,xmin(1),xmax(1),xmin(2),xmax(2))