%% Approximate Normal transformation from Cauchy-like distr.

clc; clear; close all;

%% Cauchy distribution for hitting time

% n = 100000;
% mu1 = 1;
% nu1 = 2;
% s1 = 1;
% g1 = 0.2;
% b1 = mu1 + sqrt(s1)*randn(n,1);
% v1 = nu1 + sqrt(g1)*randn(n,1);
% b2 = 1;
% v2 = -1;
% eps = 0.1;
% 
% T = -b1./v1;
% Thit = T(abs(b2 + v2.*T) < eps);
% figure;
% histogram(Thit);
% probHit = length(Thit)/n
% Tmax = 1;
% probHitTillTmax = sum(Thit < Tmax)/n
% 
% fun1 = @(T,s2,g2) sqrt(T.^2./s2 + 1/g2);
% fnc = @(T) fun1(T,s1,g1);
% normal = @(x,mu,s) 1/2 .* (1 + erf((x-mu)/(sqrt(2)*s)));
% dist = @(T) normal((nu1*T + mu1)./(sqrt(s1*g1)*fnc(T)),0,1);
% tau2 = max((eps-b2)/v2,-(eps+b2)/v2);
% tau1 = min((eps-b2)/v2,-(eps+b2)/v2);
% probHit = dist(tau2) - dist(tau1)
% probHitTillTmax = dist(Tmax) - dist(tau1)

%% Cauchy-like distribution

n = 1000000;
mu1 = -2;
mu2 = 1;
nu1 = 2;
nu2 = -2;
s1 = 1;
s2 = 0.4;
g1 = 0.2;
g2 = 0.1;
b1 = mu1 + sqrt(s1)*randn(n,1);
v1 = nu1 + sqrt(g1)*randn(n,1);
b2 = mu2 + sqrt(s2)*randn(n,1);
v2 = nu2 + sqrt(g2)*randn(n,1);
eps = 0.1;

T = -b1./v1;
Thit = T(T>0 & abs(b2 + v2.*T) < eps);
figure;
histogram(Thit);
probHit = length(Thit)/n
Tmax = 1;
probHitTillTmax = sum(Thit < Tmax)/n

% fun1 = @(T,s,g) sqrt(T.^2./s + 1/g);
% fnc = @(T) fun1(T,s1,g1);
% normal = @(x,mu,s) 1/2 .* (1 + erf((x-mu)/(sqrt(2)*s)));
% dist = @(T) normal((nu1*T + mu1)./(sqrt(s1*g1)*fnc(T)),0,1);
% tau2 = @(b2) max((eps-b2)/v2,-(eps+b2)/v2);
% tau1 = @(b2) min((eps-b2)/v2,-(eps+b2)/v2);
% % probHit = dist(tau2) - dist(tau1)
% % probHitTillTmax = dist(Tmax) - dist(tau1)
% 
% 
% gauss1 = @(x) exp(-(x-mu2).^2./(2*s2))./sqrt(2*pi*s2);
% totalFnc = @(x) gauss1(x).*(tau2(x)-tau1(x));
% probHit = integral(@(X)arrayfun(totalFnc,X),-Inf,Inf)

%% Normality experiments

% clc; clear; close all;
% n = 100000;
% mu = [1;10;5];
% s = [3;0.2;1];
% x = mu(1) + s(1) * randn(n,1);
% y = mu(2) + s(2) * randn(n,1);
% w = mu(3) + s(3) * randn(n,1);
% z = w.*x./y;
% histogram(z)

%% Racket-contact experiments

clc; clear; close all;
rx = 1;
ry = 0;
x0 = 2;
y0 = -2;
vx = -2;
vy = 3;

radiusRacket = 0.2;

% we want to end up at region D
c_x = 5;
c_y = 5 - 2/3;
radiusLand = 0.5;
a = c_y - radiusLand;
b = c_y + radiusLand;

n = 100000;

% check that the means have a soln and solve the problem at T > 0
T = x0/(rx - vx) % mean hitting time
T2 = (c_y - (y0 + vy*T)) / vy%mean landing time
assert(T == y0/(ry-vy), 'no mean solution!');
assert(T > 0, 'soln should be positive!');
assert(T2 == (c_x - x0 - vx*T)/(rx - vx), 'no mean solution!');
assert(T2 > 0, 'soln should be positive!');

% add some noise to the problem
x0 = x0 + 0.1 * randn(n,1);
y0 = y0 + 0.2 * randn(n,1);
vx = vx + 0.2 * randn(n,1);
vy = vy + 0.2 * randn(n,1);

tau = x0 ./ (rx - vx);
figure;
histogram(tau)
err = ry .* tau - y0 - vy .* tau;
figure;
histogram(err)

hit = err < radiusRacket;

tau2 = (c_x - (x0 + vx.*tau)) ./ (rx - vx);
figure;
histogram(tau2);
err2 = ry.*tau + vy.*tau2 - c_y;
figure;
histogram(err2)