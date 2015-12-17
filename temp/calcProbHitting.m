%% Calculate probability of hitting 
% for a particular ball and racket trajectory and orientation
clc; clear; close all;

n = 10;
mu = [1;-1];
d = length(mu);
xmin = [-2;-2];
xmax = [2;2];
SIGMA = [0.9, 0.4; 0.4, 0.3];

% Statistics toolbox is needed
X1 = mvnrnd(mu',SIGMA,n);
p = mvnpdf(X1,mu',SIGMA);
tic;
mvncdf(xmin',xmax',mu',SIGMA)
toc

tic;
fun = @(x,y) exp(-(1/2)*(([x;y]-mu)'*(SIGMA\([x;y]-mu)))) ./ (((2*pi)^(d/2))*(det(SIGMA)^(1/2)));
prob = quad2d(@(X,Y)arrayfun(fun,X,Y),xmin(1),xmax(1),xmin(2),xmax(2))
toc

% generate a racket trajectory
x = @(t) t;
y = @(t) t.^2 - 2*t + 5;
r = @(t) [x(t);y(t)];

dt = 0.02;
tfin = 1;
t = dt:dt:tfin;
xt = x(t);
yt = y(t);
plot(xt,yt)

racket_radius = 0.08;
T = 0.5;
xT = x(T);
yT = y(T);
% fix theta
theta = pi/4;
oT = [-sin(theta);cos(theta)];

xmin = [xT-racket_radius*abs(oT(1));xT-racket_radius*abs(oT(2))];
xmax = [xT+racket_radius*abs(oT(1));xT+racket_radius*abs(oT(2))];
probHit = quad2d(@(X,Y)arrayfun(fun,X,Y),xmin(1),xmax(1),xmin(2),xmax(2))

%% Compute the probability of hitting for a fixed racket trajectory 

% fix a racket trajectory
% fix a ball distribution (use kalman filter)
% see if simulations agree with the theoretical prob(hitting) = int(O to T)
% hitting time distribution
