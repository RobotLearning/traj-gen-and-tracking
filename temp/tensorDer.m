%% Testing tensor derivative

clc; clear; close all;

m = 2; % dim of w
n = 3; % dim of v
l = 3; % dim of x

M = @(x) [x(3)*x(1)^2, x(3)^2*x(2)^3, 2*x(1)+x(3);
          x(1)*x(2), x(3)+x(2)^2, x(1)+x(2)+x(3)];
v = rand(n,1); 
x = rand(l,1);
w = M(x) * v;

% test derivative
dx = 1e-6;
X = repmat(x,1,l) + dx * eye(l);

for i = 1:l
    W(:,i) = M(X(:,i))*v;
end

dwdx = (W - repmat(w,1,l))/dx

% test analytical derivative
dm1dx = @(x) [2*x(1)*x(3), 0, 2;
              0, 3*x(3)^2*x(2)^2, 0;
              x(1)^2, 2*x(3)*x(2)^3, 1];
dm2dx = @(x) [x(2), 0, 1;
              x(1), 2*x(2), 1;
              0, 1, 1];
dMdx = blkdiag(dm1dx(x),dm2dx(x));

dwdx2 = dMdx * repmat(v,m,1);
dwdx2 = reshape(dwdx2',l,m)'

% trying kronecker product + automatic differentiation
% x_AD = adiff(x);
% func = vectorize(M,x_AD);
% [y,dydx] = adiffget(func);
