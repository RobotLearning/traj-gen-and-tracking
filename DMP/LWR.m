% Locally weighted regression
% To learn the weights of DMPs
% Assuming fixed centers and covariances
%
% INPUTS:
%
% path - desired trajectory
% can  - canonical system class whose phase weights the regression
% force - structure which has fixed centers and covariances
% 
% OUTPUTS:
%
% force - appends a weight field to force structure

function force = LWR(path,can,alpha,beta,force)

dt = can.dt;
goal = path(end);
% TODO: interpolate over trajectory
y_des = path;
yd_des = [0,diff(path)]/dt;
ydd_des = [0,diff(yd_des)]/dt;
% calculate fd
fd = ydd_des - alpha * (beta * (goal - y_des) - yd_des);

h = force.h;
c = force.c;
% number of weights to regress 
lenw = length(force.c);
x = can.evolve();
% spatial scaling
scale = goal - path(1);

% make sure x is column vector
x = x(:);

% construct weights
w = zeros(1,lenw);
for i = 1:lenw
    phi = basis(x,h(i),c(i));
    num = x' * diag(phi) * fd(:);
    denom = x' * diag(phi) * x;
    w(i) = num / (scale * denom);
end

force.w = w;
end

% basis functions are unscaled gaussians
function out = basis(x,h,c)
out = exp(-h * (x - c).^2);
end