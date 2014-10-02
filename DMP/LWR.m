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

function dmp = LWR(path,dmp,force)

dt = dmp.can.dt;
pat = dmp.can.pattern;
[goal,scale] = dmp.setGoal(path);
alpha = dmp.alpha_g;
beta = dmp.beta_g;

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
x = dmp.can.evolve();

% make sure x is column vector
x = x(:);

% construct weights
w = zeros(1,lenw);
for i = 1:lenw
    psi = basis(x,h(i),c(i),pat);
    if strcmp(pat,'d')
        num = x' * diag(psi) * fd(:);
        denom = x' * diag(psi) * x;
    else
        num = psi' * fd(:);
        denom = sum(psi) + 1e-10;
    end
    w(i) = num / (scale * denom);
end

force.w = w;
% update dmps
dmp.setForcing(force);

end

% basis functions are unscaled gaussians
function out = basis(x,h,c,pat)

if strcmp(pat,'d')
    out = exp(-h * (x - c).^2);
else
    out = exp(h * (cos(x - c) - 1));
end

end