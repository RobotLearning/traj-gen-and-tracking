% Basic regression
% To learn the weights of DMPs
% Assuming fixed centers and covariances
%
% INPUTS:
%
% path - desired trajectory
% dmp  - includes the canonical system class whose phase 
%        weights the regression
%        includes the forcing term
%        forcing term is a structure which has fixed centers and covariances
% 
% OUTPUTS:
%
% forcing structure with the forcing weights learned

function force = Regr(path,dmp)

dt = dmp.can.dt;
pat = dmp.can.pattern;
force = dmp.FORCE;
[goal,~] = dmp.setGoal(path);
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

% regress on the weights
lent = length(fd);
Psi = zeros(lent,lenw);
for i = 1:lenw
    Psi(:,i) = basis(x,h(i),c(i),pat);
end
% scale the psi matrices
if strcmp(pat,'d')
    scale = x ./ sum(Psi,2); 
else
    scale = 1 ./ (sum(Psi,2) + 1e-10);
end
scale = repmat(scale,1,lenw);
Psi = Psi .* scale;
w = Psi \ fd(:);
force.w = w;

end

% basis functions are unscaled gaussians
function out = basis(x,h,c,pat)

if strcmp(pat,'d')
    out = exp(-h * (x - c).^2);
else
    out = exp(h * (cos(x - c) - 1));
end

end