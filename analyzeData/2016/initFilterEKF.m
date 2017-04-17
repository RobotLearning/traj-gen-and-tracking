% Initialize an Extended Kalman Filter for the ball 
% If ball is assumed to be spinning heavily, parameters are adjusted

function [filter,ball_func] = initFilterEKF(spin)

% load table parameters
loadTennisTableValues;

if ~spin.flag
    % try to compensate for effects of spin with a nonspinning ball model
    gravity = -11.06;
    Cdrag = 0.1753;
    CRT = 0.96;
    CFTY = 1.20;
    CFTX = 1.20;
end

params.Clift = Clift;
params.mu = mu; % is this accurate?
params.C = Cdrag;
params.g = gravity;
params.zTable = table_z;
params.yNet = dist_to_table - table_y;
params.table_length = table_length;
params.table_width = table_width;
params.radius = ball_radius;
% coeff of restitution-friction vector
%params.BMAT = BMAT;
params.CFTX = CFTX;
params.CFTY = CFTY;
params.CRT = CRT;
params.ALG = 'RK4';

if spin.flag
    ball_func = @(x,u,dt) discreteBallSpinModel(x,dt,params);
    dim = 12;
    C = [eye(3),zeros(3,9)];
else
    ball_func = @(x,u,dt) discreteBallFlightModel(x,dt,params);
    dim = 6;
    C = [eye(3),zeros(3)];
end

% very small but nonzero value for numerical stability
eps = 1e-3; 
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);
filter = EKF(dim,ball_func,mats);