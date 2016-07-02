% Initialize an Extended Kalman Filter for the ball 
% If ball is assumed to be spinning heavily, parameters are adjusted

function [filter,ball_func] = initFilterEKF(ballgun_with_spin)

% load table parameters
loadTennisTableValues;

if ballgun_with_spin
    CRT = 0.96;
    CFTY = 1.20;
    CFTX = 1.20;
    gravity = -11.06;
    Cdrag = 0.1753;
end

params.C = Cdrag;
params.g = gravity;
params.zTable = table_z;
params.yNet = dist_to_table - table_y;
params.table_length = table_length;
params.table_width = table_width;
params.radius = ball_radius;
% coeff of restitution-friction vector
params.BMAT = BMAT;
params.CFTX = CFTX;
params.CFTY = CFTY;
params.CRT = CRT;
params.ALG = 'RK4';

ball_func = @(x,u,dt) discreteBallFlightModel(x,dt,params);
% very small but nonzero value for numerical stability
eps = 1e-6; dim = 6; 
C = [eye(3),zeros(3)];
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);
filter = EKF(dim,ball_func,mats);