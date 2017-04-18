% Initialize an Extended Kalman Filter for the ball 
% Ball is assumed to have constant spin

function [filter,ball_func] = init_const_spin_filter(w0)

% load table parameters
loadTennisTableValues;

params.init_spin = w0;
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


ball_func = @(x,u,dt) discr_const_spin_model(x,dt,params);
dim = 6;
C = [eye(3),zeros(3)];

% very small but nonzero value for numerical stability
eps = 1e-3; 
mats.O = eps * eye(dim);
mats.C = C;
mats.M = eps * eye(3);
filter = EKF(dim,ball_func,mats);