%% Table Tennis parameters and values

% Table Variables 
if ~exist('dist_to_table')
    dist_to_table = -1.15; % -0.80;
end
table_height = 0.76;
table_length = 2.76; 
net_height   = 0.144; %0.1525
net_overhang = 0.12; %0.1525;
net_thickness = 0.01;
table_width  = 1.525; 
table_thickness = 0.02; %0.056;
net_restitution = 0.05;
table_center = 0.0;

% Table Tennis Ball Variables 
ball_radius  = 0.02;
ball_mass    = 0.0027;
ball_contact_damping = 0; 
ball_contact_spring = 0; 

% Table Tennis Racket Radius 
racket_radius = 0.076; %not 0.08

% Stand Variables 
stand_radius = 0.02;
stand_x       = -0.45;
stand_y       = -0.59;

% Floor 
floor_level   = -1.71;

% virtual hitting plane y level
vhp = -0.6;

%% Rebound and contact model parameters

% Coefficients for rebound and contact model 
% fitted for dataset realBallData1

CRT = 0.96;
CFTY = 1.20;
CFTX = 1.20;

% old coefficients from SL
% CRT = 0.86;
% CFTY = 0.72;
% CFTX = 0.68;

% nondiagonal matrix
BMAT = [1.23 0.33 0.36;
        0.10 0.65 -0.71; 
        -0.19 0.62 -0.29];


% coeff of restitution for racket
CRR = 0.78; 

%% Flight model parameters

% flight model estimated using old dataset1
% gravity = -10.2540;
% Cdrag = 0.1442; % Air drag coefficient

% estimated using new dataset1
gravity = -11.06;
Cdrag = 0.1753;

% post bounce
Cdrag_post = 0.1968;
gravity_post = -10.83;

% old coefficients from SL
% gravity = -9.801;
% Cdrag = 0.1414;

%% Values for drawing the table

% for the second demo only
%dist_to_table = dist_to_table - 0.25;

table_z = floor_level + table_height;
table_x = table_center + table_width/2;
table_y = table_length/2;

T1 = [table_center - table_x; 
    dist_to_table - table_length; 
    table_z];
T2 = [table_center + table_x;
    dist_to_table - table_length;
    table_z];
T3 = [table_center + table_x;
    dist_to_table;
    table_z];
T4 = [table_center - table_x;
    dist_to_table;
    table_z];
T5 = [table_center - table_x; 
    dist_to_table - table_length; 
    table_z - table_thickness];
T6 = [table_center + table_x;
    dist_to_table - table_length;
    table_z - table_thickness];
T7 = [table_center + table_x;
    dist_to_table;
    table_z - table_thickness];
T8 = [table_center - table_x;
    dist_to_table;
    table_z - table_thickness];
T = [T1,T2,T3,T4,T5,T6,T7,T8]';

net1 = [table_center - table_x;
        dist_to_table - table_y;
        table_z];
net2 = [table_center + table_x;
        dist_to_table - table_y;
        table_z];
net3 = [table_center + table_x;
        dist_to_table - table_y;
        table_z + net_height];
net4 = [table_center - table_x;
        dist_to_table - table_y;
        table_z + net_height];
    
loadStandValues;

net = [net1,net2,net3,net4]';

%% Ballgun values

% initialize the ball gun
ball_cannon(1) = table_center + 0.4;
ball_cannon(2) = dist_to_table - table_length - 0.2; % since dist_to_table is negative
ball_cannon(3) = floor_level + table_height + 0.15;