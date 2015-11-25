% Table Variables 
dist_to_table= -0.87; %-3.50; %-3.34;%-3.6;%-3.34%-2.74 -0.5; %0.8;  %0.2
table_height = -0.76;
table_length = 2.74; %2.76
net_height   = 0.144; %0.1525
net_overhang = 0.1525;
net_thickness = 0.01;
table_width  = 1.525; 
table_thickness = 0.056;
net_restitution = 0.05;
table_center = 0.0;

% Table Tennis Ball Variables 
ball_radius  = 0.02;
ball_mass    = 0.0027;
ball_contact_damping = 0; % filled in SimBall
ball_contact_spring = 0;  % filled in SimBall

% Table Tennis Racket Radius 
racket_radius = 0.076; %not 0.08

% Stand Variables 
stand_height  = -1.16; %-0.95;
radius_bottom = 0.1;
radius_top    = 0.02;
stand_x       = 0.02;  %-0.85
stand_y       = -0.59;
stand_z	      = 0.9;

% Floor 
floor_level   = -1.71;

% Contact parameters
CRT = 0.88; % coeff of restitution for table (i.e. for z-dir)
CFTY = 0.72; % coeff of friction on Y-dir
CFTX = 0.68; % coefficient of friction on X-dir
CRR = 0.78; % coeff of restitution for racket

% Air drag coefficients
Cdrag = 0.1414;

% gravity
gravity = -9.802;