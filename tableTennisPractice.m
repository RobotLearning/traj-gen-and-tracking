%% Table tennis practice using the table tennis class

clc; clear; close all;
initializeWAM;
wam2 = [];
OPT.draw = true; % draw the simulation
OPT.record = true; % record the simulation
OPT.train = false; % train a lookup table using optimization results
OPT.lookup = false; % use lookup table instead of optimizing online
OPT.vhp = true; % use vhp strategy
% initial ball pos and vel standard deviation
std.pos = 0.1;
std.vel = 0.1;
% measurement standard deviation 
std.camera = 0.0;
tt = TableTennis(wam,wam2,q0,std,OPT);
tt.practice(q0,5);

% Things to do for simulation
%
%{
add robot top?
can we find better sim for Barrett WAM?

can we add 4-link robot sim?

put virtual hitting plane for the 2nd method

%}