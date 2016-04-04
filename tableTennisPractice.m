%% Table tennis practice using the table tennis class

clc; clear; close all;
initializeWAM;
wam2 = [];
OPT.draw = false; % draw the simulation
OPT.record = false; % record the simulation
OPT.train = false; % train a lookup table using optimization results
OPT.lookup = false; % use lookup table instead of optimizing online
OPT.vhp = false; % use vhp strategy
% initial ball pos and vel standard deviation
STD.pos = 0.1;
STD.vel = 0.1;
% measurement standard deviation 
STD.camera = 0.0;
tt = TableTennis(wam,wam2,q0,STD,OPT);
tt.practice(q0,1);