%% Table tennis practice using the table tennis class

clc; clear; close all;
initializeWAM;
wam2 = [];
draw = true;
% initial ball pos and vel standard deviation
std.pos = 0.0;
std.vel = 0.0;
% measurement standard deviation
std.camera = 0.0;
tt = TableTennis(wam,wam2,q0,std,draw);
tt.practice(1);