%% Table tennis practice using the table tennis class

clc; clear; close all;
initializeWAM;
wam2 = [];
draw = false;
train = true; % train a lookup table using optimization results
lookup = false; % use lookup table instead of optimizing online
% initial ball pos and vel standard deviation
std.pos = 0.1;
std.vel = 0.1;
% measurement standard deviation
std.camera = 0.0;
tt = TableTennis(wam,wam2,q0,std,draw,train,lookup);
tt.practice(q0,5000);