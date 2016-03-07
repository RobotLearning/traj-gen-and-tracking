%% Table tennis practice using the table tennis class

clc; clear; close all;
initializeWAM;
wam2 = [];
draw = true;
train = false; % train a lookup table using optimization results
lookup = false; % use lookup table instead of optimizing online
% initial ball pos and vel standard deviation
std.pos = 0.1;
std.vel = 0.1;
% measurement standard deviation 
std.camera = 0.0;
tt = TableTennis(wam,wam2,q0,std,draw,train,lookup);
tt.practice(q0,1);

% Things to do for simulation
%
%{
add a proper net
add the lines for the table
add table legs
add robot top?
can we find better sim for Barrett WAM?

can we add 4-link robot sim?

put virtual hitting plane for the 2nd method

%}