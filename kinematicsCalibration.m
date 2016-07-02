%% Kinematics calibration
clc; clear; close all;

% first get some joint positions that correspond to table height
initializeWAM;
loadTennisTableValues;
Q0 = zeros(14,1);

% number of goto tasks
N = 1;
xdes = zeros(3,N); 
Qdes = zeros(7,N);
% fix the desired normal
racket.vel = zeros(3,1);
racket.normal = [-1; 0; 0]; %[0;-1/sqrt(2);1/sqrt(2)]; % 45 degrees facing up
racket.angvel = zeros(3,1);

for i = 1:N
    xdes(1,i) = 0.0; %-0.5 + (i-1)*0.5;
    xdes(2,i) = -0.80; %dist_to_table + 0.10 + racket_radius; % yanlongs table position
    xdes(3,i) = floor_level + table_height;    
    racket.pos = xdes(:,i);
    
    [qf,~] = wam.invKinTableTennis(Q0,racket);
    Qdes(:,i) = qf;
    
end
    