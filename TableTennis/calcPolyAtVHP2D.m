% Calculate the joint configuration and velocity as well as time
% at the virtual hitting point VHP using inverse kinematics
%
% Returning to the center of the table

function [qf,qfdot,timeAtVHP] = calcPolyAtVHP2D(robot,VHP,time2reach,ballDes,ballPred,ballTime,q0)

tic;
loadTennisTableValues();
dof = length(q0);

ballFull = [ballPred;ballTime];            
fprintf('Desired landing pos: y = %.3f.\n', ballDes(1));

% interpolate at VHP
vec = [2,3,4,5];
ballAtVHP = interp1(ballFull(1,:)',ballFull(vec,:)',VHP,'spline');
timeAtVHP = ballAtVHP(end);

fprintf('Expected arrival time at VHP: %.3f.\n', timeAtVHP);

ballAtVHP = [VHP;ballAtVHP(1:3)'];
ballPosAtVHP = ballAtVHP(1:2);
ballInVelAtVHP = ballAtVHP(3:4); 

% GET DESIRED OUTGOING VELOCITY OF THE BALL AT VHP            
par.fast = true;
par.g = gravity;
par.Cdrag = Cdrag;
par.CRR = CRR;
ballOutVelAtVHP = calcBallVelOut2D(ballDes,ballPosAtVHP,time2reach,par);            

% GET RACKET DESIRED VEL AND ORIENTATION AT VHP 
[racketPos,racketVel,racketNormal] = calcDesRacketState ...
               (ballPosAtVHP,ballOutVelAtVHP,ballInVelAtVHP,par);
% attach a suitable angular velocity to the racket
racketAngularVel = zeros(1);

q0dot = zeros(dof,1);
Q0 = [q0;q0dot];           

% feed to inverse kinematics to get qf
racket.pos = racketPos;
racket.normal = racketNormal;
racket.vel = racketVel;
racket.angvel = racketAngularVel;

[qf,qfdot] = robot.invKinTableTennis(Q0,racket);