% Calculate the joint configuration and velocity as well as time
% at the virtual hitting point VHP using inverse kinematics
%
% Returning to the center of the table

function [qf,qfdot,timeAtVHP] = calcPolyAtVHP(robot,VHP,time2reach,ballDes,ballPred,ballTime,q0)

tic;
loadTennisTableValues();
dof = length(q0);

ballFull = [ballPred;ballTime];
fprintf('Desired landing pos: x = %.3f, y = %.3f.\n', ballDes(1), ballDes(2));

% interpolate at VHP
vec = [1,3,4,5,6,7];
ballAtVHP = interp1(ballFull(2,:)',ballFull(vec,:)',VHP,'spline');
timeAtVHP = ballAtVHP(end);

fprintf('Expected arrival time at VHP: %.3f.\n', timeAtVHP);

ballAtVHP = [ballAtVHP(1);VHP;ballAtVHP(2:5)'];
ballPosAtVHP = ballAtVHP(1:3);
ballInVelAtVHP = ballAtVHP(4:6); 

% GET DESIRED OUTGOING VELOCITY OF THE BALL AT VHP            
par.fast = true;
par.g = gravity;
par.Cdrag = Cdrag;
par.CRR = CRR;
ballOutVelAtVHP = calcBallVelOut3D(ballDes,ballPosAtVHP,time2reach,par);            

% GET RACKET DESIRED VEL AND ORIENTATION AT VHP 
[racketPos,racketVel,racketNormal] = calcDesRacketState ...
               (ballPosAtVHP,ballOutVelAtVHP,ballInVelAtVHP,par);
% attach a suitable angular velocity to the racket
racketAngularVel = zeros(3,1);

q0dot = zeros(dof,1);
Q0 = [q0;q0dot];           

% feed to inverse kinematics to get qf
racket.pos = racketPos;
racket.normal = racketNormal;
racket.vel = racketVel;
racket.angvel = racketAngularVel;

[qf,qfdot] = robot.invKinTableTennis(Q0,racket);