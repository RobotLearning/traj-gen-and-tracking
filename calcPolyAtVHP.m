% Calculate the joint configuration and velocity as well as time
% at the virtual hitting point VHP using inverse kinematics
%
% Returning to the center of the table

function [qf,qfdot,timeAtVHP] = calcPolyAtVHP(robot,ballPred,ballTime,q0)

tic;
loadTennisTableValues();
dof = length(q0);

% define virtual hitting plane (VHP)
VHP = -0.5;
time2reach = 0.5; % time to reach desired point on opponents court
ballFull = [ballPred;ballTime];            

% land the ball on the centre of opponents court
ballDes(1) = 0.0;
ballDes(2) = dist_to_table - 3*table_y/2;
ballDes(3) = table_z + ball_radius;

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
fast = true;
ballOutVelAtVHP = calcBallVelOut3D(ballDes,ballPosAtVHP,time2reach,fast);            

% GET RACKET DESIRED VEL AND ORIENTATION AT VHP 
[racketPos,racketVel,racketNormal] = calcDesRacketState ...
               (ballPosAtVHP,ballOutVelAtVHP,ballInVelAtVHP);
% attach a suitable angular velocity to the racket
racketAngularVel = zeros(3,1);

q0dot = zeros(dof,1);
Q0 = [q0;q0dot];           

% feed to inverse kinematics to get qf
try
    % get the slide of the original racket orientation
    [~,~,o] = robot.calcRacketState(Q0);
    rotMatrix0 = quat2Rot(o);
    slide0 = rotMatrix0(1:3,2);    
    slide0 = slide0./norm(slide0,2);
    % add a comfortable slide close to original slide
    a = racketNormal;
    % project slide0 to the racket plane
    %projNormal = racketNormal*racketNormal';
    %projRacket = eye(3) - projNormal;
    %s = projRacket * slide0;    
    s = slide0 - (racketNormal'*slide0).*racketNormal;    
    assert(abs(norm(s,2)-1) < 1e-2,'slide not normalized');
    s = s./norm(s,2);
    assert(s'*a < 1e-3,'slide calculation not working!');
    n = crossProd(s,a);                
    rotMatrix = [n,s,a];
    quatRacket = rot2Quat(rotMatrix);
    % get endeffector position and quaternion
    ePos = racketPos;
    rotBack = [cos(-pi/4); -sin(-pi/4); 0; 0];
    eQuat = mult2Quat(quatRacket,rotBack);
    qf = robot.invKinematics(ePos(:),eQuat(:),q0(:));
    timeInvKin = toc;
    fprintf('InvKin took %f sec.\n',timeInvKin);
    robot.calcJacobian(qf);                
    qfdot = robot.jac \ [racketVel;racketAngularVel];
catch ME
    disp(ME.message);
    disp('InvKin problem. Not moving the robot...');                
    %disp('Virtual Hitting Point outside of workspace');
    qf = q0;
    qfdot = zeros(dof,1);
end