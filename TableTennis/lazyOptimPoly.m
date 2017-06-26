% Lazy player that as opposed to the centred player (calcOptimalPoly)
% does not fix the following:
%
% 1. Landing position
% 2. Landing time
% 3. Hitting position on the racket (no longer fixed to centre of racket)
% 4. Returning trajectory
%
% Returning trajectories are now 2nd order polynomials
% As opposed to 3rd order returning poly, q0 is not fixed
% 
% using fmincon optimizing for their parameters qf,qfdot,T

function [qf,qfdot,T] = lazyOptimPoly(tt,models,q0,Tret)

tic;
dof = length(q0);
q0dot = zeros(dof,1);
robot = tt.robot; % table tennis playing robot
% parameters are qf,qfdot,T

% constraints on joint pos and vel and time
con = robot.CON;
lb = [con.q.min;con.qd.min;0.0]; % time must be positive
ub = [con.q.max;con.qd.max;2.0];

% cost function to be minimized
fun = @(x) cost([q0;q0dot],x(1:end-1),x(end),robot,models);

% nonlinear striking inequality constraints
nonlcon = @(x) strikeConstraints(robot,models,x(1:end-1),x(end),[q0;q0dot]);

% solve with MATLAB nonlinear optimizer
options = optimoptions('fmincon','Algorithm', 'interior-point',...
                       'Display','off', ...
                       'TolX', 1e-6, ...
                       'TolCon', 1e-6, ...
                       'TolFun',1e-6, ...
                       'MaxFunEvals',2e3);

x0 = [q0;q0dot;1.0]; % initial guess
%x0 = init_soln(tt,models,q0);
[xopt,fval,exitflag,output] = fmincon(fun,x0,[],[],[],[],lb,ub,nonlcon,options);
%}
output

% release the variables
%profile viewer;
qf = xopt(1:dof);
qfdot = xopt(dof+1:end-1);
T = xopt(end);

fprintf('Fval (sum of squared acc) = %f.\n',fval);
fprintf('Calc. opt. time at T = %f.\n',T);
fprintf('Optim took %f sec.\n',toc);

end

% initialize optimization by running optimpoly
function x0 = init_soln(tabletennis,models,q0)


ballDes(1) = 0.0;
ballDes(2) = models.table.dist - 3*models.table.length/4;
ballDes(3) = models.table.height; 
ballTime = models.ball.time;
ballPred = models.ball.pred;
% land the ball on the centre of opponents court
time2return = 1.0; % time for robot to go back to q0 after hit
time2reach = 0.8; % time to reach desired point on opponents court    
% Compute traj here
racketDes = tabletennis.planRacket(ballDes,ballPred,ballTime,time2reach,q0);
[qf,qfdot,T] = optimPoly(tabletennis.robot,racketDes,ballPred,q0,time2return);

x0 = [qf(:);qfdot(:);T];

end

% G is the gradient
function [J,G] = cost(Q0,Qf,T,robot,models)

dof = length(Q0)/2;
q0 = Q0(1:dof); 
q0dot = Q0(dof+1:end); 
qf = Qf(1:dof); 
qfdot = Qf(dof+1:end); 
a3 = (2/(T^3))*(q0-qf) + (1/(T^2))*(q0dot + qfdot);
a2 = (3/(T^2))*(qf-q0) - (1/T)*(qfdot + 2*q0dot);

J = T * (3*T^2*(a3'*a3) + 3*T*(a3'*a2) + (a2'*a2));

% we assume that contact occurs
net_y = models.table.dist - models.table.length/2;
ballTime = models.ball.time;
ballPred = models.ball.pred;
table_z = models.table.height;
gravity = abs(models.gravity); % should be +9.8
contactModel = models.racket.contact;
[xf,xfd,of] = robot.calcRacketState(qf,qfdot);
nf = robot.calcRacketNormal(of);
[posBall,velBall] = interpBall(ballPred,ballTime,T);
racket.normal = nf(:);
racket.vel = xfd(:);
velBallOut = contactModel(velBall(:),racket);
% hack for now - since projectile motion is inaccurate
velBallOut = [0.90;0.90;0.83] .* velBallOut(:);
ball2Racket = posBall(:) - xf(:); 
distBall2RacketAlongNormal = nf(:)'*ball2Racket;
vecBallAlongRacket = ball2Racket - distBall2RacketAlongNormal*nf(:);
time2net = (net_y - posBall(2))/velBallOut(2);
distBall2TableZ = posBall(3) - table_z; 
time2land = max(velBallOut(3) + sqrt(velBallOut(3)^2 + 2*gravity*distBall2TableZ), ...
                velBallOut(3) - sqrt(velBallOut(3)^2 + 2*gravity*distBall2TableZ)) / gravity;

if ~isreal(time2net) || ~isreal(time2land)
    time2net = real(time2net); %-100;
    time2land = real(time2land); %-100;
end
            
xNet = posBall(1) + time2net * velBallOut(1);
zNet = posBall(3) + time2net * velBallOut(3) - 0.5*gravity*time2net^2;
xLand = posBall(1) + time2land * velBallOut(1);
yLand = posBall(2) + time2land * velBallOut(2);

zDesNet = models.table.height + models.net.height + 0.50;
yDesLand = models.table.dist - 3*models.table.length/4;

J_hit = 100 * vecBallAlongRacket'*vecBallAlongRacket;
J_net = 100 * ((xNet)^2 + (zNet - zDesNet)^2);
J_land = 100 * ((xLand)^2 + (yLand - yDesLand)^2);
J = J + J_hit + J_net + J_land;

end

% boundary conditions at the racket centre
% are relaxed so that the resulting ball lands on the table
% racket velocity no longer equal to a desired racket velocity
%
% Gc and Gceq are jacobians of the constraints respectively
function [c,ceq] = strikeConstraints(robot,models,Q1,T,Q0)

% enforce joint limits as inequality constraints throughout trajectory
dof = length(Q0)/2;
q0 = Q0(1:dof); 
q0dot = Q0(dof+1:end); 
q1 = Q1(1:dof); 
q1dot = Q1(dof+1:end); 

% inequality constraints
% makes sure ball touches the racket
[xf,xfd,of] = robot.calcRacketState(q1,q1dot);
nf = robot.calcRacketNormal(of);

ballTime = models.ball.time;
ballPred = models.ball.pred;
table_xmax = models.table.xmax;
wall_z = models.wall.height;
net_z = models.net.height;
net_y = models.table.dist - models.table.length/2;
table_ymax = net_y - models.table.length/2;
table_z = models.table.height;
ball_radius = models.ball.radius;
racket_radius = models.racket.radius;
gravity = abs(models.gravity); % should be +9.8
contactModel = models.racket.contact;
%flightModel = models.ball.flight;

[posBall,velBall] = interpBall(ballPred,ballTime,T);
racket.normal = nf(:);
racket.vel = xfd(:);
velBallOut = contactModel(velBall(:),racket);
% hack for now - since projectile motion is inaccurate
velBallOut = [0.90;0.80;0.83] .* velBallOut(:);
ball2Racket = posBall(:) - xf(:); 
distBall2RacketAlongNormal = nf(:)'*ball2Racket;
vecBallAlongRacket = ball2Racket - distBall2RacketAlongNormal*nf(:);

% using super primitive projectile motion model                 

time2net = (net_y - posBall(2))/velBallOut(2);
distBall2TableZ = posBall(3) - table_z; 
time2land = max(velBallOut(3) + sqrt(velBallOut(3)^2 + 2*gravity*distBall2TableZ), ...
                velBallOut(3) - sqrt(velBallOut(3)^2 + 2*gravity*distBall2TableZ)) / gravity;

if ~isreal(time2net) || ~isreal(time2land)
    time2net = real(time2net); %-100;
    time2land = real(time2land); %-100;
end

xNet = posBall(1) + time2net * velBallOut(1);
zNet = posBall(3) + time2net * velBallOut(3) - 0.5*gravity*time2net^2;
xLand = posBall(1) + time2land * velBallOut(1);
yLand = posBall(2) + time2land * velBallOut(2);
          
% trying to get racket close to ball 
% resulting interaction should also land the ball to the table
task_ineq = [vecBallAlongRacket'*vecBallAlongRacket - racket_radius^2;
             distBall2RacketAlongNormal - ball_radius;
             -distBall2RacketAlongNormal;
             -time2net;
             xNet - table_xmax;
             -xNet - table_xmax;
             zNet - wall_z;
             -zNet + net_z;
             -time2land;
             xLand - table_xmax;
             -xLand - table_xmax;
             yLand - net_y;
             -yLand + table_ymax];

joint_limit_ineq = calcJointViolations(robot,q0,q1,q0dot,q1dot,T,time2land);
         
c = [joint_limit_ineq;
     task_ineq];

ceq = [];
%ceq = [distBall2RacketAlongNormal - ball_radius];
end

function joint_limit_ineq = calcJointViolations(robot,q0,q1,q0dot,q1dot,T,Tland)

a3 = (2/(T^3))*(q0-q1) + (1/(T^2))*(q0dot + q1dot);
a2 = (3/(T^2))*(q1-q0) - (1/T)*(q1dot + 2*q0dot);
% compute the extrema
t1 = (-a2 + sqrt(a2.^2 - 3*a3.*q0dot))./(3*a3);
t2 = (-a2 - sqrt(a2.^2 - 3*a3.*q0dot))./(3*a3);
% discard if complex and clamp to [0,T]
t1 = min(max(~logical(imag(t1)) .* t1, 0),T);
t2 = min(max(~logical(imag(t2)) .* t2, 0),T);

% the maximum at the return trajectory occurs at Tland - endpoint
qmax_ret = q1 + q1dot * Tland/2;

% enforce both minima and maxima
qext = @(t) a3.*(t.^3) + a2.*(t.^2) + q0dot.*t + q0;
con = robot.CON;

joint_limit_ineq = [qext(t1) - con.q.max;
                     con.q.min - qext(t1);
                     qext(t2) - con.q.max;
                     con.q.min - qext(t2);
                     qmax_ret - con.q.max;
                     con.q.min - qmax_ret];
end

% interpolate to find ball at time t
function [pos,vel] = interpBall(ballPred,ballTime,t)

dim = 3;
ballInterpState = interp1(ballTime',ballPred',t);
pos = ballInterpState(1:dim);
vel = ballInterpState(dim+1:2*dim);

end