%% Solve the minimum principle with boundary conditions directly 
% mp gives 3rd degree polynomials
% q0 and q0dot are given
% using fmincon optimizing for their parameters qf,qfdot,T

function [qf,qfdot,T] = calcOptimalPoly(robot,racket,ballTime,ballPred,q0)

dof = length(q0);
q0dot = zeros(dof,1);
% parameters are qf,qfdot,T

% initialize with a VHP method
% loadTennisTableValues();
% time2reach = 0.5; % time to reach desired point on opponents court
% % land the ball on the centre of opponents court
% ballDes(1) = 0.0;
% ballDes(2) = dist_to_table - 3*table_y/2;
% ballDes(3) = table_z + ball_radius;
% VHP = -0.5;
% [qf_est,qfdot_est,timeEst] = calcPolyAtVHP(robot,VHP,time2reach,ballDes,ballPred,ballTime,q0);
% x0 = [qf_est;qfdot_est;timeEst];

timeEst = 0.5;
x0 = [q0;q0dot;timeEst];

tic;
fprintf('Initializing optimization at T = %f.\n',timeEst);

% constraints on joint pos and vel and time
con = robot.CON;
lb = [con.q.min;con.qd.min;0.0]; % time must be positive
ub = [con.q.max;con.qd.max;1.0];

% cost function to be minimized
fun = @(x) cost([q0;q0dot],x(1:end-1),x(end));

% constraint function - ball must be intersected with desirable vel.
nonlcon = @(x) calculateRacketDev(robot,racket,x(1:end-1),x(end));

% solve with a nonlinear optimizer
%options = optimoptions('fmincon','Display','off');
options = optimoptions('fmincon','Display','final-detailed');
%options = optimoptions('fmincon','Display','iter-detailed');
x = fmincon(fun,x0,[],[],[],[],lb,ub,nonlcon,options);

% release the variables
qf = x(1:dof);
qfdot = x(dof+1:end-1);
T = x(end);

fprintf('Calculated optimal time at T = %f.\n',T);
fprintf('Optimal Control for Table tennis took %f sec.\n',toc);

end

function J = cost(Q0,Qf,T)

a = calculatePolyCoeff(Q0,Qf,T);
a3 = a(1,:);
a2 = a(2,:);
J = T * (3*T^2*(a3*a3') + 3*T*(a3*a2') + (a2*a2'));

end

% solve for the first two 3rd order polynomial coeff alpha (a)
function a = calculatePolyCoeff(Q0,Qf,tf)

    dof = length(Q0)/2;
    a = zeros(4,dof);
    q0 = Q0(1:dof);
    q0dot = Q0(dof+1:end);
    qf = Qf(1:dof);
    qfdot = Qf(dof+1:end);

    M = @(t0,tf) [t0^3 t0^2 t0^1 1;
                  3*t0^2 2*t0 1 0;
                  tf^3 tf^2 tf^1 1;
                  3*tf^2 2*tf 1 0];

    for m = 1:dof
        %q0dot is zero
        Qstrike = [q0(m); q0dot(m); qf(m); qfdot(m)]; % strike
        a(:,m) = M(0,tf) \ Qstrike;
    end
end

% boundary conditions
% racket centre at the ball 
% racket velocity negative of incoming ball vel
function [c,ceq] = calculateRacketDev(robot,racket,Qf,T)

[xf,xfd,of] = robot.calcRacketState(Qf);
rot = quat2Rot(of);
nf = rot(1:3,3);
[posDes,velDes,normalDes] = calculateDesRacket(racket,T);
c = [];
ceq = [xf - posDes(:); 
       xfd - velDes(:);
       nf - normalDes(:)];


end

% interpolate to find ball at time t
function [pos,vel,orient] = calculateDesRacket(racket,t)

state = interp1(racket.time',[racket.pos',racket.vel',racket.normal'],t);
pos = state(1:3);
vel = state(4:6);
orient = state(7:9);

end