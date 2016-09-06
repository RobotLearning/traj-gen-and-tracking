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

function [qf,qfdot,T,Tland] = lazyPlayer(robot,ballPred,q0)

tic;
dof = length(q0);
q0dot = zeros(dof,1);
% parameters are qf,qfdot,T

x0 = racket.est;

% constraints on joint pos and vel and time
con = robot.CON;
lb = [con.q.min;con.qd.min;0.0]; % time must be positive
ub = [con.q.max;con.qd.max;1.0];

% cost function to be minimized
fun = @(x) cost([q0;q0dot],x(1:end-1),x(end));

% nonlinear equality constraint
% ball must be intersected with desirable vel.
nonlcon = @(x) calculateRacketDev(robot,racket,ballPred,x(1:end-1),...
                                  x(end),[q0;q0dot],Tret);

% solve with MATLAB nonlinear optimizer
%%{
%options = optimoptions('fmincon','Display','off');
options = optimoptions('fmincon','Algorithm', 'sqp',...
                       'Display','off', ...
                       'TolX', 1e-6, ...
                       'TolCon', 1e-6, ...
                       'TolFun',1e-6, ...
                       'MaxFunEvals',2e3,...
                       'GradObj','on');...
                       ... %'DerivativeCheck','on',...
                       ... %'GradConstr','on');
%options = optimoptions('fmincon','Display','iter-detailed');
%profile on;
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

% G is the gradient
function [J,G] = cost(Q0,Qf,T)

% a = calculatePolyCoeff(Q0,Qf,T);
% a3 = a(1,:)';
% a2 = a(2,:)';

dof = length(Q0)/2;
q0 = Q0(1:dof); 
q0dot = Q0(dof+1:end); 
qf = Qf(1:dof); 
qfdot = Qf(dof+1:end); 
a3 = (2/(T^3))*(q0-qf) + (1/(T^2))*(q0dot + qfdot);
a2 = (3/(T^2))*(qf-q0) - (1/T)*(qfdot + 2*q0dot);

J = T * (3*T^2*(a3'*a3) + 3*T*(a3'*a2) + (a2'*a2));

% supply gradient if algorithm supports it
if nargout > 1
    G = [(6/(T^3))*(qf-q0) - (3/T^2)*(q0dot + qfdot);
         (-3/T^2)*(qf-q0) + (1/T)*(2*qfdot + q0dot)];
    G(end+1) = (-9/(T^4))*((qf-q0)'*(qf-q0)) + ...
               (6/(T^3))*((qf-q0)'*(q0dot+qfdot)) +...
               (3/(T^2))*(q0dot'*(q0dot + qfdot)) +...
               -(1/T^2)*((qfdot+2*q0dot)'*(qfdot+2*q0dot));
end

end

% boundary conditions
% racket centre at the ball 
% racket velocity equal to desired racket velocity given in racket struct
%
% Gc and Gceq are jacobians of the constraints respectively
function [c,c_eq,Gc,Gc_eq] = calculateRacketDev(robot,racket,ballPred,Qf,T,Q0,Tret)

dim = size(racket.pos,1);
% enforce joint limits as inequality constraints throughout trajectory
dof = length(Q0)/2;
q0 = Q0(1:dof); 
q0dot = Q0(dof+1:end); 
qf = Qf(1:dof); 
qfdot = Qf(dof+1:end); 
a3 = (2/(T^3))*(q0-qf) + (1/(T^2))*(q0dot + qfdot);
a2 = (3/(T^2))*(qf-q0) - (1/T)*(qfdot + 2*q0dot);
% compute the extrema
t1 = (-a2 + sqrt(a2.^2 - 3*a3.*q0dot))./(3*a3);
t2 = (-a2 - sqrt(a2.^2 - 3*a3.*q0dot))./(3*a3);
% discard if complex and clamp to [0,T]
t1 = min(max(~logical(imag(t1)) .* t1, 0),T);
t2 = min(max(~logical(imag(t2)) .* t2, 0),T);

% do the same for the return trajectory
a3ret = (2/(Tret^3))*(qf-q0) + (1/(Tret^2))*(qfdot + q0dot);
a2ret = (3/(Tret^2))*(q0-qf) - (1/Tret)*(q0dot + 2*qfdot);
% compute the extrema
t1ret = (-a2ret + sqrt(a2ret.^2 - 3*a3ret.*qfdot))./(3*a3ret);
t2ret = (-a2ret - sqrt(a2ret.^2 - 3*a3ret.*qfdot))./(3*a3ret);
% discard if complex and clamp to [0,Treturn]
t1ret = min(max(~logical(imag(t1ret)) .* t1ret, 0),Tret);
t2ret = min(max(~logical(imag(t2ret)) .* t2ret, 0),Tret);

% enforce both minima and maxima
qext = @(t) a3.*(t.^3) + a2.*(t.^2) + q0dot.*t + q0;
qext_ret = @(t) a3ret.*(t.^3) + a2ret.*(t.^2) + qfdot.*t + qf;
con = robot.CON;

% equality constraints
[xf,xfd,of] = robot.calcRacketState(qf,qfdot);
nf = robot.calcRacketNormal(of);
[posDes,velDes,normalDes] = calculateDesRacket(racket,T);
% vecFromBallToRacket = rot'*(posDes(:) - xf); % in racket coordinates

% inequality constraints
% last one makes sure ball touches the racket
c = [qext(t1) - con.q.max;
     con.q.min - qext(t1);
     qext(t2) - con.q.max;
     con.q.min - qext(t2);
     qext_ret(t1ret) - con.q.max;
     con.q.min - qext_ret(t1ret);
     qext_ret(t2ret) - con.q.max;
     con.q.min - qext_ret(t2ret)];
     %vecFromBallToRacket(1:2)'*vecFromBallToRacket(1:2) - racket.radius^2];
     
% c = [];

% first one makes sure ball is at racket dist 
% 2cm is ball radius is not considered
% ceq =   [vecFromBallToRacket(3) - 0.02;
%         xfd - velDes(:);
%         nf - normalDes(:)];

% trying to get racket centre to ball 
c_eq = [xf - posDes(:); 
      xfd - velDes(:);
      nf - normalDes(:)];
      
      
end

% interpolate to find ball at time t
function [pos,vel,normal] = calculateDesRacket(racket,t)

dim = size(racket.pos,1);
state = interp1(racket.time',[racket.pos',racket.vel',racket.normal'],t);
pos = state(1:dim);
vel = state(dim+1:2*dim);
normal = state(2*dim+1:3*dim);

end