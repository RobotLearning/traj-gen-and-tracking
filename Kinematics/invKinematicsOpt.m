% Solve inverse kinematics with optimization 
% Using fmincon to respect joint constraints

% orientation is a quaternion
function qf = invKinematicsOpt(robot,pos,orient,qest)

con = robot.CON;
lb = con.q.min;
ub = con.q.max;

% cost function to be minimized is zero
fun = @(q) 0;

% nonlinear equality constraint
% ball must be intersected 
nonlcon = @(q) calculateEndEffDev(robot,pos,orient,q);

% solve with MATLAB nonlinear optimizer
%%{
%options = optimoptions('fmincon','Display','off');
options = optimoptions('fmincon','Algorithm', 'sqp',...
                       'Display','final-detailed', ...
                       'TolX', 1e-9, ...
                       'TolCon', 1e-9, ...
                       'MaxFunEvals', 2000, ...
                       'TolFun',1e-6);
                       %'DerivativeCheck','on',...                       
%options = optimoptions('fmincon','Display','iter-detailed');
%profile on;
[qf,~,exitflag,output] = fmincon(fun,qest,[],[],[],[],lb,ub,nonlcon,options);

% % if constraint is violated run it again
% iterMax = 10;
% iter = 0;
% while output.constrviolation > 1e-2 && iter < iterMax
%     disp('Running InvKin again');
%     qest = qest + 0.1 * randn(7,1);
%     iter = iter+1;
%     [qf,~,exitflag,output] = fmincon(fun,qest,[],[],[],[],lb,ub,nonlcon,options);
% end
% %}

end

% equality constraint
% endeff must be at a certain pos and orient
function [c,ceq] = calculateEndEffDev(robot,pos,orient,q)

%condMax = 100;
[xLink,~,~,Amats] = barrettWamKinematics(q,robot.PAR);
o = rot2Quat(Amats(6,1:3,1:3));
p = xLink(6,:);
robot.calcJacobian(q);
c = cond(robot.jac) - 100;
ceq = [p(:) - pos(:); 
       o(:) - orient(:)];
end