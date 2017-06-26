T0 = 0.5;
B0 = ballPred(:,1);

% solve numerically with Newton-Raphson / nlsq-optimization
x0 = [pi/4*ones(6,1);T0];
lb = [-Inf(6,1);0];

PAR.ball.x0 = B0(1:2);
PAR.ball.v0 = B0(3:4);
PAR.ball.g = g; % acceleration
%PAR.ball.model = @(b,v,T) funState([b;v],0,T);
PAR.ball.time = ballTime;
PAR.ball.path = ballPred;
PAR.robot.Q0 = q0;
PAR.robot.Qd0 = 0;
PAR.robot.class = rrr;
PAR.rotate = R;

tic;
%{
options = optimset('TolX',1e-12); % set TolX
fun = @(x) bvMP(x,PAR);
[x, resnorm, resval, exitflag, output, jacob] = newtonraphson(fun, x0, options);
fprintf('\nExitflag: %d, %s\n',exitflag, output.message) % display output message
%}
%%{
options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective',...
     'MaxFunEvals',50000,'MaxIter',5000,'TolFun',1e-20);
[x,resnorm,resval] = lsqnonlin(@(x) bvMP(x,PAR), x0, lb, [], options);
%fprintf('Residuals:\n');
%fprintf('%f\n', resval); 
%}
toc
T = x(end)

t = dt:dt:T;
N = length(t);
mu0 = x(1:3);
nu0 = -x(4:6) - mu0*T;
qStrike = (1/6)*mu0*t.^3 + (1/2)*nu0*t.^2 + repmat(q0,1,N);
qdStrike = (1/2)*mu0*t.^2 + nu0*t;

% calculate a returning trajectory with a polynomial
q0dot = zeros(3,1);

% GET 3RD DEGREE POLYNOMIALS
t2 = dt:dt:time2return;
qf = qStrike(:,end);
qfdot = qdStrike(:,end);
M = @(t0,tf) [t0^3 t0^2 t0^1 1;
              3*t0^2 2*t0 1 0;
              tf^3 tf^2 tf^1 1;
              3*tf^2 2*tf 1 0];
qReturn = zeros(3,length(t2));            
for m = 1:3
    %q0dot is zero
    Qreturn = [qf(m); qfdot(m); q0(m); q0dot(m)]; % return
    c = M(0,time2return) \ Qreturn;
    qReturn(m,:) = c(1)*t2.^3 + c(2)*t2.^2 + c(3)*t2 + c(4);
end
q = [qStrike,qReturn];
%}