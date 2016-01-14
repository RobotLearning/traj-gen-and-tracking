%% Here we state the boundary conditions to find p(T),x(T),T
% Can be used to solve the maximum principle bvp 
% with fsolve or lsqnonlin in MATLAB
% The last condition is the transversality condition to find T

function F = bvMP(x,PAR)

% Solve symbolically to compare
p = x(1:end-1);
dim = length(p)/2;
% divide as mu and nu
mu = p(1:dim);
nu = p(dim+1:end);
T = x(end);

b0 = PAR.ball.x0;
v0 = PAR.ball.v0;
g = PAR.ball.g;
q0 = PAR.robot.Q0;
robot = PAR.robot.class;

taskdim = length(PAR.ball.x0);

if isfield(PAR.ball, 'model') % if model is provided
    flightModel = PAR.ball.model;
    b = flightModel(b0,v0,T);
else
    assert(isfield(PAR.ball,'path'),'ball path not provided!');
    assert(isfield(PAR.ball,'time'),'ball times not provided!');
    % assume ball pos and vels are provided
    ballPath = PAR.ball.path;
    ballTime = PAR.ball.time;
    b = interp1(ballTime,ballPath',T,'linear','extrap')';
end
bT = b(1:2);
vT = b(3:4);
aT = [0;-g];

% since we're minimizing accelerations 3rd degree polynomials
mu0 = mu;
nu0 = -nu - mu0*T;
q = (1/6)*mu0*T^3 + (1/2)*nu0*T^2 + q0;
qdot = (1/2)*mu0*T^2 + nu0*T;

H = -1/2*(nu0'*nu0); %(1/2)*(mu0*T - nu0)'*(mu0*T + nu0); %0; %-1/2*(nu'*nu);

% construct M numerically
% M = zeros(taskdim,dim);
% h = 1e-3;
% QhPlus = repmat(q,1,dim) + h * eye(dim);
% QhMinus = repmat(q,1,dim) - h * eye(dim);
% for i = 1:dim
%     [~,xdotT1] = robot.getEndEffectorState(QhPlus(:,i),qdot);
%     [~,xdotT2] = robot.getEndEffectorState(QhMinus(:,i),qdot);
%     M(:,i) = (xdotT1 - xdotT2)/(2*h);
% end

%R = [0 1; -1 0];
[xT,xdotT] = robot.getEndEffectorState(q,qdot);
if isfield(PAR,'rotate')
    R = PAR.rotate;
    xT = R * xT;
    xdotT = R * xdotT;
end
J = R * robot.jac;
% M = R * M;
    
%bT = b0 + v0*T + (1/2)*[g;0]*T^2;
%vT = v0 + [g;0]*T;
vdes = -vT;

% Transversality conditions
%Dphi = [-vT, J];
% Dphi = [-vT, J, zeros(size(J));
%         aT, M, J];
%r = (eye(2*dim+1) - Dphi'*pinv(Dphi)')*[H;-mu;-nu];
%r = (eye(dim+1) - Dphi'*pinv(Dphi)')*[H;-mu];
% Final cost conditions
m = 100; % weighting for quadratic final cost


F = [xT - bT;
     nu - 0;
     %J * qdot - vdes;
     %H - mu'*(J\vT)];
     %min(T,0);
     H - vT'*m*(xT-bT);
     mu - J'*m*(xT-bT)];
     %r];
