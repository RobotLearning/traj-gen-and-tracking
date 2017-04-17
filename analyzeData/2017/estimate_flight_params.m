%% Estimate parameters Cdrag and gravity

function params = estimate_flight_params(balls3,spin)

loadTennisTableValues;

if iscell(balls3)
    num_trials = length(balls3);
else
    num_trials = 1;
    ballcell{1} = balls3;
    balls3 = ballcell;
end

for idx = 1:num_trials
    
    b3 = balls3{idx};
    % PRE BOUNCE ESTIMATION
    ballPreBounce = b3(:,2:4);
    timePreBounce = b3(:,1);
    timePreBounce = timePreBounce - timePreBounce(1);
    [~,idxBounce] = min(ballPreBounce(:,3));
    timePreBounce = timePreBounce(1:idxBounce);
    ballPreBounce = ballPreBounce(1:idxBounce,:);
    ballData = [timePreBounce,ballPreBounce];

    % Provide good initial ball estimates    
    % using polyfit on first 12 balls
    sampleSize = 12; 
    M = [ones(sampleSize,1),timePreBounce(1:sampleSize),timePreBounce(1:sampleSize).^2];
    Y = ballPreBounce(1:sampleSize,:);
    beta = M \ Y;
    ballInitPosEst = beta(1,:);
    ballInitVelEst = beta(2,:);
    b0 = [ballInitPosEst,ballInitVelEst];
    
    % data used for optimization is saved here
    data.pre.bInit(:,idx) = b0(:);
    data.pre.data{idx} = ballData;
end

ballPreBounce = [];
lenSet = zeros(1,num_trials);
dim = 6;
for i = 1:num_trials
    ballPreBounce = [ballPreBounce;data.pre.data{i}];
    lenSet(i) = size(data.pre.data{i},1);
    state0((i-1)*dim+1:i*dim) = data.pre.bInit(:,i); 
end

% Run nonlinear least squares to estimate C and g
if spin.flag
    ballFun = @(w0,C,g) predictNextBall(state0(:),C,g,ballPreBounce,lenSet,spin);
    fnc = @(x) ballFun(x(1:3),x(4),x(5));
    w0_init = [-50*2*pi;0;0]; %3000rpm topspin
    x0 = [w0_init',Cdrag,gravity];
else       
    ballFun = @(C,g) predictNextBall(state0(:),C,g,ballPreBounce,lenSet,spin);
    fnc = @(x) ballFun(x(1),x(2));                 
    x0 = [Cdrag,gravity];
end

options = optimoptions('lsqnonlin');
options.Display = 'final';
options.Algorithm = 'levenberg-marquardt';
options.MaxFunEvals = 1500;
[params,err] = lsqnonlin(fnc,x0,[],[],options);
