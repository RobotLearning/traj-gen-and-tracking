% Function that uses nonlinear least squares to estimate
% initial ball position and velocity for a table tennis ball

function [pos_est,vel_est] = estimateInitBall(time,obs)

    % load table parameters
    loadTennisTableValues;
    time = time(:);
    % using polyfit on the balls
    sampleSize = size(obs,2);
    M = [ones(sampleSize,1),time,time.^2];
    Y = obs';
    beta = pinv(M,0.01)*Y;
    ballInitPosEst = beta(1,:);
    ballInitVelEst = beta(2,:);
    ballInit = [ballInitPosEst,ballInitVelEst];
    ballData = [time(:),Y];
    
    % Run nonlinear least squares to estimate ballInit better
    x0 = ballInit(:);
    ballFun = @(b0,C,g) predictNextBall(b0,C,g,ballData,length(ballData));
    fnc = @(x) ballFun(x,Cdrag,gravity);
    options = optimoptions('lsqnonlin');
    options.Display = 'final';
    options.Algorithm = 'levenberg-marquardt';
    options.MaxFunEvals = 1500;
    [x,err] = lsqnonlin(fnc,x0,[],[],options);
    ballInit = x(1:6);
    pos_est = ballInit(1:3);
    vel_est = ballInit(4:6);

end