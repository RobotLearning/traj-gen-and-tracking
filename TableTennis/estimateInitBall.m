% Function that uses nonlinear least squares to estimate
% initial ball position and velocity for a table tennis ball

function ballInit = estimateInitBall(time,obs,spin)

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
    ballFun = @(b0,C,g) predictNextBall(b0,C,g,ballData,length(ballData),spin);
    fnc = @(x) ballFun(x,Cdrag,gravity);
    options = optimoptions('lsqnonlin');
    options.Display = 'final';
    options.Algorithm = 'levenberg-marquardt';
    options.MaxFunEvals = 1500;
    
    if spin.flag
        p0 = ballInit(1:3); 
        v0 = ballInit(4:6);
        phi0 = zeros(3,1);
        w0 = spin.est;
        x0 = [p0(:);phi0;v0(:);w0];
        [x,resnormsq] = lsqnonlin(fnc,x0,[],[],options);
        ballInit = x(1:12);        
    else
        x0 = ballInit(:);
        [x,resnormsq] = lsqnonlin(fnc,x0,[],[],options);
        ballInit = x(1:6);
    end
    fprintf('Res norm: %f\n', sqrt(resnormsq));
end