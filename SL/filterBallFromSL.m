%% Preprocess and filter incoming ball from SL

[ballObs,ballTime,lastBallPos,lastBallTime] = ...
          prefilter(ballObs,ballTime,lastBallPos,lastBallTime,table);
    
if ~isempty(ballTime)
    % get number of observations
    disp('Processing data...');
    numObs = length(ballTime);

    % if suddenly there's a jump backwards reset
    tol = 1.0;
    if (size(ballRaw,2) > 2) && (abs(ballRaw(2,end) - ballRaw(2,end-1)) > tol)
        stage = WAIT;
        ballRaw = []; 
        ballFilt = []; 
        cam = [];
        firsttime = true;  
    end        

    if firsttime
        curTime = ballTime(1);
        filter.initState([ballObs(:,1); guessBallInitVel],eps);           
        ballFilt(:,1) = ballObs(:,1);
        ballRaw(:,1) = ballObs(:,1);
        cam(:,1) = ballCam(:,1);
        firsttime = false;
        numObs = numObs - 1;
        j = 1;
    end

    % keep the observed balls
    for i = 1:numObs
        ballRaw(:,j+i) = ballObs(:,1);
        cam(:,j+i) = ballCam(:,end);
    end

    % filter up to a point
    for i = 1:numObs
        dt = ballTime(i) - curTime;
        filter.linearize(dt,0);
        filter.predict(dt,0);
        filter.update(ballObs(:,i),0);
        curTime = ballTime(i);
        ballFilt(:,j+i) = filter.x(1:3);
    end 

    j = j + numObs;
end