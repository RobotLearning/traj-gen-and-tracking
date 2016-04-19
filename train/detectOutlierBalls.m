% Estimate a rebound trajectory using polynomial fitting and 
% RANSAC, a robust estimation algorithm
% Returns the indices of the outliers

function outlierIdx = detectOutlierBalls(t,ball)

loadTennisTableValues();
% preconditioning data matrix 
t = t - t(1);
n = length(t);

% check if there is a bounce
tol = 10e-2;    
% idxBallBounce = find(ball(:,3) <= table_z + ball_radius + tol);
[~,idxBallBounce] = min(ball(:,3));

bPre = ball(1:idxBallBounce,:);
bPost = ball(idxBallBounce+1:end,:);
tPre = t(1:idxBallBounce);
matPre = [ones(idxBallBounce,1),tPre,tPre.^2];    
tPost = t(idxBallBounce+1:end);
matPost = [ones(n-idxBallBounce,1),tPost,tPost.^2];

if ~isempty(idxBallBounce)

    iter = 1000;
    num = 10; %round(n/4);
    threshDist = 0.20;
    outlierRatio = 0.1;
    inlierRatio = 1 - outlierRatio;
    bestInlierNum = 0;

    % start ransac on prebounce and postbounce trajectories
    for i = 1:iter

        % Randomly select some points
        idx = randperm(n,num);         
        idxPre = idx(idx <= idxBallBounce);
        idxPost = idx(idx > idxBallBounce);      

        ballPre = ball(idxPre,:);
        timePre = t(idxPre);
        nPre = length(idxPre);
        mat1 = [ones(nPre,1), timePre, timePre.^2];
        tol = 0.001;
        betaPreBounce = pinv(mat1,tol) * ballPre;
        
        m = betaPreBounce(:,3);
        a = m(3); 
        b = m(2);
        c = m(1) - table_z - ball_radius;
        tBounce = (-b - sqrt(b^2 - 4*a*c))/(2*a);
        M = diag([CFTX,CFTY,-CRT]);
        c = betaPreBounce(1,:);
        b = betaPreBounce(2,:);
        a = betaPreBounce(3,:);
        y = [b*tBounce + c; (2*a*tBounce + b)*M - 2*a*tBounce];
        mat = [1, tBounce; 0, 1];
        betaPostBounce = [mat \ y; a];
        
        % get the errors
        errPre = bPre - matPre * betaPreBounce;
        errPost = bPost - matPost * betaPostBounce;
        % compute inliers with less than threshold
        outlierCandIdxPre = find(sqrt(diag(errPre*errPre')) > threshDist);
        outlierCandIdxPost = find(sqrt(diag(errPost*errPost')) > threshDist);
        outlierCandIdx = [outlierCandIdxPre;idxBallBounce+outlierCandIdxPost];
        inlierNum = n - length(outlierCandIdx);
        if inlierNum > bestInlierNum %inlierNum >= round(inlierRatio*n)
            bestInlierNum = inlierNum;
            outlierIdx = outlierCandIdx;
            
        end
    end
else
    error('TODO: ball doesnt seem to bounce!');
end



end