% Estimate a rebound trajectory using polynomial fitting and 
% RANSAC, a robust estimation algorithm
% Returns the indices of the outliers

function outlierIdx = detectOutlierBalls(t,b)

loadTennisTableValues();
% preconditioning data matrix 
t = t - t(1);
n = length(t);

% check if there is a bounce
tol = 5e-2;    
idxBallBounce = find(b(:,3) <= table_z + ball_radius + tol,1);
bPre = b(1:idxBallBounce,:);
bPost = b(idxBallBounce+1:end,:);
tPre = t(1:idxBallBounce);
matPre = [ones(idxBallBounce,1),tPre,tPre.^2];    
tPost = t(idxBallBounce+1:end);
matPost = [ones(n-idxBallBounce,1),tPost,tPost.^2];

if ~isempty(idxBallBounce)

    iter = 1000;
    num = round(n/5);
    threshDist = 0.2;
    outlierRatio = 0.5;
    inlierRatio = 1 - outlierRatio;
    bestInlierNum = 0;

    % start ransac on prebounce and postbounce trajectories
    for i = 1:iter

        % Randomly select some points
        idx = randperm(n,num);         
        idxPre = idx(idx <= idxBallBounce);
        idxPost = idx(idx > idxBallBounce);      

        ballPre = b(idxPre,:);
        timePre = t(idxPre);
        nPre = length(idxPre);
        mat1 = [ones(nPre,1), timePre, timePre.^2];
        betaPreBounce = mat1 \ ballPre;
        timePost = t(idxPost);
        ballPost = b(idxPost,:);
        mat2 = [ones(num-nPre,1),timePost,timePost.^2];
        
        m = betaPreBounce(:,3);
        a = m(3);
        b = m(2);
        c = m(1) - table_z - ball_radius;
        tBounce = (-b - sqrt(b^2 - 4*a*c))/(2*a);
        
        c = betaPreBounce(1,:);
b = betaPreBounce(2,:);
a = betaPreBounce(3,:);
y = [b*tBouncePredLS + c; (2*a*tBouncePredLS + b)*M - 2*a*tBouncePredLS];
betaAfterBounce = mat \ y;
betaAfterBounce = [betaAfterBounce; a];
        
        betaPostBounce = mat2 \ ballPost;
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