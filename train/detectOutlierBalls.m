% Estimate a rebound trajectory using polynomial fitting and 
% RANSAC, a robust estimation algorithm
% Returns the indices of the outliers

function outlierIdx = detectOutlierBalls(t,ball,camId)

    % preconditioning data matrix 
    t = t - t(1);
    n = length(t);

    % check if there is a bounce
    idxBounce = findReboundIndex(ball);

    bPre = ball(1:idxBounce,:);
    bPost = ball(idxBounce+1:end,:);
    tPre = t(1:idxBounce);
    matPre = [ones(idxBounce,1),tPre,tPre.^2];    
    tPost = t(idxBounce+1:end);
    matPost = [ones(n-idxBounce,1),tPost,tPost.^2];

    if ~isempty(idxBounce)

        iter = 1000;
        num = 10; %round(n/4);
        threshDist = 0.1;
        %outlierRatio = 0.1;
        %inlierRatio = 1 - outlierRatio;
        bestInlierNum = 0;

        % start ransac on prebounce and postbounce trajectories
        for i = 1:iter

            % Randomly select some points
            idx = randperm(n,num);         
            idxPre = idx(idx <= idxBounce);
            idxPost = idx(idx > idxBounce); 

            %[betaPreBounce,betaPostBounce] = fitPoly(t,ball,idxPre,idxPost,camId);
            [betaPreBounce,betaPostBounce] = fitPolyWithGravity(t,ball,idxPre,idxPost,camId);
            
            % debug mode
            %{
            scatter3(ball(:,1),ball(:,2),ball(:,3))
            hold on;
            bnewPre = matPre * betaPreBounce;
            bnewPost = matPost * betaPostBounce;
            bnew = [bnewPre;bnewPost];
            scatter3(bnew(:,1),bnew(:,2),bnew(:,3))
            %}

            % get the errors
            errPre = bPre - matPre * betaPreBounce;
            errPost = bPost - matPost * betaPostBounce;
            % compute inliers with less than threshold
            outlierCandIdxPre = find(sqrt(diag(errPre*errPre')) > threshDist);
            outlierCandIdxPost = find(sqrt(diag(errPost*errPost')) > threshDist);
            outlierCandIdx = [outlierCandIdxPre;idxBounce+outlierCandIdxPost];
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

% function to find the index of the balldata closest to rebound in a robust way
function idx = findReboundIndex(ball,table_z,ball_radius)

    loadTennisTableValues();
    % check if there is a bounce
    tol = 5e-2;    
    idxBallBounce = [];
    while isempty(idxBallBounce)
        idxBallBounce = find(ball(:,3) <= table_z + ball_radius + tol, 1);
        tol = 1.2 * tol;
    end
    idx = idxBallBounce(1);

    
end

% Fit 2nd order polynomials in x-y-z
function [betaPre,betaPost] = fitPoly(t,ball,idxPre,idxPost,camId)

    loadTennisTableValues();
    M = diag([CFTX,CFTY,-CRT]);

    % there should be more balls after bounce
    if camId == 1
        ballPost = ball(idxPost,:);
        timePost = t(idxPost);
        nPost = length(idxPost);
        mat1 = [ones(nPost,1), timePost, timePost.^2];
        tol = 0.001;
        betaPost = pinv(mat1,tol) * ballPost;
        m = betaPost(:,3);
        a = m(3); b = m(2); c = m(1) - table_z - ball_radius;
        tBounce = (-b + sqrt(b^2 - 4*a*c))/(2*a);
        c = betaPost(1,:);
        b = betaPost(2,:);
        a = betaPost(3,:);
        y = [b*tBounce + c; (2*a*tBounce + b)/M - 2*a*tBounce];
        mat = [1, tBounce; 0, 1];
        betaPre = [mat \ y; a];
    end

    if camId == 3
        ballPre = ball(idxPre,:);
        timePre = t(idxPre);
        nPre = length(idxPre);
        mat1 = [ones(nPre,1), timePre, timePre.^2];
        tol = 0.001;
        betaPre = pinv(mat1,tol) * ballPre;            
        m = betaPre(:,3);
        a = m(3); b = m(2); c = m(1) - table_z - ball_radius;
        tBounce = (-b - sqrt(b^2 - 4*a*c))/(2*a);

        c = betaPre(1,:);
        b = betaPre(2,:);
        a = betaPre(3,:);
        y = [b*tBounce + c; (2*a*tBounce + b)*M - 2*a*tBounce];
        mat = [1, tBounce; 0, 1];
        betaPost = [mat \ y; a];
    end       

end

% Fit 2nd order polynomials in x-y
% Fit line in z with gravity constant
function [betaPre,betaPost] = fitPolyWithGravity(t,ball,idxPre,idxPost,camId)

    loadTennisTableValues();
    M = diag([CFTX,CFTY,-CRT]);

    % there should be more balls after bounce
    if camId == 1
        ballPost = ball(idxPost,:);
        timePost = t(idxPost);
        nPost = length(idxPost);
        mat1 = [ones(nPost,1), timePost, timePost.^2];
        mat2 = [ones(nPost,1), timePost];
        tol = 0.001;        
        ballPost(:,3) = ballPost(:,3) - 0.5*gravity*timePost.^2;
        betaPost(:,1:2) = pinv(mat1,tol) * ballPost(:,1:2);
        betaPost(1:2,3) = pinv(mat2,tol) * ballPost(:,3);
        betaPost(3,3) = 0.5*gravity;
        
        m = betaPost(:,3);
        a = m(3); b = m(2); c = m(1) - table_z - ball_radius;
        tBounce = (-b + sqrt(b^2 - 4*a*c))/(2*a);
        c = betaPost(1,:);
        b = betaPost(2,:);
        a = betaPost(3,:);
        y = [b*tBounce + c; (2*a*tBounce + b)/M - 2*a*tBounce];
        mat = [1, tBounce; 0, 1];
        betaPre = [mat \ y; a];
    end

    if camId == 3
        ballPre = ball(idxPre,:);
        timePre = t(idxPre);
        nPre = length(idxPre);
        mat1 = [ones(nPre,1), timePre, timePre.^2];
        mat2 = [ones(nPre,1),timePre];
        tol = 0.001;
        
        ballPre(:,3) = ballPre(:,3) - 0.5*g*timePre.^2;
        betaPre(:,1:2) = pinv(mat1,tol) * ballPre(:,1:2);
        betaPre(1:2,3) = pinv(mat2,tol) * ballPre(:,3);
        betaPre(3,3) = g;
             
        m = betaPre(:,3);
        a = m(3); b = m(2); c = m(1) - table_z - ball_radius;
        tBounce = (-b - sqrt(b^2 - 4*a*c))/(2*a);

        c = betaPre(1,:);
        b = betaPre(2,:);
        a = betaPre(3,:);
        y = [b*tBounce + c; (2*a*tBounce + b)*M - 2*a*tBounce];
        mat = [1, tBounce; 0, 1];
        betaPost = [mat \ y; a];
    end       

end