% Estimate a rebound trajectory using polynomial fitting
% Taking care of outliers using ransac

function ballTraj = fitPoly2ReboundTraj(t,b,matBounce)

    len = length(t);
    mat = [ones(len,1), t, t.^2];
    
    % check if there is a bounce
    tol = 5e-2;
    
    idxBallBounce = find(y(:,3) <= table_z + ball_radius + tol,1);
    ballPreBounce = b(1:idxBallBounce,:);
    betaPreBounce = mat \ ballPreBounce;
    %betaPreBounce = pinv(M,tol)*Y;

    % find the time till bounce

    m = betaPreBounce(:,3);
    a = m(3);
    b = m(2);
    c = m(1) - table_z - ball_radius;
    tBouncePredLS = (-b - sqrt(b^2 - 4*a*c))/(2*a);

    dt = 1/60;
    N_pred_LS = round(tBouncePredLS/dt);
    tLS = dt * (1:N_pred_LS);
    Mstar = [ones(N_pred_LS,1),tLS(:),tLS(:).^2];
    ballLS = Mstar * betaPreBounce;% + [zeros(N_pred_LS,2),0.5*gravity*tLS(:).^2];
    ballLLS = ballLS';


    % transform ballInitLS
    mat = [1, tBouncePredLS;
           0, 1];

    c = betaPreBounce(1,:);
    b = betaPreBounce(2,:);
    a = betaPreBounce(3,:);
    y = [b*tBouncePredLS + c; (2*a*tBouncePredLS + b)*matBounce - 2*a*tBouncePredLS];
    betaAfterBounce = mat \ y;
    betaAfterBounce = [betaAfterBounce; a];

    N_pred_LS = round((tPredict - tBouncePredLS)/dt);
    tLS = tBouncePredLS + dt * (1:N_pred_LS);
    MstarAfter = [ones(N_pred_LS,1),tLS(:),tLS(:).^2];
    ballLSAfter = MstarAfter * betaAfterBounce;
    ballTraj = [ballLLS,ballLSAfter'];
end