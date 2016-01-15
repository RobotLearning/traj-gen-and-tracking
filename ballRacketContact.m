function velOut = ballRacketContact(ball,robot,q,qd,robotIdx)

loadTableTennisValues();

[racketCentre,racketVel] = robot.getEndEffectorState(q,qd);
racketDir = 

tol = ball_radius;
l = min(robotIdx,size(q,2));
vecFromRacketToBall = ball(1:2) - racketCentre;
racketPlane = racket_dir(:,l);
projPlane = racketPlane*racketPlane'/(racketPlane'*racketPlane);
projOrth = eye(2) - projPlane;
distToRacketPlane = norm(projOrth * vecFromRacketToBall);
distOnRacketPlane = norm(projPlane * vecFromRacketToBall);
if distToRacketPlane < tol && distOnRacketPlane < racket_radius && ~hit            
    %disp('A hit! Well done!');
    hit = true;
    fprintf('Hit at y = %f, z = %f\n',ball(1,ballIdx),ball(2,ballIdx));
    numHits = numHits + 1;
    % Change ball velocity based on contact model
    % get ball velocity
    velIn = ball(3:4,ballIdx);
    velInAlongNormal = projOrth * velIn;
    % get racket velocity
    velRacket = (path(5:6,l+1) - path(5:6,l-1))/(2*dt);
    velRacketAlongNormal = projOrth * velRacket;
    % this is kept the same in mirror law
    velInAlongRacket = projPlane * velIn; 
    velOutAlongNormal = velRacketAlongNormal + ...
        CRR * (velRacketAlongNormal - velInAlongNormal);
    velOut = velOutAlongNormal + velInAlongRacket;
    ball(3:4,ballIdx) = velOut;
end