% Calculate desired racket positions,velocities and normals
% based on a given desired strategy,
% in this case desired landing positions desBall

function racketDes = calcRacketStrategy2D(desBall,ballPred,ballTime,...
                                        time2reach,fast)

loadTennisTableValues();     
PAR.Cdrag = Cdrag;
PAR.g = gravity;
PAR.fast = fast;
PAR.CRR = CRR;
                                                                     
for j = 1:size(ballPred,2)

    velOut = calcBallVelOut2D(desBall,ballPred(1:2,j),time2reach,PAR);              
    % Use the inverse contact model to compute racket vels and normal
    % at every point                
    [rp,rv,ro] = calcDesRacketState(ballPred(1:2,j),velOut,ballPred(3:4,j),PAR);
    racketDes.time(j) = ballTime(j);
    racketDes.pos(:,j) = rp;
    racketDes.normal(:,j) = ro;
    racketDes.vel(:,j) = rv;

end

% 8 cm radius
racketDes.radius = 0.08;