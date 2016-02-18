% Calculate desired racket positions,velocities and normals
% based on a given desired strategy,
% in this case desired landing positions desBall

function racketDes = calcRacketStrategy(desBall,ballPred,ballTime,...
                                        time2reach,fast)

for j = 1:size(ballPred,2)

    velOut = calcBallVelOut3D(desBall,ballPred(1:3,j),time2reach,fast);              
    % Use the inverse contact model to compute racket vels and normal
    % at every point                
    [rp,rv,ro] = calcDesRacketState(ballPred(1:3,j),velOut,ballPred(4:6,j));
    racketDes.time(j) = ballTime(j);
    racketDes.pos(:,j) = rp;
    racketDes.normal(:,j) = ro;
    racketDes.vel(:,j) = rv;

end