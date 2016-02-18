% Generate 3D table tennis trajectories with the VHP method
function [q,qd,qdd] = generate3DTTTwithVHP(robot,VHP,ballPred,ballTime,q0)

    loadTennisTableValues();
    dof = length(q0);

    % Check if ball will bounce twice
    tol = 1e-2;
    [M,I] = min(ballPred(3,:));
    if M < table_z + tol && ballPred(2,I) < dist_to_table
        disp('Ball will rebound twice on table! Not moving!');
        N = 100;
        q = q0 * ones(1,N);
        qd = zeros(dof,N);
        qdd = zeros(dof,N);
        return;
    end

    time2reach = 0.5; % time to reach desired point on opponents court
    % land the ball on the centre of opponents court
    ballDes(1) = 0.0;
    ballDes(2) = dist_to_table - 3*table_y/2;
    ballDes(3) = table_z + ball_radius;
    [qf,qfdot,timeAtVHP] = calcPolyAtVHP(robot,VHP,time2reach,ballDes,ballPred,ballTime,q0);

    q0dot = zeros(dof,1);
    Q0 = [q0;q0dot];  
    Qf = [qf;qfdot];            
    dt = ballTime(2) - ballTime(1);
    time2return = 0.5; % time to return to initial configuration

    % not moving in case calculation gone wrong
%     eps = 1e-2;
%     if norm(qf-q0, 2) < eps
%         disp('Not moving!');
%         totalTime = timeAtVHP + time2return;
%         N = floor(totalTime/dt);
%         q = q0 * ones(1,N);
%         qd = zeros(dof,N);
%         qdd = zeros(dof,N);
%         return;
%     end

    % GET 3RD DEGREE POLYNOMIALS            
    pStrike = generatePoly3rd(Q0,Qf,dt,timeAtVHP);
    qStrike = pStrike(1:dof,:);
    qdStrike = pStrike(dof+1:2*dof,:);
    qddStrike = pStrike(2*dof+1:end,:);

    pReturn = generatePoly3rd(Qf,Q0,dt,time2return);
    qReturn = pReturn(1:dof,:);
    qdReturn = pReturn(dof+1:2*dof,:);
    qddReturn = pReturn(2*dof+1:end,:);

    q = [qStrike,qReturn];
    qd = [qdStrike,qdReturn];
    qdd = [qddStrike,qddReturn];         

    % for debugging
    %[x,xd,o] = obj.calcRacketState([q;qd]);
    %rotMs = quat2Rot(o);
    %normals = rotMs(1:3,3,:);
end