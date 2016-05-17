% Generate 2D table tennis trajectories with the VHP method 
function [q,qd,qdd] = generate2DTTTwithVHP(robot,ballPred,ballTime,q0)

    loadTennisTableValues();
    dof = length(q0);
    par.CRR = CRR;

    % define virtual hitting plane (VHP)
    VHP = -0.6;
    time2reach = 0.5; % time to reach desired point on opponents court
    time2return = 0.5; % time to return to initial configuration          
    dt = ballTime(2)-ballTime(1);
    ballFull = [ballPred;ballTime];

    % rotate some variables for drawing in 2D simulation
    R = [0 -1; 1 0];
    % land the ball on the centre of opponents court
    ballDes(1) = dist_to_table - 3*table_y/2;
    ballDes(2) = table_z + ball_radius;            
    fprintf('Desired landing point: %f\n',ballDes(1));

    ballAtVHP = interp1(ballFull(1,:)',ballFull(2:5,:)',VHP);
    timeAtVHP = ballAtVHP(end);
    ballAtVHP = [VHP;ballAtVHP(1:end-1)'];
    ballPosAtVHP = ballAtVHP(1:2);
    ballInVelAtVHP = ballAtVHP(3:4); 

    % GET DESIRED OUTGOING VELOCITY OF THE BALL AT VHP            
    ballOutVelAtVHP = calcBallVelOut2D(ballDes,ballPosAtVHP,time2reach);            

    % GET RACKET DESIRED VEL AND ORIENTATION AT VHP 
    [racketPos,racketVel,racketOrient] = calcDesRacketState ...
                   (ballPosAtVHP,ballOutVelAtVHP,ballInVelAtVHP,par);

    % feed to inverse kinematics to get qf
    try
        normalRot = R'*racketOrient;
        phiVHP = atan2(normalRot(2),normalRot(1));
        qf = robot.invKinematics(R'*racketPos,-phiVHP);
        robot.calcJacobian(qf);
        qfdot = robot.jac \ (R'*racketVel);
    catch ME
        disp('Virtual Hitting Point outside of workspace');
        qf = q0;
        qfdot = zeros(dof,1);
    end

    q0dot = zeros(dof,1);
    Q0 = [q0;q0dot];
    Qf = [qf;qfdot];

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

end