% Generate free-time optimal table tennis trajectories

function [q,qd,qdd] = generateOptimalTTT(robot,racket,dt,q0)

    dof = length(q0);
    time2return = 1.0;

    % FOR TESTING SL CONNECTION
    %qest = [2.25;-0.38;-1.27;1.33;-1.86;-0.19;0.77];
    %qdest = [1.10;-0.53;-3.21;1.26;-0.44;-0.09;0.78];
    %timeEst = 0.8;

    %qf = qest;
    %qfdot = qdest;
    %time2hit = timeEst;
    [qf,qfdot,time2hit] = calcOptimalPoly(robot,racket,q0);
    % round time2hit to nearest dt
    time2hit = dt * ceil(time2hit/dt);

    q0dot = zeros(dof,1);
    Q0 = [q0;q0dot];  
    Qf = [qf;qfdot];

    % GET 3RD DEGREE POLYNOMIALS            
    pStrike = generatePoly3rd(Q0,Qf,dt,time2hit);
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