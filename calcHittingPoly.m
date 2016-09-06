% Generate a striking trajectory only
% Using the lazy player (optimization)
% qfdot is zero
% Based on optimization results qf,qfdot,T and given q0,q0dot
function [q,qd,qdd] = calcHittingPoly(dt,q0,q0dot,qf,T,Tland)

    dof = length(q0);
    qfdot = zeros(dof,1);
    time2hit = T;
    time2return = Tland;
    Q0 = [q0;q0dot];  
    Qf = [qf;qfdot];

    % round time2hit to nearest dt
    time2hit = dt * ceil(time2hit/dt);

    % GET 3RD DEGREE POLYNOMIAL FOR STRIKE            
    pStrike = generatePoly3rd(Q0,Qf,dt,time2hit);
    qStrike = pStrike(1:dof,:);
    qdStrike = pStrike(dof+1:2*dof,:);
    qddStrike = pStrike(2*dof+1:end,:);

    % CONTINUE WITH 2ND ORDER POLYNOMIAL TILL LANDING TIME
    pReturn = generatePoly2nd(q0,q0dot,dt,time2return);
    qReturn = pReturn(1:dof,:);
    qdReturn = pReturn(dof+1:2*dof,:);
    qddReturn = pReturn(2*dof+1:end,:);

    q = [qStrike,qReturn];
    qd = [qdStrike,qdReturn];
    qdd = [qddStrike,qddReturn];            

end