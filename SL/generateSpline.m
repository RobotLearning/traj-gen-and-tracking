% Generate a striking and a returning trajectory
% Based on optimization results qf,qfdot,T and given q0,q0dot
function [q,qd,qdd] = generateSpline(dt,q0,q0dot,qf,qfdot,T,Tret)

    dof = length(q0);
    time2hit = T;
    time2return = Tret;
    Q0 = [q0;q0dot];  
    Qf = [qf;qfdot];

    % round time2hit to nearest dt
    time2hit = dt * ceil(time2hit/dt);

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