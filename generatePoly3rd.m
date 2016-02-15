% generate 3rd degree polynomials
% output is of the form [q;qd;qdd]
%
% TODO: is it reliable? can we do this faster?
%
function p = generatePoly3rd(Q0,Qf,dt,tf)


    dof = length(Q0)/2;
    q0 = Q0(1:dof);
    q0dot = Q0(dof+1:end);
    qf = Qf(1:dof);
    qfdot = Qf(dof+1:end);
    t = dt:dt:tf;

    M = @(t0,tf) [t0^3 t0^2 t0^1 1;
                  3*t0^2 2*t0 1 0;
                  tf^3 tf^2 tf^1 1;
                  3*tf^2 2*tf 1 0];

    for m = 1:dof
        %q0dot is zero
        Qstrike = [q0(m); q0dot(m); qf(m); qfdot(m)]; % strike
        a = M(0,tf) \ Qstrike;
        qStrike(m,:) = a(1)*t.^3 + a(2)*t.^2 + a(3)*t + a(4);
        qdStrike(m,:) = 3*a(1)*t.^2 + 2*a(2)*t + a(3);
        qddStrike(m,:) = 6*a(1)*t + 2*a(2);
    end

    p = [qStrike; qdStrike; qddStrike];
end