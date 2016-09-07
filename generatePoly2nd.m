% generate2nd order polynomials
% output is of the form [q;qd;qdd]
%
% qfdot is assumed to be zero
%
function p = generatePoly2nd(q0,q0dot,dt,tf)

    dof = length(q0);
    t = dt:dt:tf;
    a = -q0dot/(2*tf); % accelerations
    qStrike = zeros(dof,length(t));
    qdStrike = zeros(dof,length(t));
    qddStrike = zeros(dof,length(t));

    for m = 1:dof
        qStrike(m,:) = a(m)*t.^2 + q0dot(m)*t + q0(m);
        qdStrike(m,:) = 2*a(m)*t + q0dot(m);
        qddStrike(m,:) = 2*a(m);
    end

    p = [qStrike; qdStrike; qddStrike];
end