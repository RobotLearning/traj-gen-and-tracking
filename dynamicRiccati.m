% Dynamic Riccati equation
% Can be iterated till it converges to the infinite-horizon matrix Pinf

function Ppre = dynamicRiccati(Pnext,Q,R,A,B)

Ppre = Q + A'*Pnext*A - A'*Pnext*B*((R + B'*Pnext*B)\((B'*Pnext)*A));

end