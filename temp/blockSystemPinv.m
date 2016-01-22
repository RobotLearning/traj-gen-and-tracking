%% Recursive plant inversion

% for LTI plant
% A,B are the system matrices 
% A is square and has full rank
% B is assumed to have full column rank

function Finv = blockSystemPinv(A,B)

n = size(A,1);
Binv = pinv(B);
Pb = B*Binv;
H = A*Pb;
Z = B;
I = eye(n);    

Hhat = I + A*Pb*A';   
%Hhatinv = inv(Hhat);
%Hhatinvsqrt = sqrtm(Hhatinv);
% weighted pseudoinverse of Z
%Zg = pinv(Hhatinvsqrt*Z)*Hhatinvsqrt;
Zbar = (Z'*(Hhat\Z))\(Z');
%Zg = (Z'*(Hhat\Z))\(Z'*Hhatinv);
Zg = Zbar / Hhat;

Zp = H' * (Hhat \ (I - Z*Zg));

M12 = Binv * Zp;
M11 = Binv  - M12 * H;
M21 = -Zg*H;
M22 = Zg;

Finv = [M11,M12;M21,M22];

end

%% Assuming B is 0
function TwoByTwoBlockInv(Ainv,Cinv,Dinv)

    

end
    

   





