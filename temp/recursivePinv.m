%% Function for constructing pseudoinverse recursively

% Proof of concept only, not a practical implementation!
% A,B,C,D are real block matrices with matching sizes

function pM = recursivePinv(A,B,C,D)

% preliminary matrices
pA = pinv(A);
K = pA * B;
H = C * pA;
Z = D - H * B;
% size of full matrix M
na1 = size(A,1);
na2 = size(A,2);
nb2 = size(B,2);
nc1 = size(C,1);
m = na1 + nc1;

% check if A is square
if na1 == na2
    Ai = zeros(na1);
    iA = zeros(na2);    
    P = zeros(na1,nb2);
    pP = zeros(nb2,na1);
else    
    % check if A has full rank
    if rank(A) == na1
        Ai = zeros(na1);
        iA = eye(na2) - pA*A;
        P = zeros(na1,nb2);
        pP = zeros(nb2,na1);
    elseif rank(A) == na2
        iA = zeros(na2);
        Ai = eye(na1) - A*pA;
        P = Ai * B;
        pP = pinv(P);
    else
        Ai = eye(na1) - A*pA;
        iA = eye(na2) - pA*A;
        P = Ai * B;
        pP = pinv(P);
    end
end

Q = C * iA;

%% Special cases where P and Q are both zero

if ~(sum(any(P))) && ~sum(any(Q))
    
    X = [-K; eye(size(K,2))];
    Y = [-H, eye(size(H,1))];
    Khat = eye(size(K,2)) + K'*K;
    Hhat = eye(size(H,1)) + H*H';   
    Hhatinvsqrt = sqrtm(inv(Hhat));
    Khatinvsqrt = sqrtm(inv(Khat));
    % weighted pseudoinverse of Z
    Zg = Khatinvsqrt*pinv(Hhatinvsqrt*Z*Khatinvsqrt)*Hhatinvsqrt;

    pM1 = eye(size(X,1)) - X*(eye(size(Zg,1)) - Zg*Z)*pinv(X);
    pM2 = [pA, zeros(size(A,2),size(Y,2)-size(A,1));
           zeros(size(X,1)-size(A,2),size(Y,2))];
    pM3 = eye(size(Y,2)) - pinv(Y)*(eye(size(Z,1)) - Z*Zg)*Y;
    pM4 = X*Zg*Y;
    
    pM = pM1*pM2*pM3 + pM4;
    
    return;
end
   

%% Switching based on rank of P, in case P != 0
% check for full column rank of P, things simplify a lot!
switch rank(P)
    
    % case where full column rank of P
    case size(P,2)
    
    U = pP;
    HZU = H + Z*U;
    
    iP = zeros(nb2);    
    W = zeros(nc1,nb2);
    Wg = zeros(nb2,nc1);
    V = zeros(nb2,nc1);
    S = zeros(nb2);    

    X1t = [-HZU, eye(nc1)];
    K1hat = eye(nc1) + (HZU*HZU');
    K1hatinvsqrt = sqrtm(inv(K1hat));
    Qg = pinv(K1hatinvsqrt*Q)*K1hatinvsqrt;
    
    pM1 = [pA - K*U, zeros(nb2,m-size(pA,2)); U, zeros(nb2,m-size(pA,2))];
    pM2 = eye(size(X1t,2)) - pinv(X1t)*(eye(nc1) - Q*Qg - W*Wg)*X1t;
    pM3 = [Qg - Qg*Z*V - K*V; V];
    pM4 = X1t;

    pM = pM1*pM2 + pM3*pM4;

    otherwise
    
    % full complexity case not fully implemented
    assert(sum(any(Q)) ~= 0,'Q zero case not implemented!');
   
    if ~sum(any(P)) % A has full row rank
        iP = eye(nb2);
        S = K'*K;
        U = pinv(K)*pA*eye(size(B,1));
    else
        iP = eye(nb2) - pP*P;
        S = iP*(K'*K)*iP;
        U = pP + iP*(S\(K'))*pA*(eye(size(B,1))-B*pP);
    end

    % secondary matrices
    HZU = H + Z*U;
    X1t = [-HZU, eye(nc1)];
    K1hat = eye(nc1) + (HZU*HZU');
    K1hatinvsqrt = sqrtm(inv(K1hat));

    Qg = pinv(K1hatinvsqrt*Q)*K1hatinvsqrt;
    T = S + iP*Z'*(Qg'*Qg)*Z*iP;
    Tinv = inv(T);
    Tinvsqrt = sqrtm(Tinv);
    W = (eye(nc1) - Q*Qg)*Z*iP;
    Wg = Tinvsqrt*pinv(K1hatinvsqrt*W*Tinvsqrt)*K1hatinvsqrt;
    V = Wg + (eye(nb2)-Wg*W)*Tinv*iP*Z'*(Qg'*Qg);

    pM1 = [pA - K*U, zeros(na2,m-size(pA,2)); U, zeros(nb2,m-size(pA,2))];
    pM2 = eye(size(X1t,2)) - pinv(X1t)*(eye(nc1) - Q*Qg - W*Wg)*X1t;
    pM3 = [Qg - Qg*Z*V - K*V; V];
    pM4 = X1t;

    pM = pM1*pM2 + pM3*pM4;

end




