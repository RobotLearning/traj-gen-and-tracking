% Discrete time finite horizon linear quadratic regulator
% lqr function gives the optimal linear feedback law
%
% Inputs:
% Q,R - Weight matrices
% A,B - Matrices of the differential equation
% Mode - Infinite horizon or finite horizon, as well as time
% varying/invariant
%
% Outputs:
% K - Feedback gain
%
% TODO: extend to include feedforward term
%       Qf can be also input 

function K = LQR(Q,R,A,B,MODE)

dimx = size(Q,1);
dimu = size(R,1);
% Qf should be set equal to Q in case Qf is not provided
Qf = Q;

% divide to cases
if MODE.N == Inf

    Pinf = Qf;
    % iterate till convergence
    delta = 1e-3;
    diff = 1;
    i = 0;
    while diff > delta
        Ppre = Pinf;
        Pinf = dynamicRiccati(Pinf,Q,R,A,B);
        matDiff = Ppre - Pinf;
        diff = norm(matDiff,'fro');
        i = i+1;
    end
    fprintf('Iterated dynamic riccati equation %d times\n', i);
    K = -(R + B'*Pinf*B)\(B'*Pinf*A);

elseif MODE.N ~= Inf && MODE.LTI == true

    % finite horizon
    N = MODE.N;

    P = zeros(dimx,dimx,N);
    P(:,:,end) = Qf;

    % calculate optimal K
    K = zeros(dimu,dimx,N);
    K(:,:,end) = -(R + B'*Qf*B)\(B'*Qf*A);
    for i = N-1:-1:1
        P(:,:,i) = dynamicRiccati(P(:,:,i+1),Q,R,A,B);
        K(:,:,i) = -(R + B'*P(:,:,i)*B)\(B'*P(:,:,i)*A);
    end

elseif MODE.N ~= Inf && MODE.LTI == false
    % Arrays A and B hold time-varying matrices

    % finite horizon
    N = MODE.N;

    P = zeros(dimx,dimx,N);
    P(:,:,end) = Qf;

    % calculate optimal K
    K = zeros(dimu,dimx,N);
    K(:,:,end) = -(R + B(:,:,end)'*Qf*B(:,:,end))\...
                  (B(:,:,end)'*Qf*A(:,:,end));
    for i = N-1:-1:1
        P(:,:,i) = dynamicRiccati(P(:,:,i+1),Q,R,A(:,:,i),B(:,:,i));
        K(:,:,i) = -(R + B(:,:,i)'*P(:,:,i)*B(:,:,i))\...
                    (B(:,:,i)'*P(:,:,i)*A(:,:,i));
    end

end

end

