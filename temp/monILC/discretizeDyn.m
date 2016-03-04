% Discretize dynamics

function [Ad,Bd] = discretizeDyn(A,B,dt)

N = size(A,3);
assert(size(A,3) == size(B,3));
assert(size(A,2) == size(A,1));
assert(size(A,2) == size(B,1));
n = size(B,1);
m = size(B,2);

Ad = zeros(n,n,N);
Bd = zeros(n,m,N);

for i = 1:N
    Mat = [A(:,:,i), B(:,:,i); zeros(m, n + m)];
    MD = expm(dt * Mat);
    Ad(:,:,i) = MD(1:n,1:n);
    Bd(:,:,i) = MD(1:n,n+1:end);
end