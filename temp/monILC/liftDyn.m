% Lift the linearized dynamics around a reference trajectory

function F = liftDyn(Ad,Bd)

N = size(Ad,3);
assert(size(Ad,3) == size(Bd,3));
assert(size(Ad,2) == size(Ad,1));
assert(size(Ad,2) == size(Bd,1));
n = size(Bd,1);
m = size(Bd,2);

F = zeros(N*n,N*m);
% construct lifted domain matrix F
for i = 1:N
    vec_x = (i-1)*n + 1:i*n;
    for j = 1:i                            
        vec_u = (j-1)*m + 1:j*m;
        mat = Bd(:,:,j);
        for k = j+1:i
            mat = Ad(:,:,k) * mat; % (Ad + Bd * K(:,:,k)) * mat;
        end
        F(vec_x,vec_u) = mat; % on diagonals only B(:,m)
    end
end                