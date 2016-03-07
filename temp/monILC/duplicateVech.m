% Produce duplication matrix for transforming n-dim block vech to vec

function D = duplicateVech(N,n,m)

% assuming A matrix is nxm

pos = ones(1,N^2*n);
vec = [];
for i = 1:N-1
    val = ((i-1)*n*(N+1))+n+1:i*N*n; %i*N*n*m+1:i*(N+1)*n*m;
    vec = [vec,val]; 
    pos(val) = 0;
end

D = diag(pos);
D = D(setdiff(1:N^2*n,vec),:);
