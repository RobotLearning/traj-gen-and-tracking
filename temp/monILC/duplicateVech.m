% Produce duplication matrix for transforming n-dim block vech to vec

function D = duplicateVech(N,m)

% assuming A matrix is nxm

pos = ones(1,N^2*m);
vec = [];
for i = 1:N-1
    val = ((i-1)*m*(N+1))+m+1:i*N*m; %i*N*n*m+1:i*(N+1)*n*m;
    vec = [vec,val]; 
    pos(val) = 0;
end

D = diag(pos);
D = D(setdiff(1:N^2*m,vec),:);
