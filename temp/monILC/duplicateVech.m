% Produce duplication matrix for transforming vech to vec

function D = duplicateVech(N,n,m)

% assuming A matrix is nxm
D = zeros(N^2 * n * m, (N*(N-1)/2 + N)*(n*m));

% assuming for now that n,m are the same
for i = 1:N
    D((i-1)*N*n + 1:i*N*n,:) = [zeros((i-1)*n, (N*(N-1)/2 + N)*(n*m)); ...
                  zeros((N-i+1)*n,(i-1)*N),eye((N-i+1)*n), zeros(N*n,(N*(N-1)/2) * (n*m))];

end