% Function for estimating the input-output matrix F

function F_est = estimatePlant(u,e,N)

    K = size(u,2);
    m = size(u,1)/N;
    n = size(e,1)/N;
    
    % estimate F
    for i = 1:K
        U((N*(i-1))+1:(N*(i-1))+N,1:m*N*N) = kron(eye(N),u(:,i)');
        E((N*(i-1))+1:(N*(i-1))+N,1:n) = reshape(e(:,i)',n,N)';
    end
    D = duplicateVech(N,m);
    M = U * D';
    
    Est = pinv(M) * E;
    EstAddZero = Est' * D;

    for i = 1:N
        F_est((i-1)*n+1:i*n,:) = EstAddZero(:,N*m*(i-1)+1:N*m*i);
    end
    
end