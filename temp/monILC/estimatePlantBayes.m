% Function for estimating the input-output matrix F 
% with Linear Bayesian Regression

function [F,G,a,b] = estimatePlantBayes(F,G,a,b,u,e,N)

    K = size(u,2);
    m = size(u,1)/N;
    n = size(e,1)/N;
    
    % get F_nom back into estimation friendly form
    for i = 1:N
        l = (i-1)*(i/2)*m*n;
        mu(l+1:l+i*m*n) = F((i-1)*m+1:i*m,1:n*i);
    end
    
    % estimate F
    for i = 1:K
        U((N*(i-1))+1:(N*(i-1))+N,1:m*N*N) = kron(eye(N),u(:,i)');
        E((N*(i-1))+1:(N*(i-1))+N,1:n) = reshape(e(:,i)',n,N)';
    end
    D = duplicateVech(N,m);
    M = U * D';
    
    % form M and E into X and y matrices
    X = kron(eye(n),M);
    y = E(:);
    [mu,G,a,b] = LBR(mu,G,a,b,X,y);
    
    % get mu back to Est form
    L = N * (N+1)/2;
    Est = reshape(mu,m,L);
    
    %Est = pinv(M) * E;
    EstAddZero = Est' * D;

    for i = 1:N
        F((i-1)*n+1:i*n,:) = EstAddZero(:,N*m*(i-1)+1:N*m*i);
    end
    
end

function [mu,Gamma,a,b] = LBR(mu,Gamma,a,b,X,y)

    n = size(X,1);
    % calculate posterior
    a = a + n/2;
    b = b + 0.5 * (y'*y + mu'*Gamma*mu);
    mu = (X'*X + Gamma)\(Gamma*mu + X'*y);
    Gamma = (X'*X + Gamma);
    b = b - 0.5 * (mu'*Gamma*mu);

end