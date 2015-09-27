function x = tls(A,b,thresh)
% Solves the linear equation Ax=b using
% truncated total least squares.
%
% Ref. Regularization by Truncated TLS
% Fierro et al, Siam J Sci Comput 1997
%
% INPUTS:
% A is the matrix
% b is the RHS vector
% thresh (optional) is the SVD threshold
%
% OUTPUTS:
% x is the estimate
%
% Modified by: Okan KOC, September 27, 2015
%

[m, n] = size(A);
if sum(size(b)-[m 1])
    error('A, b size mis-match')
end
if nargin == 2
    thresh = eps;
end

% augmented matrix
Z = [full(A');b'];
[U, S, V] = svd(Z);

% find sing val above threshold
d = diag(S);
k = sum(d < thresh);
q = n - k + 1;

% proceed to checkout
V12 = U(1:n, q:end);
V22 = U(n+1, q:end);
x = -V12 * V22' ./ norm(V22).^2;


%verbose = true;
%if verbose
    %fprintf('Truncating %d singular values\n',k);
    %fprintf('Sigma_n+1 value is %f\n', d(end));
%end

end

%{

function beta = tls(xdata,ydata)


m       = length(ydata);       %number of x,y data pairs
X       = xdata;
Y       = ydata;
n       = size(X,2);          % n is the width of X (X is m by n)
Z       = [X Y];              % Z is X augmented with Y.
[U, S, V] = svd(Z,0);           % find the SVD of Z.
VXY     = V(1:n,1+n:end);     % Take the block of V consisting of the first n rows and the n+1 to last column
VYY     = V(1+n:end,1+n:end); % Take the bottom-right block of V.
beta       = -VXY/VYY;
end
%}