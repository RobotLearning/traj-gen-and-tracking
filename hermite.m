function hp = hermite ( x, y, yp )

%HERMITE      determine the Hermite interpolating polynomial associated 
%             with a given set of interpolating points, function values
%             and derivative values
%
%     calling sequences:
%             hp = hermite ( x, y, yp )
%             hermite ( x, y, yp )
%
%     inputs:
%             x       vector containing the interpolating points
%             y       vector containing function values
%                     the i-th entry in this vector is the function
%                     value associated with the i-th entry in the 'xi'
%                     vector
%             yp      vector containing derivative values
%                     the i-th entry in this vector is the derivative
%                     value associated with the i-th entry in the 'xi'
%                     vector
%
%     output:
%             hp      vector containing coefficients of the Hermite
%                     interpolating polynomial associated with the given
%                     set of interpolating points, function values and
%                     derivative values
%
%     NOTE:
%             to evaluate the Hermite interpolating polynomial, apply the 
%             MATLAB function 'polyval' to the output from this routine 
%

n = length ( x );
z = zeros ( 1, 2*n );
f = zeros ( 1, 2*n );

z(1:2:2*n-1) = x;
z(2:2:2*n)   = x;
f(1)         = y(1);
f(3:2:2*n-1) = ( y(2:n) - y(1:n-1) ) ./ ( x(2:n) - x(1:n-1) );
f(2:2:2*n)   = yp;

for i = 3:2*n
    f(i:2*n) = ( f(i:2*n) - f(i-1:2*n-1) ) ./ ( z(i:2*n) - z(1:2*n-i+1) );
end;

hp = zeros ( 1, 2*n );
p = [1];
for i = 1:2*n
    hp = hp + f(i) * [ zeros(1,2*n-i) p ];
	p = conv ( p, [ 1 -z(i) ] );
end;