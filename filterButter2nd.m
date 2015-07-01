%% Second order Butterworth filter

function y = filterButter2nd(x,w)

    [B,A] = myButter2ndOrder(w);

    % initial conditions of the filter
    %lhs = [1+A(2), -1; A(3), 1];
    %rhs  = B(2:3) - B(1)*A(2:3);
    %zi   = lhs \ rhs(:);
    
    zi = A(2:3) \ (sum(B) - A(1));
    % set previous values to the starting value
    x1 = x(1); 
    x2 = x(1);
    y1 = zi(1)*x(1);
    y2 = zi(2)*x(1);

    for i = 1:length(x)
        [y(i),y1,y2,x1,x2] = filterNext(A,B,x(i),x1,x2,y1,y2);
    end
end

% get next filtered value
function [y,y1,y2,x1,x2] = filterNext(A,B,x,x1,x2,y1,y2)

    y = ((B(1)*x + B(2)*x1 + B(3)*x2 ) - (A(2)*y1 + A(3)*y2)) / A(1);
    x2 = x1;
    x1 = x;
    y2 = y1;
    y1 = y;

end

% construct the IIR tap weights
function [B, A] = myButter2ndOrder(W)

    V  = tan(W * pi/2);  
    Sg = V ^ 2;
    Sp = V * [-1-1i, -1+1i] / sqrt(2);
    % Bilinear transform:
    P = (1 + Sp) ./ (1 - Sp);
    G = real(Sg / prod(1 - Sp));
    % From Zeros, Poles and Gain to numerator and denominator:
    B = G * [1, 2, 1];
    A = real(poly(P));

end