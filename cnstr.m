function [L,q] = cnstr(str,d)

% constraint generating function. Each folder has a different
% implementation

C = str.C;
F = str.F;
x_con = str.x_con;

L = C*F;
q = x_con - C*d; 