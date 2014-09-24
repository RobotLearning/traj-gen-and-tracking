function dX = funDiscrete(t,X,PAR)

alpha_z = PAR.alpha_g;
beta_z = PAR.beta_g;
alpha_x = PAR.alpha_x;
tau = PAR.scale;
g = PAR.g;
y0 = PAR.X0(2);
fpar = PAR.forcing;

A = [-alpha_z/tau, -alpha_z*beta_z/tau, 0;
     1/tau, 0, 0;
     0, 0, -alpha_x/tau];

funBasis = @(x,h,c) sqExp(x,h,c);
f = forcing(X(3),fpar.w,fpar.h,fpar.c,funBasis) * (g-y0);
% forcing function acts on the accelerations
B = [alpha_z*beta_z*g/tau + f/tau; 0; 0];

dX = A*X + B;

end

% forcing function to drive nonlinear system dynamics
function f = forcing(x,w,h,c,basis)

N = length(w);
f = 0;
scale = 0;
for i = 1:N
    f = f + basis(x,h(i),c(i))*w(i)*x;
    scale = scale + basis(x,h(i),c(i));
end

f = f/scale;
end