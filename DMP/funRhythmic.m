function dX = funRhythmic(t,X,PAR)

alpha_z = PAR.time_const(1);
beta_z = PAR.time_const(2);
tau = PAR.scale;
g = PAR.goal(1);
Amp = PAR.goal(2);
fpar = PAR.forcing;

A = [-alpha_z/tau, -alpha_z*beta_z/tau, 0;
     1/tau, 0, 0;
     0, 0, 0];

funBasis = @(phi,h,c) vonMises(phi,h,c); 
f = forcing(X(3),fpar.w,fpar.h,fpar.c,funBasis) * Amp;
% forcing function acts on the accelerations
B = [alpha_z*beta_z*g/tau + f/tau; 0; 1/tau];

dX = A*X + B;

end

% forcing function to drive nonlinear system dynamics
function f = forcing(x,w,h,c,basis)

N = length(w);
f = 0;
scale = 0;
for i = 1:N
    f = f + basis(x,h(i),c(i))*w(i);
    scale = scale + basis(x,h(i),c(i));
end

f = f/scale;
end