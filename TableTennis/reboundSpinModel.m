% K is the coefficient values in x-y-z directions
function xdot = reboundSpinModel(xdot,e_t,mu,r)

    vbT = [xdot(1) - r*xdot(5);
           xdot(2) + r*xdot(4);
           0];

    nu_s = 1 - (2/5)*mu*(1+e_t)*abs(xdot(3))/norm(vbT);
    alpha = 2/5; % roll
    if nu_s > 0 % slide 
        alpha = mu * (1 + e_t) * abs(xdot(3))/norm(vbT);
    end

    Av = diag([1.2,1-alpha,-e_t]);
    Bv = [0, alpha*r, 0;
          -alpha*r, 0, 0;
          zeros(1,3)];
    Aw = [0, -3*alpha/(2*r), 0;
          3*alpha/(2*r), 0, 0;
          zeros(1,3)];
    Bw = diag([1-3*alpha/2, 1-3*alpha/2, 1]);

    M = [Av, Bv; Aw, Bw];
    xdot = M * xdot;
    
end