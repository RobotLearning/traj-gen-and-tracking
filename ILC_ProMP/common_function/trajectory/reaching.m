
function [t, x, z, a ] = reaching(xz_i, xz_f, a_i, a_f, dt, tf)

    dbg = 0; % use this to plot debugging figures
   
    t_discount = 0.0;
    
    % first generate a trajectory that is little shorter than expected
    tf_safe = tf - t_discount ;
    f  = 1/tf_safe;
    t  = 0:dt:tf_safe;
        
        
    Ax = xz_f(1)-xz_i(1);
    Az = xz_f(2)-xz_i(2);
    Aa = a_f - a_i;
    
    
    x = xz_i(1) + 0.5*Ax*(1 + cos(2*pi*f*t+pi));
    z = xz_i(2) + 0.5*Az*(1 + cos(2*pi*f*t+pi));
    a =  a_i + 0.5*Aa*(1 + cos(2*pi*f*t+pi)); % orientation of the end-effector
    
    
    % add the t_discount 
    
    tfull  = 0:dt:tf;
    x = [x  x(end)*ones(1, numel(tfull) - numel(t))];
    z = [z  z(end)*ones(1, numel(tfull) - numel(t))];
    a = [a  a(end)*ones(1, numel(tfull) - numel(t))];
    t = tfull;

    if dbg
        close all
        h = figurew('fe'); axis 'equal'
        plot(x, z, 'b.-'); hold on;

        figurew('xya')
        subplot(3,1,1); hold on;
        plot(t, x)
        subplot(3,1,2); hold on;
        plot(t, z)
        subplot(3,1,3); hold on;
        plot(t, r2d(a));
        ylabel 'end-eff (deg)'
        keyboard
    end

    
end
