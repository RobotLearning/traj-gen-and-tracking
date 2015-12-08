
function [t, x, z] = square(xz_i, Amp, dt, tf, rev)


%% DID NOT FINISH THIS IMPLEMENTATION, BUT IT IS POSSIBLE

if 0 % TO DO

    dbg = 1; % use this to plot debugging figures
    
    rev = rev+0.25; % add the final tail

    t  = 0:dt:tf;
    f  = rev*(1/tf);
    Ax = Amp(1)/2; % meters
    Az = Amp(2)/2; % meters
    
    n = length(t);
    
    x = xz_i(1) - Ax +Ax*sin(2*pi*f*t);
    z = xz_i(2) - Az +Az*cos(2*pi*f*t);

    if dbg
        h = figurew('fe'); axis 'equal'
        plot(x, z, 'b.-'); hold on;

        figurew('xy')
        subplot(2,1,1); hold on;
        plot(t, x)
        subplot(2,1,2); hold on;
        plot(t, z)
    end
    
    
    % creating smooth initial acceleration
    ini = round( (n/rev)/4 );
    x(1:ini) = x(ini)*ones(1,(length(1:ini)) );
    % creating smooth final decceleration
    z(end-ini:end) = z(end-ini)*ones(1,(length(z(end-ini:end))) );    
    
    if dbg
        subplot(2,1,1)
        plot(t, x, 'ro')
        subplot(2,1,2)    
        plot(t, z, 'ro')

        figure(h);
        plot(x, z, 'r.-')
    end
    
    
    if dbg
        figurew('vel')
        subplot(2,1,1)
        plot(t(1:end-1), diff(x)/dt, 'ro')
        subplot(2,1,2)    
        plot(t(1:end-1), diff(z)/dt, 'ro')
    end
    
end
    
    
end
