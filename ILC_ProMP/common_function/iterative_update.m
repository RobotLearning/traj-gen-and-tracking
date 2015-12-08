function [h, xyz, u0, u_ilc, accel, rmse_] = iterative_update(bioRob, h, ...
                         current_iter, ref, prm, u_init_guess)

    % compute the inverse dynamics to recover the actuator input
    % In practice we can measure the actuator input directly, so do not
    % worry that we do not have an inverse dynamics model of the real
    % robot.
    accel = diff(current_iter.qdot)./(ref.t(2)-ref.t(1));
    accel = [accel; accel(end,:)];
    u0    = bioRob.nofriction().rne(current_iter.q, current_iter.qdot, accel);
    u0(:,1) = u0(:,1)-u0(:,1); % force torques for the 1st link to be zero
    
    if exist('u_init_guess', 'var')
        % replace the inverse dynamics solution by the same initial guess
        % at all iterations
        u0 = u_init_guess;
    end    
    
    % do ILC here
    uP = (ref.q    - current_iter.q   )*diag(prm.ilc.P); % proportional learning
    uD = (ref.qdot - current_iter.qdot)*diag(prm.ilc.D); % derivative learning
    u_ilc = u0  + uP + uD;
    
    % first link
    u_ilc(:,1) = 0.*u_ilc(:,1);
    
    % plot and animate
    [h, xyz] = aux(h, ref, bioRob, current_iter, u_ilc, prm);    

    
    % compute the error metric
    rmse_ = error_metric( ref.xyz, xyz);
    

    
end




function [h, xyz] = aux(h, ref, bioRob, current_iter, u1, prm)

    for k=1:length(ref.t)
        T = bioRob.fkine( current_iter.q(k,:) );
        xyz(k,:) = T(1:3,4)';
    end
    
    figure(h.cart);
    plot(xyz(:,1), xyz(:,3), prm.style);
    plot(xyz(end,1), xyz(end,3), SBLUEBALL(5));
    
    
    figure(h.time);
    for j=1:3
        subplot(3,1,j);
        plot(ref.t, current_iter.q(:,j+1), prm.style);        
    end
    
    figure(h.ilc);
    for j=1:3
        subplot(3,1,j);
        plot(ref.t, u1(:,j+1), prm.style); 
    end
    
    
    figure(h.error);
    for j=1:3
        subplot(3,1,j);
        if isfield(h, 'colors')
            plot(ref.t, ref.q(:,j+1)-current_iter.q(:,j+1), styw(h.colors, [], 1.5) );
        else
            plot(ref.t, ref.q(:,j+1)-current_iter.q(:,j+1), prm.style);        
        end
    end
    
    
    % This can be shown in animation also
    if prm.plot.animation
        figure(h.anim);
        title('Current tracking performance');
        animate_bioRob(bioRob, current_iter.q, 10)

    end
    
    
end



