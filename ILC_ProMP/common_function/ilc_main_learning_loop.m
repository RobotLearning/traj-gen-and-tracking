

function [iter, prm]  = ilc_main_learning_loop(iter, nIter,...
                                               prm, bioRob, h, ref)

    % usefull if you continue from a given iteration
    if ~isempty(iter)
        ni = size(iter,2)+1; % current size of iter
        nStart = ni;
    else
        ni = 1;
        nStart = 1;
    end
    
    rgb_ = linspace(0.8, 0, nIter); % this is to change the line color as the number of iterations changes
    

    kk =1;
    % main learning loop is here
    for k= nStart:nStart+nIter-1

        fprintf('=========================\n');
        fprintf('Iterations remaining %g.\n', nStart+nIter-1-k);
        if isfield(prm, 'comment')            
            fprintf('%s.\n', prm.comment);
        end

        % line style
        prm.style = struct('Color', [(rgb_(kk))  (rgb_(kk)) (rgb_(kk))], ...
                           'LineStyle', '-', 'LineWidth', 1.5);

        % save the feedforward input of the current iteration
        iter(ni).uilc = prm.ilc.u;                       

        % simulate the robot without Coulomb friction (only viscous) to speed
        % up the simulation time.
        fprintf('Running simulation...\n'); tic;
            [t q qd]     = bioRob.nofriction().fdyn(ref.t, @torque_controller,...
                           ref.q(1,:), [0 0 0 0], prm);
            iter(ni).q    = q;  
            iter(ni).qdot = qd;
        fprintf('... finished in %g seconds\n\n', toc); 


        fprintf('Iterative learning control update.\n');
        fprintf('Animate current tracking performance.\n\n');
        [h, xyz, u, ilc_u, qddot, rmse_metric] = iterative_update(bioRob, h, ...
                                        iter(ni), ref, prm);

        iter(ni).xyz     = xyz;
        iter(ni).u       = u;
        iter(ni).qddot   = qddot;
        iter(ni).rmse    = rmse_metric;


        figure(h.cart_rmse);
        if isfield(h, 'colors')
            plot(k, iter(ni).rmse, sty(h.colors, 'o', 2, [], 10)); 
        else
            plot(k, iter(ni).rmse, sty('b', 'o', 2, [], 10)); 
        end

        % update the ILC for the next iteration 
        prm.ilc.u = ilc_u;

        
        ni = ni+1;
        kk = kk+1;

    end % for loop

end










