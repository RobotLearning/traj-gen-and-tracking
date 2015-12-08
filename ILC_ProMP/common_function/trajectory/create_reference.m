function ref = create_reference(bioRob, dt, refname)

    dbg = 0; % 1: plot some figures for debugging. Keep it zero.

    % ============================================
    % Create the Cartesian reference and its joint trajectories
    % ============================================
    
    % Retrieve the structure with arguments to create the trajectory
    % ============================================    
    str = get_ref_arguments(refname);

    
    % Create the reference trajectory based on its type
    % ============================================        
    if strcmp(str.type, 'p2p')
        [t, x, z, a] = reaching(str.xz_i, str.xz_f, str.a_i, str.a_f, dt, str.tf);
        [qi, T] = do_IK(bioRob, t, x, z, a);
    end
    
    if strcmp(str.type, 'ellipse')
        [t, x, z] = circ(str.xz_i, str.Amp, dt, str.tf, str.rev);
        [qi, T] = do_IK(bioRob, t, x, z);
    end

    if strcmp(str.type, 'square')
        [t, x, z] = square(str.xz_i, str.Amp, dt, str.tf, str.rev);
        [qi, T] = do_IK(bioRob, t, x, z);
    end
    
    
    ref.param = str;
    ref.param.dt = dt;
    ref.name = refname;
   
    if dbg
        dfv = [0 -1 0];
        try   close(htemp);   end
        htemp = figurew('animation_debug'); view(dfv);
        drawnow; axis 'equal'; drawnow;
        set_fig_position( [0.389 0.41 0.343 0.568]   );
        plot3(x, 0.*x, z, SRED(2));
        xlim([0      0.7]);
        ylim([-0.05   0.05]);
        zlim([-0.10  0.4]);        
        bioRob.plot(qi);
    end
        
    
    ref.t  = t;
    ref.T  = T;
    ref.q  = qi;
    ref.xyz = [x'  zeros(length(x),1) z'];

    % ============================================
    % Compute the ideal inv. dynamics command as ground truth
    % ============================================
    [ref.u, ref.qdot]  = compute_inv_dyn(bioRob, ref.q, dt );
    
end




