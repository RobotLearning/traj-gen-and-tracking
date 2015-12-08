function ref = create_circular_reference(bioRob, dt, refname)

    dbg = 0; % 1: plot some figures for debugging. Keep it zero.

  
    [xz_i, Amp, tf, rev] = get_ref_arguments(refname);
    
    [t, x, z] = circ(xz_i, Amp, dt, tf, rev);
    
    ref.param.xz_i = xz_i;
    ref.param.Amp = Amp;
    ref.param.tf  = tf;
    ref.param.dt = dt;
    ref.param.rev = rev;
    


    angle = 0*(-1+cos(2*pi*0.25*t))+0*25;
    if dbg
        figurew('wer')
        plot(angle)
    end
    
    [qi, T] = do_IK(bioRob, t, x, z);
    
   
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

    [ref.u, ref.qdot]  = compute_inv_dyn(bioRob, ref.q, dt );
    
end




