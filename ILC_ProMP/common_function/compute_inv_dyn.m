function [u0, qdot] = compute_inv_dyn(bioRob, q, dt)

    qdot = diff(q)./dt;
    qdot = [qdot ; qdot(end,:)];
    
    accel = diff(qdot)./(dt);
    accel = [accel; accel(end,:)];
    
    u0    = bioRob.nofriction().rne(q, qdot, accel);
    u0(:,1) = u0(:,1)-u0(:,1); % force torques for the 1st link to be zero    
    
    
end