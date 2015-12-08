
function [qi,T] = do_IK(bioRob, t, x, z, endEffOrient)

    % do inverse kinematics to get the joint angles (the robot is
    % controlled in joint space)
    M  = [1 1 1   0 1 0];
    q0 = [-0.0000  -63.3438  127.2693  -63.9255]; % initial guess
    
    
    if exist('endEffOrient', 'var')
        a = endEffOrient;
    else
        a = zeros(1,numel(t));
    end
    
    
    fprintf('Computing IK...\n'); tic;
    for k=1:1:length(t) 
        
        T(:,:,k) = transl(x(k), 0, z(k))*troty( a(k) )*trotz(-angle(k).*pi/180 ) ;
        %T(:,:,k) = transl(x(k), 0, z(k)) * trotz(-angle(k).*pi/180) ;
        
        qi(k,:)  = bioRob.ikineBioRobHacked(T(:,:,k), q0, M, 'pinv');
        q0 = qi(k,:); % use previous joint values as initial guess for ikine
    end
    fprintf('finished IK in %g seconds.\n', toc);
    
    if 1
        qi = wrapToPi(qi);
        qi = unwrap(qi);
    end
    
end



