function [str] = get_ref_arguments(refname)


    % common parameters for reaching task
    reachingTaskTf = 3; % seconds to reach, only for p2p task
    squareSideLength  = 0.25; %
    p2pi = [0.2 -0.1];
    p2pf = [0.3  0.1];
    p2pAnglei = d2r(60);
    p2pAnglef = d2r(0);
    
    
    
    % this setting was quite challenging
    switch refname
        
        % =================================
        % Continuous tracking with circles
        % =================================        
        case 'eA' %'up'
        str.type  = 'ellipse';            
        str.xz_i  = [0.4  0.4];    % absolute starting position
        str.Amp   = 1.5*[0.1  0.2]; % Amplitude of the ellipse in X and Z direction
        str.tf    = 10;            % Final time
        str.rev   = 2;             % Number of revolutions
    
        % make a test case that is in between these two cases
        case 'eB' %'mid'
        str.type  = 'ellipse';            
        str.xz_i  = [0.5  0.25];    % absolute starting position
        str.Amp   = 1.5*[0.15  0.15]; % Amplitude of the ellipse in X and Z direction
        str.tf    = 10;            % Final time
        str.rev   = 2;             % Number of revolutions
    
    
        % Define the reference ellipsoid trajectory here:
        case 'eC' % 'down'
        str.type  = 'ellipse';
        str.xz_i  = [0.55  0.15];    % absolute starting position
        str.Amp   = 1.5*[0.2  0.1]; % Amplitude of the ellipse in X and Z direction
        str.tf    = 10;            % Final time
        str.rev   = 2;             % Number of revolutions
        
        case 'eD' %'invUp'
        str.type  = 'ellipse';
        str.xz_i  = [0.4  -0.05];    % absolute starting position
         str.Amp   = 1.5*[0.1  0.2]; % Amplitude of the ellipse in X and Z direction
        str.tf    = 10;            % Final time
        str.rev   = 2;             % Number of revolutions
        
        case 'eE' %'invMid'
        str.type  = 'ellipse';
        str.xz_i  = [0.5  0.05];    % absolute starting position
        str.Amp   = 1.5*[0.15  0.15]; % Amplitude of the ellipse in X and Z direction
        str.tf    = 10;            % Final time
        str.rev   = 2;             % Number of revolutions
        
        % =================================
        % Reaching task (aka point-to-point)
        % =================================
        case 'rA' % southWest
        str.type  = 'p2p';
        str.xz_i  = p2pi;    % absolute starting position
        str.xz_f  = p2pf;    % absolute final position
        str.tf    = reachingTaskTf;              % Final time
        str.a_i   = p2pAnglei; % initial orientation of end-effector
        str.a_f   = p2pAnglef; % final orientation of end-effector
        
        case 'rB' % northWest
        str.type  = 'p2p';
        str.xz_i  = p2pi;    % absolute starting position
        str.xz_f  = [p2pf(1)  p2pf(2)+squareSideLength];    % absolute final position
        str.tf    = reachingTaskTf;
        str.a_i   = p2pAnglei; % initial orientation of end-effector
        str.a_f   = p2pAnglef; % final orientation of end-effector        
        
        case 'rC' % northEast
        str.type  = 'p2p';
        str.xz_i  = p2pi;    % absolute starting position
        str.xz_f  = [p2pf(1)+squareSideLength   p2pf(2)+squareSideLength];  % absolute final position
        str.tf    = reachingTaskTf;           
        str.a_i   = p2pAnglei; % initial orientation of end-effector
        str.a_f   = p2pAnglef; % final orientation of end-effector        
        
        case 'rD' % southEast
        str.type  = 'p2p';
        str.xz_i  = p2pi;    % absolute starting position
        str.xz_f  = [p2pf(1)+squareSideLength   p2pf(2)];  % absolute final position
        str.tf    = reachingTaskTf;           
        str.a_i   = p2pAnglei; % initial orientation of end-effector
        str.a_f   = p2pAnglef; % final orientation of end-effector        
        
        case 'rE' % middle
        str.type  = 'p2p';
        str.xz_i  = p2pi;    % absolute starting position
        str.xz_f  = [p2pf(1)+0.5*squareSideLength  p2pf(2)+0.5*squareSideLength];  % absolute final position
        str.tf    = reachingTaskTf;                   
        str.a_i   = p2pAnglei; % initial orientation of end-effector
        str.a_f   = p2pAnglef; % final orientation of end-effector        
        
        % =================================
        % Square task: TO DO
        % =================================
        case 'sA' % NOT IMPLEMENTED YET
        str.type  = 'square';
        str.xz_i  = [0.55  0.15];    % absolute starting position
        str.Amp   = 1.5*[0.2  0.1]; % Amplitude of the ellipse in X and Z direction
        str.tf    = 10;            % Final time
        str.rev   = 2;             % Number of revolutions
        
        otherwise
            error('up, mid, down');
    end
    
    
end


















