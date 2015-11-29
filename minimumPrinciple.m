%% Testing minimum principle for table tennis

% Using one demonstration
% 1. Estimate initial pos and vel with kalman smoother up to bounce
% 2. Predict the ball trajectory after bounce
% 3. Fix the center of the table as the goal position
%    Solve BVP problems (bvp4c) by sampling from the predicted traj.
% 4. Using the contact model estimate racket vels. and feasible
% orientations and interpolate through whole traj.
%
% 5. Using an optimal control toolbox
% State the initial position of the robot
% State the input penalties
% State the transversality conditions
%
% Transversality conditions:
% - positions of the ball should touch the racket (tubular region around
%   the path)
% - velocities of the racket are fixed
% - orientation of the racket could be fixed (although there is some
% freedom)
% - only using kinematics and jacobian
%
% 6. Compare with existing virtual hitting plane method
% - construct VHP
% - estimate pos and vel. as well as desired vel.
% - use inverse kinematics and polynomials to construct traj.
%
% 7. Evolve both with the nominal model (inverse dynamics) and compare
%   probability of returning the ball given actual ball velocities 
%   and/or control errors