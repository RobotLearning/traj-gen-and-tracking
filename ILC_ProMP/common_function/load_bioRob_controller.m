function out = load_bioRob_controller( ref, uInitGuess )
%
% INPUT
%
%   ref: a structure that was created previously and contains the inverse
%        kinematics solution for a reaching task.
%
%   uInitGuess: feedforward input commands. If this input is not given they
%                will be set to zero such that ILC will start from
%                scratch.
%
% OUTPUT
%   out: a structure with fields
%
%            ref: [1x1 struct]
%     controller: [1x1 struct]
%       feedback: [1x1 struct]
%            ilc: [1x1 struct]
%           plot: [1x1 struct]
%            ode: [1x1 struct]
%
%

    out.ref = ref;
    
    % decides if feedback is ON/OFF. Same for the feedforward.
    out.controller   = struct('feedback', 1, 'feedforward', 1);

    % Feedback PD gains that I hand tunned. Do *not* change these values!!
    % We are not interested in improving the performance of the
    % system by feedback but only by the learning feedforward controller.
    % Note that the first joint has zero gains because it should not move.
    out.feedback.P =     3*[0  5  3  2]; 
    out.feedback.D =     2*[0  2  1  1];
    
    % Feedforward ILC update gains. Do *not* change these values. Improvement
    % in performance should come from the ProMP and not by tuning these
    % gains.
    % Note that the first joint has zero gains because it should not move.
    out.ilc.P =          0.01*[0  1  1   1];  
    out.ilc.D =           1.0*[0  1  1   1];   
    

    if exist('uInitGuess', 'var')
        out.ilc.u = uInitGuess;
    else
        % ILC starts with zero torque guess
        out.ilc.u = zeros(length(ref.t), 4);
    end
    

    out.plot.animation = 1;
    
    
    % set this to fine tune ode45. Make it run faster
    out.ode.rel_tol = 0.02;
    out.ode.abs_tol = 0.025;

end