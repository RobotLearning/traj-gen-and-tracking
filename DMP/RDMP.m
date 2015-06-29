% Rhythmic Dynamic Movement Primitive
% 
% NOTES:
% Check 'Learning Motor Primitives for Robotics' for the following changes:
% 1. tau is moved to the right-hand side as opposed to original
%    convention
% 2. yd now is the second dimension, as opposed to first
%    

classdef RDMP < DMP
    
    properties
        
        % state of the dmp
        y
        % canonical system
        can
        % time constants
        alpha_g, beta_g
        % goal state
        goal
        % initial state
        y0
        % weights
        w
        % regularization constant when regressing
        lambda
    end
    
    methods
        
        %% Constructor for the rhythmic DMP
        function obj = RDMP(canonical,alpha,beta,goal,amplitude,yin)
            
            assert(strcmp(canonical.pattern, 'r'),...
                   'Please provide a rhythmic canonical system');
            obj.can = canonical;
            obj.alpha_g = alpha;
            obj.beta_g = beta;
            obj.lambda = 1e-1;
            obj.goal = [goal; amplitude];
            obj.setInitState(yin);
            obj.setWeights(zeros(1,obj.can.nbf));
            obj.resetStates();
        end
        
        %% basis functions are unscaled gaussians
        function out = basis(obj,phi)
            out = exp(obj.can.h .* (cos(phi - obj.can.c) - 1));
        end        
        
        %% one step of the DMP
        function step(obj,err)
           
            alpha = obj.alpha_g;
            beta = obj.beta_g;
            g = obj.goal(1);
            amp = obj.goal(2);
            tauStep = obj.can.tau;
            dt = obj.can.dt;

            A = [0, tauStep;
                -alpha*beta*tauStep, -alpha*tauStep];

            f = obj.forcing();% * amp;
            % forcing function acts on the accelerations
            B = [0; alpha*beta*g*tauStep + f*tauStep];
            dy = A*obj.y(1:2) + B;
            
            % integrate the position and velocity
            obj.y(1:2) = obj.y(1:2) + dt * dy;
            % acceleration
            obj.y(3) = dy(2);
            obj.can.step(err);
        end
        
        %% Forcing function to drive nonlinear system dynamics
        function f = forcing(obj)

            f = 0;
            scale = 0;
            basis = obj.basis(obj.can.x);
            for i = 1:obj.can.nbf
                f = f + basis(i)*obj.w(i);
                scale = scale + basis(i);
            end

            f = f/scale;
            end
        
    end
    
end