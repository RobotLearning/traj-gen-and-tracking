% Discrete Dynamic Movement Primitive
% 
% NOTES:
% Check 'Learning Motor Primitives for Robotics' for the following changes:
% 1. tau is moved to the right-hand side as opposed to original
%    convention
% 2. yd now is the second dimension, as opposed to first
%
classdef DDMP < DMP
    
    properties
        
        % state of the dmp
        y        
        % canonical system
        can
        % time constants
        alpha_g, beta_g
        % goal state
        goal
        % initial y,yd,ydd values
        y0
        % unscaled states
        z,zd
        % weights of the DMP
        w
        % regularization constant when regressing
        lambda
    end
    
    methods
        
        %% Constructor for discrete DMP
        function obj = DDMP(canonical,alpha,beta,goal,yin)
            
            assert(strcmp(canonical.pattern, 'd'),...
                   'Please provide a discrete canonical system');
            obj.can = canonical;
            obj.alpha_g = alpha;
            obj.beta_g = beta;
            obj.goal = goal;
            obj.lambda = 1e-1;
            assert(length(yin)==3,'please provide initial vel and acc');
            obj.setInitState(yin);
            obj.setWeights(zeros(1,obj.can.nbf));
            obj.z = 0;
            obj.zd = 0;
            % reset all states and phases
            obj.resetStates();
        end
        
                
        %% basis functions are unscaled gaussians
        function out = basis(obj,x)
            out = exp(-obj.can.h .* (x - obj.can.c).^2);
        end        
       
           
        %% one step of the DMP
        function step(obj,err)
           
            alpha = obj.alpha_g;
            beta = obj.beta_g;
            g = obj.goal;            
            amp = 1;
            %amp = g - obj.y0;
            tau = obj.can.tau;
            dt = obj.can.dt;
            
%             f = obj.forcing();
%             obj.zd = (alpha	* (beta * (g - obj.y(1)) - obj.z) + amp * f) * tau;
%         	obj.y(2) = obj.z * tau;
%             obj.y(3) = obj.zd * tau;
% 
%             obj.z = obj.z + dt * obj.zd;
%             obj.y(1) = obj.y(1) + dt * obj.y(2);

            A =  [0, 1;
                -alpha*beta*(tau^2), -alpha*tau];

            f = obj.forcing();
            % forcing function acts on the accelerations
            B = tau * tau * [0; alpha*beta*g + amp*f];

            ydot = A*obj.y(1:2) + B;
            
            % integrate the position and velocity
            obj.y(1) = obj.y(1) + dt * ydot(1);
            obj.y(2) = obj.y(2) + dt * ydot(2);
            % acceleration
            obj.y(3) = ydot(2);
            obj.can.step(err);
        end
        
        %% forcing function to drive nonlinear system dynamics
        function f = forcing(obj)
         
            f = 0;
            scale = 0;
            psi = obj.basis(obj.can.x);
            for i = 1:obj.can.nbf
                f = f + psi(i)*obj.w(i)*obj.can.x;
                scale = scale + psi(i);
            end

            f = f/(scale + obj.can.nbf*1e-10);
        end
        
    end
    
end