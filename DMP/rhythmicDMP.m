% Rhythmic Dynamic Movement Primitive
% 
% NOTES:
% Check 'Learning Motor Primitives for Robotics' for the following changes:
% 1. tau is moved to the right-hand side as opposed to original
%    convention
% 2. yd now is the second dimension, as opposed to first
%    

classdef rhythmicDMP < DMP
    
    properties
        
        % canonical system
        can
        % time constants
        alpha_g, beta_g
        % goal state
        goal
        % y, yd (or z) values
        Y
        % initial y,yd values
        Y0
        % forcing structure has weights w, widths h, and centers c
        FOR
    end
    
    methods
        
        function obj = discreteDMP(canonicalSystem,alpha,beta,goal,...
                                   amplitude,yin,force)
            
            assert(strcmp(canonicalSystem.pattern, 'r'),...
                   'Please provide a rhythmic canonical system');
            obj.can = canonicalSystem;
            obj.alpha_g = alpha;
            obj.beta_g = beta;
            obj.goal = [goal; amplitude];
            obj.Y0 = yin;
            obj.FOR = force;
            obj.resetStates();
        end
        
        function resetStates(obj)
           
            N = obj.can.N;
            obj.Y = obj.Y0;
            obj.can.reset();
        end
        
        % basis functions are unscaled gaussians
        function out = basis(obj,phi,h,c)
            out = exp(h * (cos(phi - c) - 1));
        end
        
        % evolve is the feedforward rollout function
        % TODO: apply bsxfun or arrayfun
        function [x_roll, Y_roll] = evolve(obj)
            
            Y_roll = zeros(2,obj.can.N);
            x_roll = obj.can.evolve();
            obj.can.reset();
            for i = 1:obj.can.N
                Y_roll(:,i) = obj.Y;
                obj.step(1);
            end            
        end
        
        % one step of the DMP
        function step(obj,err)
           
            alpha = obj.alpha_g;
            beta = obj.beta_g;
            g = obj.goal(1);
            amp = obj.goal(2);
            tauStep = obj.can.tau;
            yin = obj.Y0(1);
            dt = obj.can.dt;

            A = [0, tauStep;
                -alpha*beta*tauStep, -alpha*tauStep];

            f = obj.forcing() * amp;
            % forcing function acts on the accelerations
            B = [0; alpha*beta*g*tauStep + f*tauStep];

            dY = A*obj.Y + B;
            
            obj.Y = obj.Y + dt * dY;
            obj.can.step(err);
        end
        
        % forcing function to drive nonlinear system dynamics
        function f = forcing(obj)

        w = obj.FOR.w;
        h = obj.FOR.h;
        c = obj.FOR.c;
        x = obj.can.x;
        N = length(w);
        f = 0;
        scale = 0;
        for i = 1:N
            f = f + obj.basis(x,h(i),c(i))*w(i);
            scale = scale + obj.basis(x,h(i),c(i));
        end

        f = f/scale;
        end
        
        % learn weights using this function
        % basis functions are fixed
        function learnWeightsFixed(obj)
           
            % TODO:
            error('Not Implemented');
        end
        
    end
    
end