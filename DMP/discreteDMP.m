% Discrete Dynamic Movement Primitive
% 
% NOTES:
% Check 'Learning Motor Primitives for Robotics' for the following changes:
% 1. tau is moved to the right-hand side as opposed to original
%    convention
% 2. yd now is the second dimension, as opposed to first
%    
%
classdef discreteDMP < DMP
    
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
        % field containing Psi matrix
        Psi
        % forcing structure has weights w, widths h, and centers c
        FORCE
    end
    
    methods
        
        function obj = discreteDMP(canonical,alpha,beta,goal,yin,bfs)
            
            assert(strcmp(canonical.pattern, 'd'),...
                   'Please provide a discrete canonical system');
            obj.can = canonical;
            obj.alpha_g = alpha;
            obj.beta_g = beta;
            obj.goal = goal;
            obj.Y0 = yin;
            % get the last phase value
            xtr = obj.can.evolve();
            % initialize forcing function here
            obj.FORCE.w = zeros(bfs,1);
            obj.FORCE.h = ones(bfs,1) * bfs^(1.5);
            obj.FORCE.c = linspace(xtr(end),1,bfs);
            % reset all states and phases
            obj.resetStates();
        end
        
        function resetStates(obj)
           
            %N = obj.can.N;
            obj.Y = obj.Y0;
            obj.can.reset();
        end
        
        function [g,scale] = setGoal(obj,path)
            
            obj.goal = path(end);
            
            g = obj.goal;
            scale = g - obj.Y0(1);
        end
        
        function setForcing(obj,FOR)
            
            assert(isfield(FOR,'w'),'Please perform regression first');
            obj.FORCE = FOR;            
        end
        
        function setInitState(obj,y0)
            
            obj.Y0 = y0;
        end
        
        % basis functions are unscaled gaussians
        function out = basis(obj,x,h,c)
        out = exp(-h * (x - c).^2);
        end
        
        % evolve is the feedforward rollout function
        % TODO: apply bsxfun or arrayfun
        function [x_roll, Y_roll] = evolve(obj)
            
            Y_roll = zeros(2,obj.can.N);
            x_roll = obj.can.evolve();
            obj.resetStates();
            for i = 1:obj.can.N
                Y_roll(:,i) = obj.Y;
                obj.step(1);
            end            
        end
        
        % one step of the DMP
        function step(obj,err)
           
            alpha = obj.alpha_g;
            beta = obj.beta_g;
            g = obj.goal;
            tauStep = obj.can.tau;
            yin = obj.Y0(1);
            dt = obj.can.dt;

            A = [0, tauStep;
                -alpha*beta*tauStep, -alpha*tauStep];

            f = obj.forcing();
            % forcing function acts on the accelerations
            B = [0; alpha*beta*g*tauStep + f*tauStep];

            dY = A*obj.Y + B;
            
            obj.Y = obj.Y + dt * dY;
            obj.can.step(err);
        end
        
        % forcing function to drive nonlinear system dynamics
        function f = forcing(obj)

            w = obj.FORCE.w;
            h = obj.FORCE.h;
            c = obj.FORCE.c;
            x = obj.can.x;
            M = length(w);
            f = 0;
            scale = 0;
            for i = 1:M
                f = f + obj.basis(x,h(i),c(i))*w(i)*x;
                scale = scale + obj.basis(x,h(i),c(i));
            end

            f = f/scale;
        end
        
    end
    
end