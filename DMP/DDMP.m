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
            obj.y0 = yin;
            obj.w = zeros(1,obj.can.nbf);
            % reset all states and phases
            obj.resetStates();
        end
        
                
        %% basis functions are unscaled gaussians
        function out = basis(obj,x)
            out = exp(-obj.can.h .* (x - obj.can.c).^2);
        end        
        
        
        %{
        % constructs the matrix Fs, i.e. s = Fs*w + s_free
        % useful function for ILC
        function Fs = constructF(obj,t)
            
            
            alpha = obj.alpha_g;
            beta = obj.beta_g;
            g = obj.goal;
            tau = obj.can.tau;
            yin = obj.Y0(1);
            dt = obj.can.dt;
            h = obj.FORCE.h;
            c = obj.FORCE.c;            
            
            % construct As
            As = tau * [0, 1, 0; 
                 -alpha*beta, -alpha, alpha*beta; 
                        0, 0, 0];
                    
            % construct Phi
            obj.can.reset();
            N = length(t)-1;
            M = length(obj.FORCE.w);
            Phi = zeros(N,M);
            for i = 1:N
                sumphi = 0;
                for j = 1:M
                    x = obj.can.x;
                    Phi(i,j) = x * obj.basis(x,h(j),c(j));
                    sumphi = sumphi + obj.basis(x,h(j),c(j));
                end
                Phi(i,:) = Phi(i,:)/sumphi;
                obj.can.step(1);
            end
            
            % discretize As and Phi
            % we don't need the g component for relating
            % inputs w to outputs w
            As = eye(3) + dt * As;
            Phi = dt * Phi;
            
            Fs = zeros(N*2,N*M);
            % add the evolution of s0
            %sfree = zeros(N*3,1);
            %sf = s1;
            % construct Fs
            for i = 1:N
                vec_x = (i-1)*2 + 1:i*2;
                for j = 1:i        
                    vec_w = (j-1)*M + 1:j*M;
                    % put zeros in between
                    mat = [zeros(1,M); Phi(j,:)];
                    for k = j+1:i
                        mat = As(1:2,1:2) * mat;
                    end
                    Fs(vec_x,vec_w) = mat; 
                end
                %sf = As * sf;
                %sfree(vec_x) = sf;
            end

            % get the weights
            Fs = Fs * repmat(eye(M),N,1);
            obj.FORCE.Fs = Fs;
        end
        %}
           
        %% one step of the DMP
        function step(obj,err)
           
            alpha = obj.alpha_g;
            beta = obj.beta_g;
            g = obj.goal;            
            amp = 1;
            %amp = g - obj.y0;
            tauStep = obj.can.tau;
            dt = obj.can.dt;

            A = [0, tauStep;
                -alpha*beta*tauStep, -alpha*tauStep];

            f = obj.forcing();
            % forcing function acts on the accelerations
            B = [0; alpha*beta*g*tauStep + amp*f*tauStep];

            dy = A*obj.y(1:2) + B;
            
            % integrate the position and velocity
            obj.y(1:2) = obj.y(1:2) + dt * dy;
            % acceleration
            obj.y(3) = dy(2);
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