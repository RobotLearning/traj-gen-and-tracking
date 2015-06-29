% Dynamic motor primitive superclass
% rhythmic and discrete DMPs inherit from this class
%
% TODO: regression should include tau^2 and tau in the forcing function
% extraction

classdef (Abstract) DMP < handle
    
    properties (Abstract)
        
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
        % weights
        w
        % regularization constant when regressing
        lambda
    end
    
    methods (Abstract)

        % discrete and rhythmic subclasses define their own basis fnc
        basis(phi)
        
        % forcing function
        forcing(obj)
        
        % make a step - both for feedback/feedforward simulation
        step(obj)
        
    end
    
    % methods that can be implemented here in abstract class
    methods (Access = public)
                
        %% Evolve the DMP by N steps
        % feedforward simulation of DMPs
        function [x_roll, Y_roll] = evolve(obj,N)
            
            Y_roll = zeros(3,N);
            x_roll = obj.can.evolve(N);
            obj.resetStates();
            for i = 1:N
                Y_roll(:,i) = obj.y(:);
                obj.step(1);
            end            
        end
        
        %% Set and reset functions
        
        % phase of the canonical system is set to 0
        % as well as clearing Y values
        function resetStates(obj)
           
            obj.y = obj.y0;
            obj.can.reset();
        end
        
        % set the goal state from observational data
        function [g,scale] = setGoal(obj,path)
            
            if strcmp(obj.can.pattern,'d')
                obj.goal = path(end);            
                g = obj.goal;
                scale = g - obj.y0(1);
            else            
                % goal state is the center
                obj.goal(1) = (min(path) + max(path)) / 2;
                % amplitude is the difference
                obj.goal(2) = max(path) - obj.goal(1);

                g = obj.goal(1);
                scale = obj.goal(2);
            end
        end
        
         % set the weights of the forcing function after regression for instance
        function setWeights(obj,ww)
            
            assert(obj.can.nbf == length(ww),'Length of bfs do not match!');
            obj.w = ww;
        end
        
        % set the initial state of the DMP
        function setInitState(obj,yin)
            
            assert(length(yin)==3,'please provide initial vel and acc');
            obj.y0 = yin;
        end
        
        %% Regression methods here        
        
        % update the weights using regression methods
        function updateWeights(obj,path)
            
            if strcmp(obj.can.pattern,'d')
                g = path(end);
            else
                g = (min(path) + max(path)) / 2;
            end
                
            % learn the weights with locally weighted regression
            %obj.LWR(g,path);
            % or with linear regression
            obj.regression(g,path);
            
            
        end
        
        % Regression on processed actual demonstrations
        %
        % INPUTS:
        %
        % q, dq, ddq - joint positions, velocities, acc.
        %
        function regressLive(obj,q,qd,qdd,goals)

            % make sure phase is reset
            obj.resetStates();
            
            pat = obj.can.pattern;
            lenw = obj.can.nbf;
            alpha = obj.alpha_g;
            beta = obj.beta_g;
            
            % number of demonstrations
            D = length(goals);
            N = length(q)/D;
            goals = repmat(goals',N,1);
            goals = goals(:);

            % calculate fd
            fd = qdd - alpha * (beta * (goals - q) - qd);

            % number of weights to regress 
            x = obj.can.evolve(N);
            % make sure x is column vector
            x = x(:);

            % regress on the weights
            Psi = zeros(N,lenw);
            for i = 1:N
                Psi(i,:) = obj.basis(x(i));
            end
            % scale the psi matrices
            if strcmp(pat,'d')
                scale = x ./ (sum(Psi,2) + 1e-10); 
            else
                scale = 1 ./ (sum(Psi,2) + 1e-10);
            end
            scale = repmat(scale,1,lenw);
            Psi = Psi .* scale;
            
            % in case there are multiple demonstrations
            Psi = repmat(Psi,D,1);
            
            % use pinv or add lambda to smoothen inverse
            %w = pinv(Psi) * fd(:);
            %w = Psi \ fd(:);
            
            % penalize the first weights
            C = diag(ones(1,lenw))-2*diag(ones(1,lenw-1),1)+diag(ones(1,lenw-2),2);
            C(end-1:end,:) = 0;
            %D = -diag(ones(1,lenw))+diag(ones(1,lenw-1),1);
            %D(end,:) = 0;
            %lambda = 1e1;
            w = ((Psi' * Psi + obj.lambda*C) \ (Psi')) * fd(:);

            obj.setWeights(w);

        end
        
        
        % Basic regression
        % To learn the weights of DMPs
        % Assuming fixed centers and covariances
        %
        % INPUTS:
        %
        % path - desired trajectory
        % dmp  - includes the canonical system class whose phase 
        %        weights the regression
        %        includes the forcing term
        %        forcing term is a structure which has fixed centers and covariances
        % 
        % OUTPUTS:
        %
        % forcing structure with the forcing weights learned
        function regression(obj,goal,path)

            % make sure phase is reset
            obj.resetStates();
            
            dt = obj.can.dt;
            pat = obj.can.pattern;
            lenw = obj.can.nbf;
            alpha = obj.alpha_g;
            beta = obj.beta_g;

            % TODO: interpolate over trajectory
            y_des = path;
            yd_des = [0,diff(path)]/dt;
            ydd_des = [0,diff(yd_des)]/dt;
            % calculate fd
            fd = ydd_des - alpha * (beta * (goal - y_des) - yd_des);

            % number of weights to regress 
            x = obj.can.evolve(length(path));
            % make sure x is column vector
            x = x(:);

            % regress on the weights
            lent = length(fd);
            Psi = zeros(lent,lenw);
            for i = 1:lent
                Psi(i,:) = obj.basis(x(i));
            end
            % scale the psi matrices
            if strcmp(pat,'d')
                scale = x ./ (sum(Psi,2) + lenw*1e-10);
            else
                scale = 1 ./ (sum(Psi,2) + lenw*1e-10);
            end
            scale = repmat(scale,1,lenw);
            Psi = Psi .* scale;
            
            %w = pinv(Psi) * fd(:);
            %w = Psi \ fd(:);
            % add lambda to smoothen inverse
            w = ((Psi' * Psi + obj.lambda*eye(lenw)) \ (Psi')) * fd(:);
            
            obj.setWeights(w);

        end
        
        % Locally weighted regression
        % To learn the weights of DMPs
        % Assuming fixed centers and covariances
        %
        % INPUTS:
        %
        % path - desired trajectory
        % dmp  - includes the canonical system class whose phase 
        %        weights the regression
        %        includes the forcing term
        %        forcing term is a structure which has fixed centers and covariances
        % 
        % OUTPUTS:
        %
        % dmp with the forcing weights learned
        function LWR(obj,goal,path)

            % make sure phase is reset
            obj.resetStates();
            
            dt = obj.can.dt;
            pat = obj.can.pattern;
            alpha = obj.alpha_g;
            beta = obj.beta_g;

            % TODO: interpolate over trajectory
            y_des = path;
            yd_des = [0,diff(path)]/dt;
            ydd_des = [0,diff(yd_des)]/dt;
            % calculate fd
            fd = ydd_des - alpha * (beta * (goal - y_des) - yd_des);

            % number of weights to regress 
            lenw = obj.can.nbf;
            x = obj.can.evolve(length(path));

            % make sure x is column vector
            x = x(:);

            % construct weights
            w = zeros(lenw,1);
            
            lent = length(fd);
            Psi = zeros(lent,lenw);
            for i = 1:lent
                Psi(i,:) = obj.basis(x(i));
            end
            
            for i = 1:lenw
                psi = Psi(:,i);
                if strcmp(pat,'d')
                    num = x' * diag(psi) * fd(:);
                    denom = x' * diag(psi) * x;
                else
                    num = psi' * fd(:);
                    denom = sum(psi) + lenw * 1e-10;
                end
                w(i) = num / denom;
            end

            obj.setWeights(w);
        end
    end
    
end