% Dynamic motor primitive superclass
% rhythmic and discrete DMPs inherit from this class

classdef (Abstract) DMP < handle
    
    properties (Abstract)
        
        % canonical system
        can
        % time constants
        alpha_g, beta_g
        % goal state
        goal
        % y, yd (or z), ydd values
        Y
        % forcing structure has weights w, widths h, and centers c
        FORCE
    end
    
    methods (Abstract)

        % discrete and rhythmic subclasses define their own basis fnc
        basis(phi,h,c)
        
        % forcing function
        forcing(obj)
        
        % make a step - both for feedback/feedforward simulation
        step(obj)
        
        % feedforward simulation of DMPs
        evolve(obj)
        
        % phase of the canonical system is set to 0
        % as well as clearing Y values
        resetStates(obj)
        
        % set the initial state of the DMP
        setInitState(obj,Y0)
        
        % set the goal state from observational data
        setGoal(obj,path)
        
        % set the forcing function after regression for instance
        setForcing(obj,force)
        
    end
    
    % methods that can be implemented here in abstract class
    methods (Access = public)
        
        % set the weights using regression methods
        function setWeights(obj,path)
            
            [goal,~] = obj.setGoal(path);
            % learn the weights with locally weighted regression
            %force = obj.LWR(goal,path);
            % or with linear regression
            force = obj.regression(goal,path);
            
            % update dmps
            obj.setForcing(force);
        end
        
        % update the weights using regression methods
        function force = updateWeights(obj,path)
            
            if strcmp(obj.can.pattern,'d')
                g = path(end);
            else
                g = min(path) + max(path) / 2;
            end
                
            % learn the weights with locally weighted regression
            force = obj.LWR(g,path);
            % or with linear regression
            %force = obj.regression(g,path);
            
            % do not update the dmp yet!
            % obj.setForcing(force);
            
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
        function force = regression(obj,goal,path)

            dt = obj.can.dt;
            pat = obj.can.pattern;
            force = obj.FORCE;
            alpha = obj.alpha_g;
            beta = obj.beta_g;

            % TODO: interpolate over trajectory
            y_des = path;
            yd_des = [0,diff(path)]/dt;
            ydd_des = [0,diff(yd_des)]/dt;
            % calculate fd
            fd = ydd_des - alpha * (beta * (goal - y_des) - yd_des);

            h = force.h;
            c = force.c;
            % number of weights to regress 
            lenw = length(force.c);
            x = obj.can.evolve();

            % make sure x is column vector
            x = x(:);

            % regress on the weights
            lent = length(fd);
            Psi = zeros(lent,lenw);
            for i = 1:lenw
                Psi(:,i) = obj.basis(x,h(i),c(i));
            end
            % scale the psi matrices
            if strcmp(pat,'d')
                scale = x ./ sum(Psi,2); 
            else
                scale = 1 ./ (sum(Psi,2) + 1e-10);
            end
            scale = repmat(scale,1,lenw);
            obj.Psi = Psi .* scale;
            
            % TODO: use pinv or add lambda to smoothen inverse
            w = pinv(obj.Psi) * fd(:);
            %w = Psi \ fd(:);
            force.w = w;

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
        function force = LWR(obj,goal,path)

            dt = obj.can.dt;
            pat = obj.can.pattern;
            force = obj.FORCE;
            alpha = obj.alpha_g;
            beta = obj.beta_g;

            % TODO: interpolate over trajectory
            y_des = path;
            yd_des = [0,diff(path)]/dt;
            ydd_des = [0,diff(yd_des)]/dt;
            % calculate fd
            fd = ydd_des - alpha * (beta * (goal - y_des) - yd_des);

            h = force.h;
            c = force.c;
            % number of weights to regress 
            lenw = length(force.c);
            x = obj.can.evolve();

            % make sure x is column vector
            x = x(:);

            % construct weights
            w = zeros(lenw,1);
            for i = 1:lenw
                psi = obj.basis(x,h(i),c(i));
                if strcmp(pat,'d')
                    num = x' * diag(psi) * fd(:);
                    denom = x' * diag(psi) * x;
                else
                    num = psi' * fd(:);
                    denom = sum(psi) + 1e-10;
                end
                w(i) = num / denom;
            end

            force.w = w;
        end
    end
    
end