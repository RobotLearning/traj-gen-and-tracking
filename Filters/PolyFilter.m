% Class implementing an simple moving average filter
% To estimate position and velocity of a moving target

classdef PolyFilter < handle
    
    properties
          
        % Number of observations to average over
        size
        % Number of observations so far
        numobs
        % Observations to average over
        Y
        % Times of the observations
        T
        % Vandermonde matrix
        V
        % Order of the polynomial
        order
        % Prediction model
        f
        % State (position and velocity)
        x
        % current time
        time
    end
    
    % properties for robust filtering
    properties
        % Coefficients of the polynomial
        a
    end
    
    methods
        
        %% Initialize filtersize
        function obj = PolyFilter(dim,order,filter_size,predict_model)
            
            obj.order = order;
            obj.numobs = 0;
            obj.size = filter_size;
            obj.x = zeros(2*dim,1);
            obj.a = zeros(order+1,dim);
            obj.f = predict_model;
            obj.Y = [];
            obj.V = [];
            obj.time = 0.0;
            obj.T = [];
        end
 
        
        %% Initialize state
        function initState(obj,x0,eps)
            assert(length(x0) == 2*size(obj.Y,2),'dim size does not match with initialization!');
            assert(length(x0) == length(obj.x),'x0 size does not match with state!');
            obj.x = x0;
            obj.T = zeros(obj.size,1);
            obj.Y(1,:) = x0(1:size(obj.Y,2));
            obj.numobs = 1;
            obj.time = 0.0;
        end
        
        %% predict state
        % predict dt into the future
        % this step is not model based
        function xpred = predict(obj,dt)
            
            % integrating the state to get x(t+1)
            acc = obj.a(3,:);
            velPred = obj.x(4:6) + acc(:) * dt;
            xpred = obj.x(1:3) + velPred * dt;
        end
        
        %% update state from observation
        
        % delta_t is the time passed from last observation time
        function update(obj,dt,y)
            
            obj.check_rebound(dt);
            obj.time = obj.time + dt;
            if obj.numobs < obj.size 
                obj.numobs = obj.numobs + 1;
                obj.Y(obj.numobs,:) = y(:)';
                obj.T(obj.numobs,1) = obj.time;
                obj.x(1:length(y)) = y(:);
                if obj.numobs <= obj.size && obj.numobs > 1
                    obj.construct_vandermonde();
                    obj.estimate_state();
                end
            elseif obj.numobs == obj.size
                % regress over balls 
                % construct the matrix
                obj.Y = [obj.Y(2:end,:); y(:)'];
                obj.update_vandermonde();
                obj.estimate_state();
            else
                error('this shouldnt happen!');
            end
        end
        
        % construct v-matrix only once
        function construct_vandermonde(obj)
            
            obj.V = zeros(obj.numobs,obj.order+1);
            %assert(obj.numobs == obj.size, 'numobs should be equal to maxsize');
            for i = 0:1:obj.order
                obj.V(:,i+1) = (obj.T).^i;
            end
        end
        
        % update the vandermonde matrix
        function update_vandermonde(obj)
            
            assert(obj.numobs == obj.size, 'numobs should be equal to maxsize');

            vec = zeros(1,obj.order+1);
            for i = 0:1:obj.order
                vec(i+1) = obj.time^i;
            end
            obj.V = [obj.V(2:end,:);vec];
            obj.T = [obj.T(2:end); obj.time];
        end
        
        % Estimate the state (position and velocities)
        % from the last X observations 
        % after having constructed the vandermonde matrix
        function estimate_state(obj)
             obj.a = obj.V \ obj.Y;
             pos = obj.V(end,:) * obj.a;
             
             vec = zeros(1,obj.order+1);
             for i = 1:obj.order
                 vec(i+1) = factorial(i) * obj.time^(i-1);
             end
             vel = vec * obj.a;
             
             obj.x = [pos(:);vel(:)];
        end
        
        %% Robustification of filter
        % Correct for predicted bounce
        function check_rebound(obj,dt)
            
            % try to predict dt seconds into the future
            % if the ball passes the table then there
            % will be a rebound
            % in this case vandermonde matrix has to be reset            
            
            % Get first filtered position
            xpred = obj.predict(dt);
            table_z = -0.93;
            if xpred(3) < table_z
                loadTennisTableValues;
                disp('Bounce predicted');
                disp('Resetting...');
                obj.numobs = 0;
                %M = diag([CFTX; CFTY; -CRT]);
                %obj.x(4:6) = M * obj.x(4:6);
                spin = [-50*2*pi;0;0];
                vels = [obj.x(4:6); spin(:)];
                vels = reboundSpinModel(vels,CRT,mu,ball_radius);
                obj.x(4:6) = vels(1:3);
                obj.Y = [];
                obj.V = [];
                obj.T = [];
            end

        end
        
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
end