% Class implementing an Extended Kalman filter
% Assuming noise covariances are time invariant
% Assuming observation model is linear

classdef EKF < handle
    
    properties
          
        % Kalman filter covariance of the state
        P
        % state to be estimated
        x
        % function handle for state evolution model to be integrated
        % for time dt
        f
        % linearization of the function f around estimated state
        A
        % observation model g is assumed to be linear and discrete
        C
        % process noise covariance
        Q
        % observation noise covariance
        R
    end
    
    % For robust filtering we need last observation value
    properties
        y_last
    end
    
    methods
        
        %% Initialize variances
        function obj = EKF(dim,funState,mats)
            
            dimx = dim;
            obj.P = eye(dimx);
            obj.x = zeros(dimx,1);
            obj.f = funState;
            
            % initialize covariances
            obj.Q = mats.O;
            obj.R = mats.M;
            
            % initialize linearized matrices
            obj.A = [];
            obj.C = mats.C;
            
            % last observation value
            obj.y_last = zeros(size(obj.C,1),1);
            
        end
        
        %% linearize the dynamics f and g with 'TWO-SIDED-SECANT'
        % make sure you run this before predict and update
        % does not linearize observation model g (assumed linear)
        function linearize(obj,dt,u)
            
            % method for numerical differentiation
            METHOD = 'TWO-SIDED-SECANT'; 
            %METHOD = 'FIVE-POINT-STENCIL'; 
    
            switch METHOD
                case 'FIVE-POINT-STENCIL'
                h = 1e-3;
                n = length(obj.x);
                x2hPlus = repmat(obj.x,1,n) + 2*h * eye(n);
                xhPlus = repmat(obj.x,1,n) + h * eye(n);
                xhMinus = repmat(obj.x,1,n) - h * eye(n);
                x2hMinus = repmat(obj.x,1,n) - 2*h * eye(n);
                xd2hPlus = zeros(n);
                xdhPlus = zeros(n);
                xdhMinus = zeros(n);
                xd2hMinus = zeros(n);
                for i = 1:n
                    xd2hPlus(:,i) = obj.f(x2hPlus(:,i),u,dt);
                    xdhPlus(:,i) = obj.f(xhPlus(:,i),u,dt);
                    xdhMinus(:,i) = obj.f(xhMinus(:,i),u,dt);
                    xd2hMinus(:,i) = obj.f(x2hMinus(:,i),u,dt);
                end
                fder = (-xd2hPlus + 8*xdhPlus - 8*xdhMinus + xd2hMinus) / (12*h);

                case 'TWO-SIDED-SECANT'
                % Take 2n differences with h small
                h = 1e-3;
                n = length(obj.x); 
                xhPlus = repmat(obj.x,1,n) + h * eye(n);
                xhMinus = repmat(obj.x,1,n) - h * eye(n);
                fPlus = zeros(n); 
                fMinus = zeros(n);
                for i = 1:n
                    fPlus(:,i) = obj.f(xhPlus(:,i),u,dt);
                    fMinus(:,i) = obj.f(xhMinus(:,i),u,dt);
                end
                fder = (fPlus - fMinus) / (2*h);
                %dfdx = [zeros(n), eye(n); fder(n+1:2*n,:)];
            end
            obj.A = fder;
        end
        
        %% Initialize state
        function initState(obj,x0,P0)
            obj.x = x0;
            obj.P = P0;
        end
        
        % Set covariances
        function setCov(obj,Q,R)
            obj.Q = Q;
            obj.R = R;
        end
        
        %% predict state and covariance
        % predict dt into the future
        function predict(obj,dt,u)
            
            % integrating the state to get x(t+1)
            obj.x = obj.f(obj.x,u,dt);
            obj.P = obj.A * obj.P * obj.A' + obj.Q;
        end
        
        %% update Kalman filter after observation
        function update(obj,y,u)
            
            % Innovation sequence
            Inno = y(:) - obj.C * obj.x;
            % Innovation (output/residual) covariance
            Theta = obj.C * obj.P * obj.C' + obj.R;
            % Optimal Kalman gain
            K = obj.P * obj.C' * inv(Theta); %#ok
            % update state/variance
            obj.P = obj.P - K * obj.C * obj.P;
            obj.x = obj.x + K * Inno;
        end
        
        % update only if the observation is within a fixed
        % standard deviation of the filter
        %
        % considers also complete outliers where the ball is below
        % the table or above the maximum threshold zMax
        %
        % considers also the case where the ball is not updated yet
        % in this case we will also not update, nor issue warning
        function update_flag = robust_update(obj,y,u)
            
            update_flag = false;
            table_z = -0.95;
            zMin = table_z; 
            zMax = 0.5;
            validBall = y(3) >= zMin && y(3) < zMax;
            tol = 1e-3; % ball should be at least this much different in norm
            new_ball = norm(y - obj.y_last) > tol;
            obj.y_last = y(:);
            
            % standard deviation multiplier
            std_dev_mult = 3;
            % difference in measurement and predicted value
            inno = y(:) - obj.C * obj.x;
            thresh = std_dev_mult * sqrt(diag(obj.C * obj.P * obj.C'));
            
            if any(abs(inno) > thresh) || ~validBall
                % possible outlier not updating
                warning('possible outlier! not updating!');                
            elseif ~new_ball
                warning('nothing new! not updating!');
                % do not do anything
            else
                obj.update(y,u);
                update_flag = true; % update flag
            end    
        end
        
        %% Kalman smoother (batch mode)
        % class must be initialized and initState must be run before
        function [X,V] = smooth(obj,t,Y,U)
            
            N = size(Y,2);
            X = zeros(length(obj.x),N);
            X_pred = zeros(length(obj.x),N-1);
            V = zeros(length(obj.x),length(obj.x),N);
            V_pred = zeros(length(obj.x),length(obj.x),N-1);
            % forward pass
            for i = 1:N-1
                dt = t(i+1) - t(i);
                obj.linearize(dt,U(:,i));
                obj.update(Y(:,i),U(:,i));
                X(:,i) = obj.x;
                V(:,:,i) = obj.P;
                obj.predict(dt,U(:,i));
                X_pred(:,i) = obj.x;
                V_pred(:,:,i) = obj.P;
            end
            obj.linearize(dt,U(:,i));
            obj.update(Y(:,N),0);
            X(:,N) = obj.x;
            V(:,:,N) = obj.P;
            % backward pass
            for i = N-1:-1:1
                % Rauch recursion 
                H = V(:,:,i) * obj.A' * inv(V_pred(:,:,i));
                X(:,i) = X(:,i) + H * (X(:,i+1) - X_pred(:,i));
                V(:,:,i) = V(:,:,i) + H * (V(:,:,i+1) - V_pred(:,:,i)) * H';
            end
        end
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
end