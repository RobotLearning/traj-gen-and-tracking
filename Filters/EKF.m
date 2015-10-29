% Class implementing an Extended Kalman filter
% Assuming noise covariances are time invariant
% Assuming observation model is linear

classdef EKF < handle
    
    properties
          
        % Kalman filter covariance of the state
        P
        % state to be estimated
        x
        % function handle for state evolution model
        f
        % function handle for observation model
        g
        % linearization of the function f around estimated state
        A
        % linearization of the function g around estimated state
        C
        % process noise covariance
        Q
        % observation noise covariance
        R
    end
    
    methods
        
        %% Initialize variances
        function obj = EKF(dim,funState,mats)
            
            dimx = dim;
            obj.P = eye(dimx);
            obj.x = zeros(dimx,1);
            obj.f = funState;
            obj.g = []; % not used for now
            
            % initialize covariances
            obj.Q = mats.O;
            obj.R = mats.M;
            
            % initialize linearized matrices
            obj.A = [];
            obj.C = mats.C;
            
        end
        
        %% linearize the dynamics f and g with 'TWO-SIDED-SECANT'
        % make sure you run this before predict and update
        % if linObs flag is 1, do not linearize observation model g
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
        
        %% predict state and covariance
        % predict dt into the future
        function predict(obj,dt,u)
            
            % integrating the state to get x(t+1)
            obj.x = obj.f(obj.x,u,dt);
            obj.P = obj.A * obj.P * obj.A' + obj.Q;
        end
        
        %% update Kalman filter variance
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
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
end