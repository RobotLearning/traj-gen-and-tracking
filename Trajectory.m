% Class for holding inputs and states on trajectories.
% Class is used instead of a structure to enforce certain fields.
%

classdef Trajectory < handle
    
    properties
        
        % time profile
        t
        % number of discretizations
        N
        % reference trajectory (dimension = output)
        s
        % nominal u values calculated during trajectory generation
        % (optional)
        unom
        % linear feedback law calculated during trajectory generation
        % (optional)
        K
        % via points used for spline-based trajectory generation
        % (optional)
        % sp
        % particular algorithm's performance as array of structures
        PERF
        % interpolation method when upsampling
        upsampling
    end
    
    methods (Access = public)
         
        %% create a trajectory class useful for simulations
        function obj = Trajectory(t,s,unom,K)            
            obj.t = t;
            obj.N = length(t);
            obj.s = s;
            obj.unom = unom;
            obj.K = K;
            obj.upsampling = 'linear';
        end
        
        %% Add performances of learning algorithms
        % fun_cost is a cost function, not necessarily quadratic
        % controller is generally a particular ILC implementation
        function addPerformance(obj,u,y,costfun,controller)
            
            if ischar(controller)
                name = controller;
            else
                name = controller.name;
            end
            
            % downsample u,y if necessary
            n = size(y,2);
            if n ~= obj.N
                rate = n / obj.N;
                idx = rate * (1:obj.N);
                u = u(:,idx(1:end-1));
                y = y(:,idx);
            end
            
            % append performance
            i = length(obj.PERF);
            obj.PERF(i+1).name = name;
            obj.PERF(i+1).u = u;
            obj.PERF(i+1).y = y;
            obj.PERF(i+1).err = y - obj.s;
            obj.PERF(i+1).cost = costfun.fnc(y,obj.s);

            % display SSE error
            sse = sum(obj.PERF(i+1).cost);
            finCost = costfun.fnc(y(:,end),obj.s(:,end));
            rms = sqrt(sse/obj.N);
            fprintf('RMS for %s is %f \n', name, rms);
            fprintf('Final cost is %f \n', finCost);
            
            if ~ischar(controller) && ~isempty(u)
                if length(u) ~= length(controller.inp_last)
                    ds = controller.downsample;
                    idx = ds * (1:obj.N/ds);
                    u = u(:,idx(1:end-1));
                end
                controller.record(u,rms);
            end
            
        end
        
        % inverse project output trajectory s back to state space
        function sbar = projectBack(obj,C)
            
            sbar = C'*((C*C')\obj.s);
        end
        
        %% Upsample/downsample functions 
        
        function trajDown = downsample(obj,rate)
            
            if rate == 1
                trajDown = obj;
                %trajDown = Trajectory(obj.t,obj.s,obj.unom,obj.K);
                %trajDown.PERF = obj.PERF;
                return;
            end
            
            N = obj.N/rate;
            idx = rate * (1:N);
            t = obj.t(idx);
            s = obj.s(:,idx);
            unom = obj.unom(:,idx(1:end-1));
            
            if ~isempty(obj.K)
                K = obj.K(:,:,idx(1:end-1));
            else
                K = [];
            end
            
            trajDown = Trajectory(t,s,unom,K);
            
        end
        
        function trajUp = upsample(obj,rate)
            
            if rate == 1
                trajUp = Trajectory(obj.t,obj.s,obj.unom,obj.K);
                trajUp.PERF = obj.PERF;
                return;
            end
            
            N = round(obj.N * rate);
            t = linspace(obj.t(1)/rate,obj.t(end),N);
            s = interp1(obj.t,obj.s',t,obj.upsampling,'extrap')';
            unom = interp1(obj.t(1:end-1),obj.unom',t(1:end-1),obj.upsampling,'extrap')';
            
            if ~isempty(obj.K)
                % form 2d K
                K = reshape(obj.K,[],obj.N - 1);
                K = interp1(obj.t(1:end-1),K',t,obj.upsampling,'extrap')';
                K = reshape(K,[size(obj.K,1),size(obj.K,2),N]);                
            else
                K = [];
            end

            trajUp = Trajectory(t,s,unom,K);           
            
            
        end
        
    end
end