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
    end
    
    methods (Access = public)
         
        %% create a trajectory class useful for simulations
        function obj = Trajectory(t,s,unom,K)            
            obj.t = t;
            obj.N = length(t);
            obj.s = s;
            obj.unom = unom;
            obj.K = K;
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
            
            % append performance
            i = length(obj.PERF);
            obj.PERF(i+1).name = name;
            obj.PERF(i+1).u = u;
            obj.PERF(i+1).y = y;
            obj.PERF(i+1).err = y - obj.s;
            obj.PERF(i+1).cost = costfun.fnc(y,obj.s);

            % display SSE error
            sse = sum(obj.PERF(i+1).cost);
            fprintf('%d. trial: SSE for %s is %f \n', i+1, name, sse);
            
            if ~ischar(controller) && ~isempty(u)
                controller.record(u,sse);
            end
            
        end
        
        % inverse project output trajectory s back to state space
        function sbar = projectBack(obj,C)
            
            sbar = C'*((C*C')\obj.s);
        end
        
    end
end