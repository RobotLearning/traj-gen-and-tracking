% Class for holding inputs and states on trajectories.
% Class is used instead of a structure to enforce certain fields.
% TODO: do we need this class?

classdef Trajectory < handle
    
    properties
        
        % time profile
        t
        % number of discretizations
        N
        % discrete set of states on (cts) trajectory
        s      
        % nominal u values calculated during trajectory generation
        % (optional)
        unom
        % via points used for spline-based trajectory generation
        sp
        % particular algorithm's performance as array of structures
        PERF
    end
    
    methods (Access = public)
         
        function obj = Trajectory(t,sp,s,unom)            
            obj.t = t;
            obj.N = length(t);
            obj.sp = sp;
            obj.s = s;
            obj.unom = unom;
        end
        
        % fun_cost is a cost function, not necessarily quadratic
        % controller is generally a particular ILC implementation
        function addPerformance(obj,u,x,costfun,controller)
            
            if ischar(controller)
                name = controller;
            else
                name = controller.name;
            end
            
            % append performance
            i = length(obj.PERF);
            obj.PERF(i+1).name = name;
            obj.PERF(i+1).u = u;
            obj.PERF(i+1).x = x;
            obj.PERF(i+1).err = x - obj.s;
            obj.PERF(i+1).cost = costfun.fnc(x,obj.s);

            % display SSE error
            sse = sum(obj.PERF(i+1).cost);
            fprintf('Q-SSE for %s is %f \n', name, sse);
            % add to controller
            if ~ischar(controller)
                controller.record(u,sse);
            end
            
        end
        
    end
end