% Two dimensional linear dynamics model
% TODO: can be extended for Linear ND Dynamics

classdef Linear2DDynamics < Model

    properties   
        % parameters structure
        PAR
        % constraints structure
        CON
        % cost function structure (handle and weight matrix)
        COST
        % fields necessary for simulation and plotting, noise etc.
        SIM
        % cartesian coordinates of the nominal trajectory
        x, xd, xdd
        % A and B matrices
        A, B
    end
    
    methods
        
        % copies the parameter values inside the structure
        function set.PAR(obj, STR)  
            
            obj.PAR = STR;
            m1 = obj.PAR.m1;
            k1 = obj.PAR.k1;
            b1 = obj.PAR.b1;
            m2 = obj.PAR.m2;
            k2 = obj.PAR.k2;
            b2 = obj.PAR.b2;
            % create A and B matrices
            obj.A = [0, 1, 0, 0;
                     -k1/m1, -b1/m1, 0, 0;
                     0, 0, 0, 1;
                     0, 0, -k2/m2, -b2/m2];
            obj.B = [0; 1/m1; 0; 1/m2];
                     
        end
        
        % copies the constraint values inside the structure
        function set.CON(obj, STR)
            
            % TODO:
            obj.CON = STR;
            
        end 
        
        % set the simulation parameters
        function set.SIM(obj, sim)
            obj.SIM.dimx = 4;
            obj.SIM.dimu = 2;
            obj.SIM.h = sim.h;
            obj.SIM.eps = sim.eps;
            obj.SIM.eps_d = sim.eps_d;
            assert(strcmpi(sim.int,'Euler') || strcmpi(sim.int,'RK4'),...
                   'Please input Euler or RK4 as integration method');
            obj.SIM.int = sim.int;
        end
        
        % change the cost function
        function set.COST(obj, Q)
            obj.COST.Q = Q;
            obj.COST.fnc = @(x1,x2) diag((x1-x2)'*Q*(x1-x2));
            assert(length(Q) == obj.SIM.dimx);
        end
        
    end
    
    methods
        
        % constructor for convenience
        % TODO: divide into several methods?
        function obj = Linear2DDynamics(par,con,cost,sim)
            
            obj.SIM = sim;
            % set object parameter
            obj.PAR = par;
            % set object constraints
            obj.CON = con;    
            % cost function handle
            obj.COST = cost.Q;
        end
        
        % provides nominal model
        function x_dot = nominal(obj,t,x,u)

            x_dot = obj.A*x + obj.B*u;        
            
        end
        
        % provides actual model
        function x_dot = actual(obj,t,x,u)
            
            % TODO                       
            x_dot = obj.A*x + obj.B*u;
            
        end

        % using a simple inverse kinematics method
        % TODO: extend using planning to incorporate constraints
        function Traj = trajectory(obj,t,x_des)

            N = length(t);
            dimu = obj.SIM.dimu;
            h = obj.SIM.h;
            unom = zeros(dimu,N);
            
            % make sure the system is controllable/reachable
            % otherwise give an error
            
            % construct controllability Kalman matrix
            A = obj.A;
            B = obj.B;
            K = [B, A*B, A^2 * B, A^3 * B];
            assert(rank(K) == obj.SIM.dimx, 'System is not controllable!');
            
            % make a DMP that smoothens x_des
            pat = 'd';
            ax = 1;
            tau = 1;
            can = Canonical(h,ax,tau,N,pat);
            
            % create two paths
            path1 = x_des(1,:);
            goal1 = path1(end);
            path2 = x_des(2,:);
            goal2 = path2(end);

            alpha = 25;
            beta = 25/4;

            % learn the weights with locally weighted regression
            force.h = ones(100,1); % * 100^(1.5);
            force.c = linspace(0,1,100);
            force1 = LWR(path1,can,alpha,beta,force);
            force2 = LWR(path2,can,alpha,beta,force);

            % create two different DMPs
            yin1 = obj.PAR.state.init(1:2);
            dmp1 = discreteDMP(can,alpha,beta,goal1,yin1,force1);
            yin2 = obj.PAR.state.init(3:4);
            dmp2 = discreteDMP(can,alpha,beta,goal2,yin2,force2);

            [x,s1] = dmp1.evolve();
            [~,s2] = dmp2.evolve();           
            
            % calculate the nominal inputs
            
            Traj = Trajectory(t,[],s,unom);
        end
        
        % get lifted model constraints
        function [umin,umax,L,q] = lift_constraints(obj,trj,ilc)
            
            h = obj.SIM.h;
            N = trj.N - 1; 
            s = trj.s;
            u_trj = trj.unom(:,1:N);
            %dimx = obj.SIM.dimx;
            dimu = obj.SIM.dimu;
            
            % input constraints
            umin(1,:) = obj.CON.u1.min - u_trj(1,:);
            umin(2,:) = obj.CON.u2.min - u_trj(2,:);
            umax(1,:) = obj.CON.u1.max - u_trj(1,:);
            umax(2,:) = obj.CON.u2.max - u_trj(2,:);

            % arrange them in a format suitable for optimization
            umin = umin(:);
            umax = umax(:);
            
            % construct D
            D = (diag(ones(1,dimu*(N-1)),dimu) - eye(dimu*N))/h;
            D = D(1:end-dimu,:);
            
            u_dot_max = [obj.CON.u1.dot.max; obj.CON.u2.dot.max];
            u_dot_min = [obj.CON.u1.dot.min; obj.CON.u2.dot.min];
            U_dot_max = repmat(u_dot_max,N-1,1);
            U_dot_min = repmat(u_dot_min,N-1,1);
            u_star = u_trj(:);

            L1 = [D; -D];
            q1 = [U_dot_max - D*u_star;
                 -U_dot_min + D*u_star];
             
            % state constraints
            % form the constraint matrix C
            C = cell(1,N);
            m = [1 0 0; 0 1 0];
            [C{:}] = deal(m);
            C = blkdiag(C{:});
            C = [C; -C];
            
            x_min = obj.CON.state.x.min - s(1,2:end);
            x_max = obj.CON.state.x.max - s(1,2:end);
            y_min = obj.CON.state.y.min - s(2,2:end);
            y_max = obj.CON.state.y.max - s(2,2:end);

            x_con_max = [x_max(:), y_max(:)];
            x_con_max = x_con_max';
            x_con_max = x_con_max(:);
            x_con_min = [x_min(:), y_min(:)];
            x_con_min = x_con_min';
            x_con_min = x_con_min(:);
            x_con = [x_con_max; -x_con_min];
            
            F = ilc.F;
            d = ilc.filter.x;
            L2 = C*F;
            q2 = x_con - C*d; 
            
            % combine input and state constraints
            L = [L1; L2];
            q = [q1; q2];
            
            
        end
        
    end
end