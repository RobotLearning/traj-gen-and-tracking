% A simple linear dynamics model

classdef Linear < Model

    properties   
        % parameters structure
        PAR
        % constraints structure
        CON
        % cost function structure (handle and weight matrix)
        COST
        % fields necessary for simulation and plotting, noise etc.
        SIM
        % cartesian coordinates of the trajectory
        x, xd, xdd
        % A and B matrices
        A, B
        % discrete A and B matrices
        Ad, Bd
    end
    
    methods
        
        % copies the parameter values inside the structure
        function set.PAR(obj, STR)  
            
            % create A and B matrices
            obj.A = STR.A;
            obj.B = STR.B;
                     
        end
        
        % copies the constraint values inside the structure
        function set.CON(obj, STR)
            
            % TODO:
            obj.CON = STR;
            
        end 
        
        % set the simulation parameters
        function set.SIM(obj, sim)
            obj.SIM.dimx = sim.dimx;
            obj.SIM.dimu = sim.dimu;
            obj.SIM.h = sim.h;
            assert(strcmpi(sim.int,'Euler') || strcmpi(sim.int,'RK4'),...
                   'Please input Euler or RK4 as integration method');
            obj.SIM.int = sim.int;
        end
        
        % change the cost function
        function set.COST(obj, STR)
            obj.COST.Q = STR.Q;
            obj.COST.R = STR.R;
            obj.COST.fnc = @(x1,x2) diag((x1-x2)'*STR.Q*(x1-x2));
            assert(length(STR.Q) == obj.SIM.dimx);
        end
        
    end
    
    methods
        
        % constructor for convenience
        % TODO: divide into several methods?
        function obj = Linear(par,con,cost,sim)
            
            obj.SIM = sim;
            % set object parameter
            obj.PAR = par;
            % set object constraints
            obj.CON = con;    
            % cost function handle
            obj.COST = cost;
            
            % create discrete matrices
            obj.discretizeMatrices();
        end
        
        function obj = discretizeMatrices(obj)
            
            dimu = obj.SIM.dimu;
            dimx = obj.SIM.dimx;
            h = obj.SIM.h;
            % trick to get discrete time versions
            Mat = [obj.A, obj.B; zeros(dimu, dimx + dimu)];
            MD = expm(h * Mat);
            obj.Ad = MD(1:dimx,1:dimx);
            obj.Bd = MD(1:dimx,dimx+1:end);
        end
        
        % provides nominal model
        function x_dot = nominal(obj,t,x,u)

            x_dot = obj.A*x + obj.B*u;        
            
        end
        
        % provides actual model
        function x_dot = actual(obj,t,x,u)
            
            % TODO          
            %error('Not Implemented');
            x_dot = obj.A*x + obj.B*u;
            
        end
        
        function assertControllability(obj)
            
            % make sure the system is controllable/reachable
            % otherwise give an error
            
            % construct controllability Kalman matrix
            A = obj.A;
            B = obj.B;
            dim = obj.SIM.dimx;
            K = zeros(dimx,dimx*dimu);
            for i = 0:dim-1
                K(:,dimx*i:dimx*i+dimu) = A^(i)*B;
            end
            assert(rank(K) == dim, 'System is not controllable!');
        end
        
        %TODO: call nominal or actual evolution functions
        %debug testLQR first
        function [x,u,K] = lqr(obj,t,x0,s)
            
            dimx = obj.SIM.dimx;
            dimu = obj.SIM.dimu;
            N = length(t) - 1;            
                        
            % form the time varying matrices Abar and Bbar
            Abar = zeros(dimx+1,dimx+1,N);
            Bbar = zeros(dimx+1,dimu,N);
            for i = 1:N
                Abar(:,:,i) = [obj.Ad, obj.Ad*s(:,i) - s(:,i+1); ...
                               zeros(1,dimx), 0];
                Bbar(:,:,i) = [obj.Bd; 0];
            end

            MODE.N = N;
            MODE.LTI = false;
            Q = [obj.COST.Q, zeros(dimx,1); zeros(1,dimx), 0];
            R = obj.COST.R;

            K = LQR(Q,R,Abar,Bbar,MODE);

            % simulate system
            u = zeros(dimu,N);
            e0 = x0 - s(:,1);
            ebar(:,1) = [e0; 1];

            for i = 1:N
                if length(size(K)) == 3 % finite horizon
                    u(:,i) = K(:,:,i)*ebar(:,i);
                else % infinite horizon
                    u(:,i) = K*ebar(:,i);
                end

                %xtest(:,i+1) = Ad * xtest(:,i) + Bd * u(:,i);
                ebar(:,i+1) = Abar(:,:,i)*ebar(:,i) + Bbar(:,:,i)*u(:,i);
            end

            % extract the signals from error
            x = ebar(1:3,:) + s;
            
        end
        
        function s = dmpTrajectory(obj,t,numbf,goal,x_des)
            
            h = obj.SIM.h;
            dim = obj.SIM.dimx;
            N = length(t);
            % make a DMP that smoothens x_des
            pat = 'd';
            ax = 1;
            tau = 1;
            can = Canonical(h,ax,tau,N,pat);
            
            % create different DMPs, one for each dimension
            alpha = 25;
            beta = 25/4;

            % initialize the forcing functions
            force.h = ones(numbf,1) * numbf^(1.5);
            force.c = linspace(t(1),t(end),numbf);

            % initialize the dmp trajectory
            s = zeros(dimx,N);
            for i = 1:dim
                % initial states of DMPs
                yin = obj.PAR.state.init(i);
                dmp(i) = discreteDMP(can,alpha,beta,goal,yin,force);
                % learn the weights with locally weighted regression
                dmp(i) = LWR(x_des(i,:),dmp,force);
                % evolve the DMP
                [x,s(i,:)] = dmp(i).evolve();
            end

        end

        % TODO: extend using planning to incorporate constraints
        function [Traj,K] = trajectory(obj,t,x_des)

            N = length(t);
            Nu = N-1;
            dimx = obj.SIM.dimx;
            dimu = obj.SIM.dimu;
            
            % check controllability
            obj.assertControllability();
            
            % optional: make a DMP that smoothens x_des
            s = obj.dmpTrajectory(t,x_des);            
            
            % calculate the optimal feedback law
            % velocity is differences in the discrete case
            v = diff(s')';

            % form the time varying matrices Abar and Bbar
            Abar = zeros(dimx+1,dimx+1,Nu);
            Bbar = zeros(dimx+1,dimu,Nu);
            for i = 1:Nu
                Abar(:,:,i) = [obj.A, (obj.A-eye(dimx))*s(:,i) - v(:,i); ...
                               zeros(1,dimx), 0];
                Bbar(:,:,i) = [obj.B; zeros(1,dimu)];
            end

            MODE.N = N;
            MODE.LTI = false;
            % fifth dimension is for pseudovariable 1
            Q = [obj.COST.Q, zeros(dimx,1); zeros(1,dimx), 0];
            R = obj.COST.R;

            K = LQR(Q,R,Abar,Bbar,MODE);
            
            Traj = Trajectory(t,[],s,[]);
        end
        
    end
end