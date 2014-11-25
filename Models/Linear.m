% A simple linear dynamics model
% continous model always discretized

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
    end
    
    % Fields pertaining to state evolution and measurements
    properties
        % cartesian coordinates of the state
        x
        % observation of the state (output)
        y
        % A and B matrices
        A, B
        % discrete A and B matrices
        Ad, Bd
        % observation matrix
        C
    end
    
    methods
        
        % copies the parameter values inside the structure
        function set.PAR(obj, STR)  
            
            % create A and B matrices
            if ~obj.SIM.discrete
                obj.A = STR.A;
                obj.B = STR.B;                
            else
                obj.Ad = STR.Ad;
                obj.Bd = STR.Bd;
            end
            obj.C = STR.C;
                     
        end
        
        % copies the constraint values inside the structure
        function set.CON(obj, STR)
            
            % TODO:
            obj.CON = STR;
            
        end 
        
        % set the simulation parameters
        function set.SIM(obj, sim)
            
            obj.SIM.discrete = sim.discrete;
            obj.SIM.dimx = sim.dimx;
            obj.SIM.dimu = sim.dimu;
            obj.SIM.dimy = sim.dimy;
            obj.SIM.h = sim.h;
            obj.SIM.eps = sim.eps;
            assert(strcmpi(sim.int,'Euler') || strcmpi(sim.int,'RK4'),...
                   'Please input Euler or RK4 as integration method');
            obj.SIM.int = sim.int;
        end
        
        % change the cost function
        % TODO: assert Q and R matrices are of correct size
        function set.COST(obj, STR)
            obj.COST.Q = STR.Q;
            obj.COST.R = STR.R;
            obj.COST.fnc = @(x1,x2) diag((x1-x2)'*STR.Q*(x1-x2));
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
            if ~obj.SIM.discrete
                obj.SIM.discrete = true;
                obj.discretizeMatrices();
            end
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
        function x_next = nominal(obj,t,x,u)

            x_next = obj.Ad*x + obj.Bd*u;        
            
        end
        
        % provides actual model
        function x_next = actual(obj,t,x,u)
            
            % TODO          
            %error('Not Implemented');
            x_next = obj.Ad*x + obj.Bd*u;        
            
        end
        
        % TODO: isn't output controllability enough?
        function assertControllability(obj)
            
            % make sure the system is controllable/reachable
            % otherwise give an error
            
            % construct controllability Kalman matrix
            Ad = obj.Ad;
            Bd = obj.Bd;
            dimx = obj.SIM.dimx;
            dimu = obj.SIM.dimu;
            K = zeros(dimx,dimx*dimu);
            for i = 0:dimx-1
                K(:,(dimx*i+1):(dimx*i+dimu)) = (Ad^i)*Bd;
            end
            assert(rank(K) == dimx, 'System is not controllable!');
        end
        
        % Wrapper for the LQR function
        % constructs the matrices of the model required for LQR
        function K = lqr(obj,t,s)
            
            dimx = obj.SIM.dimx;
            dimu = obj.SIM.dimu;
            N = length(t) - 1;            
                        
            % form the time varying matrices Abar and Bbar
            Abar = zeros(dimx+1,dimx+1,N);
            Bbar = zeros(dimx+1,dimu,N);
            for i = 1:N
                Abar(:,:,i) = [obj.Ad, obj.Ad*s(:,i) - s(:,i+1); ...
                               zeros(1,dimx), 1];
                Bbar(:,:,i) = [obj.Bd; 0];
            end

            MODE.N = N;
            MODE.LTI = false;
            Q = [obj.COST.Q, zeros(dimx,1); zeros(1,dimx), 0];
            R = obj.COST.R;

            K = LQR(Q,R,Abar,Bbar,MODE);
            
        end
        
        % Creates a dmp trajectory
        % TODO: extend for multiple outputs
        function s = dmpTrajectory(obj,t,numbf,goal,yin,ydes)
            
            h = obj.SIM.h;
            N = length(t);
            % make a DMP that smoothens x_des
            pat = 'd';
            % scaling
            tau = 1;
            % ensure that phase decays
            exponent = 2;
            ax = exponent/(tau*t(end));            
            can = Canonical(h,ax,tau,N,pat);
            
            % create different DMPs, one for each dimension
            alpha = 25;
            beta = 25/4;

            % final phase point
            xf = exp(-tau*ax*t(end));
            
            % initialize the forcing functions
            force.h = ones(numbf,1) * numbf^(1.5);
            force.c = linspace(xf,1,numbf);

            % create the dmp trajectory
            dmp = discreteDMP(can,alpha,beta,goal,yin,force);
            % learn the weights with locally weighted regression
            dmp = LWR(ydes,dmp,force);
            % evolve the DMP
            [x,s] = dmp.evolve();         
            % add accelerations
            s = [s; 0, diff(s(2,:))/h];
            
            figure;
            plot(t,ydes,'-',t,s(1,:),'-.',t,x);
            legend('desired trajectory','state y','phase');
            title('Followed trajectory for DMP');

        end

        % create a dmp trajectory and inputs using feedback law
        function Traj = trajectory(obj,t,y0,ydes)
            
            % check controllability
            obj.assertControllability();
             
            % optional: make DMPs that smoothens x_des
            % one for each output
            goal = ydes(:,end);
            numbf = 50;
            s = obj.dmpTrajectory(t,numbf,goal,y0,ydes);

            % calculate the optimal feedback law
            K = obj.lqr(t,s);
            
            Traj = Trajectory(t,s,[],K);
        end
        
    end
end