% Simple BioRob model with 4 degrees of freedom

classdef BioRob < Robot

    properties   
        % parameters structure
        PAR
        % constraints structure
        CON
        % cost function structure (handle and weight matrix)
        COST
        % fields necessary for simulation and plotting, noise etc.
        SIM
        % observation matrix
        C
        % jacobian matrix
        jac
        % learning in joint space?
        flag_jspace
        % reference shown in joint space?
        flag_ref_jsp
    end
    
    methods
        
        % copies the parameter values inside the structure
        function set.PAR(obj, STR)  
            
            %assert(all(strcmp(fieldnames(obj.PAR), fieldnames(STR))));
            obj.PAR = STR;
            
            % set observation matrix
            obj.C = STR.C;
        end
        
        % copies the constraint values inside the structure
        function set.CON(obj, STR)
            
            % check that the input has all the fields
            %assert(all(strcmp(fieldnames(obj.CON), fieldnames(STR))));
            obj.CON = STR;
            
        end 
        
        % set the simulation parameters
        function set.SIM(obj, sim)
            obj.SIM.discrete = sim.discrete;
            obj.SIM.dimx = sim.dimx;
            obj.SIM.dimy = sim.dimy;
            obj.SIM.dimu = sim.dimu;
            obj.SIM.h = sim.h;
            obj.SIM.eps_m = sim.eps_m;
            %assert(strcmpi(sim.int,'Euler') || strcmpi(sim.int,'RK4'),...
            %       'Please input Euler or RK4 as integration method');
            obj.SIM.int = sim.int;
            obj.flag_jspace = ~sim.cartesian;
            obj.flag_ref_jsp = sim.jref;
        end
        
        % change the cost function
        function set.COST(obj, cost)
            obj.COST.Q = cost.Q;
            obj.COST.R = cost.R;
            obj.COST.fnc = @(x1,x2) diag((x1-x2)'*cost.Q*(x1-x2));
            %assert(length(Q) == obj.SIM.dimx);
        end
        
    end
    
    methods
        
        % constructor for convenience
        function obj = BioRob(par,con,cost,sim)
            
            obj.SIM = sim;            
            % set object parameter
            obj.PAR = par;
            % set object constraints
            obj.CON = con;        
            % cost function handle
            obj.COST = cost;
            % jacobian in case we need it
            obj.jac = [];
        end
        
        %% Dynamics functions here
        % provides nominal model
        function [x_dot,varargout] = nominal(obj,~,x,u,flg)
            % differential equation of the forward dynamics
            % x_dot = A(x)x + B(x)u + C(x)
            [x_dot,dfdx,dfdu] = obj.dynamics(x,u,true);
            
            % return jacobian matrices
            if flg
                varargout{1} = dfdx;
                varargout{2} = dfdu;
            end
        end
        
        % provides actual model
        function x_dot = actual(obj,t,x,u)
            
            % differential equation of the inverse dynamics
            % x_dot = A(x) + B(x)u
            
            % load actual values
            par = obj.PAR;
            par.link0 = link0;
            par.links = links;
            d = @(t) 0.0; %0.2 + 0.5 * sin(10*t);
            x_dot = bioRobDynamics(x,u,par,false) + d(t);
            
            
        end       
        
        % dynamics to get u
        function u = invDynamics(obj,q,qd,qdd)
            % inverse dynamics model taken from P.Corke toolbox
            u = bioRobInvDynamics(q,qd,qdd,obj.PAR);
            u(1,:) = 0.0; % force torques for the 1st link to be zero 
        end
        
        % dynamics to qet Qd = [qd,qdd]
        function [Qd, varargout] = dynamics(obj,Q,u,flag)
            if flag
                [Qd, dfdx, dfdu] = bioRobDynamics(Q,u,obj.PAR,flag);
                varargout{1} = dfdx;
                varargout{2} = dfdu;
            else
                Qd = bioRobDynamics(Q,u,obj.PAR,flag);
            end
        end
        
        %% Kinematics functions here
        % run kinematics using an external function
        function [x1,x2,x3,Ahmats] = kinematics(obj,q)
            
            [x1,x2,x3,Ahmats] = bioRobKinematics(q,obj.PAR);
        end
        
        % special inverse kinematics staying close to q_prev
        function q = invKinematics(obj,x,phi)

            % masking the x-z angles
            M  = [1 1 1 0 1 0];
            y = x(1,:);
            z = x(2,:);
            q0 = [-0.0000  -63.3438  127.2693  -63.9255]; % initial guess
    
            fprintf('Computing IK...\n'); tic;
            for k = 1:1:size(x,2) 
        
                T(:,:,k) = transl(y(k),0,z(k))*troty(phi(k))*trotz(-angle(k).*pi/180 );
                %T(:,:,k) = transl(x(k), 0, z(k)) * trotz(-angle(k).*pi/180) ;
        
                %q = bioRobInvKinematics(x,phi,obj.PAR);
                q(:,k) = bioRobInvKinematics(T(:,:,k), q0, M, 'pinv');
                q0 = q(:,k); % use previous joint values as initial guess for ikine
            end
            fprintf('finished IK in %g seconds.\n', toc);
    
            q = wrapToPi(q);
            q = unwrap(q);
        end            
        
        function trj = genReachingTraj(obj,tf,y,z,phi)
            
            f = 1/tf;
            dt = obj.SIM.h;
            t = 0:dt:tf;

            Ay = y(2) - y(1);
            Az = z(2) - z(1);
            Aa = phi(2) - phi(1);

            Y = y(1) + 0.5*Ay*(1 - cos(2*pi*f*t));
            Z = z(1) + 0.5*Az*(1 - cos(2*pi*f*t));
            X = [Y;Z];
            Phi = phi(1) + 0.5*Aa*(1 - cos(2*pi*f*t)); % orientation of the end-effector
            
            q = obj.invKinematics(X,Phi);
            qd = diff(q')' / dt; 
            qd(:,end+1) = qd(:,end) + qd(:,end) - qd(:,end-1);
            Q = [q;qd];
            obj.flag_ref_jsp = true;
            trj = obj.generateInputs(t,Q);
            
        end
                   

        %% Plots and animation
        % to draw the robots joints and endeffector 
        % for one posture only
        function [joints,endeff,racket] = drawPosture(obj,q)
            
            % TODO
            
        end
        
        % overwrites the model superclass function
        function h = plot_outputs(obj,trj)

            % default view angle [x y z] of the plot (this is just eye-candy)
            dfv = [0.1 -1  0.1];
            nRef = numel(trj);
            cl  = color_finder(nRef);

            trajType = trj{1}.param.type;

            h.anim = figurew('bioRob'); view(dfv);
            drawnow; axis 'equal'; drawnow;
            set_fig_position(  [0.119 0.542 0.264 0.436]  );
            % Define axes limits depending on the task
            % =================================
            ylim([-0.05   0.05]);
            if strcmp(trajType, 'ellipse')
                xlim([0       0.6]);
                zlim([-0.25   0.45]);
                txti = numel(trj{1}.xyz(:,1));
            end
            if strcmp(trajType, 'p2p')
                xlim([  -.10       0.7]);
                zlim([  -0.35   0.45]);
                txti = round(numel(trj{1}.xyz(:,1))/2);
            end    
    
            h.q = figurew('joint_position');  
            set_fig_position( [0.588 0.43 0.17 0.546]  );

            h.ilc = figurew('feedforward_motor_commands');
            set_fig_position([0.417 0.429 0.17 0.546]   );    
    
            for k = 1:nRef

                ref = trj{k};
                style = struct('Color', cl(k,:), 'LineStyle', '-', 'LineWidth', 2);
                names{k} = ref.name;
                figure(h.anim)
                plot3(ref.xyz(:,1), ref.xyz(:,2), ref.xyz(:,3), style); 
                text(ref.xyz(txti,1)+.01, ref.xyz(txti,2), ref.xyz(txti,3),...
                     [ref.name])
                figure(h.q)
                for j=1:3
                    subplot(3,1,j); hold on;
                    plot(ref.t, r2d(ref.q(:,j+1)), style);
                    xlabel('Time (s)');
                    ylabel(['Joint' num2str(j+1) ' (deg)']);        
                end  
                subplot(3,1,1);
                legend(names);

                figure(h.ilc)        
                for j=1:3
                    subplot(3,1,j); hold on;
                    plot(ref.t, ref.u(:,j+1), style);
                    xlabel('Time (s)');
                    ylabel(['u' num2str(j+1) ]);
                end
                subplot(3,1,1);
                legend(names);        
            end

            % Do this just to plot the robot in its initial configuration
            %figure(h.anim)
            %obj.animateArm(trj{1}.q(1:5,:), 1);        
    
        end
        
        % make an animation of the robot manipulator
        % TODO
        function animateArm(obj,q,num)
            
            [m,n] = size(q);

            if m >= nTraj
                t_orig  = linspace(0,1,m);
                t_small = linspace(0,1,num);
                q_small = interp1(t_orig, q, t_small);
            else
                q_small = q;
            end
    
            obj.plot(q_small);
        end
        
        % get lifted model constraints
        function [umin,umax,L,q] = lift_constraints(obj,trj,ilc)
            
            %TODO:
            
        end
        
    end
end