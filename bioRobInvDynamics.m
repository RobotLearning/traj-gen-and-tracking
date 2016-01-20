%
% Modified file from the Robotics Toolbox, Peter Corke
% 
% Compute inverse dynamics via recursive Newton-Euler formulation
%
% Recursive Newton-Euler for modified Denavit-Hartenberg notation. 
%
% See also SERIALLINK.RNE.
%
%
% Copyright (C) 1993-2011, by Peter I. Corke
%
% http://www.petercorke.com

function tau = bioRobInvDynamics(robot, Q, Qd, Qdd)

	z0 = [0;0;1];
	grav = robot.gravity;	% default gravity from the object
	fext = zeros(6, 1);
	n = 4; % number of degrees of freedom of the robot	
    np = size(Q,2);
    Q = Q'; Qd = Qd'; Qdd = Qdd';
	tau = zeros(np,n);

	for p = 1:np
		q = Q(p,:).';
		qd = Qd(p,:).';
		qdd = Qdd(p,:).';
	
		Fm = [];
		Nm = [];
		pstarm = [];
		Rm = [];
		w = zeros(3,1);
		wd = zeros(3,1);
		vd = grav(:);

	%
	% init some variables, compute the link rotation matrices
	%
		for j=1:n
			link = robot.links(j);
			Tj = link.A(q(j));
			if link.RP == 'R'
				D = link.d;
			else
				D = q(j);
			end
			alpha = link.alpha;
			pm = [link.a; -D*sin(alpha); D*cos(alpha)];	% (i-1) P i
			if j == 1
				pm = t2r(robot.base) * pm;
				Tj = robot.base * Tj;
			end
			Pm(:,j) = pm;
			Rm{j} = t2r(Tj);
			if debug>1
				Rm{j}
				Pm(:,j).'
			end
		end

	%
	%  the forward recursion
	%
		for j=1:n
			link = robot.links(j);

			R = Rm{j}.';	% transpose!!
			P = Pm(:,j);
			Pc = link.r;

			%
			% trailing underscore means new value
			%
			if link.RP == 'R'
				% revolute axis
				w_ = R*w + z0*qd(j);
				wd_ = R*wd + cross(R*w,z0*qd(j)) + z0*qdd(j);
				%v = cross(w,P) + R*v;
				vd_ = R * (cross(wd,P) + ...
					cross(w, cross(w,P)) + vd);

			else
				% prismatic axis
				w_ = R*w;
				wd_ = R*wd;
				%v = R*(z0*qd(j) + v) + cross(w,P);
				vd_ = R*(cross(wd,P) + ...
					cross(w, cross(w,P)) + vd ...
				      ) + 2*cross(R*w,z0*qd(j)) + z0*qdd(j);
			end
			% update variables
			w = w_;
			wd = wd_;
			vd = vd_;

			vdC = cross(wd,Pc).' + ...
				cross(w,cross(w,Pc)).' + vd;
			F = link.m*vdC;
			N = link.I*wd + cross(w,link.I*w);
			Fm = [Fm F];
			Nm = [Nm N];
			if debug
				fprintf('w: '); fprintf('%.3f ', w)
				fprintf('\nwd: '); fprintf('%.3f ', wd)
				fprintf('\nvd: '); fprintf('%.3f ', vd)
				fprintf('\nvdbar: '); fprintf('%.3f ', vdC)
				fprintf('\n');
			end
		end

	%
	%  the backward recursion
	%

		fext = fext(:);
		f = fext(1:3);		% force/moments on end of arm
		nn = fext(4:6);

		for j=n:-1:1
			
			%
			% order of these statements is important, since both
			% nn and f are functions of previous f.
			%
			link = robot.links(j);
			
			if j == n
				R = eye(3,3);
				P = [0;0;0];
			else
				R = Rm{j+1};
				P = Pm(:,j+1);		% i/P/(i+1)
			end
			Pc = link.r;
			
			f_ = R*f + Fm(:,j);
			nn_ = Nm(:,j) + R*nn + cross(Pc,Fm(:,j)).' + ...
				cross(P,R*f);
			
			f = f_;
			nn = nn_;

			if debug
				fprintf('f: '); fprintf('%.3f ', f)
				fprintf('\nn: '); fprintf('%.3f ', nn)
				fprintf('\n');
			end
			if link.RP == 'R'
				% revolute
				tau(p,j) = nn.'*z0 + ...
					link.G^2 * link.Jm*qdd(j) - ...
					friction(link, qd(j));
			else
				% prismatic
				tau(p,j) = f.'*z0 + ...
					link.G^2 * link.Jm*qdd(j) - ...
					friction(link, qd(j));
			end
		end
	end
