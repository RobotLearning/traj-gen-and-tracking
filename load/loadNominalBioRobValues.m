% Nominal BioRob values
%
% Apart from the masses and link lenghts the values of the robot here are
% completely random, so it does not, in fact, represent the real dynamics of
% the bioRob. As we approach the real experiment we should find the proper
% parameters of the robot. But for the moment I would leave the numbers like this.
%
% 12.05.2015 Guilherme.
%
% Edited 29.12.2015 Okan
%

deg = pi/180;
linkLength  =   [0 .307 .310 .09;];
gearRatio   = 1*[1000 100 100 100];
linkMass    =   [1.5 1.35 .53 .35]; %original biorob
linkInertia =  1/12*linkMass.*(linkLength).^2; % rod model


% this link has zero length and lots of friction.
% it is not being controlled at the moment.
% This makes the arm to reach a task on a fixed XZ plane.
links(1).m = linkMass(1);
links(1).cm = [0,0,0];
links(1).inertia = [0, 0, 0.5*linkInertia(1), 0, 0, 0];
links(1).motor.inertia = 200e-6;
links(1).motor.gearRatio = gearRatio(1);
links(1).friction.viscous = 1.48e-3;
links(1).friction.coulomb = [0.395, -0.435]; % alpha = pi/2

links(2).m = linkMass(2);
links(2).cm = [0.307/2, 0, 0];
links(2).inertia = [0, 0, linkInertia(2), 0, 0, 0];
links(2).motor.inertia = 200e-6;
links(2).motor.gearRatio = gearRatio(2);
links(2).friction.viscous = 0.817e-3;
links(2).friction.coulomb = [0.126 -0.071];

links(3).m = linkMass(3);
links(3).cm = [0.31/2, 0, 0];
links(3).inertia = [0, 0, linkInertia(3), 0, 0, 0];
links(3).motor.inertia = 200e-6;
links(3).motor.gearRatio = gearRatio(3);
links(3).friction.viscous = 1.38e-3;
links(3).friction.coulomb = [0.132, -0.105];

links(4).m = linkMass(4);
links(4).cm = [0.09/2, 0, 0];
links(4).inertia = [0, 0, linkInertia(4), 0, 0, 0];
links(4).motor.inertia = 200e-6;
links(4).motor.gearRatio = gearRatio(4);
links(4).friction.viscous = 1.2e-3;
links(4).friction.coulomb = [11.2e-3, -16.9e-3];

% some useful poses
%
some_postures.qany   = [0 -100 -10 +45] * deg;
some_postures.qz     = [0 0 0 0] * deg; % zero angles, L shaped pose
some_postures.qtop   = [0 pi/2 0 0]; % ready pose, arm up
some_postures.qdown  = [0 -pi/2 0 0 0];

%% Organize values

% Set default end effector parameters
eff(1).m = 0.0;
eff(1).mcm(1) = 0.0;
eff(1).mcm(2) = 0.0;
eff(1).mcm(3) = 0.0;
eff(1).x(1)  = 0.0;
eff(1).x(2)  = 0.0;
eff(1).x(3)  = 0;
eff(1).a(1)  = 0.0;
eff(1).a(2)  = 0.0;
eff(1).a(3)  = 0.0;

% External forces
for j = 1:3
    % I guess this is the external force to the base
    uex0.f(j) = 0.0;
    uex0.t(j) = 0.0;
    for i = 1:4
        uex(i).f(j) = 0.0;
        uex(i).t(j) = 0.0;
    end
end

% base cartesian position and orientation (quaternion)
basec.x  = [0.0,0.0,0.0];
basec.xd = [0.0,0.0,0.0];
basec.xdd = [0.0,0.0,0.0];
baseo.q = [0.0,1.0,0.0,0.0];
baseo.qd = [0.0,0.0,0.0,0.0];
baseo.qdd = [0.0,0.0,0.0,0.0];
baseo.ad = [0.0,0.0,0.0];
baseo.add = [0.0,0.0,0.0];