%% Load actual robot parameters
% These are the results of parametric regression

% definitions
ZSFE  =  0.346;              %!< z height of SAA axis above ground
ZHR  =  0.505;              %!< length of upper arm until 4.5cm before elbow link
YEB  =  0.045;              %!< elbow y offset
ZEB  =  0.045;              %!< elbow z offset
YWR  = -0.045;              %!< elbow y offset (back to forewarm)
ZWR  =  0.045;              %!< elbow z offset (back to forearm)
ZWFE  =  0.255;              %!< forearm length (minus 4.5cm)

% link 0 is the base
link0.m = 0.0;
link0.mcm(1) = 0.0;
link0.mcm(2) = 0.0;
link0.mcm(3) = 0.0;
link0.inertia(1,1) = 0.0; 
link0.inertia(1,2) = 0.0; 
link0.inertia(1,3) = 0.0; 
link0.inertia(2,2) = 0.0;   
link0.inertia(2,3) = 0.0;  
link0.inertia(3,3) = 0.0;  
% SFE joint
links(1).m = 0.00000; 
links(1).mcm(1) = -0.00000;
links(1).mcm(2) = 0.00000;
links(1).mcm(3) = -0.00000;  
links(1).inertia(1,1) = 0.00000; 
links(1).inertia(1,2) = -0.00000; 
links(1).inertia(1,3) = -0.00000; 
links(1).inertia(2,2) = 0.00000;   
links(1).inertia(2,3) = 0.00000;  
links(1).inertia(3,3) = 0.00000;  
% SAA joint
links(2).m = 0.00000;  
links(2).mcm(1) = 0.43430; 
links(2).mcm(2) = -0.00495;
links(2).mcm(3) = -0.00000;
links(2).inertia(1,1) = 0.24768;  
links(2).inertia(1,2) = 0.00364;
links(2).inertia(1,3) = 0.17270;  
links(2).inertia(2,2) = 0.53601; 
links(2).inertia(2,3) = -0.02929; 
links(2).inertia(3,3) =  0.12406;  
% HR joint    
links(3).m = 3.53923; 
links(3).mcm(1) = -0.00889;
links(3).mcm(2) = -0.02148;
links(3).mcm(3) = -1.70741;
links(3).inertia(1,1) = 0.82985;
links(3).inertia(1,2) = -0.01520;
links(3).inertia(1,3) = -0.00612;
links(3).inertia(2,2) = 0.86182;
links(3).inertia(2,3) = -0.00575;
links(3).inertia(3,3) = 0.00071;
% EB joint (elbow)
links(4).m = 1.03409;
links(4).mcm(1) = 0.14089;
links(4).mcm(2) = -0.05914;
links(4).mcm(3) = -0.00270;
links(4).inertia(1,1) = 0.01276;
links(4).inertia(1,2) = 0.00340;
links(4).inertia(1,3) = -0.00229;
links(4).inertia(2,2) = 0.02157;
links(4).inertia(2,3) = 0.00032;
links(4).inertia(3,3) = 0.03718;
% WR joint (wrist 1)
links(5).m = 2.28843;
links(5).mcm(1) = 0.00709;
links(5).mcm(2) = 0.00194; 
links(5).mcm(3) = 0.22347; 
links(5).inertia(1,1) = 0.02182;
links(5).inertia(1,2) = -0.00001;
links(5).inertia(1,3) = -0.00069;
links(5).inertia(2,2) = 0.02184;
links(5).inertia(2,3) = -0.00019;
links(5).inertia(3,3) = 0.00002;
% WFE joint (wrist 2)
links(6).m = 0.25655; 
links(6).mcm(1) = 0.03462; 
links(6).mcm(2) = -0.00415;  
links(6).mcm(3) = 0.00121; 
links(6).inertia(1,1) = 0.00167; 
links(6).inertia(1,2) = 0.00079; 
links(6).inertia(1,3) = -0.00009; 
links(6).inertia(2,2) = 0.00483; 
links(6).inertia(2,3) = -0.00084;
links(6).inertia(3,3) = 0.01101;
% WAA joint (wrist 3)
links(7).m = 0.63285; 
links(7).mcm(1) = -0.00157;  
links(7).mcm(2) = 0.00019;  
links(7).mcm(3) = 0.08286;  
links(7).inertia(1,1) = 0.01129; 
links(7).inertia(1,2) = -0.00006;
links(7).inertia(1,3) = -0.00006;  
links(7).inertia(2,2) = 0.01086; 
links(7).inertia(2,3) = 0.00001;  
links(7).inertia(3,3) = 0.00016;

% make sure inertia matrices are symmetric
for i = 1:7
    for j = 1:3
        for k = j:3
            links(i).inertia(k,j) = links(i).inertia(j,k);
        end
    end
end


% End effector parameters
eff(1).m = 0.0;
eff(1).mcm(1) = 0.0;
eff(1).mcm(2) = 0.0;
eff(1).mcm(3) = 0.0;
eff(1).x(1)  = 0.0;
eff(1).x(2)  = 0.0;
eff(1).x(3)  = 0.30; %0.06; 
eff(1).a(1)  = 0.0;
eff(1).a(2)  = 0.0;
eff(1).a(3)  = 0.0;

% External forces
for j = 1:3
    % I guess this is the external force to the base
    uex0.f(j) = 0.0;
    uex0.t(j) = 0.0;
    for i = 1:7
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