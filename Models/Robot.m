% Robot abstract class is the parent for particular
% robot arm manipulators.

classdef (Abstract) Robot < Model
        
    % Field necessary for all robotic manipulators?
    properties (Abstract)        
        
        % cartesian coordinates of the trajectory
        %x, xd, xdd
        
        % joint space coordinates
        %q, qd, qdd
    end
    
    % methods to be implemented
    methods (Abstract)
        
        % kinematics
        kinematics(q)
        % inverse kinematics
        invKinematics(x)
        % inverse dynamics to get u
        invDynamics(q,qd,qdd)
        % (direct) dynamics to qet Qd = [qd,qdd]
        dynamics(Q,u)
        % make an animation of the robot manipulator
        animateArm(x)
        
    end
    
end