% Robot abstract class is the parent for particular
% robot arm manipulators.

classdef (Abstract) Robot < Model
        
    % Field necessary for all robotic manipulators
    properties (Abstract)        
        
        % cartesian coordinates of the trajectory
        x, xd, xdd
        
        % joint space coordinates
        q, qd, qdd
    end
    
    % methods to be implemented
    methods (Abstract)
        
        % kinematics
        kinematics(q)
        % inverse kinematics
        inverseKinematics(x)
        % dynamics to get u
        dynamics(q,qd,qdd)
        % inverse dynamics to qet Qd = [qd,qdd]
        inverseDynamics(Q,u)
        % make an animation of the robot manipulator
        animateArm(x)
        
    end
    
end