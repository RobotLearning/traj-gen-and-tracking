% Calculate racket desired position, velocity and orientation 
% at a certain incoming (mean) ball position
% Using the mirror law and inverse kinematics 

function [racketPos, racketVel,racketOrient] = calcDesRacketState(ballPos,ballOutVel,ballInVel,par)

CRR = par.CRR;

% inverting the mirror law
normal = (ballOutVel - ballInVel) ...
              ./ norm(ballOutVel - ballInVel);
velOutAlongNormal = ballOutVel' * normal;
velInAlongNormal = ballInVel' * normal;
racketVelAlongNormal = (velOutAlongNormal + ... 
    CRR * velInAlongNormal) / (1 + CRR);
% racket velocity along racket plane we take to be zero
% this implies no spinning
racketVel = racketVelAlongNormal * normal;

% get desired orientation angle 
racketOrient = normal;

% racket center is equal to ballPos
racketPos = ballPos;