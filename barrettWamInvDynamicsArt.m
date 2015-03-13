% Barrett WAM inverse dynamics for control and linearization (ILC)

% Articulate Inverse Dynamics taken from SL: 
% shared/barrett/math/InvDynArt_declare.h
% shared/barrett/math/InvDynArt_math.h
% shared/barrett/math/InvDynArt_functions.h
%
% these are called from shared/barrett/src/SL_dynamics.c
%
%
% If flag is set to true, return the jacobians of the dynamics f

function u = barrettWamInvDynamicsArt(q,qd,qdd,PAR)

NDOF  =  7;

% definitions
ZSFE  =  0.346;              %!< z height of SAA axis above ground
ZHR  =  0.505;              %!< length of upper arm until 4.5cm before elbow link
YEB  =  0.045;              %!< elbow y offset
ZEB  =  0.045;              %!< elbow z offset
YWR  = -0.045;              %!< elbow y offset (back to forewarm)
ZWR  =  0.045;              %!< elbow z offset (back to forearm)
ZWFE  =  0.255;              %!< forearm length (minus 4.5cm)
FYOFF  =  0.05;              %!< y offset of rotation axis of fingers
FZOFF  =  (0.15-.015);       %!< z offset from WAA to palm (measured from robot with load cell, load cell subtracted)
F1SEGY  =  0.07;             %!< length of 1st finger segment
F1SEGZ  =  0.002;            %!< vertical offset of axis of 2nd finger segment
F2LENGTH  =  0.056;          %!< lenth of 2nd finger segment
FWIDTH  =  0.025;            %!< finger width
FT2BRATIO  =  3.11;          %!< ratio of finger movement to finger base movement
ANGLE2NDSEG  =  0.733;       %!< angle between 1st and 2nd finger segment when fingers are full extended

% the 1st to 2nd finger segment angles
rf2ndSegAngle  =  ANGLE2NDSEG;
mf2ndSegAngle  =  ANGLE2NDSEG;
lf2ndSegAngle  =  ANGLE2NDSEG;

% extract parameters
link0 = PAR.link0;
links = PAR.links;
eff = PAR.eff;
uex = PAR.uex;
uex0 = PAR.uex0;
basec = PAR.basec;
baseo = PAR.baseo;

g = 9.81;

% sine and cosine precomputation 
ss1th  =  sin(q(1));
cs1th  =  cos(q(1));
ss2th  =  sin(q(2));
cs2th  =  cos(q(2));
ss3th  =  sin(q(3));
cs3th  =  cos(q(3));
ss4th  =  sin(q(4));
cs4th  =  cos(q(4));
ss5th  =  sin(q(5));
cs5th  =  cos(q(5));
ss6th  =  sin(q(6));
cs6th  =  cos(q(6));
ss7th  =  sin(q(7));
cs7th  =  cos(q(7));

% rotation matrix sine and cosine precomputation 
rsMinusrf2ndSegAngle  =  - sin(rf2ndSegAngle);
rcMinusrf2ndSegAngle  =  cos(rf2ndSegAngle);

rsMinuslf2ndSegAngle  =  - sin(lf2ndSegAngle);
rcMinuslf2ndSegAngle  =  cos(lf2ndSegAngle);

rsmf2ndSegAngle  =  sin(mf2ndSegAngle);
rcmf2ndSegAngle  =  cos(mf2ndSegAngle);

% endeffector orientations

rseff1a1  =  sin(eff(1).a(1));
rceff1a1  =  cos(eff(1).a(1));

rseff1a2  =  sin(eff(1).a(2));
rceff1a2  =  cos(eff(1).a(2));

rseff1a3  =  sin(eff(1).a(3));
rceff1a3  =  cos(eff(1).a(3));

%% Includes the functions one by one

%wam_InvDynArtfunc1
     


%barrett_InvDynArtfunc1
     
% rotation matrices 
S00(1,1) = -1 + 2*power(baseo.q(1),2) + 2*power(baseo.q(2),2);
S00(1,2) = 2*(baseo.q(2)*baseo.q(3) + baseo.q(1)*baseo.q(4));
S00(1,3) = 2*(-(baseo.q(1)*baseo.q(3)) + baseo.q(2)*baseo.q(4));

S00(2,1) = 2*(baseo.q(2)*baseo.q(3) - baseo.q(1)*baseo.q(4));
S00(2,2) = -1 + 2*power(baseo.q(1),2) + 2*power(baseo.q(3),2);
S00(2,3) = 2*(baseo.q(1)*baseo.q(2) + baseo.q(3)*baseo.q(4));

S00(3,1) = 2*(baseo.q(1)*baseo.q(3) + baseo.q(2)*baseo.q(4));
S00(3,2) = 2*(-(baseo.q(1)*baseo.q(2)) + baseo.q(3)*baseo.q(4));
S00(3,3) = -1 + 2*power(baseo.q(1),2) + 2*power(baseo.q(4),2);


S10(1,1) = cs1th;
S10(1,2) = ss1th;

S10(2,1) = -ss1th;
S10(2,2) = cs1th;


S21(1,2) = ss2th;
S21(1,3) = cs2th;

S21(2,2) = cs2th;
S21(2,3) = -ss2th;


S32(1,2) = ss3th;
S32(1,3) = -cs3th;

S32(2,2) = cs3th;
S32(2,3) = ss3th;


S43(1,2) = ss4th;
S43(1,3) = cs4th;

S43(2,2) = cs4th;
S43(2,3) = -ss4th;


S54(1,2) = ss5th;
S54(1,3) = -cs5th;

S54(2,2) = cs5th;
S54(2,3) = ss5th;


S65(1,2) = ss6th;
S65(1,3) = cs6th;

S65(2,2) = cs6th;
S65(2,3) = -ss6th;


S76(1,2) = ss7th;
S76(1,3) = -cs7th;

S76(2,2) = cs7th;
S76(2,3) = ss7th;


S87(1,1) = rceff1a2*rceff1a3;
S87(1,2) = rceff1a3*rseff1a1*rseff1a2 + rceff1a1*rseff1a3;
S87(1,3) = -(rceff1a1*rceff1a3*rseff1a2) + rseff1a1*rseff1a3;

S87(2,1) = -(rceff1a2*rseff1a3);
S87(2,2) = rceff1a1*rceff1a3 - rseff1a1*rseff1a2*rseff1a3;
S87(2,3) = rceff1a3*rseff1a1 + rceff1a1*rseff1a2*rseff1a3;

S87(3,1) = rseff1a2;
S87(3,2) = -(rceff1a2*rseff1a1);
S87(3,3) = rceff1a1*rceff1a2;








%barrett_InvDynArtfunc2
     
% inverse rotation matrices 
Si00(1,1) = -1 + 2*power(baseo.q(1),2) + 2*power(baseo.q(2),2);
Si00(1,2) = 2*(baseo.q(2)*baseo.q(3) - baseo.q(1)*baseo.q(4));
Si00(1,3) = 2*(baseo.q(1)*baseo.q(3) + baseo.q(2)*baseo.q(4));

Si00(2,1) = 2*(baseo.q(2)*baseo.q(3) + baseo.q(1)*baseo.q(4));
Si00(2,2) = -1 + 2*power(baseo.q(1),2) + 2*power(baseo.q(3),2);
Si00(2,3) = 2*(-(baseo.q(1)*baseo.q(2)) + baseo.q(3)*baseo.q(4));

Si00(3,1) = 2*(-(baseo.q(1)*baseo.q(3)) + baseo.q(2)*baseo.q(4));
Si00(3,2) = 2*(baseo.q(1)*baseo.q(2) + baseo.q(3)*baseo.q(4));
Si00(3,3) = -1 + 2*power(baseo.q(1),2) + 2*power(baseo.q(4),2);


Si01(1,1) = cs1th;
Si01(1,2) = -ss1th;

Si01(2,1) = ss1th;
Si01(2,2) = cs1th;


Si12(2,1) = ss2th;
Si12(2,2) = cs2th;

Si12(3,1) = cs2th;
Si12(3,2) = -ss2th;


Si23(2,1) = ss3th;
Si23(2,2) = cs3th;

Si23(3,1) = -cs3th;
Si23(3,2) = ss3th;


Si34(2,1) = ss4th;
Si34(2,2) = cs4th;

Si34(3,1) = cs4th;
Si34(3,2) = -ss4th;


Si45(2,1) = ss5th;
Si45(2,2) = cs5th;

Si45(3,1) = -cs5th;
Si45(3,2) = ss5th;


Si56(2,1) = ss6th;
Si56(2,2) = cs6th;

Si56(3,1) = cs6th;
Si56(3,2) = -ss6th;


Si67(2,1) = ss7th;
Si67(2,2) = cs7th;

Si67(3,1) = -cs7th;
Si67(3,2) = ss7th;


Si78(1,1) = rceff1a2*rceff1a3;
Si78(1,2) = -(rceff1a2*rseff1a3);
Si78(1,3) = rseff1a2;

Si78(2,1) = rceff1a3*rseff1a1*rseff1a2 + rceff1a1*rseff1a3;
Si78(2,2) = rceff1a1*rceff1a3 - rseff1a1*rseff1a2*rseff1a3;
Si78(2,3) = -(rceff1a2*rseff1a1);

Si78(3,1) = -(rceff1a1*rceff1a3*rseff1a2) + rseff1a1*rseff1a3;
Si78(3,2) = rceff1a3*rseff1a1 + rceff1a1*rseff1a2*rseff1a3;
Si78(3,3) = rceff1a1*rceff1a2;








%barrett_InvDynArtfunc3
     
% rotation matrices from global to link coordinates 
SG10(1,1) = S00(1,1)*S10(1,1) + S00(2,1)*S10(1,2);
SG10(1,2) = S00(1,2)*S10(1,1) + S00(2,2)*S10(1,2);
SG10(1,3) = S00(1,3)*S10(1,1) + S00(2,3)*S10(1,2);

SG10(2,1) = S00(1,1)*S10(2,1) + S00(2,1)*S10(2,2);
SG10(2,2) = S00(1,2)*S10(2,1) + S00(2,2)*S10(2,2);
SG10(2,3) = S00(1,3)*S10(2,1) + S00(2,3)*S10(2,2);

SG10(3,1) = S00(3,1);
SG10(3,2) = S00(3,2);
SG10(3,3) = S00(3,3);


SG20(1,1) = S21(1,2)*SG10(2,1) + S21(1,3)*SG10(3,1);
SG20(1,2) = S21(1,2)*SG10(2,2) + S21(1,3)*SG10(3,2);
SG20(1,3) = S21(1,2)*SG10(2,3) + S21(1,3)*SG10(3,3);

SG20(2,1) = S21(2,2)*SG10(2,1) + S21(2,3)*SG10(3,1);
SG20(2,2) = S21(2,2)*SG10(2,2) + S21(2,3)*SG10(3,2);
SG20(2,3) = S21(2,2)*SG10(2,3) + S21(2,3)*SG10(3,3);

SG20(3,1) = -SG10(1,1);
SG20(3,2) = -SG10(1,2);
SG20(3,3) = -SG10(1,3);


SG30(1,1) = S32(1,2)*SG20(2,1) + S32(1,3)*SG20(3,1);
SG30(1,2) = S32(1,2)*SG20(2,2) + S32(1,3)*SG20(3,2);
SG30(1,3) = S32(1,2)*SG20(2,3) + S32(1,3)*SG20(3,3);

SG30(2,1) = S32(2,2)*SG20(2,1) + S32(2,3)*SG20(3,1);
SG30(2,2) = S32(2,2)*SG20(2,2) + S32(2,3)*SG20(3,2);
SG30(2,3) = S32(2,2)*SG20(2,3) + S32(2,3)*SG20(3,3);

SG30(3,1) = SG20(1,1);
SG30(3,2) = SG20(1,2);
SG30(3,3) = SG20(1,3);


SG40(1,1) = S43(1,2)*SG30(2,1) + S43(1,3)*SG30(3,1);
SG40(1,2) = S43(1,2)*SG30(2,2) + S43(1,3)*SG30(3,2);
SG40(1,3) = S43(1,2)*SG30(2,3) + S43(1,3)*SG30(3,3);

SG40(2,1) = S43(2,2)*SG30(2,1) + S43(2,3)*SG30(3,1);
SG40(2,2) = S43(2,2)*SG30(2,2) + S43(2,3)*SG30(3,2);
SG40(2,3) = S43(2,2)*SG30(2,3) + S43(2,3)*SG30(3,3);

SG40(3,1) = -SG30(1,1);
SG40(3,2) = -SG30(1,2);
SG40(3,3) = -SG30(1,3);


SG50(1,1) = S54(1,2)*SG40(2,1) + S54(1,3)*SG40(3,1);
SG50(1,2) = S54(1,2)*SG40(2,2) + S54(1,3)*SG40(3,2);
SG50(1,3) = S54(1,2)*SG40(2,3) + S54(1,3)*SG40(3,3);

SG50(2,1) = S54(2,2)*SG40(2,1) + S54(2,3)*SG40(3,1);
SG50(2,2) = S54(2,2)*SG40(2,2) + S54(2,3)*SG40(3,2);
SG50(2,3) = S54(2,2)*SG40(2,3) + S54(2,3)*SG40(3,3);

SG50(3,1) = SG40(1,1);
SG50(3,2) = SG40(1,2);
SG50(3,3) = SG40(1,3);


SG60(1,1) = S65(1,2)*SG50(2,1) + S65(1,3)*SG50(3,1);
SG60(1,2) = S65(1,2)*SG50(2,2) + S65(1,3)*SG50(3,2);
SG60(1,3) = S65(1,2)*SG50(2,3) + S65(1,3)*SG50(3,3);

SG60(2,1) = S65(2,2)*SG50(2,1) + S65(2,3)*SG50(3,1);
SG60(2,2) = S65(2,2)*SG50(2,2) + S65(2,3)*SG50(3,2);
SG60(2,3) = S65(2,2)*SG50(2,3) + S65(2,3)*SG50(3,3);

SG60(3,1) = -SG50(1,1);
SG60(3,2) = -SG50(1,2);
SG60(3,3) = -SG50(1,3);


SG70(1,1) = S76(1,2)*SG60(2,1) + S76(1,3)*SG60(3,1);
SG70(1,2) = S76(1,2)*SG60(2,2) + S76(1,3)*SG60(3,2);
SG70(1,3) = S76(1,2)*SG60(2,3) + S76(1,3)*SG60(3,3);

SG70(2,1) = S76(2,2)*SG60(2,1) + S76(2,3)*SG60(3,1);
SG70(2,2) = S76(2,2)*SG60(2,2) + S76(2,3)*SG60(3,2);
SG70(2,3) = S76(2,2)*SG60(2,3) + S76(2,3)*SG60(3,3);

SG70(3,1) = SG60(1,1);
SG70(3,2) = SG60(1,2);
SG70(3,3) = SG60(1,3);


SG80(1,1) = S87(1,1)*SG70(1,1) + S87(1,2)*SG70(2,1) + S87(1,3)*SG70(3,1);
SG80(1,2) = S87(1,1)*SG70(1,2) + S87(1,2)*SG70(2,2) + S87(1,3)*SG70(3,2);
SG80(1,3) = S87(1,1)*SG70(1,3) + S87(1,2)*SG70(2,3) + S87(1,3)*SG70(3,3);

SG80(2,1) = S87(2,1)*SG70(1,1) + S87(2,2)*SG70(2,1) + S87(2,3)*SG70(3,1);
SG80(2,2) = S87(2,1)*SG70(1,2) + S87(2,2)*SG70(2,2) + S87(2,3)*SG70(3,2);
SG80(2,3) = S87(2,1)*SG70(1,3) + S87(2,2)*SG70(2,3) + S87(2,3)*SG70(3,3);

SG80(3,1) = S87(3,1)*SG70(1,1) + S87(3,2)*SG70(2,1) + S87(3,3)*SG70(3,1);
SG80(3,2) = S87(3,1)*SG70(1,2) + S87(3,2)*SG70(2,2) + S87(3,3)*SG70(3,2);
SG80(3,3) = S87(3,1)*SG70(1,3) + S87(3,2)*SG70(2,3) + S87(3,3)*SG70(3,3);








%barrett_InvDynArtfunc4
     
% velocity vectors 
v0(1) = baseo.ad(1)*S00(1,1) + baseo.ad(2)*S00(1,2) + baseo.ad(3)*S00(1,3);
v0(2) = baseo.ad(1)*S00(2,1) + baseo.ad(2)*S00(2,2) + baseo.ad(3)*S00(2,3);
v0(3) = baseo.ad(1)*S00(3,1) + baseo.ad(2)*S00(3,2) + baseo.ad(3)*S00(3,3);
v0(4) = basec.xd(1)*S00(1,1) + basec.xd(2)*S00(1,2) + basec.xd(3)*S00(1,3);
v0(5) = basec.xd(1)*S00(2,1) + basec.xd(2)*S00(2,2) + basec.xd(3)*S00(2,3);
v0(6) = basec.xd(1)*S00(3,1) + basec.xd(2)*S00(3,2) + basec.xd(3)*S00(3,3);

v1(1) = v0(1)*S10(1,1) + v0(2)*S10(1,2);
v1(2) = v0(1)*S10(2,1) + v0(2)*S10(2,2);
v1(3) = qd(1) + v0(3);
v1(4) = ZSFE*v0(2)*S10(1,1) + v0(4)*S10(1,1) - ZSFE*v0(1)*S10(1,2) + v0(5)*S10(1,2);
v1(5) = ZSFE*v0(2)*S10(2,1) + v0(4)*S10(2,1) - ZSFE*v0(1)*S10(2,2) + v0(5)*S10(2,2);
v1(6) = v0(6);

v2(1) = v1(2)*S21(1,2) + v1(3)*S21(1,3);
v2(2) = v1(2)*S21(2,2) + v1(3)*S21(2,3);
v2(3) = qd(2) - v1(1);
v2(4) = v1(5)*S21(1,2) + v1(6)*S21(1,3);
v2(5) = v1(5)*S21(2,2) + v1(6)*S21(2,3);
v2(6) = -v1(4);

v3(1) = v2(2)*S32(1,2) + v2(3)*S32(1,3);
v3(2) = v2(2)*S32(2,2) + v2(3)*S32(2,3);
v3(3) = qd(3) + v2(1);
v3(4) = ZHR*v2(3)*S32(1,2) + v2(5)*S32(1,2) - ZHR*v2(2)*S32(1,3) + v2(6)*S32(1,3);
v3(5) = ZHR*v2(3)*S32(2,2) + v2(5)*S32(2,2) - ZHR*v2(2)*S32(2,3) + v2(6)*S32(2,3);
v3(6) = v2(4);

v4(1) = v3(2)*S43(1,2) + v3(3)*S43(1,3);
v4(2) = v3(2)*S43(2,2) + v3(3)*S43(2,3);
v4(3) = qd(4) - v3(1);
v4(4) = v3(5)*S43(1,2) + v3(6)*S43(1,3) + v3(1)*(-(ZEB*S43(1,2)) + YEB*S43(1,3));
v4(5) = v3(5)*S43(2,2) + v3(6)*S43(2,3) + v3(1)*(-(ZEB*S43(2,2)) + YEB*S43(2,3));
v4(6) = -(ZEB*v3(2)) + YEB*v3(3) - v3(4);

v5(1) = v4(2)*S54(1,2) + v4(3)*S54(1,3);
v5(2) = v4(2)*S54(2,2) + v4(3)*S54(2,3);
v5(3) = qd(5) + v4(1);
v5(4) = ZWR*v4(3)*S54(1,2) + v4(5)*S54(1,2) + YWR*v4(1)*S54(1,3) - ZWR*v4(2)*S54(1,3) + v4(6)*S54(1,3);
v5(5) = ZWR*v4(3)*S54(2,2) + v4(5)*S54(2,2) + YWR*v4(1)*S54(2,3) - ZWR*v4(2)*S54(2,3) + v4(6)*S54(2,3);
v5(6) = -(YWR*v4(3)) + v4(4);

v6(1) = v5(2)*S65(1,2) + v5(3)*S65(1,3);
v6(2) = v5(2)*S65(2,2) + v5(3)*S65(2,3);
v6(3) = qd(5) - v5(1);
v6(4) = -(ZWFE*v5(1)*S65(1,2)) + v5(5)*S65(1,2) + v5(6)*S65(1,3);
v6(5) = -(ZWFE*v5(1)*S65(2,2)) + v5(5)*S65(2,2) + v5(6)*S65(2,3);
v6(6) = -(ZWFE*v5(2)) - v5(4);

v7(1) = v6(2)*S76(1,2) + v6(3)*S76(1,3);
v7(2) = v6(2)*S76(2,2) + v6(3)*S76(2,3);
v7(3) = qd(7) + v6(1);
v7(4) = v6(5)*S76(1,2) + v6(6)*S76(1,3);
v7(5) = v6(5)*S76(2,2) + v6(6)*S76(2,3);
v7(6) = v6(4);

v8(1) = v7(1)*S87(1,1) + v7(2)*S87(1,2) + v7(3)*S87(1,3);
v8(2) = v7(1)*S87(2,1) + v7(2)*S87(2,2) + v7(3)*S87(2,3);
v8(3) = v7(1)*S87(3,1) + v7(2)*S87(3,2) + v7(3)*S87(3,3);
v8(4) = v7(4)*S87(1,1) + v7(5)*S87(1,2) + v7(3)*(-(eff(1).x(2)*S87(1,1)) + eff(1).x(1)*S87(1,2)) + v7(6)*S87(1,3) + v7(2)*(eff(1).x(3)*S87(1,1) - eff(1).x(1)*S87(1,3)) + v7(1)*(-(eff(1).x(3)*S87(1,2)) + eff(1).x(2)*S87(1,3));
v8(5) = v7(4)*S87(2,1) + v7(5)*S87(2,2) + v7(3)*(-(eff(1).x(2)*S87(2,1)) + eff(1).x(1)*S87(2,2)) + v7(6)*S87(2,3) + v7(2)*(eff(1).x(3)*S87(2,1) - eff(1).x(1)*S87(2,3)) + v7(1)*(-(eff(1).x(3)*S87(2,2)) + eff(1).x(2)*S87(2,3));
v8(6) = v7(4)*S87(3,1) + v7(5)*S87(3,2) + v7(3)*(-(eff(1).x(2)*S87(3,1)) + eff(1).x(1)*S87(3,2)) + v7(6)*S87(3,3) + v7(2)*(eff(1).x(3)*S87(3,1) - eff(1).x(1)*S87(3,3)) + v7(1)*(-(eff(1).x(3)*S87(3,2)) + eff(1).x(2)*S87(3,3));







%barrett_InvDynArtfunc5
     
% c-misc vectors 
c1(1) = qd(1)*v1(2);
c1(2) = -(qd(1)*v1(1));
c1(4) = qd(1)*v1(5);
c1(5) = -(qd(1)*v1(4));

c2(1) = qd(2)*v2(2);
c2(2) = -(qd(2)*v2(1));
c2(4) = qd(2)*v2(5);
c2(5) = -(qd(2)*v2(4));

c3(1) = qd(3)*v3(2);
c3(2) = -(qd(3)*v3(1));
c3(4) = qd(3)*v3(5);
c3(5) = -(qd(3)*v3(4));

c4(1) = qd(4)*v4(2);
c4(2) = -(qd(4)*v4(1));
c4(4) = qd(4)*v4(5);
c4(5) = -(qd(4)*v4(4));

c5(1) = qd(5)*v5(2);
c5(2) = -(qd(5)*v5(1));
c5(4) = qd(5)*v5(5);
c5(5) = -(qd(5)*v5(4));

c6(1) = qd(5)*v6(2);
c6(2) = -(qd(5)*v6(1));
c6(4) = qd(5)*v6(5);
c6(5) = -(qd(5)*v6(4));

c7(1) = qd(7)*v7(2);
c7(2) = -(qd(7)*v7(1));
c7(4) = qd(7)*v7(5);
c7(5) = -(qd(7)*v7(4));







%barrett_InvDynArtfunc6
     
% pv-misc vectors 
pv0(1) = -uex0.f(1) - link0.mcm(1)*power(v0(2),2) - link0.mcm(1)*power(v0(3),2) + v0(1)*(link0.mcm(2)*v0(2) + link0.mcm(3)*v0(3)) - link0.m*v0(3)*v0(5) + link0.m*v0(2)*v0(6) + g*link0.m*S00(1,3);
pv0(2) = -uex0.f(2) - link0.mcm(2)*power(v0(1),2) - link0.mcm(2)*power(v0(3),2) + v0(2)*(link0.mcm(1)*v0(1) + link0.mcm(3)*v0(3)) + link0.m*v0(3)*v0(4) - link0.m*v0(1)*v0(6) + g*link0.m*S00(2,3);
pv0(3) = -uex0.f(3) - link0.mcm(3)*power(v0(1),2) - link0.mcm(3)*power(v0(2),2) + (link0.mcm(1)*v0(1) + link0.mcm(2)*v0(2))*v0(3) - link0.m*v0(2)*v0(4) + link0.m*v0(1)*v0(5) + g*link0.m*S00(3,3);
pv0(4) = -uex0.t(1) + (-(link0.mcm(2)*v0(2)) - link0.mcm(3)*v0(3))*v0(4) + (link0.mcm(1)*v0(3) + link0.m*v0(5))*v0(6) + v0(5)*(link0.mcm(1)*v0(2) - link0.m*v0(6)) + v0(1)*(link0.mcm(2)*v0(5) + link0.mcm(3)*v0(6) - v0(3)*link0.inertia(1,2) + v0(2)*link0.inertia(1,3)) + v0(2)*(-(link0.mcm(1)*v0(5)) - v0(3)*link0.inertia(2,2) + v0(2)*link0.inertia(2,3)) + v0(3)*(-(link0.mcm(1)*v0(6)) - v0(3)*link0.inertia(2,3) + v0(2)*link0.inertia(3,3)) - g*link0.mcm(3)*S00(2,3) + g*link0.mcm(2)*S00(3,3);
pv0(5) = -uex0.t(2) + (-(link0.mcm(1)*v0(1)) - link0.mcm(3)*v0(3))*v0(5) + (link0.mcm(2)*v0(3) - link0.m*v0(4))*v0(6) + v0(4)*(link0.mcm(2)*v0(1) + link0.m*v0(6)) + v0(1)*(-(link0.mcm(2)*v0(4)) + v0(3)*link0.inertia(1,1) - v0(1)*link0.inertia(1,3)) + v0(2)*(link0.mcm(1)*v0(4) + link0.mcm(3)*v0(6) + v0(3)*link0.inertia(1,2) - v0(1)*link0.inertia(2,3)) + v0(3)*(-(link0.mcm(2)*v0(6)) + v0(3)*link0.inertia(1,3) - v0(1)*link0.inertia(3,3)) + g*link0.mcm(3)*S00(1,3) - g*link0.mcm(1)*S00(3,3);
pv0(6) = -uex0.t(3) + (link0.mcm(3)*v0(2) + link0.m*v0(4))*v0(5) + v0(4)*(link0.mcm(3)*v0(1) - link0.m*v0(5)) + (-(link0.mcm(1)*v0(1)) - link0.mcm(2)*v0(2))*v0(6) + v0(1)*(-(link0.mcm(3)*v0(4)) - v0(2)*link0.inertia(1,1) + v0(1)*link0.inertia(1,2)) + v0(2)*(-(link0.mcm(3)*v0(5)) - v0(2)*link0.inertia(1,2) + v0(1)*link0.inertia(2,2)) + v0(3)*(link0.mcm(1)*v0(4) + link0.mcm(2)*v0(5) - v0(2)*link0.inertia(1,3) + v0(1)*link0.inertia(2,3)) - g*link0.mcm(2)*S00(1,3) + g*link0.mcm(1)*S00(2,3);

pv1(1) = -uex(1).f(1) - links(1).mcm(1)*power(v1(2),2) - links(1).mcm(1)*power(v1(3),2) + v1(1)*(links(1).mcm(2)*v1(2) + links(1).mcm(3)*v1(3)) - links(1).m*v1(3)*v1(5) + links(1).m*v1(2)*v1(6) + g*links(1).m*SG10(1,3);
pv1(2) = -uex(1).f(2) - links(1).mcm(2)*power(v1(1),2) - links(1).mcm(2)*power(v1(3),2) + v1(2)*(links(1).mcm(1)*v1(1) + links(1).mcm(3)*v1(3)) + links(1).m*v1(3)*v1(4) - links(1).m*v1(1)*v1(6) + g*links(1).m*SG10(2,3);
pv1(3) = -uex(1).f(3) - links(1).mcm(3)*power(v1(1),2) - links(1).mcm(3)*power(v1(2),2) + (links(1).mcm(1)*v1(1) + links(1).mcm(2)*v1(2))*v1(3) - links(1).m*v1(2)*v1(4) + links(1).m*v1(1)*v1(5) + g*links(1).m*SG10(3,3);
pv1(4) = -uex(1).t(1) + (-(links(1).mcm(2)*v1(2)) - links(1).mcm(3)*v1(3))*v1(4) + (links(1).mcm(1)*v1(3) + links(1).m*v1(5))*v1(6) + v1(5)*(links(1).mcm(1)*v1(2) - links(1).m*v1(6)) + v1(1)*(links(1).mcm(2)*v1(5) + links(1).mcm(3)*v1(6) - v1(3)*links(1).inertia(1,2) + v1(2)*links(1).inertia(1,3)) + v1(2)*(-(links(1).mcm(1)*v1(5)) - v1(3)*links(1).inertia(2,2) + v1(2)*links(1).inertia(2,3)) + v1(3)*(-(links(1).mcm(1)*v1(6)) - v1(3)*links(1).inertia(2,3) + v1(2)*links(1).inertia(3,3)) - g*links(1).mcm(3)*SG10(2,3) + g*links(1).mcm(2)*SG10(3,3);
pv1(5) = -uex(1).t(2) + (-(links(1).mcm(1)*v1(1)) - links(1).mcm(3)*v1(3))*v1(5) + (links(1).mcm(2)*v1(3) - links(1).m*v1(4))*v1(6) + v1(4)*(links(1).mcm(2)*v1(1) + links(1).m*v1(6)) + v1(1)*(-(links(1).mcm(2)*v1(4)) + v1(3)*links(1).inertia(1,1) - v1(1)*links(1).inertia(1,3)) + v1(2)*(links(1).mcm(1)*v1(4) + links(1).mcm(3)*v1(6) + v1(3)*links(1).inertia(1,2) - v1(1)*links(1).inertia(2,3)) + v1(3)*(-(links(1).mcm(2)*v1(6)) + v1(3)*links(1).inertia(1,3) - v1(1)*links(1).inertia(3,3)) + g*links(1).mcm(3)*SG10(1,3) - g*links(1).mcm(1)*SG10(3,3);
pv1(6) = -uex(1).t(3) + (links(1).mcm(3)*v1(2) + links(1).m*v1(4))*v1(5) + v1(4)*(links(1).mcm(3)*v1(1) - links(1).m*v1(5)) + (-(links(1).mcm(1)*v1(1)) - links(1).mcm(2)*v1(2))*v1(6) + v1(1)*(-(links(1).mcm(3)*v1(4)) - v1(2)*links(1).inertia(1,1) + v1(1)*links(1).inertia(1,2)) + v1(2)*(-(links(1).mcm(3)*v1(5)) - v1(2)*links(1).inertia(1,2) + v1(1)*links(1).inertia(2,2)) + v1(3)*(links(1).mcm(1)*v1(4) + links(1).mcm(2)*v1(5) - v1(2)*links(1).inertia(1,3) + v1(1)*links(1).inertia(2,3)) - g*links(1).mcm(2)*SG10(1,3) + g*links(1).mcm(1)*SG10(2,3);

pv2(1) = -uex(2).f(1) - links(2).mcm(1)*power(v2(2),2) - links(2).mcm(1)*power(v2(3),2) + v2(1)*(links(2).mcm(2)*v2(2) + links(2).mcm(3)*v2(3)) - links(2).m*v2(3)*v2(5) + links(2).m*v2(2)*v2(6) + g*links(2).m*SG20(1,3);
pv2(2) = -uex(2).f(2) - links(2).mcm(2)*power(v2(1),2) - links(2).mcm(2)*power(v2(3),2) + v2(2)*(links(2).mcm(1)*v2(1) + links(2).mcm(3)*v2(3)) + links(2).m*v2(3)*v2(4) - links(2).m*v2(1)*v2(6) + g*links(2).m*SG20(2,3);
pv2(3) = -uex(2).f(3) - links(2).mcm(3)*power(v2(1),2) - links(2).mcm(3)*power(v2(2),2) + (links(2).mcm(1)*v2(1) + links(2).mcm(2)*v2(2))*v2(3) - links(2).m*v2(2)*v2(4) + links(2).m*v2(1)*v2(5) + g*links(2).m*SG20(3,3);
pv2(4) = -uex(2).t(1) + (-(links(2).mcm(2)*v2(2)) - links(2).mcm(3)*v2(3))*v2(4) + (links(2).mcm(1)*v2(3) + links(2).m*v2(5))*v2(6) + v2(5)*(links(2).mcm(1)*v2(2) - links(2).m*v2(6)) + v2(1)*(links(2).mcm(2)*v2(5) + links(2).mcm(3)*v2(6) - v2(3)*links(2).inertia(1,2) + v2(2)*links(2).inertia(1,3)) + v2(2)*(-(links(2).mcm(1)*v2(5)) - v2(3)*links(2).inertia(2,2) + v2(2)*links(2).inertia(2,3)) + v2(3)*(-(links(2).mcm(1)*v2(6)) - v2(3)*links(2).inertia(2,3) + v2(2)*links(2).inertia(3,3)) - g*links(2).mcm(3)*SG20(2,3) + g*links(2).mcm(2)*SG20(3,3);
pv2(5) = -uex(2).t(2) + (-(links(2).mcm(1)*v2(1)) - links(2).mcm(3)*v2(3))*v2(5) + (links(2).mcm(2)*v2(3) - links(2).m*v2(4))*v2(6) + v2(4)*(links(2).mcm(2)*v2(1) + links(2).m*v2(6)) + v2(1)*(-(links(2).mcm(2)*v2(4)) + v2(3)*links(2).inertia(1,1) - v2(1)*links(2).inertia(1,3)) + v2(2)*(links(2).mcm(1)*v2(4) + links(2).mcm(3)*v2(6) + v2(3)*links(2).inertia(1,2) - v2(1)*links(2).inertia(2,3)) + v2(3)*(-(links(2).mcm(2)*v2(6)) + v2(3)*links(2).inertia(1,3) - v2(1)*links(2).inertia(3,3)) + g*links(2).mcm(3)*SG20(1,3) - g*links(2).mcm(1)*SG20(3,3);
pv2(6) = -uex(2).t(3) + (links(2).mcm(3)*v2(2) + links(2).m*v2(4))*v2(5) + v2(4)*(links(2).mcm(3)*v2(1) - links(2).m*v2(5)) + (-(links(2).mcm(1)*v2(1)) - links(2).mcm(2)*v2(2))*v2(6) + v2(1)*(-(links(2).mcm(3)*v2(4)) - v2(2)*links(2).inertia(1,1) + v2(1)*links(2).inertia(1,2)) + v2(2)*(-(links(2).mcm(3)*v2(5)) - v2(2)*links(2).inertia(1,2) + v2(1)*links(2).inertia(2,2)) + v2(3)*(links(2).mcm(1)*v2(4) + links(2).mcm(2)*v2(5) - v2(2)*links(2).inertia(1,3) + v2(1)*links(2).inertia(2,3)) - g*links(2).mcm(2)*SG20(1,3) + g*links(2).mcm(1)*SG20(2,3);

pv3(1) = -uex(3).f(1) - links(3).mcm(1)*power(v3(2),2) - links(3).mcm(1)*power(v3(3),2) + v3(1)*(links(3).mcm(2)*v3(2) + links(3).mcm(3)*v3(3)) - links(3).m*v3(3)*v3(5) + links(3).m*v3(2)*v3(6) + g*links(3).m*SG30(1,3);
pv3(2) = -uex(3).f(2) - links(3).mcm(2)*power(v3(1),2) - links(3).mcm(2)*power(v3(3),2) + v3(2)*(links(3).mcm(1)*v3(1) + links(3).mcm(3)*v3(3)) + links(3).m*v3(3)*v3(4) - links(3).m*v3(1)*v3(6) + g*links(3).m*SG30(2,3);
pv3(3) = -uex(3).f(3) - links(3).mcm(3)*power(v3(1),2) - links(3).mcm(3)*power(v3(2),2) + (links(3).mcm(1)*v3(1) + links(3).mcm(2)*v3(2))*v3(3) - links(3).m*v3(2)*v3(4) + links(3).m*v3(1)*v3(5) + g*links(3).m*SG30(3,3);
pv3(4) = -uex(3).t(1) + (-(links(3).mcm(2)*v3(2)) - links(3).mcm(3)*v3(3))*v3(4) + (links(3).mcm(1)*v3(3) + links(3).m*v3(5))*v3(6) + v3(5)*(links(3).mcm(1)*v3(2) - links(3).m*v3(6)) + v3(1)*(links(3).mcm(2)*v3(5) + links(3).mcm(3)*v3(6) - v3(3)*links(3).inertia(1,2) + v3(2)*links(3).inertia(1,3)) + v3(2)*(-(links(3).mcm(1)*v3(5)) - v3(3)*links(3).inertia(2,2) + v3(2)*links(3).inertia(2,3)) + v3(3)*(-(links(3).mcm(1)*v3(6)) - v3(3)*links(3).inertia(2,3) + v3(2)*links(3).inertia(3,3)) - g*links(3).mcm(3)*SG30(2,3) + g*links(3).mcm(2)*SG30(3,3);
pv3(5) = -uex(3).t(2) + (-(links(3).mcm(1)*v3(1)) - links(3).mcm(3)*v3(3))*v3(5) + (links(3).mcm(2)*v3(3) - links(3).m*v3(4))*v3(6) + v3(4)*(links(3).mcm(2)*v3(1) + links(3).m*v3(6)) + v3(1)*(-(links(3).mcm(2)*v3(4)) + v3(3)*links(3).inertia(1,1) - v3(1)*links(3).inertia(1,3)) + v3(2)*(links(3).mcm(1)*v3(4) + links(3).mcm(3)*v3(6) + v3(3)*links(3).inertia(1,2) - v3(1)*links(3).inertia(2,3)) + v3(3)*(-(links(3).mcm(2)*v3(6)) + v3(3)*links(3).inertia(1,3) - v3(1)*links(3).inertia(3,3)) + g*links(3).mcm(3)*SG30(1,3) - g*links(3).mcm(1)*SG30(3,3);
pv3(6) = -uex(3).t(3) + (links(3).mcm(3)*v3(2) + links(3).m*v3(4))*v3(5) + v3(4)*(links(3).mcm(3)*v3(1) - links(3).m*v3(5)) + (-(links(3).mcm(1)*v3(1)) - links(3).mcm(2)*v3(2))*v3(6) + v3(1)*(-(links(3).mcm(3)*v3(4)) - v3(2)*links(3).inertia(1,1) + v3(1)*links(3).inertia(1,2)) + v3(2)*(-(links(3).mcm(3)*v3(5)) - v3(2)*links(3).inertia(1,2) + v3(1)*links(3).inertia(2,2)) + v3(3)*(links(3).mcm(1)*v3(4) + links(3).mcm(2)*v3(5) - v3(2)*links(3).inertia(1,3) + v3(1)*links(3).inertia(2,3)) - g*links(3).mcm(2)*SG30(1,3) + g*links(3).mcm(1)*SG30(2,3);

pv4(1) = -uex(4).f(1) - links(4).mcm(1)*power(v4(2),2) - links(4).mcm(1)*power(v4(3),2) + v4(1)*(links(4).mcm(2)*v4(2) + links(4).mcm(3)*v4(3)) - links(4).m*v4(3)*v4(5) + links(4).m*v4(2)*v4(6) + g*links(4).m*SG40(1,3);
pv4(2) = -uex(4).f(2) - links(4).mcm(2)*power(v4(1),2) - links(4).mcm(2)*power(v4(3),2) + v4(2)*(links(4).mcm(1)*v4(1) + links(4).mcm(3)*v4(3)) + links(4).m*v4(3)*v4(4) - links(4).m*v4(1)*v4(6) + g*links(4).m*SG40(2,3);
pv4(3) = -uex(4).f(3) - links(4).mcm(3)*power(v4(1),2) - links(4).mcm(3)*power(v4(2),2) + (links(4).mcm(1)*v4(1) + links(4).mcm(2)*v4(2))*v4(3) - links(4).m*v4(2)*v4(4) + links(4).m*v4(1)*v4(5) + g*links(4).m*SG40(3,3);
pv4(4) = -uex(4).t(1) + (-(links(4).mcm(2)*v4(2)) - links(4).mcm(3)*v4(3))*v4(4) + (links(4).mcm(1)*v4(3) + links(4).m*v4(5))*v4(6) + v4(5)*(links(4).mcm(1)*v4(2) - links(4).m*v4(6)) + v4(1)*(links(4).mcm(2)*v4(5) + links(4).mcm(3)*v4(6) - v4(3)*links(4).inertia(1,2) + v4(2)*links(4).inertia(1,3)) + v4(2)*(-(links(4).mcm(1)*v4(5)) - v4(3)*links(4).inertia(2,2) + v4(2)*links(4).inertia(2,3)) + v4(3)*(-(links(4).mcm(1)*v4(6)) - v4(3)*links(4).inertia(2,3) + v4(2)*links(4).inertia(3,3)) - g*links(4).mcm(3)*SG40(2,3) + g*links(4).mcm(2)*SG40(3,3);
pv4(5) = -uex(4).t(2) + (-(links(4).mcm(1)*v4(1)) - links(4).mcm(3)*v4(3))*v4(5) + (links(4).mcm(2)*v4(3) - links(4).m*v4(4))*v4(6) + v4(4)*(links(4).mcm(2)*v4(1) + links(4).m*v4(6)) + v4(1)*(-(links(4).mcm(2)*v4(4)) + v4(3)*links(4).inertia(1,1) - v4(1)*links(4).inertia(1,3)) + v4(2)*(links(4).mcm(1)*v4(4) + links(4).mcm(3)*v4(6) + v4(3)*links(4).inertia(1,2) - v4(1)*links(4).inertia(2,3)) + v4(3)*(-(links(4).mcm(2)*v4(6)) + v4(3)*links(4).inertia(1,3) - v4(1)*links(4).inertia(3,3)) + g*links(4).mcm(3)*SG40(1,3) - g*links(4).mcm(1)*SG40(3,3);
pv4(6) = -uex(4).t(3) + (links(4).mcm(3)*v4(2) + links(4).m*v4(4))*v4(5) + v4(4)*(links(4).mcm(3)*v4(1) - links(4).m*v4(5)) + (-(links(4).mcm(1)*v4(1)) - links(4).mcm(2)*v4(2))*v4(6) + v4(1)*(-(links(4).mcm(3)*v4(4)) - v4(2)*links(4).inertia(1,1) + v4(1)*links(4).inertia(1,2)) + v4(2)*(-(links(4).mcm(3)*v4(5)) - v4(2)*links(4).inertia(1,2) + v4(1)*links(4).inertia(2,2)) + v4(3)*(links(4).mcm(1)*v4(4) + links(4).mcm(2)*v4(5) - v4(2)*links(4).inertia(1,3) + v4(1)*links(4).inertia(2,3)) - g*links(4).mcm(2)*SG40(1,3) + g*links(4).mcm(1)*SG40(2,3);

pv5(1) = -uex(5).f(1) - links(5).mcm(1)*power(v5(2),2) - links(5).mcm(1)*power(v5(3),2) + v5(1)*(links(5).mcm(2)*v5(2) + links(5).mcm(3)*v5(3)) - links(5).m*v5(3)*v5(5) + links(5).m*v5(2)*v5(6) + g*links(5).m*SG50(1,3);
pv5(2) = -uex(5).f(2) - links(5).mcm(2)*power(v5(1),2) - links(5).mcm(2)*power(v5(3),2) + v5(2)*(links(5).mcm(1)*v5(1) + links(5).mcm(3)*v5(3)) + links(5).m*v5(3)*v5(4) - links(5).m*v5(1)*v5(6) + g*links(5).m*SG50(2,3);
pv5(3) = -uex(5).f(3) - links(5).mcm(3)*power(v5(1),2) - links(5).mcm(3)*power(v5(2),2) + (links(5).mcm(1)*v5(1) + links(5).mcm(2)*v5(2))*v5(3) - links(5).m*v5(2)*v5(4) + links(5).m*v5(1)*v5(5) + g*links(5).m*SG50(3,3);
pv5(4) = -uex(5).t(1) + (-(links(5).mcm(2)*v5(2)) - links(5).mcm(3)*v5(3))*v5(4) + (links(5).mcm(1)*v5(3) + links(5).m*v5(5))*v5(6) + v5(5)*(links(5).mcm(1)*v5(2) - links(5).m*v5(6)) + v5(1)*(links(5).mcm(2)*v5(5) + links(5).mcm(3)*v5(6) - v5(3)*links(5).inertia(1,2) + v5(2)*links(5).inertia(1,3)) + v5(2)*(-(links(5).mcm(1)*v5(5)) - v5(3)*links(5).inertia(2,2) + v5(2)*links(5).inertia(2,3)) + v5(3)*(-(links(5).mcm(1)*v5(6)) - v5(3)*links(5).inertia(2,3) + v5(2)*links(5).inertia(3,3)) - g*links(5).mcm(3)*SG50(2,3) + g*links(5).mcm(2)*SG50(3,3);
pv5(5) = -uex(5).t(2) + (-(links(5).mcm(1)*v5(1)) - links(5).mcm(3)*v5(3))*v5(5) + (links(5).mcm(2)*v5(3) - links(5).m*v5(4))*v5(6) + v5(4)*(links(5).mcm(2)*v5(1) + links(5).m*v5(6)) + v5(1)*(-(links(5).mcm(2)*v5(4)) + v5(3)*links(5).inertia(1,1) - v5(1)*links(5).inertia(1,3)) + v5(2)*(links(5).mcm(1)*v5(4) + links(5).mcm(3)*v5(6) + v5(3)*links(5).inertia(1,2) - v5(1)*links(5).inertia(2,3)) + v5(3)*(-(links(5).mcm(2)*v5(6)) + v5(3)*links(5).inertia(1,3) - v5(1)*links(5).inertia(3,3)) + g*links(5).mcm(3)*SG50(1,3) - g*links(5).mcm(1)*SG50(3,3);
pv5(6) = -uex(5).t(3) + (links(5).mcm(3)*v5(2) + links(5).m*v5(4))*v5(5) + v5(4)*(links(5).mcm(3)*v5(1) - links(5).m*v5(5)) + (-(links(5).mcm(1)*v5(1)) - links(5).mcm(2)*v5(2))*v5(6) + v5(1)*(-(links(5).mcm(3)*v5(4)) - v5(2)*links(5).inertia(1,1) + v5(1)*links(5).inertia(1,2)) + v5(2)*(-(links(5).mcm(3)*v5(5)) - v5(2)*links(5).inertia(1,2) + v5(1)*links(5).inertia(2,2)) + v5(3)*(links(5).mcm(1)*v5(4) + links(5).mcm(2)*v5(5) - v5(2)*links(5).inertia(1,3) + v5(1)*links(5).inertia(2,3)) - g*links(5).mcm(2)*SG50(1,3) + g*links(5).mcm(1)*SG50(2,3);

pv6(1) = -uex(6).f(1) - links(6).mcm(1)*power(v6(2),2) - links(6).mcm(1)*power(v6(3),2) + v6(1)*(links(6).mcm(2)*v6(2) + links(6).mcm(3)*v6(3)) - links(6).m*v6(3)*v6(5) + links(6).m*v6(2)*v6(6) + g*links(6).m*SG60(1,3);
pv6(2) = -uex(6).f(2) - links(6).mcm(2)*power(v6(1),2) - links(6).mcm(2)*power(v6(3),2) + v6(2)*(links(6).mcm(1)*v6(1) + links(6).mcm(3)*v6(3)) + links(6).m*v6(3)*v6(4) - links(6).m*v6(1)*v6(6) + g*links(6).m*SG60(2,3);
pv6(3) = -uex(6).f(3) - links(6).mcm(3)*power(v6(1),2) - links(6).mcm(3)*power(v6(2),2) + (links(6).mcm(1)*v6(1) + links(6).mcm(2)*v6(2))*v6(3) - links(6).m*v6(2)*v6(4) + links(6).m*v6(1)*v6(5) + g*links(6).m*SG60(3,3);
pv6(4) = -uex(6).t(1) + (-(links(6).mcm(2)*v6(2)) - links(6).mcm(3)*v6(3))*v6(4) + (links(6).mcm(1)*v6(3) + links(6).m*v6(5))*v6(6) + v6(5)*(links(6).mcm(1)*v6(2) - links(6).m*v6(6)) + v6(1)*(links(6).mcm(2)*v6(5) + links(6).mcm(3)*v6(6) - v6(3)*links(6).inertia(1,2) + v6(2)*links(6).inertia(1,3)) + v6(2)*(-(links(6).mcm(1)*v6(5)) - v6(3)*links(6).inertia(2,2) + v6(2)*links(6).inertia(2,3)) + v6(3)*(-(links(6).mcm(1)*v6(6)) - v6(3)*links(6).inertia(2,3) + v6(2)*links(6).inertia(3,3)) - g*links(6).mcm(3)*SG60(2,3) + g*links(6).mcm(2)*SG60(3,3);
pv6(5) = -uex(6).t(2) + (-(links(6).mcm(1)*v6(1)) - links(6).mcm(3)*v6(3))*v6(5) + (links(6).mcm(2)*v6(3) - links(6).m*v6(4))*v6(6) + v6(4)*(links(6).mcm(2)*v6(1) + links(6).m*v6(6)) + v6(1)*(-(links(6).mcm(2)*v6(4)) + v6(3)*links(6).inertia(1,1) - v6(1)*links(6).inertia(1,3)) + v6(2)*(links(6).mcm(1)*v6(4) + links(6).mcm(3)*v6(6) + v6(3)*links(6).inertia(1,2) - v6(1)*links(6).inertia(2,3)) + v6(3)*(-(links(6).mcm(2)*v6(6)) + v6(3)*links(6).inertia(1,3) - v6(1)*links(6).inertia(3,3)) + g*links(6).mcm(3)*SG60(1,3) - g*links(6).mcm(1)*SG60(3,3);
pv6(6) = -uex(6).t(3) + (links(6).mcm(3)*v6(2) + links(6).m*v6(4))*v6(5) + v6(4)*(links(6).mcm(3)*v6(1) - links(6).m*v6(5)) + (-(links(6).mcm(1)*v6(1)) - links(6).mcm(2)*v6(2))*v6(6) + v6(1)*(-(links(6).mcm(3)*v6(4)) - v6(2)*links(6).inertia(1,1) + v6(1)*links(6).inertia(1,2)) + v6(2)*(-(links(6).mcm(3)*v6(5)) - v6(2)*links(6).inertia(1,2) + v6(1)*links(6).inertia(2,2)) + v6(3)*(links(6).mcm(1)*v6(4) + links(6).mcm(2)*v6(5) - v6(2)*links(6).inertia(1,3) + v6(1)*links(6).inertia(2,3)) - g*links(6).mcm(2)*SG60(1,3) + g*links(6).mcm(1)*SG60(2,3);

pv7(1) = -uex(7).f(1) - links(7).mcm(1)*power(v7(2),2) - links(7).mcm(1)*power(v7(3),2) + v7(1)*(links(7).mcm(2)*v7(2) + links(7).mcm(3)*v7(3)) - links(7).m*v7(3)*v7(5) + links(7).m*v7(2)*v7(6) + g*links(7).m*SG70(1,3);
pv7(2) = -uex(7).f(2) - links(7).mcm(2)*power(v7(1),2) - links(7).mcm(2)*power(v7(3),2) + v7(2)*(links(7).mcm(1)*v7(1) + links(7).mcm(3)*v7(3)) + links(7).m*v7(3)*v7(4) - links(7).m*v7(1)*v7(6) + g*links(7).m*SG70(2,3);
pv7(3) = -uex(7).f(3) - links(7).mcm(3)*power(v7(1),2) - links(7).mcm(3)*power(v7(2),2) + (links(7).mcm(1)*v7(1) + links(7).mcm(2)*v7(2))*v7(3) - links(7).m*v7(2)*v7(4) + links(7).m*v7(1)*v7(5) + g*links(7).m*SG70(3,3);
pv7(4) = -uex(7).t(1) + (-(links(7).mcm(2)*v7(2)) - links(7).mcm(3)*v7(3))*v7(4) + (links(7).mcm(1)*v7(3) + links(7).m*v7(5))*v7(6) + v7(5)*(links(7).mcm(1)*v7(2) - links(7).m*v7(6)) + v7(1)*(links(7).mcm(2)*v7(5) + links(7).mcm(3)*v7(6) - v7(3)*links(7).inertia(1,2) + v7(2)*links(7).inertia(1,3)) + v7(2)*(-(links(7).mcm(1)*v7(5)) - v7(3)*links(7).inertia(2,2) + v7(2)*links(7).inertia(2,3)) + v7(3)*(-(links(7).mcm(1)*v7(6)) - v7(3)*links(7).inertia(2,3) + v7(2)*links(7).inertia(3,3)) - g*links(7).mcm(3)*SG70(2,3) + g*links(7).mcm(2)*SG70(3,3);
pv7(5) = -uex(7).t(2) + (-(links(7).mcm(1)*v7(1)) - links(7).mcm(3)*v7(3))*v7(5) + (links(7).mcm(2)*v7(3) - links(7).m*v7(4))*v7(6) + v7(4)*(links(7).mcm(2)*v7(1) + links(7).m*v7(6)) + v7(1)*(-(links(7).mcm(2)*v7(4)) + v7(3)*links(7).inertia(1,1) - v7(1)*links(7).inertia(1,3)) + v7(2)*(links(7).mcm(1)*v7(4) + links(7).mcm(3)*v7(6) + v7(3)*links(7).inertia(1,2) - v7(1)*links(7).inertia(2,3)) + v7(3)*(-(links(7).mcm(2)*v7(6)) + v7(3)*links(7).inertia(1,3) - v7(1)*links(7).inertia(3,3)) + g*links(7).mcm(3)*SG70(1,3) - g*links(7).mcm(1)*SG70(3,3);
pv7(6) = -uex(7).t(3) + (links(7).mcm(3)*v7(2) + links(7).m*v7(4))*v7(5) + v7(4)*(links(7).mcm(3)*v7(1) - links(7).m*v7(5)) + (-(links(7).mcm(1)*v7(1)) - links(7).mcm(2)*v7(2))*v7(6) + v7(1)*(-(links(7).mcm(3)*v7(4)) - v7(2)*links(7).inertia(1,1) + v7(1)*links(7).inertia(1,2)) + v7(2)*(-(links(7).mcm(3)*v7(5)) - v7(2)*links(7).inertia(1,2) + v7(1)*links(7).inertia(2,2)) + v7(3)*(links(7).mcm(1)*v7(4) + links(7).mcm(2)*v7(5) - v7(2)*links(7).inertia(1,3) + v7(1)*links(7).inertia(2,3)) - g*links(7).mcm(2)*SG70(1,3) + g*links(7).mcm(1)*SG70(2,3);

pv8(1) = -(eff(1).mcm(1)*power(v8(2),2)) - eff(1).mcm(1)*power(v8(3),2) + v8(1)*(eff(1).mcm(2)*v8(2) + eff(1).mcm(3)*v8(3)) - eff(1).m*v8(3)*v8(5) + eff(1).m*v8(2)*v8(6) + eff(1).m*g*SG80(1,3);
pv8(2) = -(eff(1).mcm(2)*power(v8(1),2)) - eff(1).mcm(2)*power(v8(3),2) + v8(2)*(eff(1).mcm(1)*v8(1) + eff(1).mcm(3)*v8(3)) + eff(1).m*v8(3)*v8(4) - eff(1).m*v8(1)*v8(6) + eff(1).m*g*SG80(2,3);
pv8(3) = -(eff(1).mcm(3)*power(v8(1),2)) - eff(1).mcm(3)*power(v8(2),2) + (eff(1).mcm(1)*v8(1) + eff(1).mcm(2)*v8(2))*v8(3) - eff(1).m*v8(2)*v8(4) + eff(1).m*v8(1)*v8(5) + eff(1).m*g*SG80(3,3);
pv8(4) = (-(eff(1).mcm(2)*v8(2)) - eff(1).mcm(3)*v8(3))*v8(4) - eff(1).mcm(1)*v8(2)*v8(5) - eff(1).mcm(1)*v8(3)*v8(6) + (eff(1).mcm(1)*v8(3) + eff(1).m*v8(5))*v8(6) + v8(5)*(eff(1).mcm(1)*v8(2) - eff(1).m*v8(6)) + v8(1)*(eff(1).mcm(2)*v8(5) + eff(1).mcm(3)*v8(6)) - g*eff(1).mcm(3)*SG80(2,3) + g*eff(1).mcm(2)*SG80(3,3);
pv8(5) = -(eff(1).mcm(2)*v8(1)*v8(4)) + (-(eff(1).mcm(1)*v8(1)) - eff(1).mcm(3)*v8(3))*v8(5) - eff(1).mcm(2)*v8(3)*v8(6) + (eff(1).mcm(2)*v8(3) - eff(1).m*v8(4))*v8(6) + v8(4)*(eff(1).mcm(2)*v8(1) + eff(1).m*v8(6)) + v8(2)*(eff(1).mcm(1)*v8(4) + eff(1).mcm(3)*v8(6)) + g*eff(1).mcm(3)*SG80(1,3) - g*eff(1).mcm(1)*SG80(3,3);
pv8(6) = -(eff(1).mcm(3)*v8(1)*v8(4)) - eff(1).mcm(3)*v8(2)*v8(5) + (eff(1).mcm(3)*v8(2) + eff(1).m*v8(4))*v8(5) + v8(4)*(eff(1).mcm(3)*v8(1) - eff(1).m*v8(5)) + v8(3)*(eff(1).mcm(1)*v8(4) + eff(1).mcm(2)*v8(5)) + (-(eff(1).mcm(1)*v8(1)) - eff(1).mcm(2)*v8(2))*v8(6) - g*eff(1).mcm(2)*SG80(1,3) + g*eff(1).mcm(1)*SG80(2,3);





% articulated body inertias and misc variables 


%barrett_InvDynArtfunc7
     
JA8(1,2) = eff(1).mcm(3);
JA8(1,3) = -eff(1).mcm(2);
JA8(1,4) = eff(1).m;

JA8(2,1) = -eff(1).mcm(3);
JA8(2,3) = eff(1).mcm(1);
JA8(2,5) = eff(1).m;

JA8(3,1) = eff(1).mcm(2);
JA8(3,2) = -eff(1).mcm(1);
JA8(3,6) = eff(1).m;

JA8(4,5) = -eff(1).mcm(3);
JA8(4,6) = eff(1).mcm(2);

JA8(5,4) = eff(1).mcm(3);
JA8(5,6) = -eff(1).mcm(1);

JA8(6,4) = -eff(1).mcm(2);
JA8(6,5) = eff(1).mcm(1);


T178(1,2) = JA8(1,2);
T178(1,3) = JA8(1,3);
T178(1,4) = JA8(1,4);

T178(2,1) = JA8(2,1);
T178(2,3) = JA8(2,3);
T178(2,5) = JA8(2,5);

T178(3,1) = JA8(3,1);
T178(3,2) = JA8(3,2);
T178(3,6) = JA8(3,6);

T178(4,5) = JA8(4,5);
T178(4,6) = JA8(4,6);

T178(5,4) = JA8(5,4);
T178(5,6) = JA8(5,6);

T178(6,4) = JA8(6,4);
T178(6,5) = JA8(6,5);


T78(1,1) = (-(eff(1).x(3)*S87(1,2)) + eff(1).x(2)*S87(1,3))*Si78(1,1)*T178(1,4) + S87(3,1)*(Si78(1,1)*T178(1,3) + Si78(1,2)*T178(2,3)) + (-(eff(1).x(3)*S87(2,2)) + eff(1).x(2)*S87(2,3))*Si78(1,2)*T178(2,5) + S87(1,1)*(Si78(1,2)*T178(2,1) + Si78(1,3)*T178(3,1)) + S87(2,1)*(Si78(1,1)*T178(1,2) + Si78(1,3)*T178(3,2)) + (-(eff(1).x(3)*S87(3,2)) + eff(1).x(2)*S87(3,3))*Si78(1,3)*T178(3,6);
T78(1,2) = (eff(1).x(3)*S87(1,1) - eff(1).x(1)*S87(1,3))*Si78(1,1)*T178(1,4) + S87(3,2)*(Si78(1,1)*T178(1,3) + Si78(1,2)*T178(2,3)) + (eff(1).x(3)*S87(2,1) - eff(1).x(1)*S87(2,3))*Si78(1,2)*T178(2,5) + S87(1,2)*(Si78(1,2)*T178(2,1) + Si78(1,3)*T178(3,1)) + S87(2,2)*(Si78(1,1)*T178(1,2) + Si78(1,3)*T178(3,2)) + (eff(1).x(3)*S87(3,1) - eff(1).x(1)*S87(3,3))*Si78(1,3)*T178(3,6);
T78(1,3) = (-(eff(1).x(2)*S87(1,1)) + eff(1).x(1)*S87(1,2))*Si78(1,1)*T178(1,4) + S87(3,3)*(Si78(1,1)*T178(1,3) + Si78(1,2)*T178(2,3)) + (-(eff(1).x(2)*S87(2,1)) + eff(1).x(1)*S87(2,2))*Si78(1,2)*T178(2,5) + S87(1,3)*(Si78(1,2)*T178(2,1) + Si78(1,3)*T178(3,1)) + S87(2,3)*(Si78(1,1)*T178(1,2) + Si78(1,3)*T178(3,2)) + (-(eff(1).x(2)*S87(3,1)) + eff(1).x(1)*S87(3,2))*Si78(1,3)*T178(3,6);
T78(1,4) = S87(1,1)*Si78(1,1)*T178(1,4) + S87(2,1)*Si78(1,2)*T178(2,5) + S87(3,1)*Si78(1,3)*T178(3,6);
T78(1,5) = S87(1,2)*Si78(1,1)*T178(1,4) + S87(2,2)*Si78(1,2)*T178(2,5) + S87(3,2)*Si78(1,3)*T178(3,6);
T78(1,6) = S87(1,3)*Si78(1,1)*T178(1,4) + S87(2,3)*Si78(1,2)*T178(2,5) + S87(3,3)*Si78(1,3)*T178(3,6);

T78(2,1) = (-(eff(1).x(3)*S87(1,2)) + eff(1).x(2)*S87(1,3))*Si78(2,1)*T178(1,4) + S87(3,1)*(Si78(2,1)*T178(1,3) + Si78(2,2)*T178(2,3)) + (-(eff(1).x(3)*S87(2,2)) + eff(1).x(2)*S87(2,3))*Si78(2,2)*T178(2,5) + S87(1,1)*(Si78(2,2)*T178(2,1) + Si78(2,3)*T178(3,1)) + S87(2,1)*(Si78(2,1)*T178(1,2) + Si78(2,3)*T178(3,2)) + (-(eff(1).x(3)*S87(3,2)) + eff(1).x(2)*S87(3,3))*Si78(2,3)*T178(3,6);
T78(2,2) = (eff(1).x(3)*S87(1,1) - eff(1).x(1)*S87(1,3))*Si78(2,1)*T178(1,4) + S87(3,2)*(Si78(2,1)*T178(1,3) + Si78(2,2)*T178(2,3)) + (eff(1).x(3)*S87(2,1) - eff(1).x(1)*S87(2,3))*Si78(2,2)*T178(2,5) + S87(1,2)*(Si78(2,2)*T178(2,1) + Si78(2,3)*T178(3,1)) + S87(2,2)*(Si78(2,1)*T178(1,2) + Si78(2,3)*T178(3,2)) + (eff(1).x(3)*S87(3,1) - eff(1).x(1)*S87(3,3))*Si78(2,3)*T178(3,6);
T78(2,3) = (-(eff(1).x(2)*S87(1,1)) + eff(1).x(1)*S87(1,2))*Si78(2,1)*T178(1,4) + S87(3,3)*(Si78(2,1)*T178(1,3) + Si78(2,2)*T178(2,3)) + (-(eff(1).x(2)*S87(2,1)) + eff(1).x(1)*S87(2,2))*Si78(2,2)*T178(2,5) + S87(1,3)*(Si78(2,2)*T178(2,1) + Si78(2,3)*T178(3,1)) + S87(2,3)*(Si78(2,1)*T178(1,2) + Si78(2,3)*T178(3,2)) + (-(eff(1).x(2)*S87(3,1)) + eff(1).x(1)*S87(3,2))*Si78(2,3)*T178(3,6);
T78(2,4) = S87(1,1)*Si78(2,1)*T178(1,4) + S87(2,1)*Si78(2,2)*T178(2,5) + S87(3,1)*Si78(2,3)*T178(3,6);
T78(2,5) = S87(1,2)*Si78(2,1)*T178(1,4) + S87(2,2)*Si78(2,2)*T178(2,5) + S87(3,2)*Si78(2,3)*T178(3,6);
T78(2,6) = S87(1,3)*Si78(2,1)*T178(1,4) + S87(2,3)*Si78(2,2)*T178(2,5) + S87(3,3)*Si78(2,3)*T178(3,6);

T78(3,1) = (-(eff(1).x(3)*S87(1,2)) + eff(1).x(2)*S87(1,3))*Si78(3,1)*T178(1,4) + S87(3,1)*(Si78(3,1)*T178(1,3) + Si78(3,2)*T178(2,3)) + (-(eff(1).x(3)*S87(2,2)) + eff(1).x(2)*S87(2,3))*Si78(3,2)*T178(2,5) + S87(1,1)*(Si78(3,2)*T178(2,1) + Si78(3,3)*T178(3,1)) + S87(2,1)*(Si78(3,1)*T178(1,2) + Si78(3,3)*T178(3,2)) + (-(eff(1).x(3)*S87(3,2)) + eff(1).x(2)*S87(3,3))*Si78(3,3)*T178(3,6);
T78(3,2) = (eff(1).x(3)*S87(1,1) - eff(1).x(1)*S87(1,3))*Si78(3,1)*T178(1,4) + S87(3,2)*(Si78(3,1)*T178(1,3) + Si78(3,2)*T178(2,3)) + (eff(1).x(3)*S87(2,1) - eff(1).x(1)*S87(2,3))*Si78(3,2)*T178(2,5) + S87(1,2)*(Si78(3,2)*T178(2,1) + Si78(3,3)*T178(3,1)) + S87(2,2)*(Si78(3,1)*T178(1,2) + Si78(3,3)*T178(3,2)) + (eff(1).x(3)*S87(3,1) - eff(1).x(1)*S87(3,3))*Si78(3,3)*T178(3,6);
T78(3,3) = (-(eff(1).x(2)*S87(1,1)) + eff(1).x(1)*S87(1,2))*Si78(3,1)*T178(1,4) + S87(3,3)*(Si78(3,1)*T178(1,3) + Si78(3,2)*T178(2,3)) + (-(eff(1).x(2)*S87(2,1)) + eff(1).x(1)*S87(2,2))*Si78(3,2)*T178(2,5) + S87(1,3)*(Si78(3,2)*T178(2,1) + Si78(3,3)*T178(3,1)) + S87(2,3)*(Si78(3,1)*T178(1,2) + Si78(3,3)*T178(3,2)) + (-(eff(1).x(2)*S87(3,1)) + eff(1).x(1)*S87(3,2))*Si78(3,3)*T178(3,6);
T78(3,4) = S87(1,1)*Si78(3,1)*T178(1,4) + S87(2,1)*Si78(3,2)*T178(2,5) + S87(3,1)*Si78(3,3)*T178(3,6);
T78(3,5) = S87(1,2)*Si78(3,1)*T178(1,4) + S87(2,2)*Si78(3,2)*T178(2,5) + S87(3,2)*Si78(3,3)*T178(3,6);
T78(3,6) = S87(1,3)*Si78(3,1)*T178(1,4) + S87(2,3)*Si78(3,2)*T178(2,5) + S87(3,3)*Si78(3,3)*T178(3,6);

T78(4,1) = S87(3,1)*((-(eff(1).x(3)*Si78(2,1)) + eff(1).x(2)*Si78(3,1))*T178(1,3) + (-(eff(1).x(3)*Si78(2,2)) + eff(1).x(2)*Si78(3,2))*T178(2,3)) + S87(1,1)*((-(eff(1).x(3)*Si78(2,2)) + eff(1).x(2)*Si78(3,2))*T178(2,1) + (-(eff(1).x(3)*Si78(2,3)) + eff(1).x(2)*Si78(3,3))*T178(3,1)) + S87(2,1)*((-(eff(1).x(3)*Si78(2,1)) + eff(1).x(2)*Si78(3,1))*T178(1,2) + (-(eff(1).x(3)*Si78(2,3)) + eff(1).x(2)*Si78(3,3))*T178(3,2)) + (-(eff(1).x(3)*S87(3,2)) + eff(1).x(2)*S87(3,3))*((-(eff(1).x(3)*Si78(2,3)) + eff(1).x(2)*Si78(3,3))*T178(3,6) + Si78(1,1)*T178(4,6) + Si78(1,2)*T178(5,6)) + (-(eff(1).x(3)*S87(1,2)) + eff(1).x(2)*S87(1,3))*((-(eff(1).x(3)*Si78(2,1)) + eff(1).x(2)*Si78(3,1))*T178(1,4) + Si78(1,2)*T178(5,4) + Si78(1,3)*T178(6,4)) + (-(eff(1).x(3)*S87(2,2)) + eff(1).x(2)*S87(2,3))*((-(eff(1).x(3)*Si78(2,2)) + eff(1).x(2)*Si78(3,2))*T178(2,5) + Si78(1,1)*T178(4,5) + Si78(1,3)*T178(6,5));
T78(4,2) = S87(3,2)*((-(eff(1).x(3)*Si78(2,1)) + eff(1).x(2)*Si78(3,1))*T178(1,3) + (-(eff(1).x(3)*Si78(2,2)) + eff(1).x(2)*Si78(3,2))*T178(2,3)) + S87(1,2)*((-(eff(1).x(3)*Si78(2,2)) + eff(1).x(2)*Si78(3,2))*T178(2,1) + (-(eff(1).x(3)*Si78(2,3)) + eff(1).x(2)*Si78(3,3))*T178(3,1)) + S87(2,2)*((-(eff(1).x(3)*Si78(2,1)) + eff(1).x(2)*Si78(3,1))*T178(1,2) + (-(eff(1).x(3)*Si78(2,3)) + eff(1).x(2)*Si78(3,3))*T178(3,2)) + (eff(1).x(3)*S87(3,1) - eff(1).x(1)*S87(3,3))*((-(eff(1).x(3)*Si78(2,3)) + eff(1).x(2)*Si78(3,3))*T178(3,6) + Si78(1,1)*T178(4,6) + Si78(1,2)*T178(5,6)) + (eff(1).x(3)*S87(1,1) - eff(1).x(1)*S87(1,3))*((-(eff(1).x(3)*Si78(2,1)) + eff(1).x(2)*Si78(3,1))*T178(1,4) + Si78(1,2)*T178(5,4) + Si78(1,3)*T178(6,4)) + (eff(1).x(3)*S87(2,1) - eff(1).x(1)*S87(2,3))*((-(eff(1).x(3)*Si78(2,2)) + eff(1).x(2)*Si78(3,2))*T178(2,5) + Si78(1,1)*T178(4,5) + Si78(1,3)*T178(6,5));
T78(4,3) = S87(3,3)*((-(eff(1).x(3)*Si78(2,1)) + eff(1).x(2)*Si78(3,1))*T178(1,3) + (-(eff(1).x(3)*Si78(2,2)) + eff(1).x(2)*Si78(3,2))*T178(2,3)) + S87(1,3)*((-(eff(1).x(3)*Si78(2,2)) + eff(1).x(2)*Si78(3,2))*T178(2,1) + (-(eff(1).x(3)*Si78(2,3)) + eff(1).x(2)*Si78(3,3))*T178(3,1)) + S87(2,3)*((-(eff(1).x(3)*Si78(2,1)) + eff(1).x(2)*Si78(3,1))*T178(1,2) + (-(eff(1).x(3)*Si78(2,3)) + eff(1).x(2)*Si78(3,3))*T178(3,2)) + (-(eff(1).x(2)*S87(3,1)) + eff(1).x(1)*S87(3,2))*((-(eff(1).x(3)*Si78(2,3)) + eff(1).x(2)*Si78(3,3))*T178(3,6) + Si78(1,1)*T178(4,6) + Si78(1,2)*T178(5,6)) + (-(eff(1).x(2)*S87(1,1)) + eff(1).x(1)*S87(1,2))*((-(eff(1).x(3)*Si78(2,1)) + eff(1).x(2)*Si78(3,1))*T178(1,4) + Si78(1,2)*T178(5,4) + Si78(1,3)*T178(6,4)) + (-(eff(1).x(2)*S87(2,1)) + eff(1).x(1)*S87(2,2))*((-(eff(1).x(3)*Si78(2,2)) + eff(1).x(2)*Si78(3,2))*T178(2,5) + Si78(1,1)*T178(4,5) + Si78(1,3)*T178(6,5));
T78(4,4) = S87(3,1)*((-(eff(1).x(3)*Si78(2,3)) + eff(1).x(2)*Si78(3,3))*T178(3,6) + Si78(1,1)*T178(4,6) + Si78(1,2)*T178(5,6)) + S87(1,1)*((-(eff(1).x(3)*Si78(2,1)) + eff(1).x(2)*Si78(3,1))*T178(1,4) + Si78(1,2)*T178(5,4) + Si78(1,3)*T178(6,4)) + S87(2,1)*((-(eff(1).x(3)*Si78(2,2)) + eff(1).x(2)*Si78(3,2))*T178(2,5) + Si78(1,1)*T178(4,5) + Si78(1,3)*T178(6,5));
T78(4,5) = S87(3,2)*((-(eff(1).x(3)*Si78(2,3)) + eff(1).x(2)*Si78(3,3))*T178(3,6) + Si78(1,1)*T178(4,6) + Si78(1,2)*T178(5,6)) + S87(1,2)*((-(eff(1).x(3)*Si78(2,1)) + eff(1).x(2)*Si78(3,1))*T178(1,4) + Si78(1,2)*T178(5,4) + Si78(1,3)*T178(6,4)) + S87(2,2)*((-(eff(1).x(3)*Si78(2,2)) + eff(1).x(2)*Si78(3,2))*T178(2,5) + Si78(1,1)*T178(4,5) + Si78(1,3)*T178(6,5));
T78(4,6) = S87(3,3)*((-(eff(1).x(3)*Si78(2,3)) + eff(1).x(2)*Si78(3,3))*T178(3,6) + Si78(1,1)*T178(4,6) + Si78(1,2)*T178(5,6)) + S87(1,3)*((-(eff(1).x(3)*Si78(2,1)) + eff(1).x(2)*Si78(3,1))*T178(1,4) + Si78(1,2)*T178(5,4) + Si78(1,3)*T178(6,4)) + S87(2,3)*((-(eff(1).x(3)*Si78(2,2)) + eff(1).x(2)*Si78(3,2))*T178(2,5) + Si78(1,1)*T178(4,5) + Si78(1,3)*T178(6,5));

T78(5,1) = S87(3,1)*((eff(1).x(3)*Si78(1,1) - eff(1).x(1)*Si78(3,1))*T178(1,3) + (eff(1).x(3)*Si78(1,2) - eff(1).x(1)*Si78(3,2))*T178(2,3)) + S87(1,1)*((eff(1).x(3)*Si78(1,2) - eff(1).x(1)*Si78(3,2))*T178(2,1) + (eff(1).x(3)*Si78(1,3) - eff(1).x(1)*Si78(3,3))*T178(3,1)) + S87(2,1)*((eff(1).x(3)*Si78(1,1) - eff(1).x(1)*Si78(3,1))*T178(1,2) + (eff(1).x(3)*Si78(1,3) - eff(1).x(1)*Si78(3,3))*T178(3,2)) + (-(eff(1).x(3)*S87(3,2)) + eff(1).x(2)*S87(3,3))*((eff(1).x(3)*Si78(1,3) - eff(1).x(1)*Si78(3,3))*T178(3,6) + Si78(2,1)*T178(4,6) + Si78(2,2)*T178(5,6)) + (-(eff(1).x(3)*S87(1,2)) + eff(1).x(2)*S87(1,3))*((eff(1).x(3)*Si78(1,1) - eff(1).x(1)*Si78(3,1))*T178(1,4) + Si78(2,2)*T178(5,4) + Si78(2,3)*T178(6,4)) + (-(eff(1).x(3)*S87(2,2)) + eff(1).x(2)*S87(2,3))*((eff(1).x(3)*Si78(1,2) - eff(1).x(1)*Si78(3,2))*T178(2,5) + Si78(2,1)*T178(4,5) + Si78(2,3)*T178(6,5));
T78(5,2) = S87(3,2)*((eff(1).x(3)*Si78(1,1) - eff(1).x(1)*Si78(3,1))*T178(1,3) + (eff(1).x(3)*Si78(1,2) - eff(1).x(1)*Si78(3,2))*T178(2,3)) + S87(1,2)*((eff(1).x(3)*Si78(1,2) - eff(1).x(1)*Si78(3,2))*T178(2,1) + (eff(1).x(3)*Si78(1,3) - eff(1).x(1)*Si78(3,3))*T178(3,1)) + S87(2,2)*((eff(1).x(3)*Si78(1,1) - eff(1).x(1)*Si78(3,1))*T178(1,2) + (eff(1).x(3)*Si78(1,3) - eff(1).x(1)*Si78(3,3))*T178(3,2)) + (eff(1).x(3)*S87(3,1) - eff(1).x(1)*S87(3,3))*((eff(1).x(3)*Si78(1,3) - eff(1).x(1)*Si78(3,3))*T178(3,6) + Si78(2,1)*T178(4,6) + Si78(2,2)*T178(5,6)) + (eff(1).x(3)*S87(1,1) - eff(1).x(1)*S87(1,3))*((eff(1).x(3)*Si78(1,1) - eff(1).x(1)*Si78(3,1))*T178(1,4) + Si78(2,2)*T178(5,4) + Si78(2,3)*T178(6,4)) + (eff(1).x(3)*S87(2,1) - eff(1).x(1)*S87(2,3))*((eff(1).x(3)*Si78(1,2) - eff(1).x(1)*Si78(3,2))*T178(2,5) + Si78(2,1)*T178(4,5) + Si78(2,3)*T178(6,5));
T78(5,3) = S87(3,3)*((eff(1).x(3)*Si78(1,1) - eff(1).x(1)*Si78(3,1))*T178(1,3) + (eff(1).x(3)*Si78(1,2) - eff(1).x(1)*Si78(3,2))*T178(2,3)) + S87(1,3)*((eff(1).x(3)*Si78(1,2) - eff(1).x(1)*Si78(3,2))*T178(2,1) + (eff(1).x(3)*Si78(1,3) - eff(1).x(1)*Si78(3,3))*T178(3,1)) + S87(2,3)*((eff(1).x(3)*Si78(1,1) - eff(1).x(1)*Si78(3,1))*T178(1,2) + (eff(1).x(3)*Si78(1,3) - eff(1).x(1)*Si78(3,3))*T178(3,2)) + (-(eff(1).x(2)*S87(3,1)) + eff(1).x(1)*S87(3,2))*((eff(1).x(3)*Si78(1,3) - eff(1).x(1)*Si78(3,3))*T178(3,6) + Si78(2,1)*T178(4,6) + Si78(2,2)*T178(5,6)) + (-(eff(1).x(2)*S87(1,1)) + eff(1).x(1)*S87(1,2))*((eff(1).x(3)*Si78(1,1) - eff(1).x(1)*Si78(3,1))*T178(1,4) + Si78(2,2)*T178(5,4) + Si78(2,3)*T178(6,4)) + (-(eff(1).x(2)*S87(2,1)) + eff(1).x(1)*S87(2,2))*((eff(1).x(3)*Si78(1,2) - eff(1).x(1)*Si78(3,2))*T178(2,5) + Si78(2,1)*T178(4,5) + Si78(2,3)*T178(6,5));
T78(5,4) = S87(3,1)*((eff(1).x(3)*Si78(1,3) - eff(1).x(1)*Si78(3,3))*T178(3,6) + Si78(2,1)*T178(4,6) + Si78(2,2)*T178(5,6)) + S87(1,1)*((eff(1).x(3)*Si78(1,1) - eff(1).x(1)*Si78(3,1))*T178(1,4) + Si78(2,2)*T178(5,4) + Si78(2,3)*T178(6,4)) + S87(2,1)*((eff(1).x(3)*Si78(1,2) - eff(1).x(1)*Si78(3,2))*T178(2,5) + Si78(2,1)*T178(4,5) + Si78(2,3)*T178(6,5));
T78(5,5) = S87(3,2)*((eff(1).x(3)*Si78(1,3) - eff(1).x(1)*Si78(3,3))*T178(3,6) + Si78(2,1)*T178(4,6) + Si78(2,2)*T178(5,6)) + S87(1,2)*((eff(1).x(3)*Si78(1,1) - eff(1).x(1)*Si78(3,1))*T178(1,4) + Si78(2,2)*T178(5,4) + Si78(2,3)*T178(6,4)) + S87(2,2)*((eff(1).x(3)*Si78(1,2) - eff(1).x(1)*Si78(3,2))*T178(2,5) + Si78(2,1)*T178(4,5) + Si78(2,3)*T178(6,5));
T78(5,6) = S87(3,3)*((eff(1).x(3)*Si78(1,3) - eff(1).x(1)*Si78(3,3))*T178(3,6) + Si78(2,1)*T178(4,6) + Si78(2,2)*T178(5,6)) + S87(1,3)*((eff(1).x(3)*Si78(1,1) - eff(1).x(1)*Si78(3,1))*T178(1,4) + Si78(2,2)*T178(5,4) + Si78(2,3)*T178(6,4)) + S87(2,3)*((eff(1).x(3)*Si78(1,2) - eff(1).x(1)*Si78(3,2))*T178(2,5) + Si78(2,1)*T178(4,5) + Si78(2,3)*T178(6,5));

T78(6,1) = S87(3,1)*((-(eff(1).x(2)*Si78(1,1)) + eff(1).x(1)*Si78(2,1))*T178(1,3) + (-(eff(1).x(2)*Si78(1,2)) + eff(1).x(1)*Si78(2,2))*T178(2,3)) + S87(1,1)*((-(eff(1).x(2)*Si78(1,2)) + eff(1).x(1)*Si78(2,2))*T178(2,1) + (-(eff(1).x(2)*Si78(1,3)) + eff(1).x(1)*Si78(2,3))*T178(3,1)) + S87(2,1)*((-(eff(1).x(2)*Si78(1,1)) + eff(1).x(1)*Si78(2,1))*T178(1,2) + (-(eff(1).x(2)*Si78(1,3)) + eff(1).x(1)*Si78(2,3))*T178(3,2)) + (-(eff(1).x(3)*S87(3,2)) + eff(1).x(2)*S87(3,3))*((-(eff(1).x(2)*Si78(1,3)) + eff(1).x(1)*Si78(2,3))*T178(3,6) + Si78(3,1)*T178(4,6) + Si78(3,2)*T178(5,6)) + (-(eff(1).x(3)*S87(1,2)) + eff(1).x(2)*S87(1,3))*((-(eff(1).x(2)*Si78(1,1)) + eff(1).x(1)*Si78(2,1))*T178(1,4) + Si78(3,2)*T178(5,4) + Si78(3,3)*T178(6,4)) + (-(eff(1).x(3)*S87(2,2)) + eff(1).x(2)*S87(2,3))*((-(eff(1).x(2)*Si78(1,2)) + eff(1).x(1)*Si78(2,2))*T178(2,5) + Si78(3,1)*T178(4,5) + Si78(3,3)*T178(6,5));
T78(6,2) = S87(3,2)*((-(eff(1).x(2)*Si78(1,1)) + eff(1).x(1)*Si78(2,1))*T178(1,3) + (-(eff(1).x(2)*Si78(1,2)) + eff(1).x(1)*Si78(2,2))*T178(2,3)) + S87(1,2)*((-(eff(1).x(2)*Si78(1,2)) + eff(1).x(1)*Si78(2,2))*T178(2,1) + (-(eff(1).x(2)*Si78(1,3)) + eff(1).x(1)*Si78(2,3))*T178(3,1)) + S87(2,2)*((-(eff(1).x(2)*Si78(1,1)) + eff(1).x(1)*Si78(2,1))*T178(1,2) + (-(eff(1).x(2)*Si78(1,3)) + eff(1).x(1)*Si78(2,3))*T178(3,2)) + (eff(1).x(3)*S87(3,1) - eff(1).x(1)*S87(3,3))*((-(eff(1).x(2)*Si78(1,3)) + eff(1).x(1)*Si78(2,3))*T178(3,6) + Si78(3,1)*T178(4,6) + Si78(3,2)*T178(5,6)) + (eff(1).x(3)*S87(1,1) - eff(1).x(1)*S87(1,3))*((-(eff(1).x(2)*Si78(1,1)) + eff(1).x(1)*Si78(2,1))*T178(1,4) + Si78(3,2)*T178(5,4) + Si78(3,3)*T178(6,4)) + (eff(1).x(3)*S87(2,1) - eff(1).x(1)*S87(2,3))*((-(eff(1).x(2)*Si78(1,2)) + eff(1).x(1)*Si78(2,2))*T178(2,5) + Si78(3,1)*T178(4,5) + Si78(3,3)*T178(6,5));
T78(6,3) = S87(3,3)*((-(eff(1).x(2)*Si78(1,1)) + eff(1).x(1)*Si78(2,1))*T178(1,3) + (-(eff(1).x(2)*Si78(1,2)) + eff(1).x(1)*Si78(2,2))*T178(2,3)) + S87(1,3)*((-(eff(1).x(2)*Si78(1,2)) + eff(1).x(1)*Si78(2,2))*T178(2,1) + (-(eff(1).x(2)*Si78(1,3)) + eff(1).x(1)*Si78(2,3))*T178(3,1)) + S87(2,3)*((-(eff(1).x(2)*Si78(1,1)) + eff(1).x(1)*Si78(2,1))*T178(1,2) + (-(eff(1).x(2)*Si78(1,3)) + eff(1).x(1)*Si78(2,3))*T178(3,2)) + (-(eff(1).x(2)*S87(3,1)) + eff(1).x(1)*S87(3,2))*((-(eff(1).x(2)*Si78(1,3)) + eff(1).x(1)*Si78(2,3))*T178(3,6) + Si78(3,1)*T178(4,6) + Si78(3,2)*T178(5,6)) + (-(eff(1).x(2)*S87(1,1)) + eff(1).x(1)*S87(1,2))*((-(eff(1).x(2)*Si78(1,1)) + eff(1).x(1)*Si78(2,1))*T178(1,4) + Si78(3,2)*T178(5,4) + Si78(3,3)*T178(6,4)) + (-(eff(1).x(2)*S87(2,1)) + eff(1).x(1)*S87(2,2))*((-(eff(1).x(2)*Si78(1,2)) + eff(1).x(1)*Si78(2,2))*T178(2,5) + Si78(3,1)*T178(4,5) + Si78(3,3)*T178(6,5));
T78(6,4) = S87(3,1)*((-(eff(1).x(2)*Si78(1,3)) + eff(1).x(1)*Si78(2,3))*T178(3,6) + Si78(3,1)*T178(4,6) + Si78(3,2)*T178(5,6)) + S87(1,1)*((-(eff(1).x(2)*Si78(1,1)) + eff(1).x(1)*Si78(2,1))*T178(1,4) + Si78(3,2)*T178(5,4) + Si78(3,3)*T178(6,4)) + S87(2,1)*((-(eff(1).x(2)*Si78(1,2)) + eff(1).x(1)*Si78(2,2))*T178(2,5) + Si78(3,1)*T178(4,5) + Si78(3,3)*T178(6,5));
T78(6,5) = S87(3,2)*((-(eff(1).x(2)*Si78(1,3)) + eff(1).x(1)*Si78(2,3))*T178(3,6) + Si78(3,1)*T178(4,6) + Si78(3,2)*T178(5,6)) + S87(1,2)*((-(eff(1).x(2)*Si78(1,1)) + eff(1).x(1)*Si78(2,1))*T178(1,4) + Si78(3,2)*T178(5,4) + Si78(3,3)*T178(6,4)) + S87(2,2)*((-(eff(1).x(2)*Si78(1,2)) + eff(1).x(1)*Si78(2,2))*T178(2,5) + Si78(3,1)*T178(4,5) + Si78(3,3)*T178(6,5));
T78(6,6) = S87(3,3)*((-(eff(1).x(2)*Si78(1,3)) + eff(1).x(1)*Si78(2,3))*T178(3,6) + Si78(3,1)*T178(4,6) + Si78(3,2)*T178(5,6)) + S87(1,3)*((-(eff(1).x(2)*Si78(1,1)) + eff(1).x(1)*Si78(2,1))*T178(1,4) + Si78(3,2)*T178(5,4) + Si78(3,3)*T178(6,4)) + S87(2,3)*((-(eff(1).x(2)*Si78(1,2)) + eff(1).x(1)*Si78(2,2))*T178(2,5) + Si78(3,1)*T178(4,5) + Si78(3,3)*T178(6,5));







%barrett_InvDynArtfunc8
     
JA7(1,1) = T78(1,1);
JA7(1,2) = links(7).mcm(3) + T78(1,2);
JA7(1,3) = -links(7).mcm(2) + T78(1,3);
JA7(1,4) = links(7).m + T78(1,4);
JA7(1,5) = T78(1,5);
JA7(1,6) = T78(1,6);

JA7(2,1) = -links(7).mcm(3) + T78(2,1);
JA7(2,2) = T78(2,2);
JA7(2,3) = links(7).mcm(1) + T78(2,3);
JA7(2,4) = T78(2,4);
JA7(2,5) = links(7).m + T78(2,5);
JA7(2,6) = T78(2,6);

JA7(3,1) = links(7).mcm(2) + T78(3,1);
JA7(3,2) = -links(7).mcm(1) + T78(3,2);
JA7(3,3) = T78(3,3);
JA7(3,4) = T78(3,4);
JA7(3,5) = T78(3,5);
JA7(3,6) = links(7).m + T78(3,6);

JA7(4,1) = links(7).inertia(1,1) + T78(4,1);
JA7(4,2) = links(7).inertia(1,2) + T78(4,2);
JA7(4,3) = links(7).inertia(1,3) + T78(4,3);
JA7(4,4) = T78(4,4);
JA7(4,5) = -links(7).mcm(3) + T78(4,5);
JA7(4,6) = links(7).mcm(2) + T78(4,6);

JA7(5,1) = links(7).inertia(1,2) + T78(5,1);
JA7(5,2) = links(7).inertia(2,2) + T78(5,2);
JA7(5,3) = links(7).inertia(2,3) + T78(5,3);
JA7(5,4) = links(7).mcm(3) + T78(5,4);
JA7(5,5) = T78(5,5);
JA7(5,6) = -links(7).mcm(1) + T78(5,6);

JA7(6,1) = links(7).inertia(1,3) + T78(6,1);
JA7(6,2) = links(7).inertia(2,3) + T78(6,2);
JA7(6,3) = links(7).inertia(3,3) + T78(6,3);
JA7(6,4) = -links(7).mcm(2) + T78(6,4);
JA7(6,5) = links(7).mcm(1) + T78(6,5);
JA7(6,6) = T78(6,6);


h7(1) = JA7(1,3);
h7(2) = JA7(2,3);
h7(3) = JA7(3,3);
h7(4) = JA7(4,3);
h7(5) = JA7(5,3);
h7(6) = JA7(6,3);

T167(1,1) = JA7(1,1);
T167(1,2) = JA7(1,2);
T167(1,3) = JA7(1,3);
T167(1,4) = JA7(1,4);
T167(1,5) = JA7(1,5);
T167(1,6) = JA7(1,6);

T167(2,1) = JA7(2,1);
T167(2,2) = JA7(2,2);
T167(2,3) = JA7(2,3);
T167(2,4) = JA7(2,4);
T167(2,5) = JA7(2,5);
T167(2,6) = JA7(2,6);

T167(3,1) = JA7(3,1);
T167(3,2) = JA7(3,2);
T167(3,3) = JA7(3,3);
T167(3,4) = JA7(3,4);
T167(3,5) = JA7(3,5);
T167(3,6) = JA7(3,6);

T167(4,1) = JA7(4,1);
T167(4,2) = JA7(4,2);
T167(4,3) = JA7(4,3);
T167(4,4) = JA7(4,4);
T167(4,5) = JA7(4,5);
T167(4,6) = JA7(4,6);

T167(5,1) = JA7(5,1);
T167(5,2) = JA7(5,2);
T167(5,3) = JA7(5,3);
T167(5,4) = JA7(5,4);
T167(5,5) = JA7(5,5);
T167(5,6) = JA7(5,6);

T167(6,1) = JA7(6,1);
T167(6,2) = JA7(6,2);
T167(6,3) = JA7(6,3);
T167(6,4) = JA7(6,4);
T167(6,5) = JA7(6,5);
T167(6,6) = JA7(6,6);


T67(1,1) = T167(3,3);
T67(1,2) = S76(1,2)*T167(3,1) + S76(2,2)*T167(3,2);
T67(1,3) = S76(1,3)*T167(3,1) + S76(2,3)*T167(3,2);
T67(1,4) = T167(3,6);
T67(1,5) = S76(1,2)*T167(3,4) + S76(2,2)*T167(3,5);
T67(1,6) = S76(1,3)*T167(3,4) + S76(2,3)*T167(3,5);

T67(2,1) = Si67(2,1)*T167(1,3) + Si67(2,2)*T167(2,3);
T67(2,2) = S76(1,2)*(Si67(2,1)*T167(1,1) + Si67(2,2)*T167(2,1)) + S76(2,2)*(Si67(2,1)*T167(1,2) + Si67(2,2)*T167(2,2));
T67(2,3) = S76(1,3)*(Si67(2,1)*T167(1,1) + Si67(2,2)*T167(2,1)) + S76(2,3)*(Si67(2,1)*T167(1,2) + Si67(2,2)*T167(2,2));
T67(2,4) = Si67(2,1)*T167(1,6) + Si67(2,2)*T167(2,6);
T67(2,5) = S76(1,2)*(Si67(2,1)*T167(1,4) + Si67(2,2)*T167(2,4)) + S76(2,2)*(Si67(2,1)*T167(1,5) + Si67(2,2)*T167(2,5));
T67(2,6) = S76(1,3)*(Si67(2,1)*T167(1,4) + Si67(2,2)*T167(2,4)) + S76(2,3)*(Si67(2,1)*T167(1,5) + Si67(2,2)*T167(2,5));

T67(3,1) = Si67(3,1)*T167(1,3) + Si67(3,2)*T167(2,3);
T67(3,2) = S76(1,2)*(Si67(3,1)*T167(1,1) + Si67(3,2)*T167(2,1)) + S76(2,2)*(Si67(3,1)*T167(1,2) + Si67(3,2)*T167(2,2));
T67(3,3) = S76(1,3)*(Si67(3,1)*T167(1,1) + Si67(3,2)*T167(2,1)) + S76(2,3)*(Si67(3,1)*T167(1,2) + Si67(3,2)*T167(2,2));
T67(3,4) = Si67(3,1)*T167(1,6) + Si67(3,2)*T167(2,6);
T67(3,5) = S76(1,2)*(Si67(3,1)*T167(1,4) + Si67(3,2)*T167(2,4)) + S76(2,2)*(Si67(3,1)*T167(1,5) + Si67(3,2)*T167(2,5));
T67(3,6) = S76(1,3)*(Si67(3,1)*T167(1,4) + Si67(3,2)*T167(2,4)) + S76(2,3)*(Si67(3,1)*T167(1,5) + Si67(3,2)*T167(2,5));

T67(4,1) = T167(6,3);
T67(4,2) = S76(1,2)*T167(6,1) + S76(2,2)*T167(6,2);
T67(4,3) = S76(1,3)*T167(6,1) + S76(2,3)*T167(6,2);
T67(4,4) = T167(6,6);
T67(4,5) = S76(1,2)*T167(6,4) + S76(2,2)*T167(6,5);
T67(4,6) = S76(1,3)*T167(6,4) + S76(2,3)*T167(6,5);

T67(5,1) = Si67(2,1)*T167(4,3) + Si67(2,2)*T167(5,3);
T67(5,2) = S76(1,2)*(Si67(2,1)*T167(4,1) + Si67(2,2)*T167(5,1)) + S76(2,2)*(Si67(2,1)*T167(4,2) + Si67(2,2)*T167(5,2));
T67(5,3) = S76(1,3)*(Si67(2,1)*T167(4,1) + Si67(2,2)*T167(5,1)) + S76(2,3)*(Si67(2,1)*T167(4,2) + Si67(2,2)*T167(5,2));
T67(5,4) = Si67(2,1)*T167(4,6) + Si67(2,2)*T167(5,6);
T67(5,5) = S76(1,2)*(Si67(2,1)*T167(4,4) + Si67(2,2)*T167(5,4)) + S76(2,2)*(Si67(2,1)*T167(4,5) + Si67(2,2)*T167(5,5));
T67(5,6) = S76(1,3)*(Si67(2,1)*T167(4,4) + Si67(2,2)*T167(5,4)) + S76(2,3)*(Si67(2,1)*T167(4,5) + Si67(2,2)*T167(5,5));

T67(6,1) = Si67(3,1)*T167(4,3) + Si67(3,2)*T167(5,3);
T67(6,2) = S76(1,2)*(Si67(3,1)*T167(4,1) + Si67(3,2)*T167(5,1)) + S76(2,2)*(Si67(3,1)*T167(4,2) + Si67(3,2)*T167(5,2));
T67(6,3) = S76(1,3)*(Si67(3,1)*T167(4,1) + Si67(3,2)*T167(5,1)) + S76(2,3)*(Si67(3,1)*T167(4,2) + Si67(3,2)*T167(5,2));
T67(6,4) = Si67(3,1)*T167(4,6) + Si67(3,2)*T167(5,6);
T67(6,5) = S76(1,2)*(Si67(3,1)*T167(4,4) + Si67(3,2)*T167(5,4)) + S76(2,2)*(Si67(3,1)*T167(4,5) + Si67(3,2)*T167(5,5));
T67(6,6) = S76(1,3)*(Si67(3,1)*T167(4,4) + Si67(3,2)*T167(5,4)) + S76(2,3)*(Si67(3,1)*T167(4,5) + Si67(3,2)*T167(5,5));







%barrett_InvDynArtfunc9
     
JA6(1,1) = T67(1,1);
JA6(1,2) = links(6).mcm(3) + T67(1,2);
JA6(1,3) = -links(6).mcm(2) + T67(1,3);
JA6(1,4) = links(6).m + T67(1,4);
JA6(1,5) = T67(1,5);
JA6(1,6) = T67(1,6);

JA6(2,1) = -links(6).mcm(3) + T67(2,1);
JA6(2,2) = T67(2,2);
JA6(2,3) = links(6).mcm(1) + T67(2,3);
JA6(2,4) = T67(2,4);
JA6(2,5) = links(6).m + T67(2,5);
JA6(2,6) = T67(2,6);

JA6(3,1) = links(6).mcm(2) + T67(3,1);
JA6(3,2) = -links(6).mcm(1) + T67(3,2);
JA6(3,3) = T67(3,3);
JA6(3,4) = T67(3,4);
JA6(3,5) = T67(3,5);
JA6(3,6) = links(6).m + T67(3,6);

JA6(4,1) = links(6).inertia(1,1) + T67(4,1);
JA6(4,2) = links(6).inertia(1,2) + T67(4,2);
JA6(4,3) = links(6).inertia(1,3) + T67(4,3);
JA6(4,4) = T67(4,4);
JA6(4,5) = -links(6).mcm(3) + T67(4,5);
JA6(4,6) = links(6).mcm(2) + T67(4,6);

JA6(5,1) = links(6).inertia(1,2) + T67(5,1);
JA6(5,2) = links(6).inertia(2,2) + T67(5,2);
JA6(5,3) = links(6).inertia(2,3) + T67(5,3);
JA6(5,4) = links(6).mcm(3) + T67(5,4);
JA6(5,5) = T67(5,5);
JA6(5,6) = -links(6).mcm(1) + T67(5,6);

JA6(6,1) = links(6).inertia(1,3) + T67(6,1);
JA6(6,2) = links(6).inertia(2,3) + T67(6,2);
JA6(6,3) = links(6).inertia(3,3) + T67(6,3);
JA6(6,4) = -links(6).mcm(2) + T67(6,4);
JA6(6,5) = links(6).mcm(1) + T67(6,5);
JA6(6,6) = T67(6,6);


h6(1) = JA6(1,3);
h6(2) = JA6(2,3);
h6(3) = JA6(3,3);
h6(4) = JA6(4,3);
h6(5) = JA6(5,3);
h6(6) = JA6(6,3);

T156(1,1) = JA6(1,1);
T156(1,2) = JA6(1,2);
T156(1,3) = JA6(1,3);
T156(1,4) = JA6(1,4);
T156(1,5) = JA6(1,5);
T156(1,6) = JA6(1,6);

T156(2,1) = JA6(2,1);
T156(2,2) = JA6(2,2);
T156(2,3) = JA6(2,3);
T156(2,4) = JA6(2,4);
T156(2,5) = JA6(2,5);
T156(2,6) = JA6(2,6);

T156(3,1) = JA6(3,1);
T156(3,2) = JA6(3,2);
T156(3,3) = JA6(3,3);
T156(3,4) = JA6(3,4);
T156(3,5) = JA6(3,5);
T156(3,6) = JA6(3,6);

T156(4,1) = JA6(4,1);
T156(4,2) = JA6(4,2);
T156(4,3) = JA6(4,3);
T156(4,4) = JA6(4,4);
T156(4,5) = JA6(4,5);
T156(4,6) = JA6(4,6);

T156(5,1) = JA6(5,1);
T156(5,2) = JA6(5,2);
T156(5,3) = JA6(5,3);
T156(5,4) = JA6(5,4);
T156(5,5) = JA6(5,5);
T156(5,6) = JA6(5,6);

T156(6,1) = JA6(6,1);
T156(6,2) = JA6(6,2);
T156(6,3) = JA6(6,3);
T156(6,4) = JA6(6,4);
T156(6,5) = JA6(6,5);
T156(6,6) = JA6(6,6);


T56(1,1) = T156(3,3) + ZWFE*S65(1,2)*T156(3,4) + ZWFE*S65(2,2)*T156(3,5);
T56(1,2) = -(S65(1,2)*T156(3,1)) - S65(2,2)*T156(3,2) + ZWFE*T156(3,6);
T56(1,3) = -(S65(1,3)*T156(3,1)) - S65(2,3)*T156(3,2);
T56(1,4) = T156(3,6);
T56(1,5) = -(S65(1,2)*T156(3,4)) - S65(2,2)*T156(3,5);
T56(1,6) = -(S65(1,3)*T156(3,4)) - S65(2,3)*T156(3,5);

T56(2,1) = -(Si56(2,1)*T156(1,3)) - Si56(2,2)*T156(2,3) - ZWFE*S65(1,2)*(Si56(2,1)*T156(1,4) + Si56(2,2)*T156(2,4)) - ZWFE*S65(2,2)*(Si56(2,1)*T156(1,5) + Si56(2,2)*T156(2,5));
T56(2,2) = S65(1,2)*(Si56(2,1)*T156(1,1) + Si56(2,2)*T156(2,1)) + S65(2,2)*(Si56(2,1)*T156(1,2) + Si56(2,2)*T156(2,2)) - ZWFE*(Si56(2,1)*T156(1,6) + Si56(2,2)*T156(2,6));
T56(2,3) = S65(1,3)*(Si56(2,1)*T156(1,1) + Si56(2,2)*T156(2,1)) + S65(2,3)*(Si56(2,1)*T156(1,2) + Si56(2,2)*T156(2,2));
T56(2,4) = -(Si56(2,1)*T156(1,6)) - Si56(2,2)*T156(2,6);
T56(2,5) = S65(1,2)*(Si56(2,1)*T156(1,4) + Si56(2,2)*T156(2,4)) + S65(2,2)*(Si56(2,1)*T156(1,5) + Si56(2,2)*T156(2,5));
T56(2,6) = S65(1,3)*(Si56(2,1)*T156(1,4) + Si56(2,2)*T156(2,4)) + S65(2,3)*(Si56(2,1)*T156(1,5) + Si56(2,2)*T156(2,5));

T56(3,1) = -(Si56(3,1)*T156(1,3)) - Si56(3,2)*T156(2,3) - ZWFE*S65(1,2)*(Si56(3,1)*T156(1,4) + Si56(3,2)*T156(2,4)) - ZWFE*S65(2,2)*(Si56(3,1)*T156(1,5) + Si56(3,2)*T156(2,5));
T56(3,2) = S65(1,2)*(Si56(3,1)*T156(1,1) + Si56(3,2)*T156(2,1)) + S65(2,2)*(Si56(3,1)*T156(1,2) + Si56(3,2)*T156(2,2)) - ZWFE*(Si56(3,1)*T156(1,6) + Si56(3,2)*T156(2,6));
T56(3,3) = S65(1,3)*(Si56(3,1)*T156(1,1) + Si56(3,2)*T156(2,1)) + S65(2,3)*(Si56(3,1)*T156(1,2) + Si56(3,2)*T156(2,2));
T56(3,4) = -(Si56(3,1)*T156(1,6)) - Si56(3,2)*T156(2,6);
T56(3,5) = S65(1,2)*(Si56(3,1)*T156(1,4) + Si56(3,2)*T156(2,4)) + S65(2,2)*(Si56(3,1)*T156(1,5) + Si56(3,2)*T156(2,5));
T56(3,6) = S65(1,3)*(Si56(3,1)*T156(1,4) + Si56(3,2)*T156(2,4)) + S65(2,3)*(Si56(3,1)*T156(1,5) + Si56(3,2)*T156(2,5));

T56(4,1) = ZWFE*Si56(2,1)*T156(1,3) + ZWFE*Si56(2,2)*T156(2,3) + T156(6,3) - ZWFE*S65(1,2)*(-(ZWFE*Si56(2,1)*T156(1,4)) - ZWFE*Si56(2,2)*T156(2,4) - T156(6,4)) - ZWFE*S65(2,2)*(-(ZWFE*Si56(2,1)*T156(1,5)) - ZWFE*Si56(2,2)*T156(2,5) - T156(6,5));
T56(4,2) = S65(1,2)*(-(ZWFE*Si56(2,1)*T156(1,1)) - ZWFE*Si56(2,2)*T156(2,1) - T156(6,1)) + S65(2,2)*(-(ZWFE*Si56(2,1)*T156(1,2)) - ZWFE*Si56(2,2)*T156(2,2) - T156(6,2)) - ZWFE*(-(ZWFE*Si56(2,1)*T156(1,6)) - ZWFE*Si56(2,2)*T156(2,6) - T156(6,6));
T56(4,3) = S65(1,3)*(-(ZWFE*Si56(2,1)*T156(1,1)) - ZWFE*Si56(2,2)*T156(2,1) - T156(6,1)) + S65(2,3)*(-(ZWFE*Si56(2,1)*T156(1,2)) - ZWFE*Si56(2,2)*T156(2,2) - T156(6,2));
T56(4,4) = ZWFE*Si56(2,1)*T156(1,6) + ZWFE*Si56(2,2)*T156(2,6) + T156(6,6);
T56(4,5) = S65(1,2)*(-(ZWFE*Si56(2,1)*T156(1,4)) - ZWFE*Si56(2,2)*T156(2,4) - T156(6,4)) + S65(2,2)*(-(ZWFE*Si56(2,1)*T156(1,5)) - ZWFE*Si56(2,2)*T156(2,5) - T156(6,5));
T56(4,6) = S65(1,3)*(-(ZWFE*Si56(2,1)*T156(1,4)) - ZWFE*Si56(2,2)*T156(2,4) - T156(6,4)) + S65(2,3)*(-(ZWFE*Si56(2,1)*T156(1,5)) - ZWFE*Si56(2,2)*T156(2,5) - T156(6,5));

T56(5,1) = ZWFE*T156(3,3) - Si56(2,1)*T156(4,3) - Si56(2,2)*T156(5,3) - ZWFE*S65(1,2)*(-(ZWFE*T156(3,4)) + Si56(2,1)*T156(4,4) + Si56(2,2)*T156(5,4)) - ZWFE*S65(2,2)*(-(ZWFE*T156(3,5)) + Si56(2,1)*T156(4,5) + Si56(2,2)*T156(5,5));
T56(5,2) = S65(1,2)*(-(ZWFE*T156(3,1)) + Si56(2,1)*T156(4,1) + Si56(2,2)*T156(5,1)) + S65(2,2)*(-(ZWFE*T156(3,2)) + Si56(2,1)*T156(4,2) + Si56(2,2)*T156(5,2)) - ZWFE*(-(ZWFE*T156(3,6)) + Si56(2,1)*T156(4,6) + Si56(2,2)*T156(5,6));
T56(5,3) = S65(1,3)*(-(ZWFE*T156(3,1)) + Si56(2,1)*T156(4,1) + Si56(2,2)*T156(5,1)) + S65(2,3)*(-(ZWFE*T156(3,2)) + Si56(2,1)*T156(4,2) + Si56(2,2)*T156(5,2));
T56(5,4) = ZWFE*T156(3,6) - Si56(2,1)*T156(4,6) - Si56(2,2)*T156(5,6);
T56(5,5) = S65(1,2)*(-(ZWFE*T156(3,4)) + Si56(2,1)*T156(4,4) + Si56(2,2)*T156(5,4)) + S65(2,2)*(-(ZWFE*T156(3,5)) + Si56(2,1)*T156(4,5) + Si56(2,2)*T156(5,5));
T56(5,6) = S65(1,3)*(-(ZWFE*T156(3,4)) + Si56(2,1)*T156(4,4) + Si56(2,2)*T156(5,4)) + S65(2,3)*(-(ZWFE*T156(3,5)) + Si56(2,1)*T156(4,5) + Si56(2,2)*T156(5,5));

T56(6,1) = -(Si56(3,1)*T156(4,3)) - Si56(3,2)*T156(5,3) - ZWFE*S65(1,2)*(Si56(3,1)*T156(4,4) + Si56(3,2)*T156(5,4)) - ZWFE*S65(2,2)*(Si56(3,1)*T156(4,5) + Si56(3,2)*T156(5,5));
T56(6,2) = S65(1,2)*(Si56(3,1)*T156(4,1) + Si56(3,2)*T156(5,1)) + S65(2,2)*(Si56(3,1)*T156(4,2) + Si56(3,2)*T156(5,2)) - ZWFE*(Si56(3,1)*T156(4,6) + Si56(3,2)*T156(5,6));
T56(6,3) = S65(1,3)*(Si56(3,1)*T156(4,1) + Si56(3,2)*T156(5,1)) + S65(2,3)*(Si56(3,1)*T156(4,2) + Si56(3,2)*T156(5,2));
T56(6,4) = -(Si56(3,1)*T156(4,6)) - Si56(3,2)*T156(5,6);
T56(6,5) = S65(1,2)*(Si56(3,1)*T156(4,4) + Si56(3,2)*T156(5,4)) + S65(2,2)*(Si56(3,1)*T156(4,5) + Si56(3,2)*T156(5,5));
T56(6,6) = S65(1,3)*(Si56(3,1)*T156(4,4) + Si56(3,2)*T156(5,4)) + S65(2,3)*(Si56(3,1)*T156(4,5) + Si56(3,2)*T156(5,5));







%barrett_InvDynArtfunc10
      
JA5(1,1) = T56(1,1);
JA5(1,2) = links(5).mcm(3) + T56(1,2);
JA5(1,3) = -links(5).mcm(2) + T56(1,3);
JA5(1,4) = links(5).m + T56(1,4);
JA5(1,5) = T56(1,5);
JA5(1,6) = T56(1,6);

JA5(2,1) = -links(5).mcm(3) + T56(2,1);
JA5(2,2) = T56(2,2);
JA5(2,3) = links(5).mcm(1) + T56(2,3);
JA5(2,4) = T56(2,4);
JA5(2,5) = links(5).m + T56(2,5);
JA5(2,6) = T56(2,6);

JA5(3,1) = links(5).mcm(2) + T56(3,1);
JA5(3,2) = -links(5).mcm(1) + T56(3,2);
JA5(3,3) = T56(3,3);
JA5(3,4) = T56(3,4);
JA5(3,5) = T56(3,5);
JA5(3,6) = links(5).m + T56(3,6);

JA5(4,1) = links(5).inertia(1,1) + T56(4,1);
JA5(4,2) = links(5).inertia(1,2) + T56(4,2);
JA5(4,3) = links(5).inertia(1,3) + T56(4,3);
JA5(4,4) = T56(4,4);
JA5(4,5) = -links(5).mcm(3) + T56(4,5);
JA5(4,6) = links(5).mcm(2) + T56(4,6);

JA5(5,1) = links(5).inertia(1,2) + T56(5,1);
JA5(5,2) = links(5).inertia(2,2) + T56(5,2);
JA5(5,3) = links(5).inertia(2,3) + T56(5,3);
JA5(5,4) = links(5).mcm(3) + T56(5,4);
JA5(5,5) = T56(5,5);
JA5(5,6) = -links(5).mcm(1) + T56(5,6);

JA5(6,1) = links(5).inertia(1,3) + T56(6,1);
JA5(6,2) = links(5).inertia(2,3) + T56(6,2);
JA5(6,3) = links(5).inertia(3,3) + T56(6,3);
JA5(6,4) = -links(5).mcm(2) + T56(6,4);
JA5(6,5) = links(5).mcm(1) + T56(6,5);
JA5(6,6) = T56(6,6);


h5(1) = JA5(1,3);
h5(2) = JA5(2,3);
h5(3) = JA5(3,3);
h5(4) = JA5(4,3);
h5(5) = JA5(5,3);
h5(6) = JA5(6,3);

T145(1,1) = JA5(1,1);
T145(1,2) = JA5(1,2);
T145(1,3) = JA5(1,3);
T145(1,4) = JA5(1,4);
T145(1,5) = JA5(1,5);
T145(1,6) = JA5(1,6);

T145(2,1) = JA5(2,1);
T145(2,2) = JA5(2,2);
T145(2,3) = JA5(2,3);
T145(2,4) = JA5(2,4);
T145(2,5) = JA5(2,5);
T145(2,6) = JA5(2,6);

T145(3,1) = JA5(3,1);
T145(3,2) = JA5(3,2);
T145(3,3) = JA5(3,3);
T145(3,4) = JA5(3,4);
T145(3,5) = JA5(3,5);
T145(3,6) = JA5(3,6);

T145(4,1) = JA5(4,1);
T145(4,2) = JA5(4,2);
T145(4,3) = JA5(4,3);
T145(4,4) = JA5(4,4);
T145(4,5) = JA5(4,5);
T145(4,6) = JA5(4,6);

T145(5,1) = JA5(5,1);
T145(5,2) = JA5(5,2);
T145(5,3) = JA5(5,3);
T145(5,4) = JA5(5,4);
T145(5,5) = JA5(5,5);
T145(5,6) = JA5(5,6);

T145(6,1) = JA5(6,1);
T145(6,2) = JA5(6,2);
T145(6,3) = JA5(6,3);
T145(6,4) = JA5(6,4);
T145(6,5) = JA5(6,5);
T145(6,6) = JA5(6,6);


T45(1,1) = T145(3,3) + YWR*S54(1,3)*T145(3,4) + YWR*S54(2,3)*T145(3,5);
T45(1,2) = S54(1,2)*T145(3,1) + S54(2,2)*T145(3,2) - ZWR*S54(1,3)*T145(3,4) - ZWR*S54(2,3)*T145(3,5);
T45(1,3) = S54(1,3)*T145(3,1) + S54(2,3)*T145(3,2) + ZWR*S54(1,2)*T145(3,4) + ZWR*S54(2,2)*T145(3,5) - YWR*T145(3,6);
T45(1,4) = T145(3,6);
T45(1,5) = S54(1,2)*T145(3,4) + S54(2,2)*T145(3,5);
T45(1,6) = S54(1,3)*T145(3,4) + S54(2,3)*T145(3,5);

T45(2,1) = Si45(2,1)*T145(1,3) + Si45(2,2)*T145(2,3) + YWR*S54(1,3)*(Si45(2,1)*T145(1,4) + Si45(2,2)*T145(2,4)) + YWR*S54(2,3)*(Si45(2,1)*T145(1,5) + Si45(2,2)*T145(2,5));
T45(2,2) = S54(1,2)*(Si45(2,1)*T145(1,1) + Si45(2,2)*T145(2,1)) + S54(2,2)*(Si45(2,1)*T145(1,2) + Si45(2,2)*T145(2,2)) - ZWR*S54(1,3)*(Si45(2,1)*T145(1,4) + Si45(2,2)*T145(2,4)) - ZWR*S54(2,3)*(Si45(2,1)*T145(1,5) + Si45(2,2)*T145(2,5));
T45(2,3) = S54(1,3)*(Si45(2,1)*T145(1,1) + Si45(2,2)*T145(2,1)) + S54(2,3)*(Si45(2,1)*T145(1,2) + Si45(2,2)*T145(2,2)) + ZWR*S54(1,2)*(Si45(2,1)*T145(1,4) + Si45(2,2)*T145(2,4)) + ZWR*S54(2,2)*(Si45(2,1)*T145(1,5) + Si45(2,2)*T145(2,5)) - YWR*(Si45(2,1)*T145(1,6) + Si45(2,2)*T145(2,6));
T45(2,4) = Si45(2,1)*T145(1,6) + Si45(2,2)*T145(2,6);
T45(2,5) = S54(1,2)*(Si45(2,1)*T145(1,4) + Si45(2,2)*T145(2,4)) + S54(2,2)*(Si45(2,1)*T145(1,5) + Si45(2,2)*T145(2,5));
T45(2,6) = S54(1,3)*(Si45(2,1)*T145(1,4) + Si45(2,2)*T145(2,4)) + S54(2,3)*(Si45(2,1)*T145(1,5) + Si45(2,2)*T145(2,5));

T45(3,1) = Si45(3,1)*T145(1,3) + Si45(3,2)*T145(2,3) + YWR*S54(1,3)*(Si45(3,1)*T145(1,4) + Si45(3,2)*T145(2,4)) + YWR*S54(2,3)*(Si45(3,1)*T145(1,5) + Si45(3,2)*T145(2,5));
T45(3,2) = S54(1,2)*(Si45(3,1)*T145(1,1) + Si45(3,2)*T145(2,1)) + S54(2,2)*(Si45(3,1)*T145(1,2) + Si45(3,2)*T145(2,2)) - ZWR*S54(1,3)*(Si45(3,1)*T145(1,4) + Si45(3,2)*T145(2,4)) - ZWR*S54(2,3)*(Si45(3,1)*T145(1,5) + Si45(3,2)*T145(2,5));
T45(3,3) = S54(1,3)*(Si45(3,1)*T145(1,1) + Si45(3,2)*T145(2,1)) + S54(2,3)*(Si45(3,1)*T145(1,2) + Si45(3,2)*T145(2,2)) + ZWR*S54(1,2)*(Si45(3,1)*T145(1,4) + Si45(3,2)*T145(2,4)) + ZWR*S54(2,2)*(Si45(3,1)*T145(1,5) + Si45(3,2)*T145(2,5)) - YWR*(Si45(3,1)*T145(1,6) + Si45(3,2)*T145(2,6));
T45(3,4) = Si45(3,1)*T145(1,6) + Si45(3,2)*T145(2,6);
T45(3,5) = S54(1,2)*(Si45(3,1)*T145(1,4) + Si45(3,2)*T145(2,4)) + S54(2,2)*(Si45(3,1)*T145(1,5) + Si45(3,2)*T145(2,5));
T45(3,6) = S54(1,3)*(Si45(3,1)*T145(1,4) + Si45(3,2)*T145(2,4)) + S54(2,3)*(Si45(3,1)*T145(1,5) + Si45(3,2)*T145(2,5));

T45(4,1) = YWR*Si45(3,1)*T145(1,3) + YWR*Si45(3,2)*T145(2,3) + T145(6,3) + YWR*S54(1,3)*(YWR*Si45(3,1)*T145(1,4) + YWR*Si45(3,2)*T145(2,4) + T145(6,4)) + YWR*S54(2,3)*(YWR*Si45(3,1)*T145(1,5) + YWR*Si45(3,2)*T145(2,5) + T145(6,5));
T45(4,2) = S54(1,2)*(YWR*Si45(3,1)*T145(1,1) + YWR*Si45(3,2)*T145(2,1) + T145(6,1)) + S54(2,2)*(YWR*Si45(3,1)*T145(1,2) + YWR*Si45(3,2)*T145(2,2) + T145(6,2)) - ZWR*S54(1,3)*(YWR*Si45(3,1)*T145(1,4) + YWR*Si45(3,2)*T145(2,4) + T145(6,4)) - ZWR*S54(2,3)*(YWR*Si45(3,1)*T145(1,5) + YWR*Si45(3,2)*T145(2,5) + T145(6,5));
T45(4,3) = S54(1,3)*(YWR*Si45(3,1)*T145(1,1) + YWR*Si45(3,2)*T145(2,1) + T145(6,1)) + S54(2,3)*(YWR*Si45(3,1)*T145(1,2) + YWR*Si45(3,2)*T145(2,2) + T145(6,2)) + ZWR*S54(1,2)*(YWR*Si45(3,1)*T145(1,4) + YWR*Si45(3,2)*T145(2,4) + T145(6,4)) + ZWR*S54(2,2)*(YWR*Si45(3,1)*T145(1,5) + YWR*Si45(3,2)*T145(2,5) + T145(6,5)) - YWR*(YWR*Si45(3,1)*T145(1,6) + YWR*Si45(3,2)*T145(2,6) + T145(6,6));
T45(4,4) = YWR*Si45(3,1)*T145(1,6) + YWR*Si45(3,2)*T145(2,6) + T145(6,6);
T45(4,5) = S54(1,2)*(YWR*Si45(3,1)*T145(1,4) + YWR*Si45(3,2)*T145(2,4) + T145(6,4)) + S54(2,2)*(YWR*Si45(3,1)*T145(1,5) + YWR*Si45(3,2)*T145(2,5) + T145(6,5));
T45(4,6) = S54(1,3)*(YWR*Si45(3,1)*T145(1,4) + YWR*Si45(3,2)*T145(2,4) + T145(6,4)) + S54(2,3)*(YWR*Si45(3,1)*T145(1,5) + YWR*Si45(3,2)*T145(2,5) + T145(6,5));

T45(5,1) = -(ZWR*Si45(3,1)*T145(1,3)) - ZWR*Si45(3,2)*T145(2,3) + Si45(2,1)*T145(4,3) + Si45(2,2)*T145(5,3) + YWR*S54(1,3)*(-(ZWR*Si45(3,1)*T145(1,4)) - ZWR*Si45(3,2)*T145(2,4) + Si45(2,1)*T145(4,4) + Si45(2,2)*T145(5,4)) + YWR*S54(2,3)*(-(ZWR*Si45(3,1)*T145(1,5)) - ZWR*Si45(3,2)*T145(2,5) + Si45(2,1)*T145(4,5) + Si45(2,2)*T145(5,5));
T45(5,2) = S54(1,2)*(-(ZWR*Si45(3,1)*T145(1,1)) - ZWR*Si45(3,2)*T145(2,1) + Si45(2,1)*T145(4,1) + Si45(2,2)*T145(5,1)) + S54(2,2)*(-(ZWR*Si45(3,1)*T145(1,2)) - ZWR*Si45(3,2)*T145(2,2) + Si45(2,1)*T145(4,2) + Si45(2,2)*T145(5,2)) - ZWR*S54(1,3)*(-(ZWR*Si45(3,1)*T145(1,4)) - ZWR*Si45(3,2)*T145(2,4) + Si45(2,1)*T145(4,4) + Si45(2,2)*T145(5,4)) - ZWR*S54(2,3)*(-(ZWR*Si45(3,1)*T145(1,5)) - ZWR*Si45(3,2)*T145(2,5) + Si45(2,1)*T145(4,5) + Si45(2,2)*T145(5,5));
T45(5,3) = S54(1,3)*(-(ZWR*Si45(3,1)*T145(1,1)) - ZWR*Si45(3,2)*T145(2,1) + Si45(2,1)*T145(4,1) + Si45(2,2)*T145(5,1)) + S54(2,3)*(-(ZWR*Si45(3,1)*T145(1,2)) - ZWR*Si45(3,2)*T145(2,2) + Si45(2,1)*T145(4,2) + Si45(2,2)*T145(5,2)) + ZWR*S54(1,2)*(-(ZWR*Si45(3,1)*T145(1,4)) - ZWR*Si45(3,2)*T145(2,4) + Si45(2,1)*T145(4,4) + Si45(2,2)*T145(5,4)) + ZWR*S54(2,2)*(-(ZWR*Si45(3,1)*T145(1,5)) - ZWR*Si45(3,2)*T145(2,5) + Si45(2,1)*T145(4,5) + Si45(2,2)*T145(5,5)) - YWR*(-(ZWR*Si45(3,1)*T145(1,6)) - ZWR*Si45(3,2)*T145(2,6) + Si45(2,1)*T145(4,6) + Si45(2,2)*T145(5,6));
T45(5,4) = -(ZWR*Si45(3,1)*T145(1,6)) - ZWR*Si45(3,2)*T145(2,6) + Si45(2,1)*T145(4,6) + Si45(2,2)*T145(5,6);
T45(5,5) = S54(1,2)*(-(ZWR*Si45(3,1)*T145(1,4)) - ZWR*Si45(3,2)*T145(2,4) + Si45(2,1)*T145(4,4) + Si45(2,2)*T145(5,4)) + S54(2,2)*(-(ZWR*Si45(3,1)*T145(1,5)) - ZWR*Si45(3,2)*T145(2,5) + Si45(2,1)*T145(4,5) + Si45(2,2)*T145(5,5));
T45(5,6) = S54(1,3)*(-(ZWR*Si45(3,1)*T145(1,4)) - ZWR*Si45(3,2)*T145(2,4) + Si45(2,1)*T145(4,4) + Si45(2,2)*T145(5,4)) + S54(2,3)*(-(ZWR*Si45(3,1)*T145(1,5)) - ZWR*Si45(3,2)*T145(2,5) + Si45(2,1)*T145(4,5) + Si45(2,2)*T145(5,5));

T45(6,1) = ZWR*Si45(2,1)*T145(1,3) + ZWR*Si45(2,2)*T145(2,3) - YWR*T145(3,3) + Si45(3,1)*T145(4,3) + Si45(3,2)*T145(5,3) + YWR*S54(1,3)*(ZWR*Si45(2,1)*T145(1,4) + ZWR*Si45(2,2)*T145(2,4) - YWR*T145(3,4) + Si45(3,1)*T145(4,4) + Si45(3,2)*T145(5,4)) + YWR*S54(2,3)*(ZWR*Si45(2,1)*T145(1,5) + ZWR*Si45(2,2)*T145(2,5) - YWR*T145(3,5) + Si45(3,1)*T145(4,5) + Si45(3,2)*T145(5,5));
T45(6,2) = S54(1,2)*(ZWR*Si45(2,1)*T145(1,1) + ZWR*Si45(2,2)*T145(2,1) - YWR*T145(3,1) + Si45(3,1)*T145(4,1) + Si45(3,2)*T145(5,1)) + S54(2,2)*(ZWR*Si45(2,1)*T145(1,2) + ZWR*Si45(2,2)*T145(2,2) - YWR*T145(3,2) + Si45(3,1)*T145(4,2) + Si45(3,2)*T145(5,2)) - ZWR*S54(1,3)*(ZWR*Si45(2,1)*T145(1,4) + ZWR*Si45(2,2)*T145(2,4) - YWR*T145(3,4) + Si45(3,1)*T145(4,4) + Si45(3,2)*T145(5,4)) - ZWR*S54(2,3)*(ZWR*Si45(2,1)*T145(1,5) + ZWR*Si45(2,2)*T145(2,5) - YWR*T145(3,5) + Si45(3,1)*T145(4,5) + Si45(3,2)*T145(5,5));
T45(6,3) = S54(1,3)*(ZWR*Si45(2,1)*T145(1,1) + ZWR*Si45(2,2)*T145(2,1) - YWR*T145(3,1) + Si45(3,1)*T145(4,1) + Si45(3,2)*T145(5,1)) + S54(2,3)*(ZWR*Si45(2,1)*T145(1,2) + ZWR*Si45(2,2)*T145(2,2) - YWR*T145(3,2) + Si45(3,1)*T145(4,2) + Si45(3,2)*T145(5,2)) + ZWR*S54(1,2)*(ZWR*Si45(2,1)*T145(1,4) + ZWR*Si45(2,2)*T145(2,4) - YWR*T145(3,4) + Si45(3,1)*T145(4,4) + Si45(3,2)*T145(5,4)) + ZWR*S54(2,2)*(ZWR*Si45(2,1)*T145(1,5) + ZWR*Si45(2,2)*T145(2,5) - YWR*T145(3,5) + Si45(3,1)*T145(4,5) + Si45(3,2)*T145(5,5)) - YWR*(ZWR*Si45(2,1)*T145(1,6) + ZWR*Si45(2,2)*T145(2,6) - YWR*T145(3,6) + Si45(3,1)*T145(4,6) + Si45(3,2)*T145(5,6));
T45(6,4) = ZWR*Si45(2,1)*T145(1,6) + ZWR*Si45(2,2)*T145(2,6) - YWR*T145(3,6) + Si45(3,1)*T145(4,6) + Si45(3,2)*T145(5,6);
T45(6,5) = S54(1,2)*(ZWR*Si45(2,1)*T145(1,4) + ZWR*Si45(2,2)*T145(2,4) - YWR*T145(3,4) + Si45(3,1)*T145(4,4) + Si45(3,2)*T145(5,4)) + S54(2,2)*(ZWR*Si45(2,1)*T145(1,5) + ZWR*Si45(2,2)*T145(2,5) - YWR*T145(3,5) + Si45(3,1)*T145(4,5) + Si45(3,2)*T145(5,5));
T45(6,6) = S54(1,3)*(ZWR*Si45(2,1)*T145(1,4) + ZWR*Si45(2,2)*T145(2,4) - YWR*T145(3,4) + Si45(3,1)*T145(4,4) + Si45(3,2)*T145(5,4)) + S54(2,3)*(ZWR*Si45(2,1)*T145(1,5) + ZWR*Si45(2,2)*T145(2,5) - YWR*T145(3,5) + Si45(3,1)*T145(4,5) + Si45(3,2)*T145(5,5));







%barrett_InvDynArtfunc11
      
JA4(1,1) = T45(1,1);
JA4(1,2) = links(4).mcm(3) + T45(1,2);
JA4(1,3) = -links(4).mcm(2) + T45(1,3);
JA4(1,4) = links(4).m + T45(1,4);
JA4(1,5) = T45(1,5);
JA4(1,6) = T45(1,6);

JA4(2,1) = -links(4).mcm(3) + T45(2,1);
JA4(2,2) = T45(2,2);
JA4(2,3) = links(4).mcm(1) + T45(2,3);
JA4(2,4) = T45(2,4);
JA4(2,5) = links(4).m + T45(2,5);
JA4(2,6) = T45(2,6);

JA4(3,1) = links(4).mcm(2) + T45(3,1);
JA4(3,2) = -links(4).mcm(1) + T45(3,2);
JA4(3,3) = T45(3,3);
JA4(3,4) = T45(3,4);
JA4(3,5) = T45(3,5);
JA4(3,6) = links(4).m + T45(3,6);

JA4(4,1) = links(4).inertia(1,1) + T45(4,1);
JA4(4,2) = links(4).inertia(1,2) + T45(4,2);
JA4(4,3) = links(4).inertia(1,3) + T45(4,3);
JA4(4,4) = T45(4,4);
JA4(4,5) = -links(4).mcm(3) + T45(4,5);
JA4(4,6) = links(4).mcm(2) + T45(4,6);

JA4(5,1) = links(4).inertia(1,2) + T45(5,1);
JA4(5,2) = links(4).inertia(2,2) + T45(5,2);
JA4(5,3) = links(4).inertia(2,3) + T45(5,3);
JA4(5,4) = links(4).mcm(3) + T45(5,4);
JA4(5,5) = T45(5,5);
JA4(5,6) = -links(4).mcm(1) + T45(5,6);

JA4(6,1) = links(4).inertia(1,3) + T45(6,1);
JA4(6,2) = links(4).inertia(2,3) + T45(6,2);
JA4(6,3) = links(4).inertia(3,3) + T45(6,3);
JA4(6,4) = -links(4).mcm(2) + T45(6,4);
JA4(6,5) = links(4).mcm(1) + T45(6,5);
JA4(6,6) = T45(6,6);


h4(1) = JA4(1,3);
h4(2) = JA4(2,3);
h4(3) = JA4(3,3);
h4(4) = JA4(4,3);
h4(5) = JA4(5,3);
h4(6) = JA4(6,3);

T134(1,1) = JA4(1,1);
T134(1,2) = JA4(1,2);
T134(1,3) = JA4(1,3);
T134(1,4) = JA4(1,4);
T134(1,5) = JA4(1,5);
T134(1,6) = JA4(1,6);

T134(2,1) = JA4(2,1);
T134(2,2) = JA4(2,2);
T134(2,3) = JA4(2,3);
T134(2,4) = JA4(2,4);
T134(2,5) = JA4(2,5);
T134(2,6) = JA4(2,6);

T134(3,1) = JA4(3,1);
T134(3,2) = JA4(3,2);
T134(3,3) = JA4(3,3);
T134(3,4) = JA4(3,4);
T134(3,5) = JA4(3,5);
T134(3,6) = JA4(3,6);

T134(4,1) = JA4(4,1);
T134(4,2) = JA4(4,2);
T134(4,3) = JA4(4,3);
T134(4,4) = JA4(4,4);
T134(4,5) = JA4(4,5);
T134(4,6) = JA4(4,6);

T134(5,1) = JA4(5,1);
T134(5,2) = JA4(5,2);
T134(5,3) = JA4(5,3);
T134(5,4) = JA4(5,4);
T134(5,5) = JA4(5,5);
T134(5,6) = JA4(5,6);

T134(6,1) = JA4(6,1);
T134(6,2) = JA4(6,2);
T134(6,3) = JA4(6,3);
T134(6,4) = JA4(6,4);
T134(6,5) = JA4(6,5);
T134(6,6) = JA4(6,6);


T34(1,1) = T134(3,3) - (-(ZEB*S43(1,2)) + YEB*S43(1,3))*T134(3,4) - (-(ZEB*S43(2,2)) + YEB*S43(2,3))*T134(3,5);
T34(1,2) = -(S43(1,2)*T134(3,1)) - S43(2,2)*T134(3,2) + ZEB*T134(3,6);
T34(1,3) = -(S43(1,3)*T134(3,1)) - S43(2,3)*T134(3,2) - YEB*T134(3,6);
T34(1,4) = T134(3,6);
T34(1,5) = -(S43(1,2)*T134(3,4)) - S43(2,2)*T134(3,5);
T34(1,6) = -(S43(1,3)*T134(3,4)) - S43(2,3)*T134(3,5);

T34(2,1) = -(Si34(2,1)*T134(1,3)) - Si34(2,2)*T134(2,3) + (-(ZEB*S43(1,2)) + YEB*S43(1,3))*(Si34(2,1)*T134(1,4) + Si34(2,2)*T134(2,4)) + (-(ZEB*S43(2,2)) + YEB*S43(2,3))*(Si34(2,1)*T134(1,5) + Si34(2,2)*T134(2,5));
T34(2,2) = S43(1,2)*(Si34(2,1)*T134(1,1) + Si34(2,2)*T134(2,1)) + S43(2,2)*(Si34(2,1)*T134(1,2) + Si34(2,2)*T134(2,2)) - ZEB*(Si34(2,1)*T134(1,6) + Si34(2,2)*T134(2,6));
T34(2,3) = S43(1,3)*(Si34(2,1)*T134(1,1) + Si34(2,2)*T134(2,1)) + S43(2,3)*(Si34(2,1)*T134(1,2) + Si34(2,2)*T134(2,2)) + YEB*(Si34(2,1)*T134(1,6) + Si34(2,2)*T134(2,6));
T34(2,4) = -(Si34(2,1)*T134(1,6)) - Si34(2,2)*T134(2,6);
T34(2,5) = S43(1,2)*(Si34(2,1)*T134(1,4) + Si34(2,2)*T134(2,4)) + S43(2,2)*(Si34(2,1)*T134(1,5) + Si34(2,2)*T134(2,5));
T34(2,6) = S43(1,3)*(Si34(2,1)*T134(1,4) + Si34(2,2)*T134(2,4)) + S43(2,3)*(Si34(2,1)*T134(1,5) + Si34(2,2)*T134(2,5));

T34(3,1) = -(Si34(3,1)*T134(1,3)) - Si34(3,2)*T134(2,3) + (-(ZEB*S43(1,2)) + YEB*S43(1,3))*(Si34(3,1)*T134(1,4) + Si34(3,2)*T134(2,4)) + (-(ZEB*S43(2,2)) + YEB*S43(2,3))*(Si34(3,1)*T134(1,5) + Si34(3,2)*T134(2,5));
T34(3,2) = S43(1,2)*(Si34(3,1)*T134(1,1) + Si34(3,2)*T134(2,1)) + S43(2,2)*(Si34(3,1)*T134(1,2) + Si34(3,2)*T134(2,2)) - ZEB*(Si34(3,1)*T134(1,6) + Si34(3,2)*T134(2,6));
T34(3,3) = S43(1,3)*(Si34(3,1)*T134(1,1) + Si34(3,2)*T134(2,1)) + S43(2,3)*(Si34(3,1)*T134(1,2) + Si34(3,2)*T134(2,2)) + YEB*(Si34(3,1)*T134(1,6) + Si34(3,2)*T134(2,6));
T34(3,4) = -(Si34(3,1)*T134(1,6)) - Si34(3,2)*T134(2,6);
T34(3,5) = S43(1,2)*(Si34(3,1)*T134(1,4) + Si34(3,2)*T134(2,4)) + S43(2,2)*(Si34(3,1)*T134(1,5) + Si34(3,2)*T134(2,5));
T34(3,6) = S43(1,3)*(Si34(3,1)*T134(1,4) + Si34(3,2)*T134(2,4)) + S43(2,3)*(Si34(3,1)*T134(1,5) + Si34(3,2)*T134(2,5));

T34(4,1) = -((-(ZEB*Si34(2,1)) + YEB*Si34(3,1))*T134(1,3)) - (-(ZEB*Si34(2,2)) + YEB*Si34(3,2))*T134(2,3) + T134(6,3) + (-(ZEB*S43(1,2)) + YEB*S43(1,3))*((-(ZEB*Si34(2,1)) + YEB*Si34(3,1))*T134(1,4) + (-(ZEB*Si34(2,2)) + YEB*Si34(3,2))*T134(2,4) - T134(6,4)) + (-(ZEB*S43(2,2)) + YEB*S43(2,3))*((-(ZEB*Si34(2,1)) + YEB*Si34(3,1))*T134(1,5) + (-(ZEB*Si34(2,2)) + YEB*Si34(3,2))*T134(2,5) - T134(6,5));
T34(4,2) = S43(1,2)*((-(ZEB*Si34(2,1)) + YEB*Si34(3,1))*T134(1,1) + (-(ZEB*Si34(2,2)) + YEB*Si34(3,2))*T134(2,1) - T134(6,1)) + S43(2,2)*((-(ZEB*Si34(2,1)) + YEB*Si34(3,1))*T134(1,2) + (-(ZEB*Si34(2,2)) + YEB*Si34(3,2))*T134(2,2) - T134(6,2)) - ZEB*((-(ZEB*Si34(2,1)) + YEB*Si34(3,1))*T134(1,6) + (-(ZEB*Si34(2,2)) + YEB*Si34(3,2))*T134(2,6) - T134(6,6));
T34(4,3) = S43(1,3)*((-(ZEB*Si34(2,1)) + YEB*Si34(3,1))*T134(1,1) + (-(ZEB*Si34(2,2)) + YEB*Si34(3,2))*T134(2,1) - T134(6,1)) + S43(2,3)*((-(ZEB*Si34(2,1)) + YEB*Si34(3,1))*T134(1,2) + (-(ZEB*Si34(2,2)) + YEB*Si34(3,2))*T134(2,2) - T134(6,2)) + YEB*((-(ZEB*Si34(2,1)) + YEB*Si34(3,1))*T134(1,6) + (-(ZEB*Si34(2,2)) + YEB*Si34(3,2))*T134(2,6) - T134(6,6));
T34(4,4) = -((-(ZEB*Si34(2,1)) + YEB*Si34(3,1))*T134(1,6)) - (-(ZEB*Si34(2,2)) + YEB*Si34(3,2))*T134(2,6) + T134(6,6);
T34(4,5) = S43(1,2)*((-(ZEB*Si34(2,1)) + YEB*Si34(3,1))*T134(1,4) + (-(ZEB*Si34(2,2)) + YEB*Si34(3,2))*T134(2,4) - T134(6,4)) + S43(2,2)*((-(ZEB*Si34(2,1)) + YEB*Si34(3,1))*T134(1,5) + (-(ZEB*Si34(2,2)) + YEB*Si34(3,2))*T134(2,5) - T134(6,5));
T34(4,6) = S43(1,3)*((-(ZEB*Si34(2,1)) + YEB*Si34(3,1))*T134(1,4) + (-(ZEB*Si34(2,2)) + YEB*Si34(3,2))*T134(2,4) - T134(6,4)) + S43(2,3)*((-(ZEB*Si34(2,1)) + YEB*Si34(3,1))*T134(1,5) + (-(ZEB*Si34(2,2)) + YEB*Si34(3,2))*T134(2,5) - T134(6,5));

T34(5,1) = ZEB*T134(3,3) - Si34(2,1)*T134(4,3) - Si34(2,2)*T134(5,3) + (-(ZEB*S43(1,2)) + YEB*S43(1,3))*(-(ZEB*T134(3,4)) + Si34(2,1)*T134(4,4) + Si34(2,2)*T134(5,4)) + (-(ZEB*S43(2,2)) + YEB*S43(2,3))*(-(ZEB*T134(3,5)) + Si34(2,1)*T134(4,5) + Si34(2,2)*T134(5,5));
T34(5,2) = S43(1,2)*(-(ZEB*T134(3,1)) + Si34(2,1)*T134(4,1) + Si34(2,2)*T134(5,1)) + S43(2,2)*(-(ZEB*T134(3,2)) + Si34(2,1)*T134(4,2) + Si34(2,2)*T134(5,2)) - ZEB*(-(ZEB*T134(3,6)) + Si34(2,1)*T134(4,6) + Si34(2,2)*T134(5,6));
T34(5,3) = S43(1,3)*(-(ZEB*T134(3,1)) + Si34(2,1)*T134(4,1) + Si34(2,2)*T134(5,1)) + S43(2,3)*(-(ZEB*T134(3,2)) + Si34(2,1)*T134(4,2) + Si34(2,2)*T134(5,2)) + YEB*(-(ZEB*T134(3,6)) + Si34(2,1)*T134(4,6) + Si34(2,2)*T134(5,6));
T34(5,4) = ZEB*T134(3,6) - Si34(2,1)*T134(4,6) - Si34(2,2)*T134(5,6);
T34(5,5) = S43(1,2)*(-(ZEB*T134(3,4)) + Si34(2,1)*T134(4,4) + Si34(2,2)*T134(5,4)) + S43(2,2)*(-(ZEB*T134(3,5)) + Si34(2,1)*T134(4,5) + Si34(2,2)*T134(5,5));
T34(5,6) = S43(1,3)*(-(ZEB*T134(3,4)) + Si34(2,1)*T134(4,4) + Si34(2,2)*T134(5,4)) + S43(2,3)*(-(ZEB*T134(3,5)) + Si34(2,1)*T134(4,5) + Si34(2,2)*T134(5,5));

T34(6,1) = -(YEB*T134(3,3)) - Si34(3,1)*T134(4,3) - Si34(3,2)*T134(5,3) + (-(ZEB*S43(1,2)) + YEB*S43(1,3))*(YEB*T134(3,4) + Si34(3,1)*T134(4,4) + Si34(3,2)*T134(5,4)) + (-(ZEB*S43(2,2)) + YEB*S43(2,3))*(YEB*T134(3,5) + Si34(3,1)*T134(4,5) + Si34(3,2)*T134(5,5));
T34(6,2) = S43(1,2)*(YEB*T134(3,1) + Si34(3,1)*T134(4,1) + Si34(3,2)*T134(5,1)) + S43(2,2)*(YEB*T134(3,2) + Si34(3,1)*T134(4,2) + Si34(3,2)*T134(5,2)) - ZEB*(YEB*T134(3,6) + Si34(3,1)*T134(4,6) + Si34(3,2)*T134(5,6));
T34(6,3) = S43(1,3)*(YEB*T134(3,1) + Si34(3,1)*T134(4,1) + Si34(3,2)*T134(5,1)) + S43(2,3)*(YEB*T134(3,2) + Si34(3,1)*T134(4,2) + Si34(3,2)*T134(5,2)) + YEB*(YEB*T134(3,6) + Si34(3,1)*T134(4,6) + Si34(3,2)*T134(5,6));
T34(6,4) = -(YEB*T134(3,6)) - Si34(3,1)*T134(4,6) - Si34(3,2)*T134(5,6);
T34(6,5) = S43(1,2)*(YEB*T134(3,4) + Si34(3,1)*T134(4,4) + Si34(3,2)*T134(5,4)) + S43(2,2)*(YEB*T134(3,5) + Si34(3,1)*T134(4,5) + Si34(3,2)*T134(5,5));
T34(6,6) = S43(1,3)*(YEB*T134(3,4) + Si34(3,1)*T134(4,4) + Si34(3,2)*T134(5,4)) + S43(2,3)*(YEB*T134(3,5) + Si34(3,1)*T134(4,5) + Si34(3,2)*T134(5,5));







%barrett_InvDynArtfunc12
      
JA3(1,1) = T34(1,1);
JA3(1,2) = links(3).mcm(3) + T34(1,2);
JA3(1,3) = -links(3).mcm(2) + T34(1,3);
JA3(1,4) = links(3).m + T34(1,4);
JA3(1,5) = T34(1,5);
JA3(1,6) = T34(1,6);

JA3(2,1) = -links(3).mcm(3) + T34(2,1);
JA3(2,2) = T34(2,2);
JA3(2,3) = links(3).mcm(1) + T34(2,3);
JA3(2,4) = T34(2,4);
JA3(2,5) = links(3).m + T34(2,5);
JA3(2,6) = T34(2,6);

JA3(3,1) = links(3).mcm(2) + T34(3,1);
JA3(3,2) = -links(3).mcm(1) + T34(3,2);
JA3(3,3) = T34(3,3);
JA3(3,4) = T34(3,4);
JA3(3,5) = T34(3,5);
JA3(3,6) = links(3).m + T34(3,6);

JA3(4,1) = links(3).inertia(1,1) + T34(4,1);
JA3(4,2) = links(3).inertia(1,2) + T34(4,2);
JA3(4,3) = links(3).inertia(1,3) + T34(4,3);
JA3(4,4) = T34(4,4);
JA3(4,5) = -links(3).mcm(3) + T34(4,5);
JA3(4,6) = links(3).mcm(2) + T34(4,6);

JA3(5,1) = links(3).inertia(1,2) + T34(5,1);
JA3(5,2) = links(3).inertia(2,2) + T34(5,2);
JA3(5,3) = links(3).inertia(2,3) + T34(5,3);
JA3(5,4) = links(3).mcm(3) + T34(5,4);
JA3(5,5) = T34(5,5);
JA3(5,6) = -links(3).mcm(1) + T34(5,6);

JA3(6,1) = links(3).inertia(1,3) + T34(6,1);
JA3(6,2) = links(3).inertia(2,3) + T34(6,2);
JA3(6,3) = links(3).inertia(3,3) + T34(6,3);
JA3(6,4) = -links(3).mcm(2) + T34(6,4);
JA3(6,5) = links(3).mcm(1) + T34(6,5);
JA3(6,6) = T34(6,6);


h3(1) = JA3(1,3);
h3(2) = JA3(2,3);
h3(3) = JA3(3,3);
h3(4) = JA3(4,3);
h3(5) = JA3(5,3);
h3(6) = JA3(6,3);

T123(1,1) = JA3(1,1);
T123(1,2) = JA3(1,2);
T123(1,3) = JA3(1,3);
T123(1,4) = JA3(1,4);
T123(1,5) = JA3(1,5);
T123(1,6) = JA3(1,6);

T123(2,1) = JA3(2,1);
T123(2,2) = JA3(2,2);
T123(2,3) = JA3(2,3);
T123(2,4) = JA3(2,4);
T123(2,5) = JA3(2,5);
T123(2,6) = JA3(2,6);

T123(3,1) = JA3(3,1);
T123(3,2) = JA3(3,2);
T123(3,3) = JA3(3,3);
T123(3,4) = JA3(3,4);
T123(3,5) = JA3(3,5);
T123(3,6) = JA3(3,6);

T123(4,1) = JA3(4,1);
T123(4,2) = JA3(4,2);
T123(4,3) = JA3(4,3);
T123(4,4) = JA3(4,4);
T123(4,5) = JA3(4,5);
T123(4,6) = JA3(4,6);

T123(5,1) = JA3(5,1);
T123(5,2) = JA3(5,2);
T123(5,3) = JA3(5,3);
T123(5,4) = JA3(5,4);
T123(5,5) = JA3(5,5);
T123(5,6) = JA3(5,6);

T123(6,1) = JA3(6,1);
T123(6,2) = JA3(6,2);
T123(6,3) = JA3(6,3);
T123(6,4) = JA3(6,4);
T123(6,5) = JA3(6,5);
T123(6,6) = JA3(6,6);


T23(1,1) = T123(3,3);
T23(1,2) = S32(1,2)*T123(3,1) + S32(2,2)*T123(3,2) - ZHR*S32(1,3)*T123(3,4) - ZHR*S32(2,3)*T123(3,5);
T23(1,3) = S32(1,3)*T123(3,1) + S32(2,3)*T123(3,2) + ZHR*S32(1,2)*T123(3,4) + ZHR*S32(2,2)*T123(3,5);
T23(1,4) = T123(3,6);
T23(1,5) = S32(1,2)*T123(3,4) + S32(2,2)*T123(3,5);
T23(1,6) = S32(1,3)*T123(3,4) + S32(2,3)*T123(3,5);

T23(2,1) = Si23(2,1)*T123(1,3) + Si23(2,2)*T123(2,3);
T23(2,2) = S32(1,2)*(Si23(2,1)*T123(1,1) + Si23(2,2)*T123(2,1)) + S32(2,2)*(Si23(2,1)*T123(1,2) + Si23(2,2)*T123(2,2)) - ZHR*S32(1,3)*(Si23(2,1)*T123(1,4) + Si23(2,2)*T123(2,4)) - ZHR*S32(2,3)*(Si23(2,1)*T123(1,5) + Si23(2,2)*T123(2,5));
T23(2,3) = S32(1,3)*(Si23(2,1)*T123(1,1) + Si23(2,2)*T123(2,1)) + S32(2,3)*(Si23(2,1)*T123(1,2) + Si23(2,2)*T123(2,2)) + ZHR*S32(1,2)*(Si23(2,1)*T123(1,4) + Si23(2,2)*T123(2,4)) + ZHR*S32(2,2)*(Si23(2,1)*T123(1,5) + Si23(2,2)*T123(2,5));
T23(2,4) = Si23(2,1)*T123(1,6) + Si23(2,2)*T123(2,6);
T23(2,5) = S32(1,2)*(Si23(2,1)*T123(1,4) + Si23(2,2)*T123(2,4)) + S32(2,2)*(Si23(2,1)*T123(1,5) + Si23(2,2)*T123(2,5));
T23(2,6) = S32(1,3)*(Si23(2,1)*T123(1,4) + Si23(2,2)*T123(2,4)) + S32(2,3)*(Si23(2,1)*T123(1,5) + Si23(2,2)*T123(2,5));

T23(3,1) = Si23(3,1)*T123(1,3) + Si23(3,2)*T123(2,3);
T23(3,2) = S32(1,2)*(Si23(3,1)*T123(1,1) + Si23(3,2)*T123(2,1)) + S32(2,2)*(Si23(3,1)*T123(1,2) + Si23(3,2)*T123(2,2)) - ZHR*S32(1,3)*(Si23(3,1)*T123(1,4) + Si23(3,2)*T123(2,4)) - ZHR*S32(2,3)*(Si23(3,1)*T123(1,5) + Si23(3,2)*T123(2,5));
T23(3,3) = S32(1,3)*(Si23(3,1)*T123(1,1) + Si23(3,2)*T123(2,1)) + S32(2,3)*(Si23(3,1)*T123(1,2) + Si23(3,2)*T123(2,2)) + ZHR*S32(1,2)*(Si23(3,1)*T123(1,4) + Si23(3,2)*T123(2,4)) + ZHR*S32(2,2)*(Si23(3,1)*T123(1,5) + Si23(3,2)*T123(2,5));
T23(3,4) = Si23(3,1)*T123(1,6) + Si23(3,2)*T123(2,6);
T23(3,5) = S32(1,2)*(Si23(3,1)*T123(1,4) + Si23(3,2)*T123(2,4)) + S32(2,2)*(Si23(3,1)*T123(1,5) + Si23(3,2)*T123(2,5));
T23(3,6) = S32(1,3)*(Si23(3,1)*T123(1,4) + Si23(3,2)*T123(2,4)) + S32(2,3)*(Si23(3,1)*T123(1,5) + Si23(3,2)*T123(2,5));

T23(4,1) = T123(6,3);
T23(4,2) = S32(1,2)*T123(6,1) + S32(2,2)*T123(6,2) - ZHR*S32(1,3)*T123(6,4) - ZHR*S32(2,3)*T123(6,5);
T23(4,3) = S32(1,3)*T123(6,1) + S32(2,3)*T123(6,2) + ZHR*S32(1,2)*T123(6,4) + ZHR*S32(2,2)*T123(6,5);
T23(4,4) = T123(6,6);
T23(4,5) = S32(1,2)*T123(6,4) + S32(2,2)*T123(6,5);
T23(4,6) = S32(1,3)*T123(6,4) + S32(2,3)*T123(6,5);

T23(5,1) = -(ZHR*Si23(3,1)*T123(1,3)) - ZHR*Si23(3,2)*T123(2,3) + Si23(2,1)*T123(4,3) + Si23(2,2)*T123(5,3);
T23(5,2) = S32(1,2)*(-(ZHR*Si23(3,1)*T123(1,1)) - ZHR*Si23(3,2)*T123(2,1) + Si23(2,1)*T123(4,1) + Si23(2,2)*T123(5,1)) + S32(2,2)*(-(ZHR*Si23(3,1)*T123(1,2)) - ZHR*Si23(3,2)*T123(2,2) + Si23(2,1)*T123(4,2) + Si23(2,2)*T123(5,2)) - ZHR*S32(1,3)*(-(ZHR*Si23(3,1)*T123(1,4)) - ZHR*Si23(3,2)*T123(2,4) + Si23(2,1)*T123(4,4) + Si23(2,2)*T123(5,4)) - ZHR*S32(2,3)*(-(ZHR*Si23(3,1)*T123(1,5)) - ZHR*Si23(3,2)*T123(2,5) + Si23(2,1)*T123(4,5) + Si23(2,2)*T123(5,5));
T23(5,3) = S32(1,3)*(-(ZHR*Si23(3,1)*T123(1,1)) - ZHR*Si23(3,2)*T123(2,1) + Si23(2,1)*T123(4,1) + Si23(2,2)*T123(5,1)) + S32(2,3)*(-(ZHR*Si23(3,1)*T123(1,2)) - ZHR*Si23(3,2)*T123(2,2) + Si23(2,1)*T123(4,2) + Si23(2,2)*T123(5,2)) + ZHR*S32(1,2)*(-(ZHR*Si23(3,1)*T123(1,4)) - ZHR*Si23(3,2)*T123(2,4) + Si23(2,1)*T123(4,4) + Si23(2,2)*T123(5,4)) + ZHR*S32(2,2)*(-(ZHR*Si23(3,1)*T123(1,5)) - ZHR*Si23(3,2)*T123(2,5) + Si23(2,1)*T123(4,5) + Si23(2,2)*T123(5,5));
T23(5,4) = -(ZHR*Si23(3,1)*T123(1,6)) - ZHR*Si23(3,2)*T123(2,6) + Si23(2,1)*T123(4,6) + Si23(2,2)*T123(5,6);
T23(5,5) = S32(1,2)*(-(ZHR*Si23(3,1)*T123(1,4)) - ZHR*Si23(3,2)*T123(2,4) + Si23(2,1)*T123(4,4) + Si23(2,2)*T123(5,4)) + S32(2,2)*(-(ZHR*Si23(3,1)*T123(1,5)) - ZHR*Si23(3,2)*T123(2,5) + Si23(2,1)*T123(4,5) + Si23(2,2)*T123(5,5));
T23(5,6) = S32(1,3)*(-(ZHR*Si23(3,1)*T123(1,4)) - ZHR*Si23(3,2)*T123(2,4) + Si23(2,1)*T123(4,4) + Si23(2,2)*T123(5,4)) + S32(2,3)*(-(ZHR*Si23(3,1)*T123(1,5)) - ZHR*Si23(3,2)*T123(2,5) + Si23(2,1)*T123(4,5) + Si23(2,2)*T123(5,5));

T23(6,1) = ZHR*Si23(2,1)*T123(1,3) + ZHR*Si23(2,2)*T123(2,3) + Si23(3,1)*T123(4,3) + Si23(3,2)*T123(5,3);
T23(6,2) = S32(1,2)*(ZHR*Si23(2,1)*T123(1,1) + ZHR*Si23(2,2)*T123(2,1) + Si23(3,1)*T123(4,1) + Si23(3,2)*T123(5,1)) + S32(2,2)*(ZHR*Si23(2,1)*T123(1,2) + ZHR*Si23(2,2)*T123(2,2) + Si23(3,1)*T123(4,2) + Si23(3,2)*T123(5,2)) - ZHR*S32(1,3)*(ZHR*Si23(2,1)*T123(1,4) + ZHR*Si23(2,2)*T123(2,4) + Si23(3,1)*T123(4,4) + Si23(3,2)*T123(5,4)) - ZHR*S32(2,3)*(ZHR*Si23(2,1)*T123(1,5) + ZHR*Si23(2,2)*T123(2,5) + Si23(3,1)*T123(4,5) + Si23(3,2)*T123(5,5));
T23(6,3) = S32(1,3)*(ZHR*Si23(2,1)*T123(1,1) + ZHR*Si23(2,2)*T123(2,1) + Si23(3,1)*T123(4,1) + Si23(3,2)*T123(5,1)) + S32(2,3)*(ZHR*Si23(2,1)*T123(1,2) + ZHR*Si23(2,2)*T123(2,2) + Si23(3,1)*T123(4,2) + Si23(3,2)*T123(5,2)) + ZHR*S32(1,2)*(ZHR*Si23(2,1)*T123(1,4) + ZHR*Si23(2,2)*T123(2,4) + Si23(3,1)*T123(4,4) + Si23(3,2)*T123(5,4)) + ZHR*S32(2,2)*(ZHR*Si23(2,1)*T123(1,5) + ZHR*Si23(2,2)*T123(2,5) + Si23(3,1)*T123(4,5) + Si23(3,2)*T123(5,5));
T23(6,4) = ZHR*Si23(2,1)*T123(1,6) + ZHR*Si23(2,2)*T123(2,6) + Si23(3,1)*T123(4,6) + Si23(3,2)*T123(5,6);
T23(6,5) = S32(1,2)*(ZHR*Si23(2,1)*T123(1,4) + ZHR*Si23(2,2)*T123(2,4) + Si23(3,1)*T123(4,4) + Si23(3,2)*T123(5,4)) + S32(2,2)*(ZHR*Si23(2,1)*T123(1,5) + ZHR*Si23(2,2)*T123(2,5) + Si23(3,1)*T123(4,5) + Si23(3,2)*T123(5,5));
T23(6,6) = S32(1,3)*(ZHR*Si23(2,1)*T123(1,4) + ZHR*Si23(2,2)*T123(2,4) + Si23(3,1)*T123(4,4) + Si23(3,2)*T123(5,4)) + S32(2,3)*(ZHR*Si23(2,1)*T123(1,5) + ZHR*Si23(2,2)*T123(2,5) + Si23(3,1)*T123(4,5) + Si23(3,2)*T123(5,5));







%barrett_InvDynArtfunc13
      
JA2(1,1) = T23(1,1);
JA2(1,2) = links(2).mcm(3) + T23(1,2);
JA2(1,3) = -links(2).mcm(2) + T23(1,3);
JA2(1,4) = links(2).m + T23(1,4);
JA2(1,5) = T23(1,5);
JA2(1,6) = T23(1,6);

JA2(2,1) = -links(2).mcm(3) + T23(2,1);
JA2(2,2) = T23(2,2);
JA2(2,3) = links(2).mcm(1) + T23(2,3);
JA2(2,4) = T23(2,4);
JA2(2,5) = links(2).m + T23(2,5);
JA2(2,6) = T23(2,6);

JA2(3,1) = links(2).mcm(2) + T23(3,1);
JA2(3,2) = -links(2).mcm(1) + T23(3,2);
JA2(3,3) = T23(3,3);
JA2(3,4) = T23(3,4);
JA2(3,5) = T23(3,5);
JA2(3,6) = links(2).m + T23(3,6);

JA2(4,1) = links(2).inertia(1,1) + T23(4,1);
JA2(4,2) = links(2).inertia(1,2) + T23(4,2);
JA2(4,3) = links(2).inertia(1,3) + T23(4,3);
JA2(4,4) = T23(4,4);
JA2(4,5) = -links(2).mcm(3) + T23(4,5);
JA2(4,6) = links(2).mcm(2) + T23(4,6);

JA2(5,1) = links(2).inertia(1,2) + T23(5,1);
JA2(5,2) = links(2).inertia(2,2) + T23(5,2);
JA2(5,3) = links(2).inertia(2,3) + T23(5,3);
JA2(5,4) = links(2).mcm(3) + T23(5,4);
JA2(5,5) = T23(5,5);
JA2(5,6) = -links(2).mcm(1) + T23(5,6);

JA2(6,1) = links(2).inertia(1,3) + T23(6,1);
JA2(6,2) = links(2).inertia(2,3) + T23(6,2);
JA2(6,3) = links(2).inertia(3,3) + T23(6,3);
JA2(6,4) = -links(2).mcm(2) + T23(6,4);
JA2(6,5) = links(2).mcm(1) + T23(6,5);
JA2(6,6) = T23(6,6);


h2(1) = JA2(1,3);
h2(2) = JA2(2,3);
h2(3) = JA2(3,3);
h2(4) = JA2(4,3);
h2(5) = JA2(5,3);
h2(6) = JA2(6,3);

T112(1,1) = JA2(1,1);
T112(1,2) = JA2(1,2);
T112(1,3) = JA2(1,3);
T112(1,4) = JA2(1,4);
T112(1,5) = JA2(1,5);
T112(1,6) = JA2(1,6);

T112(2,1) = JA2(2,1);
T112(2,2) = JA2(2,2);
T112(2,3) = JA2(2,3);
T112(2,4) = JA2(2,4);
T112(2,5) = JA2(2,5);
T112(2,6) = JA2(2,6);

T112(3,1) = JA2(3,1);
T112(3,2) = JA2(3,2);
T112(3,3) = JA2(3,3);
T112(3,4) = JA2(3,4);
T112(3,5) = JA2(3,5);
T112(3,6) = JA2(3,6);

T112(4,1) = JA2(4,1);
T112(4,2) = JA2(4,2);
T112(4,3) = JA2(4,3);
T112(4,4) = JA2(4,4);
T112(4,5) = JA2(4,5);
T112(4,6) = JA2(4,6);

T112(5,1) = JA2(5,1);
T112(5,2) = JA2(5,2);
T112(5,3) = JA2(5,3);
T112(5,4) = JA2(5,4);
T112(5,5) = JA2(5,5);
T112(5,6) = JA2(5,6);

T112(6,1) = JA2(6,1);
T112(6,2) = JA2(6,2);
T112(6,3) = JA2(6,3);
T112(6,4) = JA2(6,4);
T112(6,5) = JA2(6,5);
T112(6,6) = JA2(6,6);


T12(1,1) = T112(3,3);
T12(1,2) = -(S21(1,2)*T112(3,1)) - S21(2,2)*T112(3,2);
T12(1,3) = -(S21(1,3)*T112(3,1)) - S21(2,3)*T112(3,2);
T12(1,4) = T112(3,6);
T12(1,5) = -(S21(1,2)*T112(3,4)) - S21(2,2)*T112(3,5);
T12(1,6) = -(S21(1,3)*T112(3,4)) - S21(2,3)*T112(3,5);

T12(2,1) = -(Si12(2,1)*T112(1,3)) - Si12(2,2)*T112(2,3);
T12(2,2) = S21(1,2)*(Si12(2,1)*T112(1,1) + Si12(2,2)*T112(2,1)) + S21(2,2)*(Si12(2,1)*T112(1,2) + Si12(2,2)*T112(2,2));
T12(2,3) = S21(1,3)*(Si12(2,1)*T112(1,1) + Si12(2,2)*T112(2,1)) + S21(2,3)*(Si12(2,1)*T112(1,2) + Si12(2,2)*T112(2,2));
T12(2,4) = -(Si12(2,1)*T112(1,6)) - Si12(2,2)*T112(2,6);
T12(2,5) = S21(1,2)*(Si12(2,1)*T112(1,4) + Si12(2,2)*T112(2,4)) + S21(2,2)*(Si12(2,1)*T112(1,5) + Si12(2,2)*T112(2,5));
T12(2,6) = S21(1,3)*(Si12(2,1)*T112(1,4) + Si12(2,2)*T112(2,4)) + S21(2,3)*(Si12(2,1)*T112(1,5) + Si12(2,2)*T112(2,5));

T12(3,1) = -(Si12(3,1)*T112(1,3)) - Si12(3,2)*T112(2,3);
T12(3,2) = S21(1,2)*(Si12(3,1)*T112(1,1) + Si12(3,2)*T112(2,1)) + S21(2,2)*(Si12(3,1)*T112(1,2) + Si12(3,2)*T112(2,2));
T12(3,3) = S21(1,3)*(Si12(3,1)*T112(1,1) + Si12(3,2)*T112(2,1)) + S21(2,3)*(Si12(3,1)*T112(1,2) + Si12(3,2)*T112(2,2));
T12(3,4) = -(Si12(3,1)*T112(1,6)) - Si12(3,2)*T112(2,6);
T12(3,5) = S21(1,2)*(Si12(3,1)*T112(1,4) + Si12(3,2)*T112(2,4)) + S21(2,2)*(Si12(3,1)*T112(1,5) + Si12(3,2)*T112(2,5));
T12(3,6) = S21(1,3)*(Si12(3,1)*T112(1,4) + Si12(3,2)*T112(2,4)) + S21(2,3)*(Si12(3,1)*T112(1,5) + Si12(3,2)*T112(2,5));

T12(4,1) = T112(6,3);
T12(4,2) = -(S21(1,2)*T112(6,1)) - S21(2,2)*T112(6,2);
T12(4,3) = -(S21(1,3)*T112(6,1)) - S21(2,3)*T112(6,2);
T12(4,4) = T112(6,6);
T12(4,5) = -(S21(1,2)*T112(6,4)) - S21(2,2)*T112(6,5);
T12(4,6) = -(S21(1,3)*T112(6,4)) - S21(2,3)*T112(6,5);

T12(5,1) = -(Si12(2,1)*T112(4,3)) - Si12(2,2)*T112(5,3);
T12(5,2) = S21(1,2)*(Si12(2,1)*T112(4,1) + Si12(2,2)*T112(5,1)) + S21(2,2)*(Si12(2,1)*T112(4,2) + Si12(2,2)*T112(5,2));
T12(5,3) = S21(1,3)*(Si12(2,1)*T112(4,1) + Si12(2,2)*T112(5,1)) + S21(2,3)*(Si12(2,1)*T112(4,2) + Si12(2,2)*T112(5,2));
T12(5,4) = -(Si12(2,1)*T112(4,6)) - Si12(2,2)*T112(5,6);
T12(5,5) = S21(1,2)*(Si12(2,1)*T112(4,4) + Si12(2,2)*T112(5,4)) + S21(2,2)*(Si12(2,1)*T112(4,5) + Si12(2,2)*T112(5,5));
T12(5,6) = S21(1,3)*(Si12(2,1)*T112(4,4) + Si12(2,2)*T112(5,4)) + S21(2,3)*(Si12(2,1)*T112(4,5) + Si12(2,2)*T112(5,5));

T12(6,1) = -(Si12(3,1)*T112(4,3)) - Si12(3,2)*T112(5,3);
T12(6,2) = S21(1,2)*(Si12(3,1)*T112(4,1) + Si12(3,2)*T112(5,1)) + S21(2,2)*(Si12(3,1)*T112(4,2) + Si12(3,2)*T112(5,2));
T12(6,3) = S21(1,3)*(Si12(3,1)*T112(4,1) + Si12(3,2)*T112(5,1)) + S21(2,3)*(Si12(3,1)*T112(4,2) + Si12(3,2)*T112(5,2));
T12(6,4) = -(Si12(3,1)*T112(4,6)) - Si12(3,2)*T112(5,6);
T12(6,5) = S21(1,2)*(Si12(3,1)*T112(4,4) + Si12(3,2)*T112(5,4)) + S21(2,2)*(Si12(3,1)*T112(4,5) + Si12(3,2)*T112(5,5));
T12(6,6) = S21(1,3)*(Si12(3,1)*T112(4,4) + Si12(3,2)*T112(5,4)) + S21(2,3)*(Si12(3,1)*T112(4,5) + Si12(3,2)*T112(5,5));







%barrett_InvDynArtfunc14
      
JA1(1,1) = T12(1,1);
JA1(1,2) = links(1).mcm(3) + T12(1,2);
JA1(1,3) = -links(1).mcm(2) + T12(1,3);
JA1(1,4) = links(1).m + T12(1,4);
JA1(1,5) = T12(1,5);
JA1(1,6) = T12(1,6);

JA1(2,1) = -links(1).mcm(3) + T12(2,1);
JA1(2,2) = T12(2,2);
JA1(2,3) = links(1).mcm(1) + T12(2,3);
JA1(2,4) = T12(2,4);
JA1(2,5) = links(1).m + T12(2,5);
JA1(2,6) = T12(2,6);

JA1(3,1) = links(1).mcm(2) + T12(3,1);
JA1(3,2) = -links(1).mcm(1) + T12(3,2);
JA1(3,3) = T12(3,3);
JA1(3,4) = T12(3,4);
JA1(3,5) = T12(3,5);
JA1(3,6) = links(1).m + T12(3,6);

JA1(4,1) = links(1).inertia(1,1) + T12(4,1);
JA1(4,2) = links(1).inertia(1,2) + T12(4,2);
JA1(4,3) = links(1).inertia(1,3) + T12(4,3);
JA1(4,4) = T12(4,4);
JA1(4,5) = -links(1).mcm(3) + T12(4,5);
JA1(4,6) = links(1).mcm(2) + T12(4,6);

JA1(5,1) = links(1).inertia(1,2) + T12(5,1);
JA1(5,2) = links(1).inertia(2,2) + T12(5,2);
JA1(5,3) = links(1).inertia(2,3) + T12(5,3);
JA1(5,4) = links(1).mcm(3) + T12(5,4);
JA1(5,5) = T12(5,5);
JA1(5,6) = -links(1).mcm(1) + T12(5,6);

JA1(6,1) = links(1).inertia(1,3) + T12(6,1);
JA1(6,2) = links(1).inertia(2,3) + T12(6,2);
JA1(6,3) = links(1).inertia(3,3) + T12(6,3);
JA1(6,4) = -links(1).mcm(2) + T12(6,4);
JA1(6,5) = links(1).mcm(1) + T12(6,5);
JA1(6,6) = T12(6,6);


h1(1) = JA1(1,3);
h1(2) = JA1(2,3);
h1(3) = JA1(3,3);
h1(4) = JA1(4,3);
h1(5) = JA1(5,3);
h1(6) = JA1(6,3);

T101(1,1) = JA1(1,1);
T101(1,2) = JA1(1,2);
T101(1,3) = JA1(1,3);
T101(1,4) = JA1(1,4);
T101(1,5) = JA1(1,5);
T101(1,6) = JA1(1,6);

T101(2,1) = JA1(2,1);
T101(2,2) = JA1(2,2);
T101(2,3) = JA1(2,3);
T101(2,4) = JA1(2,4);
T101(2,5) = JA1(2,5);
T101(2,6) = JA1(2,6);

T101(3,1) = JA1(3,1);
T101(3,2) = JA1(3,2);
T101(3,3) = JA1(3,3);
T101(3,4) = JA1(3,4);
T101(3,5) = JA1(3,5);
T101(3,6) = JA1(3,6);

T101(4,1) = JA1(4,1);
T101(4,2) = JA1(4,2);
T101(4,3) = JA1(4,3);
T101(4,4) = JA1(4,4);
T101(4,5) = JA1(4,5);
T101(4,6) = JA1(4,6);

T101(5,1) = JA1(5,1);
T101(5,2) = JA1(5,2);
T101(5,3) = JA1(5,3);
T101(5,4) = JA1(5,4);
T101(5,5) = JA1(5,5);
T101(5,6) = JA1(5,6);

T101(6,1) = JA1(6,1);
T101(6,2) = JA1(6,2);
T101(6,3) = JA1(6,3);
T101(6,4) = JA1(6,4);
T101(6,5) = JA1(6,5);
T101(6,6) = JA1(6,6);


T01(1,1) = S10(1,1)*(Si01(1,1)*T101(1,1) + Si01(1,2)*T101(2,1)) + S10(2,1)*(Si01(1,1)*T101(1,2) + Si01(1,2)*T101(2,2)) - ZSFE*S10(1,2)*(Si01(1,1)*T101(1,4) + Si01(1,2)*T101(2,4)) - ZSFE*S10(2,2)*(Si01(1,1)*T101(1,5) + Si01(1,2)*T101(2,5));
T01(1,2) = S10(1,2)*(Si01(1,1)*T101(1,1) + Si01(1,2)*T101(2,1)) + S10(2,2)*(Si01(1,1)*T101(1,2) + Si01(1,2)*T101(2,2)) + ZSFE*S10(1,1)*(Si01(1,1)*T101(1,4) + Si01(1,2)*T101(2,4)) + ZSFE*S10(2,1)*(Si01(1,1)*T101(1,5) + Si01(1,2)*T101(2,5));
T01(1,3) = Si01(1,1)*T101(1,3) + Si01(1,2)*T101(2,3);
T01(1,4) = S10(1,1)*(Si01(1,1)*T101(1,4) + Si01(1,2)*T101(2,4)) + S10(2,1)*(Si01(1,1)*T101(1,5) + Si01(1,2)*T101(2,5));
T01(1,5) = S10(1,2)*(Si01(1,1)*T101(1,4) + Si01(1,2)*T101(2,4)) + S10(2,2)*(Si01(1,1)*T101(1,5) + Si01(1,2)*T101(2,5));
T01(1,6) = Si01(1,1)*T101(1,6) + Si01(1,2)*T101(2,6);

T01(2,1) = S10(1,1)*(Si01(2,1)*T101(1,1) + Si01(2,2)*T101(2,1)) + S10(2,1)*(Si01(2,1)*T101(1,2) + Si01(2,2)*T101(2,2)) - ZSFE*S10(1,2)*(Si01(2,1)*T101(1,4) + Si01(2,2)*T101(2,4)) - ZSFE*S10(2,2)*(Si01(2,1)*T101(1,5) + Si01(2,2)*T101(2,5));
T01(2,2) = S10(1,2)*(Si01(2,1)*T101(1,1) + Si01(2,2)*T101(2,1)) + S10(2,2)*(Si01(2,1)*T101(1,2) + Si01(2,2)*T101(2,2)) + ZSFE*S10(1,1)*(Si01(2,1)*T101(1,4) + Si01(2,2)*T101(2,4)) + ZSFE*S10(2,1)*(Si01(2,1)*T101(1,5) + Si01(2,2)*T101(2,5));
T01(2,3) = Si01(2,1)*T101(1,3) + Si01(2,2)*T101(2,3);
T01(2,4) = S10(1,1)*(Si01(2,1)*T101(1,4) + Si01(2,2)*T101(2,4)) + S10(2,1)*(Si01(2,1)*T101(1,5) + Si01(2,2)*T101(2,5));
T01(2,5) = S10(1,2)*(Si01(2,1)*T101(1,4) + Si01(2,2)*T101(2,4)) + S10(2,2)*(Si01(2,1)*T101(1,5) + Si01(2,2)*T101(2,5));
T01(2,6) = Si01(2,1)*T101(1,6) + Si01(2,2)*T101(2,6);

T01(3,1) = S10(1,1)*T101(3,1) + S10(2,1)*T101(3,2) - ZSFE*S10(1,2)*T101(3,4) - ZSFE*S10(2,2)*T101(3,5);
T01(3,2) = S10(1,2)*T101(3,1) + S10(2,2)*T101(3,2) + ZSFE*S10(1,1)*T101(3,4) + ZSFE*S10(2,1)*T101(3,5);
T01(3,3) = T101(3,3);
T01(3,4) = S10(1,1)*T101(3,4) + S10(2,1)*T101(3,5);
T01(3,5) = S10(1,2)*T101(3,4) + S10(2,2)*T101(3,5);
T01(3,6) = T101(3,6);

T01(4,1) = S10(1,1)*(-(ZSFE*Si01(2,1)*T101(1,1)) - ZSFE*Si01(2,2)*T101(2,1) + Si01(1,1)*T101(4,1) + Si01(1,2)*T101(5,1)) + S10(2,1)*(-(ZSFE*Si01(2,1)*T101(1,2)) - ZSFE*Si01(2,2)*T101(2,2) + Si01(1,1)*T101(4,2) + Si01(1,2)*T101(5,2)) - ZSFE*S10(1,2)*(-(ZSFE*Si01(2,1)*T101(1,4)) - ZSFE*Si01(2,2)*T101(2,4) + Si01(1,1)*T101(4,4) + Si01(1,2)*T101(5,4)) - ZSFE*S10(2,2)*(-(ZSFE*Si01(2,1)*T101(1,5)) - ZSFE*Si01(2,2)*T101(2,5) + Si01(1,1)*T101(4,5) + Si01(1,2)*T101(5,5));
T01(4,2) = S10(1,2)*(-(ZSFE*Si01(2,1)*T101(1,1)) - ZSFE*Si01(2,2)*T101(2,1) + Si01(1,1)*T101(4,1) + Si01(1,2)*T101(5,1)) + S10(2,2)*(-(ZSFE*Si01(2,1)*T101(1,2)) - ZSFE*Si01(2,2)*T101(2,2) + Si01(1,1)*T101(4,2) + Si01(1,2)*T101(5,2)) + ZSFE*S10(1,1)*(-(ZSFE*Si01(2,1)*T101(1,4)) - ZSFE*Si01(2,2)*T101(2,4) + Si01(1,1)*T101(4,4) + Si01(1,2)*T101(5,4)) + ZSFE*S10(2,1)*(-(ZSFE*Si01(2,1)*T101(1,5)) - ZSFE*Si01(2,2)*T101(2,5) + Si01(1,1)*T101(4,5) + Si01(1,2)*T101(5,5));
T01(4,3) = -(ZSFE*Si01(2,1)*T101(1,3)) - ZSFE*Si01(2,2)*T101(2,3) + Si01(1,1)*T101(4,3) + Si01(1,2)*T101(5,3);
T01(4,4) = S10(1,1)*(-(ZSFE*Si01(2,1)*T101(1,4)) - ZSFE*Si01(2,2)*T101(2,4) + Si01(1,1)*T101(4,4) + Si01(1,2)*T101(5,4)) + S10(2,1)*(-(ZSFE*Si01(2,1)*T101(1,5)) - ZSFE*Si01(2,2)*T101(2,5) + Si01(1,1)*T101(4,5) + Si01(1,2)*T101(5,5));
T01(4,5) = S10(1,2)*(-(ZSFE*Si01(2,1)*T101(1,4)) - ZSFE*Si01(2,2)*T101(2,4) + Si01(1,1)*T101(4,4) + Si01(1,2)*T101(5,4)) + S10(2,2)*(-(ZSFE*Si01(2,1)*T101(1,5)) - ZSFE*Si01(2,2)*T101(2,5) + Si01(1,1)*T101(4,5) + Si01(1,2)*T101(5,5));
T01(4,6) = -(ZSFE*Si01(2,1)*T101(1,6)) - ZSFE*Si01(2,2)*T101(2,6) + Si01(1,1)*T101(4,6) + Si01(1,2)*T101(5,6);

T01(5,1) = S10(1,1)*(ZSFE*Si01(1,1)*T101(1,1) + ZSFE*Si01(1,2)*T101(2,1) + Si01(2,1)*T101(4,1) + Si01(2,2)*T101(5,1)) + S10(2,1)*(ZSFE*Si01(1,1)*T101(1,2) + ZSFE*Si01(1,2)*T101(2,2) + Si01(2,1)*T101(4,2) + Si01(2,2)*T101(5,2)) - ZSFE*S10(1,2)*(ZSFE*Si01(1,1)*T101(1,4) + ZSFE*Si01(1,2)*T101(2,4) + Si01(2,1)*T101(4,4) + Si01(2,2)*T101(5,4)) - ZSFE*S10(2,2)*(ZSFE*Si01(1,1)*T101(1,5) + ZSFE*Si01(1,2)*T101(2,5) + Si01(2,1)*T101(4,5) + Si01(2,2)*T101(5,5));
T01(5,2) = S10(1,2)*(ZSFE*Si01(1,1)*T101(1,1) + ZSFE*Si01(1,2)*T101(2,1) + Si01(2,1)*T101(4,1) + Si01(2,2)*T101(5,1)) + S10(2,2)*(ZSFE*Si01(1,1)*T101(1,2) + ZSFE*Si01(1,2)*T101(2,2) + Si01(2,1)*T101(4,2) + Si01(2,2)*T101(5,2)) + ZSFE*S10(1,1)*(ZSFE*Si01(1,1)*T101(1,4) + ZSFE*Si01(1,2)*T101(2,4) + Si01(2,1)*T101(4,4) + Si01(2,2)*T101(5,4)) + ZSFE*S10(2,1)*(ZSFE*Si01(1,1)*T101(1,5) + ZSFE*Si01(1,2)*T101(2,5) + Si01(2,1)*T101(4,5) + Si01(2,2)*T101(5,5));
T01(5,3) = ZSFE*Si01(1,1)*T101(1,3) + ZSFE*Si01(1,2)*T101(2,3) + Si01(2,1)*T101(4,3) + Si01(2,2)*T101(5,3);
T01(5,4) = S10(1,1)*(ZSFE*Si01(1,1)*T101(1,4) + ZSFE*Si01(1,2)*T101(2,4) + Si01(2,1)*T101(4,4) + Si01(2,2)*T101(5,4)) + S10(2,1)*(ZSFE*Si01(1,1)*T101(1,5) + ZSFE*Si01(1,2)*T101(2,5) + Si01(2,1)*T101(4,5) + Si01(2,2)*T101(5,5));
T01(5,5) = S10(1,2)*(ZSFE*Si01(1,1)*T101(1,4) + ZSFE*Si01(1,2)*T101(2,4) + Si01(2,1)*T101(4,4) + Si01(2,2)*T101(5,4)) + S10(2,2)*(ZSFE*Si01(1,1)*T101(1,5) + ZSFE*Si01(1,2)*T101(2,5) + Si01(2,1)*T101(4,5) + Si01(2,2)*T101(5,5));
T01(5,6) = ZSFE*Si01(1,1)*T101(1,6) + ZSFE*Si01(1,2)*T101(2,6) + Si01(2,1)*T101(4,6) + Si01(2,2)*T101(5,6);

T01(6,1) = S10(1,1)*T101(6,1) + S10(2,1)*T101(6,2) - ZSFE*S10(1,2)*T101(6,4) - ZSFE*S10(2,2)*T101(6,5);
T01(6,2) = S10(1,2)*T101(6,1) + S10(2,2)*T101(6,2) + ZSFE*S10(1,1)*T101(6,4) + ZSFE*S10(2,1)*T101(6,5);
T01(6,3) = T101(6,3);
T01(6,4) = S10(1,1)*T101(6,4) + S10(2,1)*T101(6,5);
T01(6,5) = S10(1,2)*T101(6,4) + S10(2,2)*T101(6,5);
T01(6,6) = T101(6,6);







%barrett_InvDynArtfunc15
      
JA0(1,1) = T01(1,1);
JA0(1,2) = link0.mcm(3) + T01(1,2);
JA0(1,3) = -link0.mcm(2) + T01(1,3);
JA0(1,4) = link0.m + T01(1,4);
JA0(1,5) = T01(1,5);
JA0(1,6) = T01(1,6);

JA0(2,1) = -link0.mcm(3) + T01(2,1);
JA0(2,2) = T01(2,2);
JA0(2,3) = link0.mcm(1) + T01(2,3);
JA0(2,4) = T01(2,4);
JA0(2,5) = link0.m + T01(2,5);
JA0(2,6) = T01(2,6);

JA0(3,1) = link0.mcm(2) + T01(3,1);
JA0(3,2) = -link0.mcm(1) + T01(3,2);
JA0(3,3) = T01(3,3);
JA0(3,4) = T01(3,4);
JA0(3,5) = T01(3,5);
JA0(3,6) = link0.m + T01(3,6);

JA0(4,1) = link0.inertia(1,1) + T01(4,1);
JA0(4,2) = link0.inertia(1,2) + T01(4,2);
JA0(4,3) = link0.inertia(1,3) + T01(4,3);
JA0(4,4) = T01(4,4);
JA0(4,5) = -link0.mcm(3) + T01(4,5);
JA0(4,6) = link0.mcm(2) + T01(4,6);

JA0(5,1) = link0.inertia(1,2) + T01(5,1);
JA0(5,2) = link0.inertia(2,2) + T01(5,2);
JA0(5,3) = link0.inertia(2,3) + T01(5,3);
JA0(5,4) = link0.mcm(3) + T01(5,4);
JA0(5,5) = T01(5,5);
JA0(5,6) = -link0.mcm(1) + T01(5,6);

JA0(6,1) = link0.inertia(1,3) + T01(6,1);
JA0(6,2) = link0.inertia(2,3) + T01(6,2);
JA0(6,3) = link0.inertia(3,3) + T01(6,3);
JA0(6,4) = -link0.mcm(2) + T01(6,4);
JA0(6,5) = link0.mcm(1) + T01(6,5);
JA0(6,6) = T01(6,6);


%barrett_InvDynArtfunc16
      
% bias forces 
p8(1) = pv8(1);
p8(2) = pv8(2);
p8(3) = pv8(3);
p8(4) = pv8(4);
p8(5) = pv8(5);
p8(6) = pv8(6);

pmm8(1) = p8(1);
pmm8(2) = p8(2);
pmm8(3) = p8(3);
pmm8(4) = p8(4);
pmm8(5) = p8(5);
pmm8(6) = p8(6);

pm8(1) = pmm8(1)*Si78(1,1) + pmm8(2)*Si78(1,2) + pmm8(3)*Si78(1,3);
pm8(2) = pmm8(1)*Si78(2,1) + pmm8(2)*Si78(2,2) + pmm8(3)*Si78(2,3);
pm8(3) = pmm8(1)*Si78(3,1) + pmm8(2)*Si78(3,2) + pmm8(3)*Si78(3,3);
pm8(4) = pmm8(4)*Si78(1,1) + pmm8(5)*Si78(1,2) + pmm8(6)*Si78(1,3) + pmm8(1)*(-(eff(1).x(3)*Si78(2,1)) + eff(1).x(2)*Si78(3,1)) + pmm8(2)*(-(eff(1).x(3)*Si78(2,2)) + eff(1).x(2)*Si78(3,2)) + pmm8(3)*(-(eff(1).x(3)*Si78(2,3)) + eff(1).x(2)*Si78(3,3));
pm8(5) = pmm8(4)*Si78(2,1) + pmm8(5)*Si78(2,2) + pmm8(6)*Si78(2,3) + pmm8(1)*(eff(1).x(3)*Si78(1,1) - eff(1).x(1)*Si78(3,1)) + pmm8(2)*(eff(1).x(3)*Si78(1,2) - eff(1).x(1)*Si78(3,2)) + pmm8(3)*(eff(1).x(3)*Si78(1,3) - eff(1).x(1)*Si78(3,3));
pm8(6) = pmm8(1)*(-(eff(1).x(2)*Si78(1,1)) + eff(1).x(1)*Si78(2,1)) + pmm8(2)*(-(eff(1).x(2)*Si78(1,2)) + eff(1).x(1)*Si78(2,2)) + pmm8(3)*(-(eff(1).x(2)*Si78(1,3)) + eff(1).x(1)*Si78(2,3)) + pmm8(4)*Si78(3,1) + pmm8(5)*Si78(3,2) + pmm8(6)*Si78(3,3);

p7(1) = pm8(1) + pv7(1);
p7(2) = pm8(2) + pv7(2);
p7(3) = pm8(3) + pv7(3);
p7(4) = pm8(4) + pv7(4);
p7(5) = pm8(5) + pv7(5);
p7(6) = pm8(6) + pv7(6);

pmm7(1) = qdd(7)*h7(1) + p7(1) + c7(1)*JA7(1,1) + c7(2)*JA7(1,2) + c7(4)*JA7(1,4) + c7(5)*JA7(1,5);
pmm7(2) = qdd(7)*h7(2) + p7(2) + c7(1)*JA7(2,1) + c7(2)*JA7(2,2) + c7(4)*JA7(2,4) + c7(5)*JA7(2,5);
pmm7(3) = qdd(7)*h7(3) + p7(3) + c7(1)*JA7(3,1) + c7(2)*JA7(3,2) + c7(4)*JA7(3,4) + c7(5)*JA7(3,5);
pmm7(4) = qdd(7)*h7(4) + p7(4) + c7(1)*JA7(4,1) + c7(2)*JA7(4,2) + c7(4)*JA7(4,4) + c7(5)*JA7(4,5);
pmm7(5) = qdd(7)*h7(5) + p7(5) + c7(1)*JA7(5,1) + c7(2)*JA7(5,2) + c7(4)*JA7(5,4) + c7(5)*JA7(5,5);
pmm7(6) = qdd(7)*h7(6) + p7(6) + c7(1)*JA7(6,1) + c7(2)*JA7(6,2) + c7(4)*JA7(6,4) + c7(5)*JA7(6,5);

pm7(1) = pmm7(3);
pm7(2) = pmm7(1)*Si67(2,1) + pmm7(2)*Si67(2,2);
pm7(3) = pmm7(1)*Si67(3,1) + pmm7(2)*Si67(3,2);
pm7(4) = pmm7(6);
pm7(5) = pmm7(4)*Si67(2,1) + pmm7(5)*Si67(2,2);
pm7(6) = pmm7(4)*Si67(3,1) + pmm7(5)*Si67(3,2);

p6(1) = pm7(1) + pv6(1);
p6(2) = pm7(2) + pv6(2);
p6(3) = pm7(3) + pv6(3);
p6(4) = pm7(4) + pv6(4);
p6(5) = pm7(5) + pv6(5);
p6(6) = pm7(6) + pv6(6);

pmm6(1) = qdd(6)*h6(1) + p6(1) + c6(1)*JA6(1,1) + c6(2)*JA6(1,2) + c6(4)*JA6(1,4) + c6(5)*JA6(1,5);
pmm6(2) = qdd(6)*h6(2) + p6(2) + c6(1)*JA6(2,1) + c6(2)*JA6(2,2) + c6(4)*JA6(2,4) + c6(5)*JA6(2,5);
pmm6(3) = qdd(6)*h6(3) + p6(3) + c6(1)*JA6(3,1) + c6(2)*JA6(3,2) + c6(4)*JA6(3,4) + c6(5)*JA6(3,5);
pmm6(4) = qdd(6)*h6(4) + p6(4) + c6(1)*JA6(4,1) + c6(2)*JA6(4,2) + c6(4)*JA6(4,4) + c6(5)*JA6(4,5);
pmm6(5) = qdd(6)*h6(5) + p6(5) + c6(1)*JA6(5,1) + c6(2)*JA6(5,2) + c6(4)*JA6(5,4) + c6(5)*JA6(5,5);
pmm6(6) = qdd(6)*h6(6) + p6(6) + c6(1)*JA6(6,1) + c6(2)*JA6(6,2) + c6(4)*JA6(6,4) + c6(5)*JA6(6,5);

pm6(1) = -pmm6(3);
pm6(2) = pmm6(1)*Si56(2,1) + pmm6(2)*Si56(2,2);
pm6(3) = pmm6(1)*Si56(3,1) + pmm6(2)*Si56(3,2);
pm6(4) = -pmm6(6) - ZWFE*pmm6(1)*Si56(2,1) - ZWFE*pmm6(2)*Si56(2,2);
pm6(5) = -(ZWFE*pmm6(3)) + pmm6(4)*Si56(2,1) + pmm6(5)*Si56(2,2);
pm6(6) = pmm6(4)*Si56(3,1) + pmm6(5)*Si56(3,2);

p5(1) = pm6(1) + pv5(1);
p5(2) = pm6(2) + pv5(2);
p5(3) = pm6(3) + pv5(3);
p5(4) = pm6(4) + pv5(4);
p5(5) = pm6(5) + pv5(5);
p5(6) = pm6(6) + pv5(6);

pmm5(1) = qdd(5)*h5(1) + p5(1) + c5(1)*JA5(1,1) + c5(2)*JA5(1,2) + c5(4)*JA5(1,4) + c5(5)*JA5(1,5);
pmm5(2) = qdd(5)*h5(2) + p5(2) + c5(1)*JA5(2,1) + c5(2)*JA5(2,2) + c5(4)*JA5(2,4) + c5(5)*JA5(2,5);
pmm5(3) = qdd(5)*h5(3) + p5(3) + c5(1)*JA5(3,1) + c5(2)*JA5(3,2) + c5(4)*JA5(3,4) + c5(5)*JA5(3,5);
pmm5(4) = qdd(5)*h5(4) + p5(4) + c5(1)*JA5(4,1) + c5(2)*JA5(4,2) + c5(4)*JA5(4,4) + c5(5)*JA5(4,5);
pmm5(5) = qdd(5)*h5(5) + p5(5) + c5(1)*JA5(5,1) + c5(2)*JA5(5,2) + c5(4)*JA5(5,4) + c5(5)*JA5(5,5);
pmm5(6) = qdd(5)*h5(6) + p5(6) + c5(1)*JA5(6,1) + c5(2)*JA5(6,2) + c5(4)*JA5(6,4) + c5(5)*JA5(6,5);

pm5(1) = pmm5(3);
pm5(2) = pmm5(1)*Si45(2,1) + pmm5(2)*Si45(2,2);
pm5(3) = pmm5(1)*Si45(3,1) + pmm5(2)*Si45(3,2);
pm5(4) = pmm5(6) + YWR*pmm5(1)*Si45(3,1) + YWR*pmm5(2)*Si45(3,2);
pm5(5) = pmm5(4)*Si45(2,1) + pmm5(5)*Si45(2,2) - ZWR*pmm5(1)*Si45(3,1) - ZWR*pmm5(2)*Si45(3,2);
pm5(6) = -(YWR*pmm5(3)) + ZWR*pmm5(1)*Si45(2,1) + ZWR*pmm5(2)*Si45(2,2) + pmm5(4)*Si45(3,1) + pmm5(5)*Si45(3,2);

p4(1) = pm5(1) + pv4(1);
p4(2) = pm5(2) + pv4(2);
p4(3) = pm5(3) + pv4(3);
p4(4) = pm5(4) + pv4(4);
p4(5) = pm5(5) + pv4(5);
p4(6) = pm5(6) + pv4(6);

pmm4(1) = qdd(4)*h4(1) + p4(1) + c4(1)*JA4(1,1) + c4(2)*JA4(1,2) + c4(4)*JA4(1,4) + c4(5)*JA4(1,5);
pmm4(2) = qdd(4)*h4(2) + p4(2) + c4(1)*JA4(2,1) + c4(2)*JA4(2,2) + c4(4)*JA4(2,4) + c4(5)*JA4(2,5);
pmm4(3) = qdd(4)*h4(3) + p4(3) + c4(1)*JA4(3,1) + c4(2)*JA4(3,2) + c4(4)*JA4(3,4) + c4(5)*JA4(3,5);
pmm4(4) = qdd(4)*h4(4) + p4(4) + c4(1)*JA4(4,1) + c4(2)*JA4(4,2) + c4(4)*JA4(4,4) + c4(5)*JA4(4,5);
pmm4(5) = qdd(4)*h4(5) + p4(5) + c4(1)*JA4(5,1) + c4(2)*JA4(5,2) + c4(4)*JA4(5,4) + c4(5)*JA4(5,5);
pmm4(6) = qdd(4)*h4(6) + p4(6) + c4(1)*JA4(6,1) + c4(2)*JA4(6,2) + c4(4)*JA4(6,4) + c4(5)*JA4(6,5);

pm4(1) = -pmm4(3);
pm4(2) = pmm4(1)*Si34(2,1) + pmm4(2)*Si34(2,2);
pm4(3) = pmm4(1)*Si34(3,1) + pmm4(2)*Si34(3,2);
pm4(4) = -pmm4(6) + pmm4(1)*(-(ZEB*Si34(2,1)) + YEB*Si34(3,1)) + pmm4(2)*(-(ZEB*Si34(2,2)) + YEB*Si34(3,2));
pm4(5) = -(ZEB*pmm4(3)) + pmm4(4)*Si34(2,1) + pmm4(5)*Si34(2,2);
pm4(6) = YEB*pmm4(3) + pmm4(4)*Si34(3,1) + pmm4(5)*Si34(3,2);

p3(1) = pm4(1) + pv3(1);
p3(2) = pm4(2) + pv3(2);
p3(3) = pm4(3) + pv3(3);
p3(4) = pm4(4) + pv3(4);
p3(5) = pm4(5) + pv3(5);
p3(6) = pm4(6) + pv3(6);

pmm3(1) = qdd(3)*h3(1) + p3(1) + c3(1)*JA3(1,1) + c3(2)*JA3(1,2) + c3(4)*JA3(1,4) + c3(5)*JA3(1,5);
pmm3(2) = qdd(3)*h3(2) + p3(2) + c3(1)*JA3(2,1) + c3(2)*JA3(2,2) + c3(4)*JA3(2,4) + c3(5)*JA3(2,5);
pmm3(3) = qdd(3)*h3(3) + p3(3) + c3(1)*JA3(3,1) + c3(2)*JA3(3,2) + c3(4)*JA3(3,4) + c3(5)*JA3(3,5);
pmm3(4) = qdd(3)*h3(4) + p3(4) + c3(1)*JA3(4,1) + c3(2)*JA3(4,2) + c3(4)*JA3(4,4) + c3(5)*JA3(4,5);
pmm3(5) = qdd(3)*h3(5) + p3(5) + c3(1)*JA3(5,1) + c3(2)*JA3(5,2) + c3(4)*JA3(5,4) + c3(5)*JA3(5,5);
pmm3(6) = qdd(3)*h3(6) + p3(6) + c3(1)*JA3(6,1) + c3(2)*JA3(6,2) + c3(4)*JA3(6,4) + c3(5)*JA3(6,5);

pm3(1) = pmm3(3);
pm3(2) = pmm3(1)*Si23(2,1) + pmm3(2)*Si23(2,2);
pm3(3) = pmm3(1)*Si23(3,1) + pmm3(2)*Si23(3,2);
pm3(4) = pmm3(6);
pm3(5) = pmm3(4)*Si23(2,1) + pmm3(5)*Si23(2,2) - ZHR*pmm3(1)*Si23(3,1) - ZHR*pmm3(2)*Si23(3,2);
pm3(6) = ZHR*pmm3(1)*Si23(2,1) + ZHR*pmm3(2)*Si23(2,2) + pmm3(4)*Si23(3,1) + pmm3(5)*Si23(3,2);

p2(1) = pm3(1) + pv2(1);
p2(2) = pm3(2) + pv2(2);
p2(3) = pm3(3) + pv2(3);
p2(4) = pm3(4) + pv2(4);
p2(5) = pm3(5) + pv2(5);
p2(6) = pm3(6) + pv2(6);

pmm2(1) = qdd(2)*h2(1) + p2(1) + c2(1)*JA2(1,1) + c2(2)*JA2(1,2) + c2(4)*JA2(1,4) + c2(5)*JA2(1,5);
pmm2(2) = qdd(2)*h2(2) + p2(2) + c2(1)*JA2(2,1) + c2(2)*JA2(2,2) + c2(4)*JA2(2,4) + c2(5)*JA2(2,5);
pmm2(3) = qdd(2)*h2(3) + p2(3) + c2(1)*JA2(3,1) + c2(2)*JA2(3,2) + c2(4)*JA2(3,4) + c2(5)*JA2(3,5);
pmm2(4) = qdd(2)*h2(4) + p2(4) + c2(1)*JA2(4,1) + c2(2)*JA2(4,2) + c2(4)*JA2(4,4) + c2(5)*JA2(4,5);
pmm2(5) = qdd(2)*h2(5) + p2(5) + c2(1)*JA2(5,1) + c2(2)*JA2(5,2) + c2(4)*JA2(5,4) + c2(5)*JA2(5,5);
pmm2(6) = qdd(2)*h2(6) + p2(6) + c2(1)*JA2(6,1) + c2(2)*JA2(6,2) + c2(4)*JA2(6,4) + c2(5)*JA2(6,5);

pm2(1) = -pmm2(3);
pm2(2) = pmm2(1)*Si12(2,1) + pmm2(2)*Si12(2,2);
pm2(3) = pmm2(1)*Si12(3,1) + pmm2(2)*Si12(3,2);
pm2(4) = -pmm2(6);
pm2(5) = pmm2(4)*Si12(2,1) + pmm2(5)*Si12(2,2);
pm2(6) = pmm2(4)*Si12(3,1) + pmm2(5)*Si12(3,2);

p1(1) = pm2(1) + pv1(1);
p1(2) = pm2(2) + pv1(2);
p1(3) = pm2(3) + pv1(3);
p1(4) = pm2(4) + pv1(4);
p1(5) = pm2(5) + pv1(5);
p1(6) = pm2(6) + pv1(6);

pmm1(1) = qdd(1)*h1(1) + p1(1) + c1(1)*JA1(1,1) + c1(2)*JA1(1,2) + c1(4)*JA1(1,4) + c1(5)*JA1(1,5);
pmm1(2) = qdd(1)*h1(2) + p1(2) + c1(1)*JA1(2,1) + c1(2)*JA1(2,2) + c1(4)*JA1(2,4) + c1(5)*JA1(2,5);
pmm1(3) = qdd(1)*h1(3) + p1(3) + c1(1)*JA1(3,1) + c1(2)*JA1(3,2) + c1(4)*JA1(3,4) + c1(5)*JA1(3,5);
pmm1(4) = qdd(1)*h1(4) + p1(4) + c1(1)*JA1(4,1) + c1(2)*JA1(4,2) + c1(4)*JA1(4,4) + c1(5)*JA1(4,5);
pmm1(5) = qdd(1)*h1(5) + p1(5) + c1(1)*JA1(5,1) + c1(2)*JA1(5,2) + c1(4)*JA1(5,4) + c1(5)*JA1(5,5);
pmm1(6) = qdd(1)*h1(6) + p1(6) + c1(1)*JA1(6,1) + c1(2)*JA1(6,2) + c1(4)*JA1(6,4) + c1(5)*JA1(6,5);

pm1(1) = pmm1(1)*Si01(1,1) + pmm1(2)*Si01(1,2);
pm1(2) = pmm1(1)*Si01(2,1) + pmm1(2)*Si01(2,2);
pm1(3) = pmm1(3);
pm1(4) = pmm1(4)*Si01(1,1) + pmm1(5)*Si01(1,2) - ZSFE*pmm1(1)*Si01(2,1) - ZSFE*pmm1(2)*Si01(2,2);
pm1(5) = ZSFE*pmm1(1)*Si01(1,1) + ZSFE*pmm1(2)*Si01(1,2) + pmm1(4)*Si01(2,1) + pmm1(5)*Si01(2,2);
pm1(6) = pmm1(6);

p0(1) = pm1(1) + pv0(1);
p0(2) = pm1(2) + pv0(2);
p0(3) = pm1(3) + pv0(3);
p0(4) = pm1(4) + pv0(4);
p0(5) = pm1(5) + pv0(5);
p0(6) = pm1(6) + pv0(6);

%% Acceleration vectors, base acceleration, and joint torques 
a1(1) = c1(1);
a1(2) = c1(2);
a1(3) = qdd(1);
a1(4) = c1(4);
a1(5) = c1(5);

a2(1) = c2(1) + a1(2)*S21(1,2) + a1(3)*S21(1,3);
a2(2) = c2(2) + a1(2)*S21(2,2) + a1(3)*S21(2,3);
a2(3) = qdd(2) - a1(1);
a2(4) = c2(4) + a1(5)*S21(1,2);
a2(5) = c2(5) + a1(5)*S21(2,2);
a2(6) = -a1(4);

a3(1) = c3(1) + a2(2)*S32(1,2) + a2(3)*S32(1,3);
a3(2) = c3(2) + a2(2)*S32(2,2) + a2(3)*S32(2,3);
a3(3) = qdd(3) + a2(1);
a3(4) = c3(4) + ZHR*a2(3)*S32(1,2) + a2(5)*S32(1,2) - ZHR*a2(2)*S32(1,3) + a2(6)*S32(1,3);
a3(5) = c3(5) + ZHR*a2(3)*S32(2,2) + a2(5)*S32(2,2) - ZHR*a2(2)*S32(2,3) + a2(6)*S32(2,3);
a3(6) = a2(4);

a4(1) = c4(1) + a3(2)*S43(1,2) + a3(3)*S43(1,3);
a4(2) = c4(2) + a3(2)*S43(2,2) + a3(3)*S43(2,3);
a4(3) = qdd(4) - a3(1);
a4(4) = c4(4) + a3(5)*S43(1,2) + a3(6)*S43(1,3) + a3(1)*(-(ZEB*S43(1,2)) + YEB*S43(1,3));
a4(5) = c4(5) + a3(5)*S43(2,2) + a3(6)*S43(2,3) + a3(1)*(-(ZEB*S43(2,2)) + YEB*S43(2,3));
a4(6) = -(ZEB*a3(2)) + YEB*a3(3) - a3(4);

a5(1) = c5(1) + a4(2)*S54(1,2) + a4(3)*S54(1,3);
a5(2) = c5(2) + a4(2)*S54(2,2) + a4(3)*S54(2,3);
a5(3) = qdd(5) + a4(1);
a5(4) = c5(4) + ZWR*a4(3)*S54(1,2) + a4(5)*S54(1,2) + YWR*a4(1)*S54(1,3) - ZWR*a4(2)*S54(1,3) + a4(6)*S54(1,3);
a5(5) = c5(5) + ZWR*a4(3)*S54(2,2) + a4(5)*S54(2,2) + YWR*a4(1)*S54(2,3) - ZWR*a4(2)*S54(2,3) + a4(6)*S54(2,3);
a5(6) = -(YWR*a4(3)) + a4(4);

a6(1) = c6(1) + a5(2)*S65(1,2) + a5(3)*S65(1,3);
a6(2) = c6(2) + a5(2)*S65(2,2) + a5(3)*S65(2,3);
a6(3) = qdd(6) - a5(1);
a6(4) = c6(4) - ZWFE*a5(1)*S65(1,2) + a5(5)*S65(1,2) + a5(6)*S65(1,3);
a6(5) = c6(5) - ZWFE*a5(1)*S65(2,2) + a5(5)*S65(2,2) + a5(6)*S65(2,3);
a6(6) = -(ZWFE*a5(2)) - a5(4);

a7(1) = c7(1) + a6(2)*S76(1,2) + a6(3)*S76(1,3);
a7(2) = c7(2) + a6(2)*S76(2,2) + a6(3)*S76(2,3);
a7(3) = qdd(7) + a6(1);
a7(4) = c7(4) + a6(5)*S76(1,2) + a6(6)*S76(1,3);
a7(5) = c7(5) + a6(5)*S76(2,2) + a6(6)*S76(2,3);
a7(6) = a6(4);

a8(1) = a7(1)*S87(1,1) + a7(2)*S87(1,2) + a7(3)*S87(1,3);
a8(2) = a7(1)*S87(2,1) + a7(2)*S87(2,2) + a7(3)*S87(2,3);
a8(3) = a7(1)*S87(3,1) + a7(2)*S87(3,2) + a7(3)*S87(3,3);
a8(4) = a7(4)*S87(1,1) + a7(5)*S87(1,2) + a7(3)*(-(eff(1).x(2)*S87(1,1)) + eff(1).x(1)*S87(1,2)) + a7(6)*S87(1,3) + a7(2)*(eff(1).x(3)*S87(1,1) - eff(1).x(1)*S87(1,3)) + a7(1)*(-(eff(1).x(3)*S87(1,2)) + eff(1).x(2)*S87(1,3));
a8(5) = a7(4)*S87(2,1) + a7(5)*S87(2,2) + a7(3)*(-(eff(1).x(2)*S87(2,1)) + eff(1).x(1)*S87(2,2)) + a7(6)*S87(2,3) + a7(2)*(eff(1).x(3)*S87(2,1) - eff(1).x(1)*S87(2,3)) + a7(1)*(-(eff(1).x(3)*S87(2,2)) + eff(1).x(2)*S87(2,3));
a8(6) = a7(4)*S87(3,1) + a7(5)*S87(3,2) + a7(3)*(-(eff(1).x(2)*S87(3,1)) + eff(1).x(1)*S87(3,2)) + a7(6)*S87(3,3) + a7(2)*(eff(1).x(3)*S87(3,1) - eff(1).x(1)*S87(3,3)) + a7(1)*(-(eff(1).x(3)*S87(3,2)) + eff(1).x(2)*S87(3,3));

% inverse dynamics torques
u(1) = p1(6) + a1(1)*JA1(6,1) + a1(2)*JA1(6,2) + a1(3)*JA1(6,3) + a1(4)*JA1(6,4) + a1(5)*JA1(6,5);
u(2) = p2(6) + a2(1)*JA2(6,1) + a2(2)*JA2(6,2) + a2(3)*JA2(6,3) + a2(4)*JA2(6,4) + a2(5)*JA2(6,5) + a2(6)*JA2(6,6);
u(3) = p3(6) + a3(1)*JA3(6,1) + a3(2)*JA3(6,2) + a3(3)*JA3(6,3) + a3(4)*JA3(6,4) + a3(5)*JA3(6,5) + a3(6)*JA3(6,6);
u(4) = p4(6) + a4(1)*JA4(6,1) + a4(2)*JA4(6,2) + a4(3)*JA4(6,3) + a4(4)*JA4(6,4) + a4(5)*JA4(6,5) + a4(6)*JA4(6,6);
u(5) = p5(6) + a5(1)*JA5(6,1) + a5(2)*JA5(6,2) + a5(3)*JA5(6,3) + a5(4)*JA5(6,4) + a5(5)*JA5(6,5) + a5(6)*JA5(6,6);
u(6) = p6(6) + a6(1)*JA6(6,1) + a6(2)*JA6(6,2) + a6(3)*JA6(6,3) + a6(4)*JA6(6,4) + a6(5)*JA6(6,5) + a6(6)*JA6(6,6);
u(7) = p7(6) + a7(1)*JA7(6,1) + a7(2)*JA7(6,2) + a7(3)*JA7(6,3) + a7(4)*JA7(6,4) + a7(5)*JA7(6,5) + a7(6)*JA7(6,6);


