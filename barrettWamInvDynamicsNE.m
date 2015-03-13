% Barrett WAM inverse dynamics for control and linearization (ILC)

% Newton-Euler based Inverse Dynamics taken from SL: 
% shared/barrett/math/InvDynNE_declare.h
% shared/barrett/math/InvDynNE_math.h
% shared/barrett/math/InvDynNE_functions.h
%
% these are called from shared/barrett/src/SL_dynamics.c
%
%
% If flag is set to true, return the jacobians of the dynamics f

function u  =  barrettWamInvDynamicsNE(q,qd,qdd,PAR)

NDOF   =   7;

% definitions
ZSFE   =   0.346;              %!< z height of SAA axis above ground
ZHR   =   0.505;              %!< length of upper arm until 4.5cm before elbow link
YEB   =   0.045;              %!< elbow y offset
ZEB   =   0.045;              %!< elbow z offset
YWR   =  -0.045;              %!< elbow y offset (back to forewarm)
ZWR   =   0.045;              %!< elbow z offset (back to forearm)
ZWFE   =   0.255;              %!< forearm length (minus 4.5cm)

% extract parameters
link0  =  PAR.link0;
links  =  PAR.links;
eff  =  PAR.eff;
uex  =  PAR.uex;
uex0  =  PAR.uex0;
basec  =  PAR.basec;
baseo  =  PAR.baseo;

% warning! uex_des is what desired state structure in SL carries
uex_des = zeros(1,7);

g  =  9.81;

% sine and cosine precomputation 
ss1th   =   sin(q(1));
cs1th   =   cos(q(1));
ss2th   =   sin(q(2));
cs2th   =   cos(q(2));
ss3th   =   sin(q(3));
cs3th   =   cos(q(3));
ss4th   =   sin(q(4));
cs4th   =   cos(q(4));
ss5th   =   sin(q(5));
cs5th   =   cos(q(5));
ss6th   =   sin(q(6));
cs6th   =   cos(q(6));
ss7th   =   sin(q(7));
cs7th   =   cos(q(7));

% endeffector orientations

rseff1a1   =   sin(eff(1).a(1));
rceff1a1   =   cos(eff(1).a(1));

rseff1a2   =   sin(eff(1).a(2));
rceff1a2   =   cos(eff(1).a(2));

rseff1a3   =   sin(eff(1).a(3));
rceff1a3   =   cos(eff(1).a(3));

%% Includes the functions one by one

%barrett_InvDynNEfunc1(void)
     
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

%barrett_InvDynNEfunc2(void)
     
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


%barrett_InvDynNEfunc3(void)
     
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

%barrett_InvDynNEfunc4(void)
     
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
v6(3) = qd(6) - v5(1);
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

%barrett_InvDynNEfunc5(void)
     
% acceleration vectors 
a0(1) = baseo.add(1)*S00(1,1) + baseo.add(2)*S00(1,2) + baseo.add(3)*S00(1,3);
a0(2) = baseo.add(1)*S00(2,1) + baseo.add(2)*S00(2,2) + baseo.add(3)*S00(2,3);
a0(3) = baseo.add(1)*S00(3,1) + baseo.add(2)*S00(3,2) + baseo.add(3)*S00(3,3);
a0(4) = basec.xdd(1)*S00(1,1) + basec.xdd(2)*S00(1,2) + (gravity + basec.xdd(3))*S00(1,3);
a0(5) = basec.xdd(1)*S00(2,1) + basec.xdd(2)*S00(2,2) + (gravity + basec.xdd(3))*S00(2,3);
a0(6) = basec.xdd(1)*S00(3,1) + basec.xdd(2)*S00(3,2) + (gravity + basec.xdd(3))*S00(3,3);

a1(1) = qd(1)*v1(2) + a0(1)*S10(1,1) + a0(2)*S10(1,2);
a1(2) = -(qd(1)*v1(1)) + a0(1)*S10(2,1) + a0(2)*S10(2,2);
a1(3) = qdd(1) + a0(3);
a1(4) = qd(1)*v1(5) + ZSFE*a0(2)*S10(1,1) + a0(4)*S10(1,1) - ZSFE*a0(1)*S10(1,2) + a0(5)*S10(1,2);
a1(5) = -(qd(1)*v1(4)) + ZSFE*a0(2)*S10(2,1) + a0(4)*S10(2,1) - ZSFE*a0(1)*S10(2,2) + a0(5)*S10(2,2);
a1(6) = a0(6);

a2(1) = qd(2)*v2(2) + a1(2)*S21(1,2) + a1(3)*S21(1,3);
a2(2) = -(qd(2)*v2(1)) + a1(2)*S21(2,2) + a1(3)*S21(2,3);
a2(3) = qdd(2) - a1(1);
a2(4) = qd(2)*v2(5) + a1(5)*S21(1,2) + a1(6)*S21(1,3);
a2(5) = -(qd(2)*v2(4)) + a1(5)*S21(2,2) + a1(6)*S21(2,3);
a2(6) = -a1(4);

a3(1) = qd(3)*v3(2) + a2(2)*S32(1,2) + a2(3)*S32(1,3);
a3(2) = -(qd(3)*v3(1)) + a2(2)*S32(2,2) + a2(3)*S32(2,3);
a3(3) = qdd(3) + a2(1);
a3(4) = qd(3)*v3(5) + ZHR*a2(3)*S32(1,2) + a2(5)*S32(1,2) - ZHR*a2(2)*S32(1,3) + a2(6)*S32(1,3);
a3(5) = -(qd(3)*v3(4)) + ZHR*a2(3)*S32(2,2) + a2(5)*S32(2,2) - ZHR*a2(2)*S32(2,3) + a2(6)*S32(2,3);
a3(6) = a2(4);

a4(1) = qd(4)*v4(2) + a3(2)*S43(1,2) + a3(3)*S43(1,3);
a4(2) = -(qd(4)*v4(1)) + a3(2)*S43(2,2) + a3(3)*S43(2,3);
a4(3) = qdd(4) - a3(1);
a4(4) = qd(4)*v4(5) + a3(5)*S43(1,2) + a3(6)*S43(1,3) + a3(1)*(-(ZEB*S43(1,2)) + YEB*S43(1,3));
a4(5) = -(qd(4)*v4(4)) + a3(5)*S43(2,2) + a3(6)*S43(2,3) + a3(1)*(-(ZEB*S43(2,2)) + YEB*S43(2,3));
a4(6) = -(ZEB*a3(2)) + YEB*a3(3) - a3(4);

a5(1) = qd(5)*v5(2) + a4(2)*S54(1,2) + a4(3)*S54(1,3);
a5(2) = -(qd(5)*v5(1)) + a4(2)*S54(2,2) + a4(3)*S54(2,3);
a5(3) = qdd(5) + a4(1);
a5(4) = qd(5)*v5(5) + ZWR*a4(3)*S54(1,2) + a4(5)*S54(1,2) + YWR*a4(1)*S54(1,3) - ZWR*a4(2)*S54(1,3) + a4(6)*S54(1,3);
a5(5) = -(qd(5)*v5(4)) + ZWR*a4(3)*S54(2,2) + a4(5)*S54(2,2) + YWR*a4(1)*S54(2,3) - ZWR*a4(2)*S54(2,3) + a4(6)*S54(2,3);
a5(6) = -(YWR*a4(3)) + a4(4);

a6(1) = qd(6)*v6(2) + a5(2)*S65(1,2) + a5(3)*S65(1,3);
a6(2) = -(qd(6)*v6(1)) + a5(2)*S65(2,2) + a5(3)*S65(2,3);
a6(3) = qdd(6) - a5(1);
a6(4) = qd(6)*v6(5) - ZWFE*a5(1)*S65(1,2) + a5(5)*S65(1,2) + a5(6)*S65(1,3);
a6(5) = -(qd(6)*v6(4)) - ZWFE*a5(1)*S65(2,2) + a5(5)*S65(2,2) + a5(6)*S65(2,3);
a6(6) = -(ZWFE*a5(2)) - a5(4);

a7(1) = qd(7)*v7(2) + a6(2)*S76(1,2) + a6(3)*S76(1,3);
a7(2) = -(qd(7)*v7(1)) + a6(2)*S76(2,2) + a6(3)*S76(2,3);
a7(3) = qdd(7) + a6(1);
a7(4) = qd(7)*v7(5) + a6(5)*S76(1,2) + a6(6)*S76(1,3);
a7(5) = -(qd(7)*v7(4)) + a6(5)*S76(2,2) + a6(6)*S76(2,3);
a7(6) = a6(4);

a8(1) = a7(1)*S87(1,1) + a7(2)*S87(1,2) + a7(3)*S87(1,3);
a8(2) = a7(1)*S87(2,1) + a7(2)*S87(2,2) + a7(3)*S87(2,3);
a8(3) = a7(1)*S87(3,1) + a7(2)*S87(3,2) + a7(3)*S87(3,3);
a8(4) = a7(4)*S87(1,1) + a7(5)*S87(1,2) + a7(3)*(-(eff(1).x(2)*S87(1,1)) + eff(1).x(1)*S87(1,2)) + a7(6)*S87(1,3) + a7(2)*(eff(1).x(3)*S87(1,1) - eff(1).x(1)*S87(1,3)) + a7(1)*(-(eff(1).x(3)*S87(1,2)) + eff(1).x(2)*S87(1,3));
a8(5) = a7(4)*S87(2,1) + a7(5)*S87(2,2) + a7(3)*(-(eff(1).x(2)*S87(2,1)) + eff(1).x(1)*S87(2,2)) + a7(6)*S87(2,3) + a7(2)*(eff(1).x(3)*S87(2,1) - eff(1).x(1)*S87(2,3)) + a7(1)*(-(eff(1).x(3)*S87(2,2)) + eff(1).x(2)*S87(2,3));
a8(6) = a7(4)*S87(3,1) + a7(5)*S87(3,2) + a7(3)*(-(eff(1).x(2)*S87(3,1)) + eff(1).x(1)*S87(3,2)) + a7(6)*S87(3,3) + a7(2)*(eff(1).x(3)*S87(3,1) - eff(1).x(1)*S87(3,3)) + a7(1)*(-(eff(1).x(3)*S87(3,2)) + eff(1).x(2)*S87(3,3));

%barrett_InvDynNEfunc6(void)
     
% net forces and external forces for each joint 
fnet0(1) = link0.m*a0(4) - a0(3)*link0.mcm(2) + a0(2)*link0.mcm(3) - link0.mcm(1)*power(v0(2),2) - link0.mcm(1)*power(v0(3),2) + v0(1)*(link0.mcm(2)*v0(2) + link0.mcm(3)*v0(3)) - link0.m*v0(3)*v0(5) + link0.m*v0(2)*v0(6);
fnet0(2) = link0.m*a0(5) + a0(3)*link0.mcm(1) - a0(1)*link0.mcm(3) - link0.mcm(2)*power(v0(1),2) - link0.mcm(2)*power(v0(3),2) + v0(2)*(link0.mcm(1)*v0(1) + link0.mcm(3)*v0(3)) + link0.m*v0(3)*v0(4) - link0.m*v0(1)*v0(6);
fnet0(3) = link0.m*a0(6) - a0(2)*link0.mcm(1) + a0(1)*link0.mcm(2) - link0.mcm(3)*power(v0(1),2) - link0.mcm(3)*power(v0(2),2) + (link0.mcm(1)*v0(1) + link0.mcm(2)*v0(2))*v0(3) - link0.m*v0(2)*v0(4) + link0.m*v0(1)*v0(5);
fnet0(4) = a0(6)*link0.mcm(2) - a0(5)*link0.mcm(3) + (-(link0.mcm(2)*v0(2)) - link0.mcm(3)*v0(3))*v0(4) + (link0.mcm(1)*v0(3) + link0.m*v0(5))*v0(6) + v0(5)*(link0.mcm(1)*v0(2) - link0.m*v0(6)) + a0(1)*link0.inertia(1,1) + a0(2)*link0.inertia(1,2) + a0(3)*link0.inertia(1,3) + v0(1)*(link0.mcm(2)*v0(5) + link0.mcm(3)*v0(6) - v0(3)*link0.inertia(1,2) + v0(2)*link0.inertia(1,3)) + v0(2)*(-(link0.mcm(1)*v0(5)) - v0(3)*link0.inertia(2,2) + v0(2)*link0.inertia(2,3)) + v0(3)*(-(link0.mcm(1)*v0(6)) - v0(3)*link0.inertia(2,3) + v0(2)*link0.inertia(3,3));
fnet0(5) = -(a0(6)*link0.mcm(1)) + a0(4)*link0.mcm(3) + (-(link0.mcm(1)*v0(1)) - link0.mcm(3)*v0(3))*v0(5) + (link0.mcm(2)*v0(3) - link0.m*v0(4))*v0(6) + v0(4)*(link0.mcm(2)*v0(1) + link0.m*v0(6)) + a0(1)*link0.inertia(1,2) + v0(1)*(-(link0.mcm(2)*v0(4)) + v0(3)*link0.inertia(1,1) - v0(1)*link0.inertia(1,3)) + a0(2)*link0.inertia(2,2) + a0(3)*link0.inertia(2,3) + v0(2)*(link0.mcm(1)*v0(4) + link0.mcm(3)*v0(6) + v0(3)*link0.inertia(1,2) - v0(1)*link0.inertia(2,3)) + v0(3)*(-(link0.mcm(2)*v0(6)) + v0(3)*link0.inertia(1,3) - v0(1)*link0.inertia(3,3));
fnet0(6) = a0(5)*link0.mcm(1) - a0(4)*link0.mcm(2) + (link0.mcm(3)*v0(2) + link0.m*v0(4))*v0(5) + v0(4)*(link0.mcm(3)*v0(1) - link0.m*v0(5)) + (-(link0.mcm(1)*v0(1)) - link0.mcm(2)*v0(2))*v0(6) + v0(1)*(-(link0.mcm(3)*v0(4)) - v0(2)*link0.inertia(1,1) + v0(1)*link0.inertia(1,2)) + a0(1)*link0.inertia(1,3) + v0(2)*(-(link0.mcm(3)*v0(5)) - v0(2)*link0.inertia(1,2) + v0(1)*link0.inertia(2,2)) + a0(2)*link0.inertia(2,3) + v0(3)*(link0.mcm(1)*v0(4) + link0.mcm(2)*v0(5) - v0(2)*link0.inertia(1,3) + v0(1)*link0.inertia(2,3)) + a0(3)*link0.inertia(3,3);

fnet1(1) = links(1).m*a1(4) - a1(3)*links(1).mcm(2) + a1(2)*links(1).mcm(3) - links(1).mcm(1)*power(v1(2),2) - links(1).mcm(1)*power(v1(3),2) + v1(1)*(links(1).mcm(2)*v1(2) + links(1).mcm(3)*v1(3)) - links(1).m*v1(3)*v1(5) + links(1).m*v1(2)*v1(6);
fnet1(2) = links(1).m*a1(5) + a1(3)*links(1).mcm(1) - a1(1)*links(1).mcm(3) - links(1).mcm(2)*power(v1(1),2) - links(1).mcm(2)*power(v1(3),2) + v1(2)*(links(1).mcm(1)*v1(1) + links(1).mcm(3)*v1(3)) + links(1).m*v1(3)*v1(4) - links(1).m*v1(1)*v1(6);
fnet1(3) = links(1).m*a1(6) - a1(2)*links(1).mcm(1) + a1(1)*links(1).mcm(2) - links(1).mcm(3)*power(v1(1),2) - links(1).mcm(3)*power(v1(2),2) + (links(1).mcm(1)*v1(1) + links(1).mcm(2)*v1(2))*v1(3) - links(1).m*v1(2)*v1(4) + links(1).m*v1(1)*v1(5);
fnet1(4) = a1(6)*links(1).mcm(2) - a1(5)*links(1).mcm(3) + (-(links(1).mcm(2)*v1(2)) - links(1).mcm(3)*v1(3))*v1(4) + (links(1).mcm(1)*v1(3) + links(1).m*v1(5))*v1(6) + v1(5)*(links(1).mcm(1)*v1(2) - links(1).m*v1(6)) + a1(1)*links(1).inertia(1,1) + a1(2)*links(1).inertia(1,2) + a1(3)*links(1).inertia(1,3) + v1(1)*(links(1).mcm(2)*v1(5) + links(1).mcm(3)*v1(6) - v1(3)*links(1).inertia(1,2) + v1(2)*links(1).inertia(1,3)) + v1(2)*(-(links(1).mcm(1)*v1(5)) - v1(3)*links(1).inertia(2,2) + v1(2)*links(1).inertia(2,3)) + v1(3)*(-(links(1).mcm(1)*v1(6)) - v1(3)*links(1).inertia(2,3) + v1(2)*links(1).inertia(3,3));
fnet1(5) = -(a1(6)*links(1).mcm(1)) + a1(4)*links(1).mcm(3) + (-(links(1).mcm(1)*v1(1)) - links(1).mcm(3)*v1(3))*v1(5) + (links(1).mcm(2)*v1(3) - links(1).m*v1(4))*v1(6) + v1(4)*(links(1).mcm(2)*v1(1) + links(1).m*v1(6)) + a1(1)*links(1).inertia(1,2) + v1(1)*(-(links(1).mcm(2)*v1(4)) + v1(3)*links(1).inertia(1,1) - v1(1)*links(1).inertia(1,3)) + a1(2)*links(1).inertia(2,2) + a1(3)*links(1).inertia(2,3) + v1(2)*(links(1).mcm(1)*v1(4) + links(1).mcm(3)*v1(6) + v1(3)*links(1).inertia(1,2) - v1(1)*links(1).inertia(2,3)) + v1(3)*(-(links(1).mcm(2)*v1(6)) + v1(3)*links(1).inertia(1,3) - v1(1)*links(1).inertia(3,3));
fnet1(6) = a1(5)*links(1).mcm(1) - a1(4)*links(1).mcm(2) + (links(1).mcm(3)*v1(2) + links(1).m*v1(4))*v1(5) + v1(4)*(links(1).mcm(3)*v1(1) - links(1).m*v1(5)) + (-(links(1).mcm(1)*v1(1)) - links(1).mcm(2)*v1(2))*v1(6) + v1(1)*(-(links(1).mcm(3)*v1(4)) - v1(2)*links(1).inertia(1,1) + v1(1)*links(1).inertia(1,2)) + a1(1)*links(1).inertia(1,3) + v1(2)*(-(links(1).mcm(3)*v1(5)) - v1(2)*links(1).inertia(1,2) + v1(1)*links(1).inertia(2,2)) + a1(2)*links(1).inertia(2,3) + v1(3)*(links(1).mcm(1)*v1(4) + links(1).mcm(2)*v1(5) - v1(2)*links(1).inertia(1,3) + v1(1)*links(1).inertia(2,3)) + a1(3)*links(1).inertia(3,3);

fnet2(1) = links(2).m*a2(4) - a2(3)*links(2).mcm(2) + a2(2)*links(2).mcm(3) - links(2).mcm(1)*power(v2(2),2) - links(2).mcm(1)*power(v2(3),2) + v2(1)*(links(2).mcm(2)*v2(2) + links(2).mcm(3)*v2(3)) - links(2).m*v2(3)*v2(5) + links(2).m*v2(2)*v2(6);
fnet2(2) = links(2).m*a2(5) + a2(3)*links(2).mcm(1) - a2(1)*links(2).mcm(3) - links(2).mcm(2)*power(v2(1),2) - links(2).mcm(2)*power(v2(3),2) + v2(2)*(links(2).mcm(1)*v2(1) + links(2).mcm(3)*v2(3)) + links(2).m*v2(3)*v2(4) - links(2).m*v2(1)*v2(6);
fnet2(3) = links(2).m*a2(6) - a2(2)*links(2).mcm(1) + a2(1)*links(2).mcm(2) - links(2).mcm(3)*power(v2(1),2) - links(2).mcm(3)*power(v2(2),2) + (links(2).mcm(1)*v2(1) + links(2).mcm(2)*v2(2))*v2(3) - links(2).m*v2(2)*v2(4) + links(2).m*v2(1)*v2(5);
fnet2(4) = a2(6)*links(2).mcm(2) - a2(5)*links(2).mcm(3) + (-(links(2).mcm(2)*v2(2)) - links(2).mcm(3)*v2(3))*v2(4) + (links(2).mcm(1)*v2(3) + links(2).m*v2(5))*v2(6) + v2(5)*(links(2).mcm(1)*v2(2) - links(2).m*v2(6)) + a2(1)*links(2).inertia(1,1) + a2(2)*links(2).inertia(1,2) + a2(3)*links(2).inertia(1,3) + v2(1)*(links(2).mcm(2)*v2(5) + links(2).mcm(3)*v2(6) - v2(3)*links(2).inertia(1,2) + v2(2)*links(2).inertia(1,3)) + v2(2)*(-(links(2).mcm(1)*v2(5)) - v2(3)*links(2).inertia(2,2) + v2(2)*links(2).inertia(2,3)) + v2(3)*(-(links(2).mcm(1)*v2(6)) - v2(3)*links(2).inertia(2,3) + v2(2)*links(2).inertia(3,3));
fnet2(5) = -(a2(6)*links(2).mcm(1)) + a2(4)*links(2).mcm(3) + (-(links(2).mcm(1)*v2(1)) - links(2).mcm(3)*v2(3))*v2(5) + (links(2).mcm(2)*v2(3) - links(2).m*v2(4))*v2(6) + v2(4)*(links(2).mcm(2)*v2(1) + links(2).m*v2(6)) + a2(1)*links(2).inertia(1,2) + v2(1)*(-(links(2).mcm(2)*v2(4)) + v2(3)*links(2).inertia(1,1) - v2(1)*links(2).inertia(1,3)) + a2(2)*links(2).inertia(2,2) + a2(3)*links(2).inertia(2,3) + v2(2)*(links(2).mcm(1)*v2(4) + links(2).mcm(3)*v2(6) + v2(3)*links(2).inertia(1,2) - v2(1)*links(2).inertia(2,3)) + v2(3)*(-(links(2).mcm(2)*v2(6)) + v2(3)*links(2).inertia(1,3) - v2(1)*links(2).inertia(3,3));
fnet2(6) = a2(5)*links(2).mcm(1) - a2(4)*links(2).mcm(2) + (links(2).mcm(3)*v2(2) + links(2).m*v2(4))*v2(5) + v2(4)*(links(2).mcm(3)*v2(1) - links(2).m*v2(5)) + (-(links(2).mcm(1)*v2(1)) - links(2).mcm(2)*v2(2))*v2(6) + v2(1)*(-(links(2).mcm(3)*v2(4)) - v2(2)*links(2).inertia(1,1) + v2(1)*links(2).inertia(1,2)) + a2(1)*links(2).inertia(1,3) + v2(2)*(-(links(2).mcm(3)*v2(5)) - v2(2)*links(2).inertia(1,2) + v2(1)*links(2).inertia(2,2)) + a2(2)*links(2).inertia(2,3) + v2(3)*(links(2).mcm(1)*v2(4) + links(2).mcm(2)*v2(5) - v2(2)*links(2).inertia(1,3) + v2(1)*links(2).inertia(2,3)) + a2(3)*links(2).inertia(3,3);

fnet3(1) = links(3).m*a3(4) - a3(3)*links(3).mcm(2) + a3(2)*links(3).mcm(3) - links(3).mcm(1)*power(v3(2),2) - links(3).mcm(1)*power(v3(3),2) + v3(1)*(links(3).mcm(2)*v3(2) + links(3).mcm(3)*v3(3)) - links(3).m*v3(3)*v3(5) + links(3).m*v3(2)*v3(6);
fnet3(2) = links(3).m*a3(5) + a3(3)*links(3).mcm(1) - a3(1)*links(3).mcm(3) - links(3).mcm(2)*power(v3(1),2) - links(3).mcm(2)*power(v3(3),2) + v3(2)*(links(3).mcm(1)*v3(1) + links(3).mcm(3)*v3(3)) + links(3).m*v3(3)*v3(4) - links(3).m*v3(1)*v3(6);
fnet3(3) = links(3).m*a3(6) - a3(2)*links(3).mcm(1) + a3(1)*links(3).mcm(2) - links(3).mcm(3)*power(v3(1),2) - links(3).mcm(3)*power(v3(2),2) + (links(3).mcm(1)*v3(1) + links(3).mcm(2)*v3(2))*v3(3) - links(3).m*v3(2)*v3(4) + links(3).m*v3(1)*v3(5);
fnet3(4) = a3(6)*links(3).mcm(2) - a3(5)*links(3).mcm(3) + (-(links(3).mcm(2)*v3(2)) - links(3).mcm(3)*v3(3))*v3(4) + (links(3).mcm(1)*v3(3) + links(3).m*v3(5))*v3(6) + v3(5)*(links(3).mcm(1)*v3(2) - links(3).m*v3(6)) + a3(1)*links(3).inertia(1,1) + a3(2)*links(3).inertia(1,2) + a3(3)*links(3).inertia(1,3) + v3(1)*(links(3).mcm(2)*v3(5) + links(3).mcm(3)*v3(6) - v3(3)*links(3).inertia(1,2) + v3(2)*links(3).inertia(1,3)) + v3(2)*(-(links(3).mcm(1)*v3(5)) - v3(3)*links(3).inertia(2,2) + v3(2)*links(3).inertia(2,3)) + v3(3)*(-(links(3).mcm(1)*v3(6)) - v3(3)*links(3).inertia(2,3) + v3(2)*links(3).inertia(3,3));
fnet3(5) = -(a3(6)*links(3).mcm(1)) + a3(4)*links(3).mcm(3) + (-(links(3).mcm(1)*v3(1)) - links(3).mcm(3)*v3(3))*v3(5) + (links(3).mcm(2)*v3(3) - links(3).m*v3(4))*v3(6) + v3(4)*(links(3).mcm(2)*v3(1) + links(3).m*v3(6)) + a3(1)*links(3).inertia(1,2) + v3(1)*(-(links(3).mcm(2)*v3(4)) + v3(3)*links(3).inertia(1,1) - v3(1)*links(3).inertia(1,3)) + a3(2)*links(3).inertia(2,2) + a3(3)*links(3).inertia(2,3) + v3(2)*(links(3).mcm(1)*v3(4) + links(3).mcm(3)*v3(6) + v3(3)*links(3).inertia(1,2) - v3(1)*links(3).inertia(2,3)) + v3(3)*(-(links(3).mcm(2)*v3(6)) + v3(3)*links(3).inertia(1,3) - v3(1)*links(3).inertia(3,3));
fnet3(6) = a3(5)*links(3).mcm(1) - a3(4)*links(3).mcm(2) + (links(3).mcm(3)*v3(2) + links(3).m*v3(4))*v3(5) + v3(4)*(links(3).mcm(3)*v3(1) - links(3).m*v3(5)) + (-(links(3).mcm(1)*v3(1)) - links(3).mcm(2)*v3(2))*v3(6) + v3(1)*(-(links(3).mcm(3)*v3(4)) - v3(2)*links(3).inertia(1,1) + v3(1)*links(3).inertia(1,2)) + a3(1)*links(3).inertia(1,3) + v3(2)*(-(links(3).mcm(3)*v3(5)) - v3(2)*links(3).inertia(1,2) + v3(1)*links(3).inertia(2,2)) + a3(2)*links(3).inertia(2,3) + v3(3)*(links(3).mcm(1)*v3(4) + links(3).mcm(2)*v3(5) - v3(2)*links(3).inertia(1,3) + v3(1)*links(3).inertia(2,3)) + a3(3)*links(3).inertia(3,3);

fnet4(1) = links(4).m*a4(4) - a4(3)*links(4).mcm(2) + a4(2)*links(4).mcm(3) - links(4).mcm(1)*power(v4(2),2) - links(4).mcm(1)*power(v4(3),2) + v4(1)*(links(4).mcm(2)*v4(2) + links(4).mcm(3)*v4(3)) - links(4).m*v4(3)*v4(5) + links(4).m*v4(2)*v4(6);
fnet4(2) = links(4).m*a4(5) + a4(3)*links(4).mcm(1) - a4(1)*links(4).mcm(3) - links(4).mcm(2)*power(v4(1),2) - links(4).mcm(2)*power(v4(3),2) + v4(2)*(links(4).mcm(1)*v4(1) + links(4).mcm(3)*v4(3)) + links(4).m*v4(3)*v4(4) - links(4).m*v4(1)*v4(6);
fnet4(3) = links(4).m*a4(6) - a4(2)*links(4).mcm(1) + a4(1)*links(4).mcm(2) - links(4).mcm(3)*power(v4(1),2) - links(4).mcm(3)*power(v4(2),2) + (links(4).mcm(1)*v4(1) + links(4).mcm(2)*v4(2))*v4(3) - links(4).m*v4(2)*v4(4) + links(4).m*v4(1)*v4(5);
fnet4(4) = a4(6)*links(4).mcm(2) - a4(5)*links(4).mcm(3) + (-(links(4).mcm(2)*v4(2)) - links(4).mcm(3)*v4(3))*v4(4) + (links(4).mcm(1)*v4(3) + links(4).m*v4(5))*v4(6) + v4(5)*(links(4).mcm(1)*v4(2) - links(4).m*v4(6)) + a4(1)*links(4).inertia(1,1) + a4(2)*links(4).inertia(1,2) + a4(3)*links(4).inertia(1,3) + v4(1)*(links(4).mcm(2)*v4(5) + links(4).mcm(3)*v4(6) - v4(3)*links(4).inertia(1,2) + v4(2)*links(4).inertia(1,3)) + v4(2)*(-(links(4).mcm(1)*v4(5)) - v4(3)*links(4).inertia(2,2) + v4(2)*links(4).inertia(2,3)) + v4(3)*(-(links(4).mcm(1)*v4(6)) - v4(3)*links(4).inertia(2,3) + v4(2)*links(4).inertia(3,3));
fnet4(5) = -(a4(6)*links(4).mcm(1)) + a4(4)*links(4).mcm(3) + (-(links(4).mcm(1)*v4(1)) - links(4).mcm(3)*v4(3))*v4(5) + (links(4).mcm(2)*v4(3) - links(4).m*v4(4))*v4(6) + v4(4)*(links(4).mcm(2)*v4(1) + links(4).m*v4(6)) + a4(1)*links(4).inertia(1,2) + v4(1)*(-(links(4).mcm(2)*v4(4)) + v4(3)*links(4).inertia(1,1) - v4(1)*links(4).inertia(1,3)) + a4(2)*links(4).inertia(2,2) + a4(3)*links(4).inertia(2,3) + v4(2)*(links(4).mcm(1)*v4(4) + links(4).mcm(3)*v4(6) + v4(3)*links(4).inertia(1,2) - v4(1)*links(4).inertia(2,3)) + v4(3)*(-(links(4).mcm(2)*v4(6)) + v4(3)*links(4).inertia(1,3) - v4(1)*links(4).inertia(3,3));
fnet4(6) = a4(5)*links(4).mcm(1) - a4(4)*links(4).mcm(2) + (links(4).mcm(3)*v4(2) + links(4).m*v4(4))*v4(5) + v4(4)*(links(4).mcm(3)*v4(1) - links(4).m*v4(5)) + (-(links(4).mcm(1)*v4(1)) - links(4).mcm(2)*v4(2))*v4(6) + v4(1)*(-(links(4).mcm(3)*v4(4)) - v4(2)*links(4).inertia(1,1) + v4(1)*links(4).inertia(1,2)) + a4(1)*links(4).inertia(1,3) + v4(2)*(-(links(4).mcm(3)*v4(5)) - v4(2)*links(4).inertia(1,2) + v4(1)*links(4).inertia(2,2)) + a4(2)*links(4).inertia(2,3) + v4(3)*(links(4).mcm(1)*v4(4) + links(4).mcm(2)*v4(5) - v4(2)*links(4).inertia(1,3) + v4(1)*links(4).inertia(2,3)) + a4(3)*links(4).inertia(3,3);

fnet5(1) = links(5).m*a5(4) - a5(3)*links(5).mcm(2) + a5(2)*links(5).mcm(3) - links(5).mcm(1)*power(v5(2),2) - links(5).mcm(1)*power(v5(3),2) + v5(1)*(links(5).mcm(2)*v5(2) + links(5).mcm(3)*v5(3)) - links(5).m*v5(3)*v5(5) + links(5).m*v5(2)*v5(6);
fnet5(2) = links(5).m*a5(5) + a5(3)*links(5).mcm(1) - a5(1)*links(5).mcm(3) - links(5).mcm(2)*power(v5(1),2) - links(5).mcm(2)*power(v5(3),2) + v5(2)*(links(5).mcm(1)*v5(1) + links(5).mcm(3)*v5(3)) + links(5).m*v5(3)*v5(4) - links(5).m*v5(1)*v5(6);
fnet5(3) = links(5).m*a5(6) - a5(2)*links(5).mcm(1) + a5(1)*links(5).mcm(2) - links(5).mcm(3)*power(v5(1),2) - links(5).mcm(3)*power(v5(2),2) + (links(5).mcm(1)*v5(1) + links(5).mcm(2)*v5(2))*v5(3) - links(5).m*v5(2)*v5(4) + links(5).m*v5(1)*v5(5);
fnet5(4) = a5(6)*links(5).mcm(2) - a5(5)*links(5).mcm(3) + (-(links(5).mcm(2)*v5(2)) - links(5).mcm(3)*v5(3))*v5(4) + (links(5).mcm(1)*v5(3) + links(5).m*v5(5))*v5(6) + v5(5)*(links(5).mcm(1)*v5(2) - links(5).m*v5(6)) + a5(1)*links(5).inertia(1,1) + a5(2)*links(5).inertia(1,2) + a5(3)*links(5).inertia(1,3) + v5(1)*(links(5).mcm(2)*v5(5) + links(5).mcm(3)*v5(6) - v5(3)*links(5).inertia(1,2) + v5(2)*links(5).inertia(1,3)) + v5(2)*(-(links(5).mcm(1)*v5(5)) - v5(3)*links(5).inertia(2,2) + v5(2)*links(5).inertia(2,3)) + v5(3)*(-(links(5).mcm(1)*v5(6)) - v5(3)*links(5).inertia(2,3) + v5(2)*links(5).inertia(3,3));
fnet5(5) = -(a5(6)*links(5).mcm(1)) + a5(4)*links(5).mcm(3) + (-(links(5).mcm(1)*v5(1)) - links(5).mcm(3)*v5(3))*v5(5) + (links(5).mcm(2)*v5(3) - links(5).m*v5(4))*v5(6) + v5(4)*(links(5).mcm(2)*v5(1) + links(5).m*v5(6)) + a5(1)*links(5).inertia(1,2) + v5(1)*(-(links(5).mcm(2)*v5(4)) + v5(3)*links(5).inertia(1,1) - v5(1)*links(5).inertia(1,3)) + a5(2)*links(5).inertia(2,2) + a5(3)*links(5).inertia(2,3) + v5(2)*(links(5).mcm(1)*v5(4) + links(5).mcm(3)*v5(6) + v5(3)*links(5).inertia(1,2) - v5(1)*links(5).inertia(2,3)) + v5(3)*(-(links(5).mcm(2)*v5(6)) + v5(3)*links(5).inertia(1,3) - v5(1)*links(5).inertia(3,3));
fnet5(6) = a5(5)*links(5).mcm(1) - a5(4)*links(5).mcm(2) + (links(5).mcm(3)*v5(2) + links(5).m*v5(4))*v5(5) + v5(4)*(links(5).mcm(3)*v5(1) - links(5).m*v5(5)) + (-(links(5).mcm(1)*v5(1)) - links(5).mcm(2)*v5(2))*v5(6) + v5(1)*(-(links(5).mcm(3)*v5(4)) - v5(2)*links(5).inertia(1,1) + v5(1)*links(5).inertia(1,2)) + a5(1)*links(5).inertia(1,3) + v5(2)*(-(links(5).mcm(3)*v5(5)) - v5(2)*links(5).inertia(1,2) + v5(1)*links(5).inertia(2,2)) + a5(2)*links(5).inertia(2,3) + v5(3)*(links(5).mcm(1)*v5(4) + links(5).mcm(2)*v5(5) - v5(2)*links(5).inertia(1,3) + v5(1)*links(5).inertia(2,3)) + a5(3)*links(5).inertia(3,3);

fnet6(1) = links(6).m*a6(4) - a6(3)*links(6).mcm(2) + a6(2)*links(6).mcm(3) - links(6).mcm(1)*power(v6(2),2) - links(6).mcm(1)*power(v6(3),2) + v6(1)*(links(6).mcm(2)*v6(2) + links(6).mcm(3)*v6(3)) - links(6).m*v6(3)*v6(5) + links(6).m*v6(2)*v6(6);
fnet6(2) = links(6).m*a6(5) + a6(3)*links(6).mcm(1) - a6(1)*links(6).mcm(3) - links(6).mcm(2)*power(v6(1),2) - links(6).mcm(2)*power(v6(3),2) + v6(2)*(links(6).mcm(1)*v6(1) + links(6).mcm(3)*v6(3)) + links(6).m*v6(3)*v6(4) - links(6).m*v6(1)*v6(6);
fnet6(3) = links(6).m*a6(6) - a6(2)*links(6).mcm(1) + a6(1)*links(6).mcm(2) - links(6).mcm(3)*power(v6(1),2) - links(6).mcm(3)*power(v6(2),2) + (links(6).mcm(1)*v6(1) + links(6).mcm(2)*v6(2))*v6(3) - links(6).m*v6(2)*v6(4) + links(6).m*v6(1)*v6(5);
fnet6(4) = a6(6)*links(6).mcm(2) - a6(5)*links(6).mcm(3) + (-(links(6).mcm(2)*v6(2)) - links(6).mcm(3)*v6(3))*v6(4) + (links(6).mcm(1)*v6(3) + links(6).m*v6(5))*v6(6) + v6(5)*(links(6).mcm(1)*v6(2) - links(6).m*v6(6)) + a6(1)*links(6).inertia(1,1) + a6(2)*links(6).inertia(1,2) + a6(3)*links(6).inertia(1,3) + v6(1)*(links(6).mcm(2)*v6(5) + links(6).mcm(3)*v6(6) - v6(3)*links(6).inertia(1,2) + v6(2)*links(6).inertia(1,3)) + v6(2)*(-(links(6).mcm(1)*v6(5)) - v6(3)*links(6).inertia(2,2) + v6(2)*links(6).inertia(2,3)) + v6(3)*(-(links(6).mcm(1)*v6(6)) - v6(3)*links(6).inertia(2,3) + v6(2)*links(6).inertia(3,3));
fnet6(5) = -(a6(6)*links(6).mcm(1)) + a6(4)*links(6).mcm(3) + (-(links(6).mcm(1)*v6(1)) - links(6).mcm(3)*v6(3))*v6(5) + (links(6).mcm(2)*v6(3) - links(6).m*v6(4))*v6(6) + v6(4)*(links(6).mcm(2)*v6(1) + links(6).m*v6(6)) + a6(1)*links(6).inertia(1,2) + v6(1)*(-(links(6).mcm(2)*v6(4)) + v6(3)*links(6).inertia(1,1) - v6(1)*links(6).inertia(1,3)) + a6(2)*links(6).inertia(2,2) + a6(3)*links(6).inertia(2,3) + v6(2)*(links(6).mcm(1)*v6(4) + links(6).mcm(3)*v6(6) + v6(3)*links(6).inertia(1,2) - v6(1)*links(6).inertia(2,3)) + v6(3)*(-(links(6).mcm(2)*v6(6)) + v6(3)*links(6).inertia(1,3) - v6(1)*links(6).inertia(3,3));
fnet6(6) = a6(5)*links(6).mcm(1) - a6(4)*links(6).mcm(2) + (links(6).mcm(3)*v6(2) + links(6).m*v6(4))*v6(5) + v6(4)*(links(6).mcm(3)*v6(1) - links(6).m*v6(5)) + (-(links(6).mcm(1)*v6(1)) - links(6).mcm(2)*v6(2))*v6(6) + v6(1)*(-(links(6).mcm(3)*v6(4)) - v6(2)*links(6).inertia(1,1) + v6(1)*links(6).inertia(1,2)) + a6(1)*links(6).inertia(1,3) + v6(2)*(-(links(6).mcm(3)*v6(5)) - v6(2)*links(6).inertia(1,2) + v6(1)*links(6).inertia(2,2)) + a6(2)*links(6).inertia(2,3) + v6(3)*(links(6).mcm(1)*v6(4) + links(6).mcm(2)*v6(5) - v6(2)*links(6).inertia(1,3) + v6(1)*links(6).inertia(2,3)) + a6(3)*links(6).inertia(3,3);

fnet7(1) = links(7).m*a7(4) - a7(3)*links(7).mcm(2) + a7(2)*links(7).mcm(3) - links(7).mcm(1)*power(v7(2),2) - links(7).mcm(1)*power(v7(3),2) + v7(1)*(links(7).mcm(2)*v7(2) + links(7).mcm(3)*v7(3)) - links(7).m*v7(3)*v7(5) + links(7).m*v7(2)*v7(6);
fnet7(2) = links(7).m*a7(5) + a7(3)*links(7).mcm(1) - a7(1)*links(7).mcm(3) - links(7).mcm(2)*power(v7(1),2) - links(7).mcm(2)*power(v7(3),2) + v7(2)*(links(7).mcm(1)*v7(1) + links(7).mcm(3)*v7(3)) + links(7).m*v7(3)*v7(4) - links(7).m*v7(1)*v7(6);
fnet7(3) = links(7).m*a7(6) - a7(2)*links(7).mcm(1) + a7(1)*links(7).mcm(2) - links(7).mcm(3)*power(v7(1),2) - links(7).mcm(3)*power(v7(2),2) + (links(7).mcm(1)*v7(1) + links(7).mcm(2)*v7(2))*v7(3) - links(7).m*v7(2)*v7(4) + links(7).m*v7(1)*v7(5);
fnet7(4) = a7(6)*links(7).mcm(2) - a7(5)*links(7).mcm(3) + (-(links(7).mcm(2)*v7(2)) - links(7).mcm(3)*v7(3))*v7(4) + (links(7).mcm(1)*v7(3) + links(7).m*v7(5))*v7(6) + v7(5)*(links(7).mcm(1)*v7(2) - links(7).m*v7(6)) + a7(1)*links(7).inertia(1,1) + a7(2)*links(7).inertia(1,2) + a7(3)*links(7).inertia(1,3) + v7(1)*(links(7).mcm(2)*v7(5) + links(7).mcm(3)*v7(6) - v7(3)*links(7).inertia(1,2) + v7(2)*links(7).inertia(1,3)) + v7(2)*(-(links(7).mcm(1)*v7(5)) - v7(3)*links(7).inertia(2,2) + v7(2)*links(7).inertia(2,3)) + v7(3)*(-(links(7).mcm(1)*v7(6)) - v7(3)*links(7).inertia(2,3) + v7(2)*links(7).inertia(3,3));
fnet7(5) = -(a7(6)*links(7).mcm(1)) + a7(4)*links(7).mcm(3) + (-(links(7).mcm(1)*v7(1)) - links(7).mcm(3)*v7(3))*v7(5) + (links(7).mcm(2)*v7(3) - links(7).m*v7(4))*v7(6) + v7(4)*(links(7).mcm(2)*v7(1) + links(7).m*v7(6)) + a7(1)*links(7).inertia(1,2) + v7(1)*(-(links(7).mcm(2)*v7(4)) + v7(3)*links(7).inertia(1,1) - v7(1)*links(7).inertia(1,3)) + a7(2)*links(7).inertia(2,2) + a7(3)*links(7).inertia(2,3) + v7(2)*(links(7).mcm(1)*v7(4) + links(7).mcm(3)*v7(6) + v7(3)*links(7).inertia(1,2) - v7(1)*links(7).inertia(2,3)) + v7(3)*(-(links(7).mcm(2)*v7(6)) + v7(3)*links(7).inertia(1,3) - v7(1)*links(7).inertia(3,3));
fnet7(6) = a7(5)*links(7).mcm(1) - a7(4)*links(7).mcm(2) + (links(7).mcm(3)*v7(2) + links(7).m*v7(4))*v7(5) + v7(4)*(links(7).mcm(3)*v7(1) - links(7).m*v7(5)) + (-(links(7).mcm(1)*v7(1)) - links(7).mcm(2)*v7(2))*v7(6) + v7(1)*(-(links(7).mcm(3)*v7(4)) - v7(2)*links(7).inertia(1,1) + v7(1)*links(7).inertia(1,2)) + a7(1)*links(7).inertia(1,3) + v7(2)*(-(links(7).mcm(3)*v7(5)) - v7(2)*links(7).inertia(1,2) + v7(1)*links(7).inertia(2,2)) + a7(2)*links(7).inertia(2,3) + v7(3)*(links(7).mcm(1)*v7(4) + links(7).mcm(2)*v7(5) - v7(2)*links(7).inertia(1,3) + v7(1)*links(7).inertia(2,3)) + a7(3)*links(7).inertia(3,3);

fnet8(1) = eff(1).m*a8(4) - a8(3)*eff(1).mcm(2) + a8(2)*eff(1).mcm(3) - eff(1).mcm(1)*power(v8(2),2) - eff(1).mcm(1)*power(v8(3),2) + v8(1)*(eff(1).mcm(2)*v8(2) + eff(1).mcm(3)*v8(3)) - eff(1).m*v8(3)*v8(5) + eff(1).m*v8(2)*v8(6);
fnet8(2) = eff(1).m*a8(5) + a8(3)*eff(1).mcm(1) - a8(1)*eff(1).mcm(3) - eff(1).mcm(2)*power(v8(1),2) - eff(1).mcm(2)*power(v8(3),2) + v8(2)*(eff(1).mcm(1)*v8(1) + eff(1).mcm(3)*v8(3)) + eff(1).m*v8(3)*v8(4) - eff(1).m*v8(1)*v8(6);
fnet8(3) = eff(1).m*a8(6) - a8(2)*eff(1).mcm(1) + a8(1)*eff(1).mcm(2) - eff(1).mcm(3)*power(v8(1),2) - eff(1).mcm(3)*power(v8(2),2) + (eff(1).mcm(1)*v8(1) + eff(1).mcm(2)*v8(2))*v8(3) - eff(1).m*v8(2)*v8(4) + eff(1).m*v8(1)*v8(5);
fnet8(4) = a8(6)*eff(1).mcm(2) - a8(5)*eff(1).mcm(3) + (-(eff(1).mcm(2)*v8(2)) - eff(1).mcm(3)*v8(3))*v8(4) - eff(1).mcm(1)*v8(2)*v8(5) - eff(1).mcm(1)*v8(3)*v8(6) + (eff(1).mcm(1)*v8(3) + eff(1).m*v8(5))*v8(6) + v8(5)*(eff(1).mcm(1)*v8(2) - eff(1).m*v8(6)) + v8(1)*(eff(1).mcm(2)*v8(5) + eff(1).mcm(3)*v8(6));
fnet8(5) = -(a8(6)*eff(1).mcm(1)) + a8(4)*eff(1).mcm(3) - eff(1).mcm(2)*v8(1)*v8(4) + (-(eff(1).mcm(1)*v8(1)) - eff(1).mcm(3)*v8(3))*v8(5) - eff(1).mcm(2)*v8(3)*v8(6) + (eff(1).mcm(2)*v8(3) - eff(1).m*v8(4))*v8(6) + v8(4)*(eff(1).mcm(2)*v8(1) + eff(1).m*v8(6)) + v8(2)*(eff(1).mcm(1)*v8(4) + eff(1).mcm(3)*v8(6));
fnet8(6) = a8(5)*eff(1).mcm(1) - a8(4)*eff(1).mcm(2) - eff(1).mcm(3)*v8(1)*v8(4) - eff(1).mcm(3)*v8(2)*v8(5) + (eff(1).mcm(3)*v8(2) + eff(1).m*v8(4))*v8(5) + v8(4)*(eff(1).mcm(3)*v8(1) - eff(1).m*v8(5)) + v8(3)*(eff(1).mcm(1)*v8(4) + eff(1).mcm(2)*v8(5)) + (-(eff(1).mcm(1)*v8(1)) - eff(1).mcm(2)*v8(2))*v8(6);


fex0(1) = -(uex0.f(1)*S00(1,1)) - uex0.f(2)*S00(1,2) - uex0.f(3)*S00(1,3);
fex0(2) = -(uex0.f(1)*S00(2,1)) - uex0.f(2)*S00(2,2) - uex0.f(3)*S00(2,3);
fex0(3) = -(uex0.f(1)*S00(3,1)) - uex0.f(2)*S00(3,2) - uex0.f(3)*S00(3,3);
fex0(4) = -(uex0.t(1)*S00(1,1)) - uex0.t(2)*S00(1,2) - uex0.t(3)*S00(1,3);
fex0(5) = -(uex0.t(1)*S00(2,1)) - uex0.t(2)*S00(2,2) - uex0.t(3)*S00(2,3);
fex0(6) = -(uex0.t(1)*S00(3,1)) - uex0.t(2)*S00(3,2) - uex0.t(3)*S00(3,3);

fex1(1) = -(uex(1).f(1)*SG10(1,1)) - uex(1).f(2)*SG10(1,2) - uex(1).f(3)*SG10(1,3);
fex1(2) = -(uex(1).f(1)*SG10(2,1)) - uex(1).f(2)*SG10(2,2) - uex(1).f(3)*SG10(2,3);
fex1(3) = -(uex(1).f(1)*SG10(3,1)) - uex(1).f(2)*SG10(3,2) - uex(1).f(3)*SG10(3,3);
fex1(4) = -(uex(1).t(1)*SG10(1,1)) - uex(1).t(2)*SG10(1,2) - uex(1).t(3)*SG10(1,3);
fex1(5) = -(uex(1).t(1)*SG10(2,1)) - uex(1).t(2)*SG10(2,2) - uex(1).t(3)*SG10(2,3);
fex1(6) = -(uex(1).t(1)*SG10(3,1)) - uex(1).t(2)*SG10(3,2) - uex(1).t(3)*SG10(3,3);

fex2(1) = -(uex(2).f(1)*SG20(1,1)) - uex(2).f(2)*SG20(1,2) - uex(2).f(3)*SG20(1,3);
fex2(2) = -(uex(2).f(1)*SG20(2,1)) - uex(2).f(2)*SG20(2,2) - uex(2).f(3)*SG20(2,3);
fex2(3) = -(uex(2).f(1)*SG20(3,1)) - uex(2).f(2)*SG20(3,2) - uex(2).f(3)*SG20(3,3);
fex2(4) = -(uex(2).t(1)*SG20(1,1)) - uex(2).t(2)*SG20(1,2) - uex(2).t(3)*SG20(1,3);
fex2(5) = -(uex(2).t(1)*SG20(2,1)) - uex(2).t(2)*SG20(2,2) - uex(2).t(3)*SG20(2,3);
fex2(6) = -(uex(2).t(1)*SG20(3,1)) - uex(2).t(2)*SG20(3,2) - uex(2).t(3)*SG20(3,3);

fex3(1) = -(uex(3).f(1)*SG30(1,1)) - uex(3).f(2)*SG30(1,2) - uex(3).f(3)*SG30(1,3);
fex3(2) = -(uex(3).f(1)*SG30(2,1)) - uex(3).f(2)*SG30(2,2) - uex(3).f(3)*SG30(2,3);
fex3(3) = -(uex(3).f(1)*SG30(3,1)) - uex(3).f(2)*SG30(3,2) - uex(3).f(3)*SG30(3,3);
fex3(4) = -(uex(3).t(1)*SG30(1,1)) - uex(3).t(2)*SG30(1,2) - uex(3).t(3)*SG30(1,3);
fex3(5) = -(uex(3).t(1)*SG30(2,1)) - uex(3).t(2)*SG30(2,2) - uex(3).t(3)*SG30(2,3);
fex3(6) = -(uex(3).t(1)*SG30(3,1)) - uex(3).t(2)*SG30(3,2) - uex(3).t(3)*SG30(3,3);

fex4(1) = -(uex(4).f(1)*SG40(1,1)) - uex(4).f(2)*SG40(1,2) - uex(4).f(3)*SG40(1,3);
fex4(2) = -(uex(4).f(1)*SG40(2,1)) - uex(4).f(2)*SG40(2,2) - uex(4).f(3)*SG40(2,3);
fex4(3) = -(uex(4).f(1)*SG40(3,1)) - uex(4).f(2)*SG40(3,2) - uex(4).f(3)*SG40(3,3);
fex4(4) = -(uex(4).t(1)*SG40(1,1)) - uex(4).t(2)*SG40(1,2) - uex(4).t(3)*SG40(1,3);
fex4(5) = -(uex(4).t(1)*SG40(2,1)) - uex(4).t(2)*SG40(2,2) - uex(4).t(3)*SG40(2,3);
fex4(6) = -(uex(4).t(1)*SG40(3,1)) - uex(4).t(2)*SG40(3,2) - uex(4).t(3)*SG40(3,3);

fex5(1) = -(uex(5).f(1)*SG50(1,1)) - uex(5).f(2)*SG50(1,2) - uex(5).f(3)*SG50(1,3);
fex5(2) = -(uex(5).f(1)*SG50(2,1)) - uex(5).f(2)*SG50(2,2) - uex(5).f(3)*SG50(2,3);
fex5(3) = -(uex(5).f(1)*SG50(3,1)) - uex(5).f(2)*SG50(3,2) - uex(5).f(3)*SG50(3,3);
fex5(4) = -(uex(5).t(1)*SG50(1,1)) - uex(5).t(2)*SG50(1,2) - uex(5).t(3)*SG50(1,3);
fex5(5) = -(uex(5).t(1)*SG50(2,1)) - uex(5).t(2)*SG50(2,2) - uex(5).t(3)*SG50(2,3);
fex5(6) = -(uex(5).t(1)*SG50(3,1)) - uex(5).t(2)*SG50(3,2) - uex(5).t(3)*SG50(3,3);

fex6(1) = -(uex(6).f(1)*SG60(1,1)) - uex(6).f(2)*SG60(1,2) - uex(6).f(3)*SG60(1,3);
fex6(2) = -(uex(6).f(1)*SG60(2,1)) - uex(6).f(2)*SG60(2,2) - uex(6).f(3)*SG60(2,3);
fex6(3) = -(uex(6).f(1)*SG60(3,1)) - uex(6).f(2)*SG60(3,2) - uex(6).f(3)*SG60(3,3);
fex6(4) = -(uex(6).t(1)*SG60(1,1)) - uex(6).t(2)*SG60(1,2) - uex(6).t(3)*SG60(1,3);
fex6(5) = -(uex(6).t(1)*SG60(2,1)) - uex(6).t(2)*SG60(2,2) - uex(6).t(3)*SG60(2,3);
fex6(6) = -(uex(6).t(1)*SG60(3,1)) - uex(6).t(2)*SG60(3,2) - uex(6).t(3)*SG60(3,3);

fex7(1) = -(uex(7).f(1)*SG70(1,1)) - uex(7).f(2)*SG70(1,2) - uex(7).f(3)*SG70(1,3);
fex7(2) = -(uex(7).f(1)*SG70(2,1)) - uex(7).f(2)*SG70(2,2) - uex(7).f(3)*SG70(2,3);
fex7(3) = -(uex(7).f(1)*SG70(3,1)) - uex(7).f(2)*SG70(3,2) - uex(7).f(3)*SG70(3,3);
fex7(4) = -(uex(7).t(1)*SG70(1,1)) - uex(7).t(2)*SG70(1,2) - uex(7).t(3)*SG70(1,3);
fex7(5) = -(uex(7).t(1)*SG70(2,1)) - uex(7).t(2)*SG70(2,2) - uex(7).t(3)*SG70(2,3);
fex7(6) = -(uex(7).t(1)*SG70(3,1)) - uex(7).t(2)*SG70(3,2) - uex(7).t(3)*SG70(3,3);

%barrett_InvDynNEfunc7(void)
     
% total forces and external forces for each joint 
f8(1) = fnet8(1);
f8(2) = fnet8(2);
f8(3) = fnet8(3);
f8(4) = fnet8(4);
f8(5) = fnet8(5);
f8(6) = fnet8(6);

f7(1) = fnet7(1) + f8(1)*Si78(1,1) + f8(2)*Si78(1,2) + f8(3)*Si78(1,3);
f7(2) = fnet7(2) + f8(1)*Si78(2,1) + f8(2)*Si78(2,2) + f8(3)*Si78(2,3);
f7(3) = fnet7(3) + f8(1)*Si78(3,1) + f8(2)*Si78(3,2) + f8(3)*Si78(3,3);
f7(4) = fnet7(4) + f8(4)*Si78(1,1) + f8(5)*Si78(1,2) + f8(6)*Si78(1,3) + f8(1)*(-(eff(1).x(3)*Si78(2,1)) + eff(1).x(2)*Si78(3,1)) + f8(2)*(-(eff(1).x(3)*Si78(2,2)) + eff(1).x(2)*Si78(3,2)) + f8(3)*(-(eff(1).x(3)*Si78(2,3)) + eff(1).x(2)*Si78(3,3));
f7(5) = fnet7(5) + f8(4)*Si78(2,1) + f8(5)*Si78(2,2) + f8(6)*Si78(2,3) + f8(1)*(eff(1).x(3)*Si78(1,1) - eff(1).x(1)*Si78(3,1)) + f8(2)*(eff(1).x(3)*Si78(1,2) - eff(1).x(1)*Si78(3,2)) + f8(3)*(eff(1).x(3)*Si78(1,3) - eff(1).x(1)*Si78(3,3));
f7(6) = fnet7(6) + f8(1)*(-(eff(1).x(2)*Si78(1,1)) + eff(1).x(1)*Si78(2,1)) + f8(2)*(-(eff(1).x(2)*Si78(1,2)) + eff(1).x(1)*Si78(2,2)) + f8(3)*(-(eff(1).x(2)*Si78(1,3)) + eff(1).x(1)*Si78(2,3)) + f8(4)*Si78(3,1) + f8(5)*Si78(3,2) + f8(6)*Si78(3,3);

f6(1) = f7(3) + fnet6(1);
f6(2) = fnet6(2) + f7(1)*Si67(2,1) + f7(2)*Si67(2,2);
f6(3) = fnet6(3) + f7(1)*Si67(3,1) + f7(2)*Si67(3,2);
f6(4) = f7(6) + fnet6(4);
f6(5) = fnet6(5) + f7(4)*Si67(2,1) + f7(5)*Si67(2,2);
f6(6) = fnet6(6) + f7(4)*Si67(3,1) + f7(5)*Si67(3,2);

f5(1) = -f6(3) + fnet5(1);
f5(2) = fnet5(2) + f6(1)*Si56(2,1) + f6(2)*Si56(2,2);
f5(3) = fnet5(3) + f6(1)*Si56(3,1) + f6(2)*Si56(3,2);
f5(4) = -f6(6) + fnet5(4) - ZWFE*f6(1)*Si56(2,1) - ZWFE*f6(2)*Si56(2,2);
f5(5) = -(ZWFE*f6(3)) + fnet5(5) + f6(4)*Si56(2,1) + f6(5)*Si56(2,2);
f5(6) = fnet5(6) + f6(4)*Si56(3,1) + f6(5)*Si56(3,2);

f4(1) = f5(3) + fnet4(1);
f4(2) = fnet4(2) + f5(1)*Si45(2,1) + f5(2)*Si45(2,2);
f4(3) = fnet4(3) + f5(1)*Si45(3,1) + f5(2)*Si45(3,2);
f4(4) = f5(6) + fnet4(4) + YWR*f5(1)*Si45(3,1) + YWR*f5(2)*Si45(3,2);
f4(5) = fnet4(5) + f5(4)*Si45(2,1) + f5(5)*Si45(2,2) - ZWR*f5(1)*Si45(3,1) - ZWR*f5(2)*Si45(3,2);
f4(6) = -(YWR*f5(3)) + fnet4(6) + ZWR*f5(1)*Si45(2,1) + ZWR*f5(2)*Si45(2,2) + f5(4)*Si45(3,1) + f5(5)*Si45(3,2);

f3(1) = -f4(3) + fnet3(1);
f3(2) = fnet3(2) + f4(1)*Si34(2,1) + f4(2)*Si34(2,2);
f3(3) = fnet3(3) + f4(1)*Si34(3,1) + f4(2)*Si34(3,2);
f3(4) = -f4(6) + fnet3(4) + f4(1)*(-(ZEB*Si34(2,1)) + YEB*Si34(3,1)) + f4(2)*(-(ZEB*Si34(2,2)) + YEB*Si34(3,2));
f3(5) = -(ZEB*f4(3)) + fnet3(5) + f4(4)*Si34(2,1) + f4(5)*Si34(2,2);
f3(6) = YEB*f4(3) + fnet3(6) + f4(4)*Si34(3,1) + f4(5)*Si34(3,2);

f2(1) = f3(3) + fnet2(1);
f2(2) = fnet2(2) + f3(1)*Si23(2,1) + f3(2)*Si23(2,2);
f2(3) = fnet2(3) + f3(1)*Si23(3,1) + f3(2)*Si23(3,2);
f2(4) = f3(6) + fnet2(4);
f2(5) = fnet2(5) + f3(4)*Si23(2,1) + f3(5)*Si23(2,2) - ZHR*f3(1)*Si23(3,1) - ZHR*f3(2)*Si23(3,2);
f2(6) = fnet2(6) + ZHR*f3(1)*Si23(2,1) + ZHR*f3(2)*Si23(2,2) + f3(4)*Si23(3,1) + f3(5)*Si23(3,2);

f1(1) = -f2(3) + fnet1(1);
f1(2) = fnet1(2) + f2(1)*Si12(2,1) + f2(2)*Si12(2,2);
f1(3) = fnet1(3) + f2(1)*Si12(3,1) + f2(2)*Si12(3,2);
f1(4) = -f2(6) + fnet1(4);
f1(5) = fnet1(5) + f2(4)*Si12(2,1) + f2(5)*Si12(2,2);
f1(6) = fnet1(6) + f2(4)*Si12(3,1) + f2(5)*Si12(3,2);

f0(1) = fnet0(1) + f1(1)*Si01(1,1) + f1(2)*Si01(1,2);
f0(2) = fnet0(2) + f1(1)*Si01(2,1) + f1(2)*Si01(2,2);
f0(3) = f1(3) + fnet0(3);
f0(4) = fnet0(4) + f1(4)*Si01(1,1) + f1(5)*Si01(1,2) - ZSFE*f1(1)*Si01(2,1) - ZSFE*f1(2)*Si01(2,2);
f0(5) = fnet0(5) + ZSFE*f1(1)*Si01(1,1) + ZSFE*f1(2)*Si01(1,2) + f1(4)*Si01(2,1) + f1(5)*Si01(2,2);
f0(6) = f1(6) + fnet0(6);


fext7(1) = fex7(1);
fext7(2) = fex7(2);
fext7(3) = fex7(3);
fext7(4) = fex7(4);
fext7(5) = fex7(5);
fext7(6) = fex7(6);

fext6(1) = fex6(1) + fext7(3);
fext6(2) = fex6(2) + fext7(1)*Si67(2,1) + fext7(2)*Si67(2,2);
fext6(3) = fex6(3) + fext7(1)*Si67(3,1) + fext7(2)*Si67(3,2);
fext6(4) = fex6(4) + fext7(6);
fext6(5) = fex6(5) + fext7(4)*Si67(2,1) + fext7(5)*Si67(2,2);
fext6(6) = fex6(6) + fext7(4)*Si67(3,1) + fext7(5)*Si67(3,2);

fext5(1) = fex5(1) - fext6(3);
fext5(2) = fex5(2) + fext6(1)*Si56(2,1) + fext6(2)*Si56(2,2);
fext5(3) = fex5(3) + fext6(1)*Si56(3,1) + fext6(2)*Si56(3,2);
fext5(4) = fex5(4) - fext6(6) - ZWFE*fext6(1)*Si56(2,1) - ZWFE*fext6(2)*Si56(2,2);
fext5(5) = fex5(5) - ZWFE*fext6(3) + fext6(4)*Si56(2,1) + fext6(5)*Si56(2,2);
fext5(6) = fex5(6) + fext6(4)*Si56(3,1) + fext6(5)*Si56(3,2);

fext4(1) = fex4(1) + fext5(3);
fext4(2) = fex4(2) + fext5(1)*Si45(2,1) + fext5(2)*Si45(2,2);
fext4(3) = fex4(3) + fext5(1)*Si45(3,1) + fext5(2)*Si45(3,2);
fext4(4) = fex4(4) + fext5(6) + YWR*fext5(1)*Si45(3,1) + YWR*fext5(2)*Si45(3,2);
fext4(5) = fex4(5) + fext5(4)*Si45(2,1) + fext5(5)*Si45(2,2) - ZWR*fext5(1)*Si45(3,1) - ZWR*fext5(2)*Si45(3,2);
fext4(6) = fex4(6) - YWR*fext5(3) + ZWR*fext5(1)*Si45(2,1) + ZWR*fext5(2)*Si45(2,2) + fext5(4)*Si45(3,1) + fext5(5)*Si45(3,2);

fext3(1) = fex3(1) - fext4(3);
fext3(2) = fex3(2) + fext4(1)*Si34(2,1) + fext4(2)*Si34(2,2);
fext3(3) = fex3(3) + fext4(1)*Si34(3,1) + fext4(2)*Si34(3,2);
fext3(4) = fex3(4) - fext4(6) + fext4(1)*(-(ZEB*Si34(2,1)) + YEB*Si34(3,1)) + fext4(2)*(-(ZEB*Si34(2,2)) + YEB*Si34(3,2));
fext3(5) = fex3(5) - ZEB*fext4(3) + fext4(4)*Si34(2,1) + fext4(5)*Si34(2,2);
fext3(6) = fex3(6) + YEB*fext4(3) + fext4(4)*Si34(3,1) + fext4(5)*Si34(3,2);

fext2(1) = fex2(1) + fext3(3);
fext2(2) = fex2(2) + fext3(1)*Si23(2,1) + fext3(2)*Si23(2,2);
fext2(3) = fex2(3) + fext3(1)*Si23(3,1) + fext3(2)*Si23(3,2);
fext2(4) = fex2(4) + fext3(6);
fext2(5) = fex2(5) + fext3(4)*Si23(2,1) + fext3(5)*Si23(2,2) - ZHR*fext3(1)*Si23(3,1) - ZHR*fext3(2)*Si23(3,2);
fext2(6) = fex2(6) + ZHR*fext3(1)*Si23(2,1) + ZHR*fext3(2)*Si23(2,2) + fext3(4)*Si23(3,1) + fext3(5)*Si23(3,2);

fext1(1) = fex1(1) - fext2(3);
fext1(2) = fex1(2) + fext2(1)*Si12(2,1) + fext2(2)*Si12(2,2);
fext1(3) = fex1(3) + fext2(1)*Si12(3,1) + fext2(2)*Si12(3,2);
fext1(4) = fex1(4) - fext2(6);
fext1(5) = fex1(5) + fext2(4)*Si12(2,1) + fext2(5)*Si12(2,2);
fext1(6) = fex1(6) + fext2(4)*Si12(3,1) + fext2(5)*Si12(3,2);

fext0(1) = fex0(1) + fext1(1)*Si01(1,1) + fext1(2)*Si01(1,2);
fext0(2) = fex0(2) + fext1(1)*Si01(2,1) + fext1(2)*Si01(2,2);
fext0(3) = fex0(3) + fext1(3);
fext0(4) = fex0(4) + fext1(4)*Si01(1,1) + fext1(5)*Si01(1,2) - ZSFE*fext1(1)*Si01(2,1) - ZSFE*fext1(2)*Si01(2,2);
fext0(5) = fex0(5) + ZSFE*fext1(1)*Si01(1,1) + ZSFE*fext1(2)*Si01(1,2) + fext1(4)*Si01(2,1) + fext1(5)*Si01(2,2);
fext0(6) = fex0(6) + fext1(6);


%% Acceleration vectors, joint torques 

% force/torque of base in world coordinates 
fbase(1) = f0(1)*Si00(1,1) + f0(2)*Si00(1,2) + f0(3)*Si00(1,3);
fbase(2) = f0(1)*Si00(2,1) + f0(2)*Si00(2,2) + f0(3)*Si00(2,3);
fbase(3) = f0(1)*Si00(3,1) + f0(2)*Si00(3,2) + f0(3)*Si00(3,3);
fbase(4) = f0(4)*Si00(1,1) + f0(5)*Si00(1,2) + f0(6)*Si00(1,3);
fbase(5) = f0(4)*Si00(2,1) + f0(5)*Si00(2,2) + f0(6)*Si00(2,3);
fbase(6) = f0(4)*Si00(3,1) + f0(5)*Si00(3,2) + f0(6)*Si00(3,3);

% inverse dynamics torques 
u(1) = -uex_des(1) + f1(6) + fext1(6);
u(2) = -uex_des(2) + f2(6) + fext2(6);
u(3) = -uex_des(3) + f3(6) + fext3(6);
u(4) = -uex_des(4) + f4(6) + fext4(6);
u(5) = -uex_des(5) + f5(6) + fext5(6);
u(6) = -uex_des(6) + f6(6) + fext6(6);
u(7) = -uex_des(7) + f7(6) + fext7(6);

% torques due to external forces 
qext(1) = -uex_des(1) + fext1(6);
qext(2) = -uex_des(2) + fext2(6);
qext(3) = -uex_des(3) + fext3(6);
qext(4) = -uex_des(4) + fext4(6);
qext(5) = -uex_des(5) + fext5(6);
qext(6) = -uex_des(6) + fext6(6);
qext(7) = -uex_des(7) + fext7(6);



