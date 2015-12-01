% Barrett WAM forward kinematics
% used to show trajectories in cartesian space
%
% Function taken from SL: 
% shared/barrett/math/LInfo_declare.h
% shared/barrett/math/LInfo_math.h
%
% these are called from the kinematics method of Barrett WAM class
%
% \param[out]    Xaxis   : array of rotation axes (z)
% \param[out]    Xorigin : array of coord.sys. origin vectors
% \param[out]    Xlink   : array of link position
% \param[out]    Amats   : homogeneous transformation matrices of each link
%

function [Xlink,Xorigin,Xaxis,Ahmat] = barrettWamKinematics(Q,PAR)

NDOF = 7;
% system states are X   =   (q(1),...q(7));
q = Q(1:NDOF);

% definitions
ZSFE = 0.346;              %!< z height of SAA axis above ground
ZHR = 0.505;              %!< length of upper arm until 4.5cm before elbow link
YEB = 0.045;              %!< elbow y offset
ZEB = 0.045;              %!< elbow z offset
YWR = -0.045;              %!< elbow y offset (back to forewarm)
ZWR = 0.045;              %!< elbow z offset (back to forearm)
ZWFE = 0.255;              %!< forearm length (minus 4.5cm)

% extract parameters
link0 = PAR.link0;
links = PAR.links;
eff = PAR.eff;
uex = PAR.uex;
uex0 = PAR.uex0;
basec = PAR.basec;
baseo = PAR.baseo;

% sine and cosine precomputation 
ss1th = sin(q(1));
cs1th = cos(q(1));
ss2th = sin(q(2));
cs2th = cos(q(2));
ss3th = sin(q(3));
cs3th = cos(q(3));
ss4th = sin(q(4));
cs4th = cos(q(4));
ss5th = sin(q(5));
cs5th = cos(q(5));
ss6th = sin(q(6));
cs6th = cos(q(6));
ss7th = sin(q(7));
cs7th = cos(q(7));

% endeffector orientations

rseff1a1 = sin(eff(1).a(1));
rceff1a1 = cos(eff(1).a(1));
rseff1a2 = sin(eff(1).a(2));
rceff1a2 = cos(eff(1).a(2));
rseff1a3 = sin(eff(1).a(3));
rceff1a3 = cos(eff(1).a(3));

%% Calculations are done here

% inverse homogeneous rotation matrices 
Hi00(1,1) = -1 + 2*power(baseo.q(1),2) + 2*power(baseo.q(2),2);
Hi00(1,2) = 2*(baseo.q(2)*baseo.q(3) - baseo.q(1)*baseo.q(4));
Hi00(1,3) = 2*(baseo.q(1)*baseo.q(3) + baseo.q(2)*baseo.q(4));
Hi00(1,4) = basec.x(1);
Hi00(2,1) = 2*(baseo.q(2)*baseo.q(3) + baseo.q(1)*baseo.q(4));
Hi00(2,2) = -1 + 2*power(baseo.q(1),2) + 2*power(baseo.q(3),2);
Hi00(2,3) = 2*(-(baseo.q(1)*baseo.q(2)) + baseo.q(3)*baseo.q(4));
Hi00(2,4) = basec.x(2);
Hi00(3,1) = 2*(-(baseo.q(1)*baseo.q(3)) + baseo.q(2)*baseo.q(4));
Hi00(3,2) = 2*(baseo.q(1)*baseo.q(2) + baseo.q(3)*baseo.q(4));
Hi00(3,3) = -1 + 2*power(baseo.q(1),2) + 2*power(baseo.q(4),2);
Hi00(3,4) = basec.x(3);
Hi01(1,1) = cs1th;
Hi01(1,2) = -ss1th;
Hi01(2,1) = ss1th;
Hi01(2,2) = cs1th;
Hi01(3,4) = ZSFE;
Hi12(2,1) = ss2th;
Hi12(2,2) = cs2th;
Hi12(3,1) = cs2th;
Hi12(3,2) = -ss2th;
Hi23(1,4) = ZHR;
Hi23(2,1) = ss3th;
Hi23(2,2) = cs3th;
Hi23(3,1) = -cs3th;
Hi23(3,2) = ss3th;
Hi34(2,1) = ss4th;
Hi34(2,2) = cs4th;
Hi34(2,4) = YEB;
Hi34(3,1) = cs4th;
Hi34(3,2) = -ss4th;
Hi34(3,4) = ZEB;
Hi45(1,4) = ZWR;
Hi45(2,1) = ss5th;
Hi45(2,2) = cs5th;
Hi45(2,4) = YWR;
Hi45(3,1) = -cs5th;
Hi45(3,2) = ss5th;
Hi56(2,1) = ss6th;
Hi56(2,2) = cs6th;
Hi56(3,1) = cs6th;
Hi56(3,2) = -ss6th;
Hi56(3,4) = ZWFE;
Hi67(2,1) = ss7th;
Hi67(2,2) = cs7th;
Hi67(3,1) = -cs7th;
Hi67(3,2) = ss7th;
Hi78(1,1) = rceff1a2*rceff1a3;
Hi78(1,2) = -(rceff1a2*rseff1a3);
Hi78(1,3) = rseff1a2;
Hi78(1,4) = eff(1).x(1);
Hi78(2,1) = rceff1a3*rseff1a1*rseff1a2 + rceff1a1*rseff1a3;
Hi78(2,2) = rceff1a1*rceff1a3 - rseff1a1*rseff1a2*rseff1a3;
Hi78(2,3) = -(rceff1a2*rseff1a1);
Hi78(2,4) = eff(1).x(2);
Hi78(3,1) = -(rceff1a1*rceff1a3*rseff1a2) + rseff1a1*rseff1a3;
Hi78(3,2) = rceff1a3*rseff1a1 + rceff1a1*rseff1a2*rseff1a3;
Hi78(3,3) = rceff1a1*rceff1a2;
Hi78(3,4) = eff(1).x(3);

% per link inverse homogeneous rotation matrices 
Ai01(1,1) = Hi00(1,1)*Hi01(1,1) + Hi00(1,2)*Hi01(2,1);
Ai01(1,2) = Hi00(1,1)*Hi01(1,2) + Hi00(1,2)*Hi01(2,2);
Ai01(1,3) = Hi00(1,3);
Ai01(1,4) = Hi00(1,4) + Hi00(1,3)*Hi01(3,4);
Ai01(2,1) = Hi00(2,1)*Hi01(1,1) + Hi00(2,2)*Hi01(2,1);
Ai01(2,2) = Hi00(2,1)*Hi01(1,2) + Hi00(2,2)*Hi01(2,2);
Ai01(2,3) = Hi00(2,3);
Ai01(2,4) = Hi00(2,4) + Hi00(2,3)*Hi01(3,4);
Ai01(3,1) = Hi00(3,1)*Hi01(1,1) + Hi00(3,2)*Hi01(2,1);
Ai01(3,2) = Hi00(3,1)*Hi01(1,2) + Hi00(3,2)*Hi01(2,2);
Ai01(3,3) = Hi00(3,3);
Ai01(3,4) = Hi00(3,4) + Hi00(3,3)*Hi01(3,4);
Ai02(1,1) = Ai01(1,2)*Hi12(2,1) + Ai01(1,3)*Hi12(3,1);
Ai02(1,2) = Ai01(1,2)*Hi12(2,2) + Ai01(1,3)*Hi12(3,2);
Ai02(1,3) = -Ai01(1,1);
Ai02(1,4) = Ai01(1,4);
Ai02(2,1) = Ai01(2,2)*Hi12(2,1) + Ai01(2,3)*Hi12(3,1);
Ai02(2,2) = Ai01(2,2)*Hi12(2,2) + Ai01(2,3)*Hi12(3,2);
Ai02(2,3) = -Ai01(2,1);
Ai02(2,4) = Ai01(2,4);
Ai02(3,1) = Ai01(3,2)*Hi12(2,1) + Ai01(3,3)*Hi12(3,1);
Ai02(3,2) = Ai01(3,2)*Hi12(2,2) + Ai01(3,3)*Hi12(3,2);
Ai02(3,3) = -Ai01(3,1);
Ai02(3,4) = Ai01(3,4);
Ai03(1,1) = Ai02(1,2)*Hi23(2,1) + Ai02(1,3)*Hi23(3,1);
Ai03(1,2) = Ai02(1,2)*Hi23(2,2) + Ai02(1,3)*Hi23(3,2);
Ai03(1,3) = Ai02(1,1);
Ai03(1,4) = Ai02(1,4) + Ai02(1,1)*Hi23(1,4);
Ai03(2,1) = Ai02(2,2)*Hi23(2,1) + Ai02(2,3)*Hi23(3,1);
Ai03(2,2) = Ai02(2,2)*Hi23(2,2) + Ai02(2,3)*Hi23(3,2);
Ai03(2,3) = Ai02(2,1);
Ai03(2,4) = Ai02(2,4) + Ai02(2,1)*Hi23(1,4);
Ai03(3,1) = Ai02(3,2)*Hi23(2,1) + Ai02(3,3)*Hi23(3,1);
Ai03(3,2) = Ai02(3,2)*Hi23(2,2) + Ai02(3,3)*Hi23(3,2);
Ai03(3,3) = Ai02(3,1);
Ai03(3,4) = Ai02(3,4) + Ai02(3,1)*Hi23(1,4);
Ai04(1,1) = Ai03(1,2)*Hi34(2,1) + Ai03(1,3)*Hi34(3,1);
Ai04(1,2) = Ai03(1,2)*Hi34(2,2) + Ai03(1,3)*Hi34(3,2);
Ai04(1,3) = -Ai03(1,1);
Ai04(1,4) = Ai03(1,4) + Ai03(1,2)*Hi34(2,4) + Ai03(1,3)*Hi34(3,4);
Ai04(2,1) = Ai03(2,2)*Hi34(2,1) + Ai03(2,3)*Hi34(3,1);
Ai04(2,2) = Ai03(2,2)*Hi34(2,2) + Ai03(2,3)*Hi34(3,2);
Ai04(2,3) = -Ai03(2,1);
Ai04(2,4) = Ai03(2,4) + Ai03(2,2)*Hi34(2,4) + Ai03(2,3)*Hi34(3,4);
Ai04(3,1) = Ai03(3,2)*Hi34(2,1) + Ai03(3,3)*Hi34(3,1);
Ai04(3,2) = Ai03(3,2)*Hi34(2,2) + Ai03(3,3)*Hi34(3,2);
Ai04(3,3) = -Ai03(3,1);
Ai04(3,4) = Ai03(3,4) + Ai03(3,2)*Hi34(2,4) + Ai03(3,3)*Hi34(3,4);
Ai05(1,1) = Ai04(1,2)*Hi45(2,1) + Ai04(1,3)*Hi45(3,1);
Ai05(1,2) = Ai04(1,2)*Hi45(2,2) + Ai04(1,3)*Hi45(3,2);
Ai05(1,3) = Ai04(1,1);
Ai05(1,4) = Ai04(1,4) + Ai04(1,1)*Hi45(1,4) + Ai04(1,2)*Hi45(2,4);
Ai05(2,1) = Ai04(2,2)*Hi45(2,1) + Ai04(2,3)*Hi45(3,1);
Ai05(2,2) = Ai04(2,2)*Hi45(2,2) + Ai04(2,3)*Hi45(3,2);
Ai05(2,3) = Ai04(2,1);
Ai05(2,4) = Ai04(2,4) + Ai04(2,1)*Hi45(1,4) + Ai04(2,2)*Hi45(2,4);
Ai05(3,1) = Ai04(3,2)*Hi45(2,1) + Ai04(3,3)*Hi45(3,1);
Ai05(3,2) = Ai04(3,2)*Hi45(2,2) + Ai04(3,3)*Hi45(3,2);
Ai05(3,3) = Ai04(3,1);
Ai05(3,4) = Ai04(3,4) + Ai04(3,1)*Hi45(1,4) + Ai04(3,2)*Hi45(2,4);
Ai06(1,1) = Ai05(1,2)*Hi56(2,1) + Ai05(1,3)*Hi56(3,1);
Ai06(1,2) = Ai05(1,2)*Hi56(2,2) + Ai05(1,3)*Hi56(3,2);
Ai06(1,3) = -Ai05(1,1);
Ai06(1,4) = Ai05(1,4) + Ai05(1,3)*Hi56(3,4);
Ai06(2,1) = Ai05(2,2)*Hi56(2,1) + Ai05(2,3)*Hi56(3,1);
Ai06(2,2) = Ai05(2,2)*Hi56(2,2) + Ai05(2,3)*Hi56(3,2);
Ai06(2,3) = -Ai05(2,1);
Ai06(2,4) = Ai05(2,4) + Ai05(2,3)*Hi56(3,4);
Ai06(3,1) = Ai05(3,2)*Hi56(2,1) + Ai05(3,3)*Hi56(3,1);
Ai06(3,2) = Ai05(3,2)*Hi56(2,2) + Ai05(3,3)*Hi56(3,2);
Ai06(3,3) = -Ai05(3,1);
Ai06(3,4) = Ai05(3,4) + Ai05(3,3)*Hi56(3,4);
Ai07(1,1) = Ai06(1,2)*Hi67(2,1) + Ai06(1,3)*Hi67(3,1);
Ai07(1,2) = Ai06(1,2)*Hi67(2,2) + Ai06(1,3)*Hi67(3,2);
Ai07(1,3) = Ai06(1,1);
Ai07(1,4) = Ai06(1,4);
Ai07(2,1) = Ai06(2,2)*Hi67(2,1) + Ai06(2,3)*Hi67(3,1);
Ai07(2,2) = Ai06(2,2)*Hi67(2,2) + Ai06(2,3)*Hi67(3,2);
Ai07(2,3) = Ai06(2,1);
Ai07(2,4) = Ai06(2,4);
Ai07(3,1) = Ai06(3,2)*Hi67(2,1) + Ai06(3,3)*Hi67(3,1);
Ai07(3,2) = Ai06(3,2)*Hi67(2,2) + Ai06(3,3)*Hi67(3,2);
Ai07(3,3) = Ai06(3,1);
Ai07(3,4) = Ai06(3,4);
Ai08(1,1) = Ai07(1,1)*Hi78(1,1) + Ai07(1,2)*Hi78(2,1) + Ai07(1,3)*Hi78(3,1);
Ai08(1,2) = Ai07(1,1)*Hi78(1,2) + Ai07(1,2)*Hi78(2,2) + Ai07(1,3)*Hi78(3,2);
Ai08(1,3) = Ai07(1,1)*Hi78(1,3) + Ai07(1,2)*Hi78(2,3) + Ai07(1,3)*Hi78(3,3);
Ai08(1,4) = Ai07(1,4) + Ai07(1,1)*Hi78(1,4) + Ai07(1,2)*Hi78(2,4) + Ai07(1,3)*Hi78(3,4);
Ai08(2,1) = Ai07(2,1)*Hi78(1,1) + Ai07(2,2)*Hi78(2,1) + Ai07(2,3)*Hi78(3,1);
Ai08(2,2) = Ai07(2,1)*Hi78(1,2) + Ai07(2,2)*Hi78(2,2) + Ai07(2,3)*Hi78(3,2);
Ai08(2,3) = Ai07(2,1)*Hi78(1,3) + Ai07(2,2)*Hi78(2,3) + Ai07(2,3)*Hi78(3,3);
Ai08(2,4) = Ai07(2,4) + Ai07(2,1)*Hi78(1,4) + Ai07(2,2)*Hi78(2,4) + Ai07(2,3)*Hi78(3,4);
Ai08(3,1) = Ai07(3,1)*Hi78(1,1) + Ai07(3,2)*Hi78(2,1) + Ai07(3,3)*Hi78(3,1);
Ai08(3,2) = Ai07(3,1)*Hi78(1,2) + Ai07(3,2)*Hi78(2,2) + Ai07(3,3)*Hi78(3,2);
Ai08(3,3) = Ai07(3,1)*Hi78(1,3) + Ai07(3,2)*Hi78(2,3) + Ai07(3,3)*Hi78(3,3);
Ai08(3,4) = Ai07(3,4) + Ai07(3,1)*Hi78(1,4) + Ai07(3,2)*Hi78(2,4) + Ai07(3,3)*Hi78(3,4);

% Ahmat0(1:3,1:4) = Hi00;
% Ahmat0(4,4) = 1;
Ahmat(1,1:3,1:4) = Ai02;
Ahmat(1,4,4) = 1;
Ahmat(2,1:3,1:4) = Ai03;
Ahmat(2,4,4) = 1;
Ahmat(3,1:3,1:4) = Ai04;
Ahmat(3,4,4) = 1;
Ahmat(4,1:3,1:4) = Ai05;
Ahmat(4,4,4) = 1;
Ahmat(5,1:3,1:4) = Ai07;
Ahmat(5,4,4) = 1;
Ahmat(6,1:3,1:4) = Ai08;
Ahmat(6,4,4) = 1;

% joint ID: 1 
Xorigin(1,:) = Ai01(:,4)';
Xaxis(1,:)   = Ai01(:,3)';
% joint ID: 2 
Xorigin(2,:) = Ai02(:,4)';
Xaxis(2,:) = Ai02(:,3)';
% joint ID: 3 
Xorigin(3,:) = Ai03(:,4)';
Xaxis(3,:) = Ai03(:,3)';
% joint ID: 4 
Xorigin(4,:) = Ai04(:,4)';
Xaxis(4,:) = Ai04(:,3)';
% joint ID: 5 
Xorigin(5,:) = Ai05(:,4)';
Xaxis(5,:) = Ai05(:,3)';
% joint ID: 6 
Xorigin(6,:) = Ai06(:,4)';
Xaxis(6,:) = Ai06(:,3)';
% joint ID: 7 
Xorigin(7,:) = Ai07(:,4)';
Xaxis(7,:) = Ai07(:,3)';

% NOTES:
Xlink(1,:) = Xorigin(1,:);
for i = 2:6
    Xlink(i,:) = Ahmat(i,1:3,4)';
end

% Xmcog0(1) = link0.mcm(1)*Hi00(1,1) + link0.mcm(2)*Hi00(1,2) + link0.mcm(3)*Hi00(1,3) + link0.m*Hi00(1,4);
% Xmcog0(2) = link0.mcm(1)*Hi00(2,1) + link0.mcm(2)*Hi00(2,2) + link0.mcm(3)*Hi00(2,3) + link0.m*Hi00(2,4);
% Xmcog0(3) = link0.mcm(1)*Hi00(3,1) + link0.mcm(2)*Hi00(3,2) + link0.mcm(3)*Hi00(3,3) + link0.m*Hi00(3,4);
% Xmcog(1,1) = links(1).mcm(1)*Ai01(1,1) + links(1).mcm(2)*Ai01(1,2) + links(1).mcm(3)*Ai01(1,3) + links(1).m*Ai01(1,4);
% Xmcog(1,2) = links(1).mcm(1)*Ai01(2,1) + links(1).mcm(2)*Ai01(2,2) + links(1).mcm(3)*Ai01(2,3) + links(1).m*Ai01(2,4);
% Xmcog(1,3) = links(1).mcm(1)*Ai01(3,1) + links(1).mcm(2)*Ai01(3,2) + links(1).mcm(3)*Ai01(3,3) + links(1).m*Ai01(3,4);
% Xmcog(2,1) = links(2).mcm(1)*Ai02(1,1) + links(2).mcm(2)*Ai02(1,2) + links(2).mcm(3)*Ai02(1,3) + links(2).m*Ai02(1,4);
% Xmcog(2,2) = links(2).mcm(1)*Ai02(2,1) + links(2).mcm(2)*Ai02(2,2) + links(2).mcm(3)*Ai02(2,3) + links(2).m*Ai02(2,4);
% Xmcog(2,3) = links(2).mcm(1)*Ai02(3,1) + links(2).mcm(2)*Ai02(3,2) + links(2).mcm(3)*Ai02(3,3) + links(2).m*Ai02(3,4);
% Xmcog(3,1) = links(3).mcm(1)*Ai03(1,1) + links(3).mcm(2)*Ai03(1,2) + links(3).mcm(3)*Ai03(1,3) + links(3).m*Ai03(1,4);
% Xmcog(3,2) = links(3).mcm(1)*Ai03(2,1) + links(3).mcm(2)*Ai03(2,2) + links(3).mcm(3)*Ai03(2,3) + links(3).m*Ai03(2,4);
% Xmcog(3,3) = links(3).mcm(1)*Ai03(3,1) + links(3).mcm(2)*Ai03(3,2) + links(3).mcm(3)*Ai03(3,3) + links(3).m*Ai03(3,4);
% Xmcog(4,1) = links(4).mcm(1)*Ai04(1,1) + links(4).mcm(2)*Ai04(1,2) + links(4).mcm(3)*Ai04(1,3) + links(4).m*Ai04(1,4);
% Xmcog(4,2) = links(4).mcm(1)*Ai04(2,1) + links(4).mcm(2)*Ai04(2,2) + links(4).mcm(3)*Ai04(2,3) + links(4).m*Ai04(2,4);
% Xmcog(4,3) = links(4).mcm(1)*Ai04(3,1) + links(4).mcm(2)*Ai04(3,2) + links(4).mcm(3)*Ai04(3,3) + links(4).m*Ai04(3,4);
% Xmcog(5,1) = links(5).mcm(1)*Ai05(1,1) + links(5).mcm(2)*Ai05(1,2) + links(5).mcm(3)*Ai05(1,3) + links(5).m*Ai05(1,4);
% Xmcog(5,2) = links(5).mcm(1)*Ai05(2,1) + links(5).mcm(2)*Ai05(2,2) + links(5).mcm(3)*Ai05(2,3) + links(5).m*Ai05(2,4);
% Xmcog(5,3) = links(5).mcm(1)*Ai05(3,1) + links(5).mcm(2)*Ai05(3,2) + links(5).mcm(3)*Ai05(3,3) + links(5).m*Ai05(3,4);
% Xmcog(6,1) = links(6).mcm(1)*Ai06(1,1) + links(6).mcm(2)*Ai06(1,2) + links(6).mcm(3)*Ai06(1,3) + links(6).m*Ai06(1,4);
% Xmcog(6,2) = links(6).mcm(1)*Ai06(2,1) + links(6).mcm(2)*Ai06(2,2) + links(6).mcm(3)*Ai06(2,3) + links(6).m*Ai06(2,4);
% Xmcog(6,3) = links(6).mcm(1)*Ai06(3,1) + links(6).mcm(2)*Ai06(3,2) + links(6).mcm(3)*Ai06(3,3) + links(6).m*Ai06(3,4);
% Xmcog(7,1) = links(7).mcm(1)*Ai07(1,1) + links(7).mcm(2)*Ai07(1,2) + links(7).mcm(3)*Ai07(1,3) + links(7).m*Ai07(1,4);
% Xmcog(7,2) = links(7).mcm(1)*Ai07(2,1) + links(7).mcm(2)*Ai07(2,2) + links(7).mcm(3)*Ai07(2,3) + links(7).m*Ai07(2,4);
% Xmcog(7,3) = links(7).mcm(1)*Ai07(3,1) + links(7).mcm(2)*Ai07(3,2) + links(7).mcm(3)*Ai07(3,3) + links(7).m*Ai07(3,4);
% 
