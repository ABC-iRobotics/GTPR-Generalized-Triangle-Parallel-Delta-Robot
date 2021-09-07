
clear; clc;
close all;




%% Classic Delta:
aA = 0.5;   aB = 0.5;   aC = 0.5; % [m] Base Triangle Parameters 
eA = 0.28;  eB = 0.28;  eC = 0.28;  % [m] Work Triangle Parameters
LA = 0.57;  LB = 0.57;  LC = 0.57; % [m] Upper Arm lengths
lA = 0.9;   lB = 0.9;   lC = 0.9; % [m] Lower arm lengths
alphaAA = 0;    alphaAB = deg2rad(120);     alphaAC = deg2rad(240); % [deg] Angles between the arms, as seen on the base/work triangle

%% GTPR Parameters:
% aA = 0.4;   aB = 0.45;   aC = 0.5; % [m] Base Triangle Parameters 
% eA = 0.15;  eB = 0.175;  eC = 0.2;  % [m] Work Triangle Parameters
% LA = 0.5;  LB = 0.57;  LC = 0.64; % [m] Upper Arm lengths
% lA = 0.85;  lB = 0.9;   lC = 0.95; % [m] Lower arm lengths
% alphaAA = 0;    alphaAB = deg2rad(100);     alphaAC = deg2rad(200); % [deg] Angles between the arms, as seen on the base/work triangle

%% GTPR Parameters2:
% aA = 0.5;   aB = 0.6;   aC = 0.5; % [m] Base Triangle Parameters 
% eA = 0.15;  eB = 0.175;  eC = 0.2;  % [m] Work Triangle Parameters
% LA = 0.4;  LB = 0.5;  LC = 0.64; % [m] Upper Arm lengths
% lA = 0.85;  lB = 0.9;   lC = 1; % [m] Lower arm lengths
% alphaAA = 0;    alphaAB = deg2rad(100);     alphaAC = deg2rad(200); % [deg] Angles between the arms, as seen on the base/work triangle

%% GTPR Parameters3:
% aA = 0.4;   aB = 0.5;   aC = 0.6; % [m] Base Triangle Parameters 
% eA = 0.08;  eB = 0.18;  eC = 0.28;  % [m] Work Triangle Parameters
% LA = 0.47;  LB = 0.77;  LC = 0.57; % [m] Upper Arm lengths
% lA = 1.09;  lB = 0.9;   lC = 0.7; % [m] Lower arm lengths
% alphaAA = 0;    alphaAB = deg2rad(100);     alphaAC = deg2rad(200); % [deg] Angles between the arms, as seen on the base/work triangle

%% Modeling Parameters:
% Dimension Parameters
RLA = 0.01; RLB = 0.01; RLC = 0.01; %[m] radius of upper arm
RlA = 0.005; RlB = 0.005; RlC = 0.005; %[m] radius of lower arm
m_unaccounted = 1e-16; % [kg] Mass which should not be taken into consideration
lCrossLink = 0.05; %[m] Cross link length
mCrossLink = m_unaccounted; % Cross Link Weight
maxBaseRadius = max([aA, aB, aC]);
maxTCPRadius = max([eA, eB, eC]);

hTCP = 0.002; %[m] work triangle height
r_endEffector = 0.01; %[m] end effector radius
l_endEffector = 0.01; %[m] end effector length

% Masses and Inertia:
density_steel = 8050; % [kg/m^3]
density_aluminium = 2700; % [kg/m^3]

density_upperArm = density_steel;

% density_LowerArm = density_steel;
% density_workTriangle = density_steel;

density_LowerArm = density_aluminium;
density_workTriangle = density_aluminium;

m_unaccounted = 1e-16;% vrep can only handle this as the smallest mass value.

m_endEffector = m_unaccounted; %[kg] end effector mass

% volumes:
VLA = pi*(RLA^2)*LA;
VLB = pi*(RLB^2)*LB;
VLC = pi*(RLC^2)*LC;

VlA = pi*(RlA^2)*lA;
VlB = pi*(RlB^2)*lB;
VlC = pi*(RlC^2)*lC;

mLA = density_upperArm*VLA;
mLB = density_upperArm*VLB;
mLC = density_upperArm*VLC;

mlA = density_LowerArm*VlA;
mlB = density_LowerArm*VlB;
mlC = density_LowerArm*VlC;

ILA = zeros(3);
ILB = zeros(3);
ILC = zeros(3);

IlA = zeros(3);
IlB = zeros(3);
IlC = zeros(3);

ILA(1,1) = (mLA/12) * (3*RLA^2+LA^2);
ILA(2,2) = (mLA/12) * (3*RLA^2+LA^2);
ILA(3,3) = (mLA*RLA^2)/2;

ILB(1,1) = (mLB/12) * (3*RLB^2+LB^2);
ILB(2,2) = (mLB/12) * (3*RLB^2+LB^2);
ILB(3,3) = (mLB*RLB^2)/2;

ILC(1,1) = (mLC/12) * (3*RLC^2+LC^2);
ILC(2,2) = (mLC/12) * (3*RLC^2+LC^2);
ILC(3,3) = (mLC*RLC^2)/2;

IlA(1,1) = (mlA/12) * (3*RlA^2+lA^2);
IlA(2,2) = (mlA/12) * (3*RlA^2+lA^2);
IlA(3,3) = (mlA*RlA^2)/2;

IlB(1,1) = (mlB/12) * (3*RlB^2+lB^2);
IlB(2,2) = (mlB/12) * (3*RlB^2+lB^2);
IlB(3,3) = (mlB*RlB^2)/2;

IlC(1,1) = (mlC/12) * (3*RlC^2+lC^2);
IlC(2,2) = (mlC/12) * (3*RlC^2+lC^2);
IlC(3,3) = (mlC*RlC^2)/2;

IlA2 = IlA; IlB2 = IlB; IlC2 = IlC;

% IlA2 = zeros(3); IlB2 = zeros(3); IlC2 = zeros(3);

% mLA = 0.5; mLB = 0.5; mLC = 0.5; % Upper Arm Weights
% mlA = 0.3; mlB = 0.3; mlC = 0.3; % Lower arm rod weights (Note: 1 rod Mass, parallelogram has 2 rods)
% mTCP = 0.8; % Work Triangle Mass

% ILA = diag([0.0839583, 0.0839583, 0.00125]); ILB = ILA; ILC = ILA;
% IlA = diag([0.0839583, 0.0839583, 0.00125]); IlB = IlA; IlC = IlA;
% IlA2 = IlA; IlB2 = IlA2; IlC2 = IlA2;
% IlA2 = zeros(3); IlB2 = IlA2; IlC2 = IlA2


% Environment Parameters
g = -9.81; % Gravity along Z axis


%% GTPR Class Initialization:
gtpr = GTPR(aA,aB,aC,LA,LB,LC,lA,lB,lC,eA,eB,eC,alphaAB,alphaAC);

gtpr.directGeometry([0;0;0]);

% work triangle calculation parameters:
a = norm(gtpr.A3-gtpr.B3);
b = norm(gtpr.A3-gtpr.C3);
c = norm(gtpr.C3-gtpr.B3);

VTCP = (1/4)*(hTCP)*sqrt(-a^4 + 2*(a*b)^2 + 2*(a*c)^2 - b^4 + 2*(b*c)^2 - c^4);
mTCP = density_workTriangle * VTCP;


%set weights and inertias:
gtpr.setDynamicParameters(mLA, mLB, mLC, mlA, mlB, mlC, mTCP, ILA, ILB, ILC, IlA2, IlB2, IlC2, g);

%% Initial Parameters:
q = deg2rad([60;60;60]); % for classic delta
% q = deg2rad([45;45;45]); % for classic delta

t_step = 0.001; % Simulation Step Time

v_start = 0; % Initial Velocity
a_start = 0; % Initial Acceleration



% Path Parameters:
gtpr.directGeometry(q);

% Circular Path:
% Path Parameters:
NoOfCircles = 3;
q0_center = q;
radi = 0.05;
normalToPlane = [0;0;1];

phi0 = 0;
omega0 = 0; % initial angular velocity
omega_max = 15; % [rad/s]
beta = 100; % [rad/s^2]
[t, tau_des_c, ddq_des_c, dq_des_c, q_des_c, CirclePath_c, v, a, J11_c, J12_c, J13_c, J21_c, J22_c, J23_c, J31_c, J32_c, J33_c] = gtpr.CircularMotion_Dynamic(NoOfCircles, q0_center, radi, normalToPlane, phi0, omega0, omega_max, beta, t_step);

path_des_c = CirclePath_c;

%% Simulation
t_end = t(end);
q_init = q_des_c(:,1);
dq_init = dq_des_c(:,1);

tau_des.time = t;

tau_des.signals.values = tau_des_c';
tau_des.signals.dimensions = 3;

q_des.time = t;
q_des.signals.values = q_des_c(:,1:length(t))';
q_des_signals.dimensions = 3;

dq_des.time = t;
dq_des.signals.values = dq_des_c(:,1:length(t))';
dq_des_signals.dimensions = 3;

ddq_des.time = t;
ddq_des.signals.values = ddq_des_c(:,1:length(t))';
ddq_des_signals.dimensions = 3;

path_des.time = t;
path_des.signals.values = path_des_c';
path_des_signals.dimensions = 3;

a_des.time = t;
a_des.signals.values = a';
a_des_signals.dimensions = 3;

Fext_des.time = t;
Fext_des.signals.values = zeros(length(t),3);
Fext_des.signals.dimensions = 3;

v_des.time = t;
v_des.signals.values = v';
v_des.signals.dimensions = 3;

J11.time = t;
J11.signals.values = J11_c';
J11.signals.dimensions = 1;

J12.time = t;
J12.signals.values = J12_c';
J12.signals.dimensions = 1;

J13.time = t;
J13.signals.values = J13_c';
J13.signals.dimensions = 1;

J21.time = t;
J21.signals.values = J21_c';
J21.signals.dimensions = 1;

J22.time = t;
J22.signals.values = J22_c';
J22.signals.dimensions = 1;

J23.time = t;
J23.signals.values = J23_c';
J23.signals.dimensions = 1;

J31.time = t;
J31.signals.values = J31_c';
J31.signals.dimensions= 1;

J32.time = t;
J32.signals.values = J32_c';
J32.signals.dimensions = 1;

J33.time = t;
J33.signals.values = J33_c';
J33.signals.dimensions = 1;

KpA = 0;
KpB = 0;
KpC = 0;

KiA = 0;
KiB = 0;
KiC = 0;

KvA = 0;
KvB = 0;
KvC = 0;


% Open loop case:
path_des_c = CirclePath_c;
sim("GTPR_OpenLoop");

TCP_Pos_OpenLoop = TCP_Pos;
TCP_Vel_OpenLoop = TCP_Vel;
TCP_Accel_OpenLoop = TCP_Accel;

sm_q_OpenLoop = sm_q;
sm_dq_OpenLoop = sm_dq;

sm_tau_OpenLoop = sm_tau;

rmse_TCP_Pos_OpenLoop = rms(TCP_Pos_OpenLoop.Data - CirclePath_c');

mrmse_TCP_Pos_OpenLoop = mean(rmse_TCP_Pos_OpenLoop);

% Closed loop case:
KpA = 1000000;
KpB = 1000000;
KpC = 1000000;

KiA = 10000;
KiB = 10000;
KiC = 10000;

KvA = 1000;
KvB = 1000;
KvC = 1000;

% sim("GTPR_ClosedLoop");
sim("GTPR_ClosedLoop_PID");

TCP_Pos_ClosedLoop = TCP_Pos;
TCP_Vel_ClosedLoop = TCP_Vel;
TCP_Accel_ClosedLoop = TCP_Accel;

sm_q_ClosedLoop = sm_q;
sm_dq_ClosedLoop = sm_dq;

sm_tau_ClosedLoop = sm_tau;

rmse_TCP_Pos_ClosedLoop = rms(TCP_Pos_ClosedLoop.Data - CirclePath_c');

mrmse_TCP_Pos_OpenLoop
mrmse_TCP_Pos_ClosedLoop = mean(rmse_TCP_Pos_ClosedLoop)



figure(1);
clf;
grid on;
axis equal;
xlabel("X");
ylabel("Y");
zlabel("Z");
hold on;

plot3(path_des_c(1,1), path_des_c(2,1), path_des_c(3,1), 'k*');
plot3(path_des_c(1,end), path_des_c(2,end), path_des_c(3,end), 'ko');
plot3(path_des_c(1,:), path_des_c(2,:), path_des_c(3,:), 'k:');

plot3(TCP_Pos_OpenLoop.Data(:,1), TCP_Pos_OpenLoop.Data(:,2), TCP_Pos_OpenLoop.Data(:,3), 'r');
plot3(TCP_Pos_OpenLoop.Data(1,1), TCP_Pos_OpenLoop.Data(1,2), TCP_Pos_OpenLoop.Data(1,3), 'r*');
plot3(TCP_Pos_OpenLoop.Data(end,1), TCP_Pos_OpenLoop.Data(end,2), TCP_Pos_OpenLoop.Data(end,3), 'ro');

plot3(TCP_Pos_ClosedLoop.Data(:,1), TCP_Pos_ClosedLoop.Data(:,2), TCP_Pos_ClosedLoop.Data(:,3), 'b');
plot3(TCP_Pos_ClosedLoop.Data(1,1), TCP_Pos_ClosedLoop.Data(1,2), TCP_Pos_ClosedLoop.Data(1,3), 'b*');
plot3(TCP_Pos_ClosedLoop.Data(end,1), TCP_Pos_ClosedLoop.Data(end,2), TCP_Pos_ClosedLoop.Data(end,3), 'bo');

view(30,30);


% figure(2);
% clf;
% grid on;
% xlabel("t");
% ylabel("q [rad]");
% hold on;
% 
% plot(t, q_des_c(1,1:length(t)), 'r');
% plot(t, q_des_c(2,1:length(t)), 'g');
% plot(t, q_des_c(3,1:length(t)), 'b');
% 
% plot(sm_q.Time, sm_q.Data(:,1), 'r--');
% plot(sm_q.Time, sm_q.Data(:,2), 'g--');
% plot(sm_q.Time, sm_q.Data(:,3), 'b--');
% 
% title("Actuator Angles");
% legend("q1_calc", "q2_calc", "q3_calc", "q1_sm", "q2_sm", "q3_sm");
% 
% 
% 
% figure(3);
% clf;
% grid on;
% xlabel("t");
% ylabel("dq [rad/s]");
% hold on;
% 
% plot(t, dq_des_c(1,1:length(t)), 'r');
% plot(t, dq_des_c(2,1:length(t)), 'g');
% plot(t, dq_des_c(3,1:length(t)), 'b');
% 
% plot(sm_dq.Time, sm_dq.Data(:,1), 'r--');
% plot(sm_dq.Time, sm_dq.Data(:,2), 'g--');
% plot(sm_dq.Time, sm_dq.Data(:,3), 'b--');
% 
% title("Actuator Angular Velocity");
% legend("dq1_calc", "dq2_calc", "dq3_calc", "dq1_sm", "dq2_sm", "dq3_sm");
% 
% 
% figure(4);
% clf;
% grid on;
% xlabel("t");
% ylabel("TCP Position [m]");
% hold on;
% 
% plot(t, path_des_c(1,1:length(t)), 'r');
% plot(t, path_des_c(2,1:length(t)), 'g');
% plot(t, path_des_c(3,1:length(t)), 'b');
% 
% plot(TCP_Pos.Time, TCP_Pos.Data(:,1), 'r--');
% plot(TCP_Pos.Time, TCP_Pos.Data(:,2), 'g--');
% plot(TCP_Pos.Time, TCP_Pos.Data(:,3), 'b--');
% 
% title("TCP Position");
% legend("X_calc", "Y_calc", "Z_calc", "X_sm", "Y_sm", "Z_sm");
% 
% 

figure(6);
clf;
grid on;
xlabel("t");
ylabel("\tau [Nm]");
hold on;

plot(t, sm_tau_OpenLoop.Data(:,1), 'r');
plot(t, sm_tau_OpenLoop.Data(:,2), 'g');
plot(t, sm_tau_OpenLoop.Data(:,3), 'b');

plot(sm_tau.Time, sm_tau_ClosedLoop.Data(:,1), 'r--');
plot(sm_tau.Time, sm_tau_ClosedLoop.Data(:,2), 'g--');
plot(sm_tau.Time, sm_tau_ClosedLoop.Data(:,3), 'b--');

title("Actuator Torques");
legend("\tau1_calc", "\tau2_calc", "\tau3_calc", "\tau1_sm", "\tau2_sm", "\tau3_sm");

beep();

% 
% figure(5);
% clf;
% grid on;
% xlabel("t");
% ylabel("TCP Velocity [m/s]");
% hold on;
% 
% plot(t, v_des_c(1,:), 'r');
% plot(t, v_des_c(2,:), 'g');
% plot(t, v_des_c(3,:), 'b');
% 
% plot(t, TCP_Vel.Data(:,1), 'r--');
% plot(t, TCP_Vel.Data(:,2), 'g--');
% plot(t, TCP_Vel.Data(:,3), 'b--');
% 
% title("TCP Velocity");
% legend("dX_calc", "dY_calc", "dZ_calc", "dX_sm", "dY_sm", "dZ_sm");
% 
