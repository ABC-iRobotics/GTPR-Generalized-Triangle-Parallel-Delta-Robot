% Based on that stuff


% In notations s stands for symbolic e stenads for equation form of the
% variable
clear;
clear livescript

syms g x y z dx dy dz ddx ddy ddz q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3 
syms eA eB eC LA LB LC lA lB lC aA aB aC alphaAA alphaAB alphaAC 
syms mLA mLB mLC mlA mlB mlC mTCP
syms Fextx Fexty Fextz pi

syms A3x A3y A3z
syms B3x B3y B3z
syms C3x C3y C3z

syms J11 J12 J13
syms J21 J22 J23
syms J31 J32 J33

assumeAlso(alphaAA,'real');

Rz = @(theta)([cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0,0,1]);

syms ILAxx ILAxy ILAxz ILAyy ILAyz ILAzz
syms ILBxx ILBxy ILBxz ILByy ILByz ILBzz
syms ILCxx ILCxy ILCxz ILCyy ILCyz ILCzz

syms IlAxx IlAxy IlAxz IlAyy IlAyz IlAzz
syms IlBxx IlBxy IlBxz IlByy IlByz IlBzz
syms IlCxx IlCxy IlCxz IlCyy IlCyz IlCzz

syms tau1 tau2 tau3
syms lambdaA lambdaB lambdaC % constraints

assumeAlso(aA > 0);
assumeAlso(aB > 0);
assumeAlso(aC > 0);
assumeAlso(LA > 0);
assumeAlso(LB > 0);
assumeAlso(LC > 0);
assumeAlso(lA > 0);
assumeAlso(lB > 0);
assumeAlso(lC > 0);
assumeAlso(eA > 0);
assumeAlso(eB > 0);
assumeAlso(eC > 0);
assumeAlso(mlA > 0);
assumeAlso(mlB > 0);
assumeAlso(mlC > 0);
assumeAlso(mLA > 0);
assumeAlso(mLB > 0);
assumeAlso(mLC > 0);
assumeAlso(mlA > 0);
assumeAlso(mlB > 0);
assumeAlso(mTCP > 0);
assumeAlso(alphaAB > 0);
assumeAlso(alphaAC > 0);
assumeAlso(z ~= 0);
assumeAlso(g ~= 0);

Fext = [Fextx; Fexty; Fextz];

ILA = [ILAxx, ILAxy, ILAxz;
       ILAxy, ILAyy, ILAyz;
       ILAxz, ILAyz, ILAzz];
   
ILB = [ILBxx, ILBxy, ILBxz;
       ILBxy, ILByy, ILByz;
       ILBxz, ILByz, ILBzz];

ILC = [ILCxx, ILCxy, ILCxz;
       ILCxy, ILCyy, ILCyz;
       ILCxz, ILCyz, ILCzz];
   
IlA = [IlAxx, IlAxy, IlAxz;
       IlAxy, IlAyy, IlAyz;
       IlAxz, IlAyz, IlAzz];
   
IlB = [IlBxx, IlBxy, IlBxz;
       IlBxy, IlByy, IlByz;
       IlBxz, IlByz, IlBzz];

IlC = [IlCxx, IlCxy, IlCxz;
       IlCxy, IlCyy, IlCyz;
       IlCxz, IlCyz, IlCzz];

J = [J11 J12 J13;
     J21 J22 J23;
     J31 J32 J33];
   
TCP = [x;y;z];
v = [dx;dy;dz];

q = [q1;q2;q3];
dq = [dq1;dq2;dq3];
ddq = [ddq1;ddq2;ddq3];


names = whos;
for i = 1:length(names)
    if ~strcmp(names(i).name, 'Rz')
        cmd = strcat("assumeAlso(", names(i).name,", 'real');");
        eval(cmd);
    end
end
clear i names cmd;

% alphaAA = 0;

% In world frame:
A1 = Rz(alphaAA)'*[0;-aA;0];
% A2 = Rz(alphaAA)*[0;aA;0]
B1 = Rz(alphaAB)'*[0;-aB;0];
C1 = Rz(alphaAC)'*[0;-aC;0];

A2 = A1 + LA * Rz(alphaAA)'*-[0;cos(q1);sin(q1)];
B2 = B1 + LB * Rz(alphaAB)'*-[0;cos(q2);sin(q2)];
C2 = C1 + LC * Rz(alphaAC)'*-[0;cos(q3);sin(q3)];

A3 = TCP - Rz(alphaAA)'*[0;eA;0];
B3 = TCP - Rz(alphaAB)'*[0;eB;0] ;
C3 = TCP - Rz(alphaAC)'*[0;eC;0];


S12A = Rz(pi/2) * ((A1)/norm(A1));
S12B = Rz(pi/2) * ((B1)/norm(B1));
S12C = Rz(pi/2) * ((C1)/norm(C1));

% Rotational Velocities (upper arm):
wLA = formula(simplify(dq1*S12A(1:3)));
wLB = formula(simplify(dq2*S12B(1:3)));
wLC = formula(simplify(dq3*S12C(1:3)));

% Linear velocities: (upper arm):
vLA = formula(simplify(cross(wLA, A2-A1)));
vLB = formula(simplify(cross(wLB, B2-B1)));
vLC = formula(simplify(cross(wLC, C2-C1)));

% Direction of lower limbs:
elA = formula(simplify((A3-A2)/lA));
elB = formula(simplify((B3-B2)/lB));
elC = formula(simplify((C3-C2)/lC));

% Rotational Velocities (lower limbs):
wlA = formula(simplify((1/lA)*(cross(elA, vLA-v))));
wlB = formula(simplify((1/lB)*(cross(elB, vLB-v))));
wlC = formula(simplify((1/lC)*(cross(elC, vLC-v))));

% Linear velocities (lower limbs) at CoMs:
vlA = formula(simplify((1/2)*(vLA+v)));
vlB = formula(simplify((1/2)*(vLB+v)));
vlC = formula(simplify((1/2)*(vLC+v)));

%%
% ENERGIES:

% Kinetic Energies
% From Roational Inertia:
    % upper arms
TLAw = formula(simplify((1/2)*wLA'*ILA*wLA));
TLBw = formula(simplify((1/2)*wLB'*ILB*wLB));
TLCw = formula(simplify((1/2)*wLC'*ILC*wLC));


% lower arms:
TlAw = formula(simplify((1/2)*wlA'*IlA*wlA));
TlBw = formula(simplify((1/2)*wlB'*IlB*wlB));
TlCw = formula(simplify((1/2)*wlC'*IlC*wlC));

% From linear velocities:
TlAv = formula(simplify((1/2)*vlA'*mlA*vlA));
TlBv = formula(simplify((1/2)*vlB'*mlB*vlB));
TlCv = formula(simplify((1/2)*vlC'*mlC*vlC));

TTCP = formula(simplify((1/2)*v'*mTCP*v));

% Potential Energy:
CoMLA = A1 + ((A2-A1)/2);
CoMLB = B1 + ((B2-B1)/2);
CoMLC = C1 + ((C2-C1)/2);

% VLA = formula(simplify(mLA*-g*A2(3)));
% VLB = formula(simplify(mLB*-g*B2(3)));
% VLC = formula(simplify(mLC*-g*C2(3)));

VLA = formula(simplify(mLA*-g*CoMLA(3)));
VLB = formula(simplify(mLB*-g*CoMLB(3)));
VLC = formula(simplify(mLC*-g*CoMLC(3)));

CoMlA = A2 + ((A3-A2)/2);
CoMlB = B2 + ((B3-B2)/2);
CoMlC = C2 + ((C3-C2)/2);

VlA = formula(simplify(mlA*-g*CoMlA(3)));
VlB = formula(simplify(mlB*-g*CoMlB(3)));
VlC = formula(simplify(mlC*-g*CoMlC(3)));

VTCP = formula(simplify(mTCP*-g*z));


T = formula(TTCP + TLAw + TLBw + TLCw + 2*TlAw + 2*TlAv + 2*TlBw + 2*TlBv + 2*TlCw + 2*TlCv);
V = formula(simplify(VLA + VLB + VLC + 2*VlA + 2*VlB + 2*VlC + VTCP));


L = formula(T-V);

% Constraint Equations:
GammaA = simplify((norm(A2-A3)^2)-lA^2);% = 0
GammaB = simplify((norm(B2-B3)^2)-lB^2);% = 0
GammaC = simplify((norm(C2-C3)^2)-lC^2);% = 0

dGammaAdx = formula(simplify(diff(GammaA, x)));
dGammaBdx = formula(simplify(diff(GammaB, x)));
dGammaCdx = formula(simplify(diff(GammaC, x)));

dGammaAdy = formula(simplify(diff(GammaA, y)));
dGammaBdy = formula(simplify(diff(GammaB, y)));
dGammaCdy = formula(simplify(diff(GammaC, y)));

dGammaAdz = formula(simplify(diff(GammaA, z)));
dGammaBdz = formula(simplify(diff(GammaB, z)));
dGammaCdz = formula(simplify(diff(GammaC, z)));

dGammaAdq1 = formula(simplify(diff(GammaA, q1)));
dGammaBdq1 = formula(simplify(diff(GammaB, q1)));
dGammaCdq1 = formula(simplify(diff(GammaC, q1)));

dGammaAdq2 = formula(simplify(diff(GammaA, q2)));
dGammaBdq2 = formula(simplify(diff(GammaB, q2)));
dGammaCdq2 = formula(simplify(diff(GammaC, q2)));

dGammaAdq3 = formula(simplify(diff(GammaA, q3)));
dGammaBdq3 = formula(simplify(diff(GammaB, q3)));
dGammaCdq3 = formula(simplify(diff(GammaC, q3)));

dGammasXYZ = [dGammaAdx, dGammaBdx, dGammaCdx;
              dGammaAdy, dGammaBdy, dGammaCdy;
              dGammaAdz, dGammaBdz, dGammaCdz];

dGammasq = [dGammaAdq1, dGammaBdq1, dGammaCdq1;
              dGammaAdq2, dGammaBdq2, dGammaCdq2;
              dGammaAdq3, dGammaBdq3, dGammaCdq3];

          

% Derivate the lagrangian according to the generalized cordinates
dLdx = diff(L, x);
dLdy = diff(L, y);
dLdz = diff(L, z);

dLdq1 = diff(L, q1);
dLdq2 = diff(L, q2);
dLdq3 = diff(L, q3);


% Derivate the lagrangian according to the generalized cordinates rate of change
dLddx = diff(L, dx);
dLddy = diff(L, dy);
dLddz = diff(L, dz);

dLddq1 = diff(L, dq1);
dLddq2 = diff(L, dq2);
dLddq3 = diff(L, dq3);

% derivate the above equations (rate of change) w respect to time


% dx
ddLddxdt = subs(dLddx, 'x', str2sym('x(t)'));
ddLddxdt = subs(ddLddxdt, 'y', str2sym('y(t)'));
ddLddxdt = subs(ddLddxdt, 'z', str2sym('z(t)'));

ddLddxdt = subs(ddLddxdt, 'dx', str2sym('diff(x(t),t)'));
ddLddxdt = subs(ddLddxdt, 'dy', str2sym('diff(y(t),t)'));
ddLddxdt = subs(ddLddxdt, 'dz', str2sym('diff(z(t),t)'));

ddLddxdt = subs(ddLddxdt, 'q1', str2sym('q1(t)'));
ddLddxdt = subs(ddLddxdt, 'q2', str2sym('q2(t)'));
ddLddxdt = subs(ddLddxdt, 'q3', str2sym('q3(t)'));

ddLddxdt = subs(ddLddxdt, 'dq1', str2sym('diff(q1(t),t)'));
ddLddxdt = subs(ddLddxdt, 'dq2', str2sym('diff(q2(t),t)'));
ddLddxdt = subs(ddLddxdt, 'dq3', str2sym('diff(q3(t),t)'));

syms t
assumeAlso(t,'real');
ddLddxdt = diff(ddLddxdt, t);
clear t;

for i = 1:2
    ddLddxdt = subs(ddLddxdt, str2sym('diff(x(t), t, t)'), ddx);
    ddLddxdt = subs(ddLddxdt, str2sym('diff(y(t), t, t)'), ddy);
    ddLddxdt = subs(ddLddxdt, str2sym('diff(z(t), t, t)'), ddz);

    ddLddxdt = subs(ddLddxdt, str2sym('diff(x(t), t)'), dx);
    ddLddxdt = subs(ddLddxdt, str2sym('diff(y(t), t)'), dy);
    ddLddxdt = subs(ddLddxdt, str2sym('diff(z(t), t)'), dz);

    ddLddxdt = subs(ddLddxdt, str2sym('x(t)'), x);
    ddLddxdt = subs(ddLddxdt, str2sym('y(t)'), y);
    ddLddxdt = subs(ddLddxdt, str2sym('z(t)'), z);

    ddLddxdt = subs(ddLddxdt, str2sym('diff(q1(t), t, t)'), ddq1);
    ddLddxdt = subs(ddLddxdt, str2sym('diff(q2(t), t, t)'), ddq2);
    ddLddxdt = subs(ddLddxdt, str2sym('diff(q3(t), t, t)'), ddq3);

    ddLddxdt = subs(ddLddxdt, str2sym('diff(q1(t), t)'), dq1);
    ddLddxdt = subs(ddLddxdt, str2sym('diff(q2(t), t)'), dq2);
    ddLddxdt = subs(ddLddxdt, str2sym('diff(q3(t), t)'), dq3);

    ddLddxdt = subs(ddLddxdt, str2sym('q1(t)'), q1);
    ddLddxdt = subs(ddLddxdt, str2sym('q2(t)'), q2);
    ddLddxdt = subs(ddLddxdt, str2sym('q3(t)'), q3);
end

% dy
ddLddydt = subs(dLddy, 'dx', str2sym('diff(x(t),t)'));
ddLddydt = subs(ddLddydt, 'dy', str2sym('diff(y(t),t)'));
ddLddydt = subs(ddLddydt, 'dz', str2sym('diff(z(t),t)'));
	 
ddLddydt = subs(ddLddydt, 'x', str2sym('x(t)'));
ddLddydt = subs(ddLddydt, 'y', str2sym('y(t)'));
ddLddydt = subs(ddLddydt, 'z', str2sym('z(t)'));
	 
ddLddydt = subs(ddLddydt, 'dq1', str2sym('diff(q1(t),t)'));
ddLddydt = subs(ddLddydt, 'dq2', str2sym('diff(q2(t),t)'));
ddLddydt = subs(ddLddydt, 'dq3', str2sym('diff(q3(t),t)'));
	 
ddLddydt = subs(ddLddydt, 'q1', str2sym('q1(t)'));
ddLddydt = subs(ddLddydt, 'q2', str2sym('q2(t)'));
ddLddydt = subs(ddLddydt, 'q3', str2sym('q3(t)'));

syms t
assumeAlso(t,'real');
ddLddydt = diff(ddLddydt, t);
clear t;

for i = 1:2
    ddLddydt = subs(ddLddydt, str2sym('diff(x(t), t, t)'), ddx);
    ddLddydt = subs(ddLddydt, str2sym('diff(y(t), t, t)'), ddy);
    ddLddydt = subs(ddLddydt, str2sym('diff(z(t), t, t)'), ddz);

    ddLddydt = subs(ddLddydt, str2sym('diff(x(t), t)'), dx);
    ddLddydt = subs(ddLddydt, str2sym('diff(y(t), t)'), dy);
    ddLddydt = subs(ddLddydt, str2sym('diff(z(t), t)'), dz);

    ddLddydt = subs(ddLddydt, str2sym('x(t)'), x);
    ddLddydt = subs(ddLddydt, str2sym('y(t)'), y);
    ddLddydt = subs(ddLddydt, str2sym('z(t)'), z);

    ddLddydt = subs(ddLddydt, str2sym('diff(q1(t), t, t)'), ddq1);
    ddLddydt = subs(ddLddydt, str2sym('diff(q2(t), t, t)'), ddq2);
    ddLddydt = subs(ddLddydt, str2sym('diff(q3(t), t, t)'), ddq3);

    ddLddydt = subs(ddLddydt, str2sym('diff(q1(t), t)'), dq1);
    ddLddydt = subs(ddLddydt, str2sym('diff(q2(t), t)'), dq2);
    ddLddydt = subs(ddLddydt, str2sym('diff(q3(t), t)'), dq3);

    ddLddydt = subs(ddLddydt, str2sym('q1(t)'), q1);
    ddLddydt = subs(ddLddydt, str2sym('q2(t)'), q2);
    ddLddydt = subs(ddLddydt, str2sym('q3(t)'), q3);
end

% dz
ddLddzdt = subs(dLddz, 'dx', str2sym('diff(x(t),t)'));
ddLddzdt = subs(ddLddzdt, 'dy', str2sym('diff(y(t),t)'));
ddLddzdt = subs(ddLddzdt, 'dz', str2sym('diff(z(t),t)'));
	 
ddLddzdt = subs(ddLddzdt, 'x', str2sym('x(t)'));
ddLddzdt = subs(ddLddzdt, 'y', str2sym('y(t)'));
ddLddzdt = subs(ddLddzdt, 'z', str2sym('z(t)'));
	 
ddLddzdt = subs(ddLddzdt, 'dq1', str2sym('diff(q1(t),t)'));
ddLddzdt = subs(ddLddzdt, 'dq2', str2sym('diff(q2(t),t)'));
ddLddzdt = subs(ddLddzdt, 'dq3', str2sym('diff(q3(t),t)'));
	 
ddLddzdt = subs(ddLddzdt, 'q1', str2sym('q1(t)'));
ddLddzdt = subs(ddLddzdt, 'q2', str2sym('q2(t)'));
ddLddzdt = subs(ddLddzdt, 'q3', str2sym('q3(t)'));

syms t
assumeAlso(t,'real');
ddLddzdt = diff(ddLddzdt, t);
clear t;

for i = 1:2
    ddLddzdt = subs(ddLddzdt, str2sym('diff(x(t), t, t)'), ddx);
    ddLddzdt = subs(ddLddzdt, str2sym('diff(y(t), t, t)'), ddy);
    ddLddzdt = subs(ddLddzdt, str2sym('diff(z(t), t, t)'), ddz);

    ddLddzdt = subs(ddLddzdt, str2sym('diff(x(t), t)'), dx);
    ddLddzdt = subs(ddLddzdt, str2sym('diff(y(t), t)'), dy);
    ddLddzdt = subs(ddLddzdt, str2sym('diff(z(t), t)'), dz);

    ddLddzdt = subs(ddLddzdt, str2sym('x(t)'), x);
    ddLddzdt = subs(ddLddzdt, str2sym('y(t)'), y);
    ddLddzdt = subs(ddLddzdt, str2sym('z(t)'), z);

    ddLddzdt = subs(ddLddzdt, str2sym('diff(q1(t), t, t)'), ddq1);
    ddLddzdt = subs(ddLddzdt, str2sym('diff(q2(t), t, t)'), ddq2);
    ddLddzdt = subs(ddLddzdt, str2sym('diff(q3(t), t, t)'), ddq3);

    ddLddzdt = subs(ddLddzdt, str2sym('diff(q1(t), t)'), dq1);
    ddLddzdt = subs(ddLddzdt, str2sym('diff(q2(t), t)'), dq2);
    ddLddzdt = subs(ddLddzdt, str2sym('diff(q3(t), t)'), dq3);

    ddLddzdt = subs(ddLddzdt, str2sym('q1(t)'), q1);
    ddLddzdt = subs(ddLddzdt, str2sym('q2(t)'), q2);
    ddLddzdt = subs(ddLddzdt, str2sym('q3(t)'), q3);
end

% dq1
ddLddq1dt = subs(dLddq1, 'dx', str2sym('diff(x(t),t)'));
ddLddq1dt = subs(ddLddq1dt, 'dy', str2sym('diff(y(t),t)'));
ddLddq1dt = subs(ddLddq1dt, 'dz', str2sym('diff(z(t),t)'));
	 
ddLddq1dt = subs(ddLddq1dt, 'x', str2sym('x(t)'));
ddLddq1dt = subs(ddLddq1dt, 'y', str2sym('y(t)'));
ddLddq1dt = subs(ddLddq1dt, 'z', str2sym('z(t)'));
	 
ddLddq1dt = subs(ddLddq1dt, 'dq1', str2sym('diff(q1(t),t)'));
ddLddq1dt = subs(ddLddq1dt, 'dq2', str2sym('diff(q2(t),t)'));
ddLddq1dt = subs(ddLddq1dt, 'dq3', str2sym('diff(q3(t),t)'));
	 
ddLddq1dt = subs(ddLddq1dt, 'q1', str2sym('q1(t)'));
ddLddq1dt = subs(ddLddq1dt, 'q2', str2sym('q2(t)'));
ddLddq1dt = subs(ddLddq1dt, 'q3', str2sym('q3(t)'));

syms t
assumeAlso(t,'real');
ddLddq1dt = diff(ddLddq1dt, t);
clear t;

for i = 1:2
    ddLddq1dt = subs(ddLddq1dt, str2sym('diff(x(t), t, t)'), ddx);
    ddLddq1dt = subs(ddLddq1dt, str2sym('diff(y(t), t, t)'), ddy);
    ddLddq1dt = subs(ddLddq1dt, str2sym('diff(z(t), t, t)'), ddz);

    ddLddq1dt = subs(ddLddq1dt, str2sym('diff(x(t), t)'), dx);
    ddLddq1dt = subs(ddLddq1dt, str2sym('diff(y(t), t)'), dy);
    ddLddq1dt = subs(ddLddq1dt, str2sym('diff(z(t), t)'), dz);

    ddLddq1dt = subs(ddLddq1dt, str2sym('x(t)'), x);
    ddLddq1dt = subs(ddLddq1dt, str2sym('y(t)'), y);
    ddLddq1dt = subs(ddLddq1dt, str2sym('z(t)'), z);

    ddLddq1dt = subs(ddLddq1dt, str2sym('diff(q1(t), t, t)'), ddq1);
    ddLddq1dt = subs(ddLddq1dt, str2sym('diff(q2(t), t, t)'), ddq2);
    ddLddq1dt = subs(ddLddq1dt, str2sym('diff(q3(t), t, t)'), ddq3);

    ddLddq1dt = subs(ddLddq1dt, str2sym('diff(q1(t), t)'), dq1);
    ddLddq1dt = subs(ddLddq1dt, str2sym('diff(q2(t), t)'), dq2);
    ddLddq1dt = subs(ddLddq1dt, str2sym('diff(q3(t), t)'), dq3);

    ddLddq1dt = subs(ddLddq1dt, str2sym('q1(t)'), q1);
    ddLddq1dt = subs(ddLddq1dt, str2sym('q2(t)'), q2);
    ddLddq1dt = subs(ddLddq1dt, str2sym('q3(t)'), q3);
end

% dq2
ddLddq2dt = subs(dLddq2, 'dx', str2sym('diff(x(t),t)'));
ddLddq2dt = subs(ddLddq2dt, 'dy', str2sym('diff(y(t),t)'));
ddLddq2dt = subs(ddLddq2dt, 'dz', str2sym('diff(z(t),t)'));
	 
ddLddq2dt = subs(ddLddq2dt, 'x', str2sym('x(t)'));
ddLddq2dt = subs(ddLddq2dt, 'y', str2sym('y(t)'));
ddLddq2dt = subs(ddLddq2dt, 'z', str2sym('z(t)'));
	 
ddLddq2dt = subs(ddLddq2dt, 'dq1', str2sym('diff(q1(t),t)'));
ddLddq2dt = subs(ddLddq2dt, 'dq2', str2sym('diff(q2(t),t)'));
ddLddq2dt = subs(ddLddq2dt, 'dq3', str2sym('diff(q3(t),t)'));
	 
ddLddq2dt = subs(ddLddq2dt, 'q1', str2sym('q1(t)'));
ddLddq2dt = subs(ddLddq2dt, 'q2', str2sym('q2(t)'));
ddLddq2dt = subs(ddLddq2dt, 'q3', str2sym('q3(t)'));

syms t
assumeAlso(t,'real');
ddLddq2dt = diff(ddLddq2dt, t);
clear t;

for i = 1:2
    ddLddq2dt = subs(ddLddq2dt, str2sym('diff(x(t), t, t)'), ddx);
    ddLddq2dt = subs(ddLddq2dt, str2sym('diff(y(t), t, t)'), ddy);
    ddLddq2dt = subs(ddLddq2dt, str2sym('diff(z(t), t, t)'), ddz);

    ddLddq2dt = subs(ddLddq2dt, str2sym('diff(x(t), t)'), dx);
    ddLddq2dt = subs(ddLddq2dt, str2sym('diff(y(t), t)'), dy);
    ddLddq2dt = subs(ddLddq2dt, str2sym('diff(z(t), t)'), dz);

    ddLddq2dt = subs(ddLddq2dt, str2sym('x(t)'), x);
    ddLddq2dt = subs(ddLddq2dt, str2sym('y(t)'), y);
    ddLddq2dt = subs(ddLddq2dt, str2sym('z(t)'), z);

    ddLddq2dt = subs(ddLddq2dt, str2sym('diff(q1(t), t, t)'), ddq1);
    ddLddq2dt = subs(ddLddq2dt, str2sym('diff(q2(t), t, t)'), ddq2);
    ddLddq2dt = subs(ddLddq2dt, str2sym('diff(q3(t), t, t)'), ddq3);

    ddLddq2dt = subs(ddLddq2dt, str2sym('diff(q1(t), t)'), dq1);
    ddLddq2dt = subs(ddLddq2dt, str2sym('diff(q2(t), t)'), dq2);
    ddLddq2dt = subs(ddLddq2dt, str2sym('diff(q3(t), t)'), dq3);

    ddLddq2dt = subs(ddLddq2dt, str2sym('q1(t)'), q1);
    ddLddq2dt = subs(ddLddq2dt, str2sym('q2(t)'), q2);
    ddLddq2dt = subs(ddLddq2dt, str2sym('q3(t)'), q3);
end

% dq3
ddLddq3dt = subs(dLddq3, 'dx', str2sym('diff(x(t),t)'));
ddLddq3dt = subs(ddLddq3dt, 'dy', str2sym('diff(y(t),t)'));
ddLddq3dt = subs(ddLddq3dt, 'dz', str2sym('diff(z(t),t)'));
	 
ddLddq3dt = subs(ddLddq3dt, 'x', str2sym('x(t)'));
ddLddq3dt = subs(ddLddq3dt, 'y', str2sym('y(t)'));
ddLddq3dt = subs(ddLddq3dt, 'z', str2sym('z(t)'));
	 
ddLddq3dt = subs(ddLddq3dt, 'dq1', str2sym('diff(q1(t),t)'));
ddLddq3dt = subs(ddLddq3dt, 'dq2', str2sym('diff(q2(t),t)'));
ddLddq3dt = subs(ddLddq3dt, 'dq3', str2sym('diff(q3(t),t)'));
	 
ddLddq3dt = subs(ddLddq3dt, 'q1', str2sym('q1(t)'));
ddLddq3dt = subs(ddLddq3dt, 'q2', str2sym('q2(t)'));
ddLddq3dt = subs(ddLddq3dt, 'q3', str2sym('q3(t)'));

syms t
assumeAlso(t,'real');
ddLddq3dt = diff(ddLddq3dt, t);
clear t;

for i = 1:2
    ddLddq3dt = subs(ddLddq3dt, str2sym('diff(x(t), t, t)'), ddx);
    ddLddq3dt = subs(ddLddq3dt, str2sym('diff(y(t), t, t)'), ddy);
    ddLddq3dt = subs(ddLddq3dt, str2sym('diff(z(t), t, t)'), ddz);

    ddLddq3dt = subs(ddLddq3dt, str2sym('diff(x(t), t)'), dx);
    ddLddq3dt = subs(ddLddq3dt, str2sym('diff(y(t), t)'), dy);
    ddLddq3dt = subs(ddLddq3dt, str2sym('diff(z(t), t)'), dz);

    ddLddq3dt = subs(ddLddq3dt, str2sym('x(t)'), x);
    ddLddq3dt = subs(ddLddq3dt, str2sym('y(t)'), y);
    ddLddq3dt = subs(ddLddq3dt, str2sym('z(t)'), z);

    ddLddq3dt = subs(ddLddq3dt, str2sym('diff(q1(t), t, t)'), ddq1);
    ddLddq3dt = subs(ddLddq3dt, str2sym('diff(q2(t), t, t)'), ddq2);
    ddLddq3dt = subs(ddLddq3dt, str2sym('diff(q3(t), t, t)'), ddq3);

    ddLddq3dt = subs(ddLddq3dt, str2sym('diff(q1(t), t)'), dq1);
    ddLddq3dt = subs(ddLddq3dt, str2sym('diff(q2(t), t)'), dq2);
    ddLddq3dt = subs(ddLddq3dt, str2sym('diff(q3(t), t)'), dq3);

    ddLddq3dt = subs(ddLddq3dt, str2sym('q1(t)'), q1);
    ddLddq3dt = subs(ddLddq3dt, str2sym('q2(t)'), q2);
    ddLddq3dt = subs(ddLddq3dt, str2sym('q3(t)'), q3);
end

% J jacobian of the system:

% Mxyz1 = Lagrange(L, [x dx ddx y dy ddy z dz ddz]);

Mxyz = [ddLddxdt - dLdx;
        ddLddydt - dLdy;
        ddLddzdt - dLdz];

% Lambdas Part of the equation:
% Mxyz = Lagrange(L, [x dx ddx y dy ddy z dz ddz]);


% lambdas = Mxyz\(dGammasXYZ-Fext); % inv(dGammasXYZ)*(Mxyz-Fext);

% Mddq1 = Lagrange(L, [q1 dq1 ddq1 q2 dq2 ddq2 q3 dq3 ddq3]);

Mddq = [ddLddq1dt - dLdq1;
        ddLddq2dt - dLdq2;
        ddLddq3dt - dLdq3];
    
lambdas = [lambdaA; 
           lambdaB; 
           lambdaC];
    
dGammasq = simplify(dGammasq);

tau = Mddq - dGammasq*lambdas;



% Separating into matrices:

[M, f] = equationsToMatrix(tau, [ddq1;ddq2; ddq3]);
% M = simplify(M);
% % f = simplify(f);
% [C,f2] = myequationsToMatrix(f, [dq1; dq2; dq3]);
% 
% [taug, c] = myequationsToMatrix(f2, [q1; q2; q3]);

save("GTPR_Dynamics_equationsComplete");

fname = "DynamicEquations";

dGammasdxyzfs = func2str(matlabFunction(dGammasXYZ));
Mxyzfs = func2str(matlabFunction(Mxyz));
Mddqfs = func2str(matlabFunction(Mddq));
taufs = func2str(matlabFunction(tau));
dGammasqfs = func2str(matlabFunction(dGammasq));

Mfs = func2str(matlabFunction(M));
ffs = func2str(matlabFunction(f));

fid = fopen(strcat(fname,".m"), 'w+');
fprintf(fid, "obj.dGammasdxyzf = " + dGammasdxyzfs + ";\n" + ...
                "obj.Mxyzf = " + Mxyzfs + ";\n" + ...
                "obj.Mddqf = " + Mddqfs + ";\n" + ...
                "obj.dGammasdqf = " + dGammasqfs + ";\n" + ...
                "obj.tauf = " + taufs + ";\n" + ...
                "obj.Mf = " + Mfs + ";\n" + ...
                "obj.ff = " + ffs + ";\n");
fclose(fid);

beep();