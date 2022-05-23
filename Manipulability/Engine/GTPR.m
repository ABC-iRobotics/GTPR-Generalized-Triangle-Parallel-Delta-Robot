classdef GTPR < handle
    %GTPR Class for Generalize-d Triangle Parallel Robot Calculations
    
    properties (Access = public)
        % Robot Parameters
        
        a1 % Base triangvle "radius"
        a2
        a3
        L1 % Upper arm lengths
        L2
        L3
        l1 % Lower arm lengths
        l2
        l3
        e1 % Work triangle "radius"
        e2
        e3
        alpha12 % Angles enclosed by the arms
        alpha13
        
        % Characteristic Points:
        % Letters denote each serial arms,
        % the numbers different characteristic points
        % v - represents the virtual center points of spheres
        
        O = [0;0;0];

        A1
        A2
        A3 
        
        B1
        B2
        B3
        B1v
        B2v
        B3v
        
        C1
        C2
        C3
        
        W
        
        B1L
        B1R
        B2L
        B2R
        B3L
        B3R
        C1L
        C1R
        C2L
        C2R
        C3L
        C3R
        
        v
        a
        
        theta
        dtheta
        ddtheta
        
        
        % Masses of the system:
        mL1 % Upper arms masse
        mL2
        mL3
        ml1 % Lower arms masse
        ml2
        ml3
        mW % Mass of tcp
        
        % Inertias of the system:
        IL_1 % Upper arms inertia
        IL_2
        IL_3
        Il_1 % Lower arms inertia 
        Il_2
        Il_3
        
        % initial inertia matrices: the inertia matrix axis is aligned with
        % the gtpr-s z axis
        IL0_1 % Upper arms inertia
        IL0_2
        IL0_3
        Il0_1 % Lower arms inertia 
        Il0_2
        Il0_3
        
        g
        
        
        % values of inline functions:
        el1
        el2
        el3
        tau
        dGammaAdtheta1 
        dGammaBdtheta2 
        dGammaCdtheta3
        Mxyz
        dGammasdxyz
        dGammasdtheta
        lambda
        Mddtheta 
        
        % angular velocities relating to 23, and 34 screws (screws located at A2)
        w23_1
        w23_2
        w23_3
        w34_1
        w34_2
        w34_3
        
        % Screws:
        S01_1
		S12_1
		S23_1
		S34_1
		S45_1
		S56_1
		
		S01_2
		S12_2
		S23_2
		S34_2
		S45_2
		S56_2
		
		S01_3
		S12_3
		S23_3
		S34_3
		S45_3
		S56_3
		
        % Screw of upper arm
        SL1
        SL2
        SL3
        
        % Screw of lower arm
        Sl1
        Sl2
        Sl3
        
        ScrewConvention
        
        J_A % The kinematic matrix of the active joints
        J_B % The kinematic matrix of the reciprocal screw
        J_C % Complementary matrix of acceleration
        J % The Jacobian Matrix for velocity
		Jinv % the Inverse Jacobian matrix of velocity
        
        
        J_arm1
        J_arm2
        J_arm3
        phi1
        phi2
        phi3
        
        
        phi
        beta
        d
        e_perp
        m_perp
        p
        XiL
        SlXf
        
        RL1
        RL2
        RL3
        Rl1
        Rl2
        Rl3
        
        Velocity_EllipsePts
        
        eigVecs
        eigValues
        
        mu1
        mu2
        mu3
        
        zeroTolerance % in case of imaginary solution: if the imaginary part is below this value, it is considered to be real
        
        CoM
        
        %dynamic Equations:
        dGammasdxyzf = @(L1,L2,L3,a1,a2,a3,alpha11,alpha12,alpha13,e1,e2,e3,theta1,theta2,theta3,x,y,z)reshape([x.*2.0+a1.*sin(alpha11).*2.0-e1.*sin(alpha11).*2.0+L1.*cos(theta1).*sin(alpha11).*2.0,y.*2.0+a1.*cos(alpha11).*2.0-e1.*cos(alpha11).*2.0+L1.*cos(alpha11).*cos(theta1).*2.0,z.*2.0+L1.*sin(theta1).*2.0,x.*2.0+a2.*sin(alpha12).*2.0-e2.*sin(alpha12).*2.0+L2.*cos(theta2).*sin(alpha12).*2.0,y.*2.0+a2.*cos(alpha12).*2.0-e2.*cos(alpha12).*2.0+L2.*cos(alpha12).*cos(theta2).*2.0,z.*2.0+L2.*sin(theta2).*2.0,x.*2.0+a3.*sin(alpha13).*2.0-e3.*sin(alpha13).*2.0+L3.*cos(theta3).*sin(alpha13).*2.0,y.*2.0+a3.*cos(alpha13).*2.0-e3.*cos(alpha13).*2.0+L3.*cos(alpha13).*cos(theta3).*2.0,z.*2.0+L3.*sin(theta3).*2.0],[3,3]);
        Mxyzf = @(IlAxx,IlAxy,IlAxz,IlAyy,IlAyz,IlAzz,IlBxx,IlBxy,IlBxz,IlByy,IlByz,IlBzz,IlCxx,IlCxy,IlCxz,IlCyy,IlCyz,IlCzz,LA,LB,LC,aA,aB,aC,alphaAA,alphaAB,alphaAC,ddq1,ddq2,ddq3,ddx,ddy,ddz,dq1,dq2,dq3,dx,dy,dz,eA,eB,eC,g,lA,lB,lC,mTCP,mlA,mlB,mlC,pi,q1,q2,q3,x,y,z)[(mlA.*(ddx./2.0+(LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0))./2.0+(LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0))./2.0))./2.0+mlA.*(ddx./4.0+(LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0))./4.0+(LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0))./4.0)+(mlB.*(ddx./2.0+(LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0))./2.0+(LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0))./2.0))./2.0+mlB.*(ddx./4.0+(LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0))./4.0+(LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0))./4.0)+(mlC.*(ddx./2.0+(LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0))./2.0+(LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0))./2.0))./2.0+mlC.*(ddx./4.0+(LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0))./4.0+(LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0))./4.0)+ddx.*mTCP+(((IlAxz.*1.0./lA.^2.*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./2.0-(IlAxy.*1.0./lA.^2.*(z+LA.*sin(q1)))./2.0).*(((dz+LA.*dq1.*cos(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA).*2.0)./lA+(((IlBxz.*1.0./lB.^2.*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./2.0-(IlBxy.*1.0./lB.^2.*(z+LB.*sin(q2)))./2.0).*(((dz+LB.*dq2.*cos(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB).*2.0)./lB+(((IlCxz.*1.0./lC.^2.*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./2.0-(IlCxy.*1.0./lC.^2.*(z+LC.*sin(q3)))./2.0).*(((dz+LC.*dq3.*cos(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC).*2.0)./lC-(((IlAxz.*1.0./lA.^2.*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./2.0-(IlAxy.*1.0./lA.^2.*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./2.0).*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA).*2.0)./lA-(((IlBxz.*1.0./lB.^2.*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./2.0-(IlBxy.*1.0./lB.^2.*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./2.0).*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB).*2.0)./lB-(((IlCxz.*1.0./lC.^2.*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./2.0-(IlCxy.*1.0./lC.^2.*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./2.0).*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC).*2.0)./lC-1.0./lA.^2.*(dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*((IlAzz.*(((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*(-1.0./2.0))./lA+(IlAxz.*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA))./(lA.*2.0)+(IlAyz.*(((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA))./(lA.*2.0)).*2.0-1.0./lB.^2.*(dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*((IlBzz.*(((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*(-1.0./2.0))./lB+(IlBxz.*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB))./(lB.*2.0)+(IlByz.*(((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB))./(lB.*2.0)).*2.0-1.0./lC.^2.*(dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*((IlCzz.*(((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*(-1.0./2.0))./lC+(IlCxz.*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC))./(lC.*2.0)+(IlCyz.*(((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC))./(lC.*2.0)).*2.0-((((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*((IlAyz.*1.0./lA.^2.*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./2.0-(IlAyy.*1.0./lA.^2.*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./2.0).*2.0)./lA-((((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*((IlByz.*1.0./lB.^2.*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./2.0-(IlByy.*1.0./lB.^2.*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./2.0).*2.0)./lB-((((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*((IlCyz.*1.0./lC.^2.*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./2.0-(IlCyy.*1.0./lC.^2.*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./2.0).*2.0)./lC+1.0./lA.^2.*(dz+LA.*dq1.*cos(q1)).*((IlAyz.*(((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*(-1.0./2.0))./lA+(IlAxy.*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA))./(lA.*2.0)+(IlAyy.*(((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA))./(lA.*2.0)).*2.0+1.0./lB.^2.*(dz+LB.*dq2.*cos(q2)).*((IlByz.*(((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*(-1.0./2.0))./lB+(IlBxy.*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB))./(lB.*2.0)+(IlByy.*(((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB))./(lB.*2.0)).*2.0+1.0./lC.^2.*(dz+LC.*dq3.*cos(q3)).*((IlCyz.*(((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*(-1.0./2.0))./lC+(IlCxy.*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC))./(lC.*2.0)+(IlCyy.*(((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC))./(lC.*2.0)).*2.0-(((IlAyz.*1.0./lA.^2.*(dz+LA.*dq1.*cos(q1)))./2.0-(IlAzz.*1.0./lA.^2.*(dy-LA.*dq1.*cos(alphaAA).*sin(q1)))./2.0).*(((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*2.0)./lA-(((IlByz.*1.0./lB.^2.*(dz+LB.*dq2.*cos(q2)))./2.0-(IlBzz.*1.0./lB.^2.*(dy-LB.*dq2.*cos(alphaAB).*sin(q2)))./2.0).*(((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*2.0)./lB-(((IlCyz.*1.0./lC.^2.*(dz+LC.*dq3.*cos(q3)))./2.0-(IlCzz.*1.0./lC.^2.*(dy-LC.*dq3.*cos(alphaAC).*sin(q3)))./2.0).*(((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*2.0)./lC-1.0./lA.^2.*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*((IlAzz.*(((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*(-1.0./2.0))./lA+(IlAxz.*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA))./(lA.*2.0)+(IlAyz.*(((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA))./(lA.*2.0)).*2.0-1.0./lB.^2.*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*((IlBzz.*(((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*(-1.0./2.0))./lB+(IlBxz.*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB))./(lB.*2.0)+(IlByz.*(((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB))./(lB.*2.0)).*2.0-1.0./lC.^2.*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*((IlCzz.*(((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*(-1.0./2.0))./lC+(IlCxz.*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC))./(lC.*2.0)+(IlCyz.*(((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC))./(lC.*2.0)).*2.0-(((IlAzz.*1.0./lA.^2.*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./2.0-(IlAyz.*1.0./lA.^2.*(z+LA.*sin(q1)))./2.0).*(-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA+((ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA+((-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA).*2.0)./lA-(((IlBzz.*1.0./lB.^2.*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./2.0-(IlByz.*1.0./lB.^2.*(z+LB.*sin(q2)))./2.0).*(-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB+((ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB+((-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB).*2.0)./lB-(((IlCzz.*1.0./lC.^2.*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./2.0-(IlCyz.*1.0./lC.^2.*(z+LC.*sin(q3)))./2.0).*(-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC+((ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC+((-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC).*2.0)./lC+(((IlAyz.*1.0./lA.^2.*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./2.0-(IlAyy.*1.0./lA.^2.*(z+LA.*sin(q1)))./2.0).*(((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA-((dz+LA.*dq1.*cos(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA).*2.0)./lA+(((IlByz.*1.0./lB.^2.*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./2.0-(IlByy.*1.0./lB.^2.*(z+LB.*sin(q2)))./2.0).*(((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB-((dz+LB.*dq2.*cos(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB).*2.0)./lB+(((IlCyz.*1.0./lC.^2.*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./2.0-(IlCyy.*1.0./lC.^2.*(z+LC.*sin(q3)))./2.0).*(((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC-((dz+LC.*dq3.*cos(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC).*2.0)./lC+1.0./lA.^2.*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*((IlAyz.*(((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*(-1.0./2.0))./lA+(IlAxy.*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA))./(lA.*2.0)+(IlAyy.*(((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA))./(lA.*2.0)).*2.0+1.0./lB.^2.*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*((IlByz.*(((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*(-1.0./2.0))./lB+(IlBxy.*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB))./(lB.*2.0)+(IlByy.*(((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB))./(lB.*2.0)).*2.0+1.0./lC.^2.*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*((IlCyz.*(((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*(-1.0./2.0))./lC+(IlCxy.*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC))./(lC.*2.0)+(IlCyy.*(((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC))./(lC.*2.0)).*2.0+(((IlAxy.*1.0./lA.^2.*(dz+LA.*dq1.*cos(q1)))./2.0-(IlAxz.*1.0./lA.^2.*(dy-LA.*dq1.*cos(alphaAA).*sin(q1)))./2.0).*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA).*2.0)./lA+(((IlBxy.*1.0./lB.^2.*(dz+LB.*dq2.*cos(q2)))./2.0-(IlBxz.*1.0./lB.^2.*(dy-LB.*dq2.*cos(alphaAB).*sin(q2)))./2.0).*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB).*2.0)./lB+(((IlCxy.*1.0./lC.^2.*(dz+LC.*dq3.*cos(q3)))./2.0-(IlCxz.*1.0./lC.^2.*(dy-LC.*dq3.*cos(alphaAC).*sin(q3)))./2.0).*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC).*2.0)./lC+((((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*((IlAyy.*1.0./lA.^2.*(dz+LA.*dq1.*cos(q1)))./2.0-(IlAyz.*1.0./lA.^2.*(dy-LA.*dq1.*cos(alphaAA).*sin(q1)))./2.0).*2.0)./lA+((((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*((IlByy.*1.0./lB.^2.*(dz+LB.*dq2.*cos(q2)))./2.0-(IlByz.*1.0./lB.^2.*(dy-LB.*dq2.*cos(alphaAB).*sin(q2)))./2.0).*2.0)./lB+((((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*((IlCyy.*1.0./lC.^2.*(dz+LC.*dq3.*cos(q3)))./2.0-(IlCyz.*1.0./lC.^2.*(dy-LC.*dq3.*cos(alphaAC).*sin(q3)))./2.0).*2.0)./lC+1.0./lA.^2.*((IlAxz.*(((dz+LA.*dq1.*cos(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA))./(lA.*2.0)-(IlAzz.*(-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA+((ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA+((-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA))./(lA.*2.0)+(IlAyz.*(((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA-((dz+LA.*dq1.*cos(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA))./(lA.*2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)).*2.0+1.0./lB.^2.*((IlBxz.*(((dz+LB.*dq2.*cos(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB))./(lB.*2.0)-(IlBzz.*(-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB+((ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB+((-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB))./(lB.*2.0)+(IlByz.*(((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB-((dz+LB.*dq2.*cos(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB))./(lB.*2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)).*2.0+1.0./lC.^2.*((IlCxz.*(((dz+LC.*dq3.*cos(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC))./(lC.*2.0)-(IlCzz.*(-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC+((ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC+((-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC))./(lC.*2.0)+(IlCyz.*(((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC-((dz+LC.*dq3.*cos(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC))./(lC.*2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)).*2.0+((((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*((IlAzz.*1.0./lA.^2.*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./2.0-(IlAyz.*1.0./lA.^2.*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./2.0).*2.0)./lA+((((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*((IlBzz.*1.0./lB.^2.*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./2.0-(IlByz.*1.0./lB.^2.*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./2.0).*2.0)./lB+((((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*((IlCzz.*1.0./lC.^2.*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./2.0-(IlCyz.*1.0./lC.^2.*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./2.0).*2.0)./lC-1.0./lA.^2.*(z+LA.*sin(q1)).*((IlAxy.*(((dz+LA.*dq1.*cos(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA))./(lA.*2.0)-(IlAyz.*(-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA+((ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA+((-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA))./(lA.*2.0)+(IlAyy.*(((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA-((dz+LA.*dq1.*cos(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA))./(lA.*2.0)).*2.0-1.0./lB.^2.*(z+LB.*sin(q2)).*((IlBxy.*(((dz+LB.*dq2.*cos(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB))./(lB.*2.0)-(IlByz.*(-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB+((ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB+((-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB))./(lB.*2.0)+(IlByy.*(((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB-((dz+LB.*dq2.*cos(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB))./(lB.*2.0)).*2.0-1.0./lC.^2.*(z+LC.*sin(q3)).*((IlCxy.*(((dz+LC.*dq3.*cos(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC))./(lC.*2.0)-(IlCyz.*(-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC+((ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC+((-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC))./(lC.*2.0)+(IlCyy.*(((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC-((dz+LC.*dq3.*cos(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC))./(lC.*2.0)).*2.0;mlA.*(ddy.*(-1.0./2.0)+(LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1))./2.0+(LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1))./2.0).*(-1.0./2.0)-mlA.*(ddy.*(-1.0./4.0)+(LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1))./4.0+(LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1))./4.0)-(mlB.*(ddy.*(-1.0./2.0)+(LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2))./2.0+(LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2))./2.0))./2.0-mlB.*(ddy.*(-1.0./4.0)+(LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2))./4.0+(LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2))./4.0)-(mlC.*(ddy.*(-1.0./2.0)+(LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3))./2.0+(LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3))./2.0))./2.0-mlC.*(ddy.*(-1.0./4.0)+(LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3))./4.0+(LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3))./4.0)+ddy.*mTCP+(((IlCxz.*1.0./lC.^2.*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./2.0-(IlCxx.*1.0./lC.^2.*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./2.0).*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC).*2.0)./lC+((((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*((IlAyz.*1.0./lA.^2.*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./2.0-(IlAxy.*1.0./lA.^2.*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./2.0).*2.0)./lA+((((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*((IlByz.*1.0./lB.^2.*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./2.0-(IlBxy.*1.0./lB.^2.*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./2.0).*2.0)./lB+((((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*((IlCyz.*1.0./lC.^2.*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./2.0-(IlCxy.*1.0./lC.^2.*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./2.0).*2.0)./lC+1.0./lA.^2.*(dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*((IlAzz.*(((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*(-1.0./2.0))./lA+(IlAxz.*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA))./(lA.*2.0)+(IlAyz.*(((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA))./(lA.*2.0)).*2.0+1.0./lB.^2.*(dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*((IlBzz.*(((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*(-1.0./2.0))./lB+(IlBxz.*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB))./(lB.*2.0)+(IlByz.*(((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB))./(lB.*2.0)).*2.0+1.0./lC.^2.*(dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*((IlCzz.*(((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*(-1.0./2.0))./lC+(IlCxz.*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC))./(lC.*2.0)+(IlCyz.*(((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC))./(lC.*2.0)).*2.0-(((IlAxz.*1.0./lA.^2.*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./2.0-(IlAxx.*1.0./lA.^2.*(z+LA.*sin(q1)))./2.0).*(((dz+LA.*dq1.*cos(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA).*2.0)./lA-(((IlBxz.*1.0./lB.^2.*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./2.0-(IlBxx.*1.0./lB.^2.*(z+LB.*sin(q2)))./2.0).*(((dz+LB.*dq2.*cos(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB).*2.0)./lB-(((IlCxz.*1.0./lC.^2.*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./2.0-(IlCxx.*1.0./lC.^2.*(z+LC.*sin(q3)))./2.0).*(((dz+LC.*dq3.*cos(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC).*2.0)./lC-1.0./lA.^2.*(dz+LA.*dq1.*cos(q1)).*((IlAxz.*(((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*(-1.0./2.0))./lA+(IlAxx.*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA))./(lA.*2.0)+(IlAxy.*(((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA))./(lA.*2.0)).*2.0-1.0./lB.^2.*(dz+LB.*dq2.*cos(q2)).*((IlBxz.*(((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*(-1.0./2.0))./lB+(IlBxx.*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB))./(lB.*2.0)+(IlBxy.*(((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB))./(lB.*2.0)).*2.0-1.0./lC.^2.*(dz+LC.*dq3.*cos(q3)).*((IlCxz.*(((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*(-1.0./2.0))./lC+(IlCxx.*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC))./(lC.*2.0)+(IlCxy.*(((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC))./(lC.*2.0)).*2.0+1.0./lA.^2.*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*((IlAzz.*(((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*(-1.0./2.0))./lA+(IlAxz.*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA))./(lA.*2.0)+(IlAyz.*(((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA))./(lA.*2.0)).*2.0+1.0./lB.^2.*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*((IlBzz.*(((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*(-1.0./2.0))./lB+(IlBxz.*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB))./(lB.*2.0)+(IlByz.*(((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB))./(lB.*2.0)).*2.0+1.0./lC.^2.*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*((IlCzz.*(((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*(-1.0./2.0))./lC+(IlCxz.*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC))./(lC.*2.0)+(IlCyz.*(((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC))./(lC.*2.0)).*2.0-(((IlAzz.*1.0./lA.^2.*(dx-LA.*dq1.*sin(alphaAA).*sin(q1)))./2.0-(IlAxz.*1.0./lA.^2.*(dz+LA.*dq1.*cos(q1)))./2.0).*(((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*2.0)./lA-(((IlBzz.*1.0./lB.^2.*(dx-LB.*dq2.*sin(alphaAB).*sin(q2)))./2.0-(IlBxz.*1.0./lB.^2.*(dz+LB.*dq2.*cos(q2)))./2.0).*(((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*2.0)./lB-(((IlCzz.*1.0./lC.^2.*(dx-LC.*dq3.*sin(alphaAC).*sin(q3)))./2.0-(IlCxz.*1.0./lC.^2.*(dz+LC.*dq3.*cos(q3)))./2.0).*(((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*2.0)./lC-1.0./lA.^2.*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*((IlAxz.*(((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*(-1.0./2.0))./lA+(IlAxx.*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA))./(lA.*2.0)+(IlAxy.*(((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA))./(lA.*2.0)).*2.0-1.0./lB.^2.*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*((IlBxz.*(((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*(-1.0./2.0))./lB+(IlBxx.*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB))./(lB.*2.0)+(IlBxy.*(((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB))./(lB.*2.0)).*2.0-1.0./lC.^2.*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*((IlCxz.*(((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*(-1.0./2.0))./lC+(IlCxx.*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC))./(lC.*2.0)+(IlCxy.*(((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC))./(lC.*2.0)).*2.0+(((IlAzz.*1.0./lA.^2.*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./2.0-(IlAxz.*1.0./lA.^2.*(z+LA.*sin(q1)))./2.0).*(-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA+((ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA+((-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA).*2.0)./lA+(((IlBzz.*1.0./lB.^2.*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./2.0-(IlBxz.*1.0./lB.^2.*(z+LB.*sin(q2)))./2.0).*(-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB+((ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB+((-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB).*2.0)./lB+(((IlCzz.*1.0./lC.^2.*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./2.0-(IlCxz.*1.0./lC.^2.*(z+LC.*sin(q3)))./2.0).*(-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC+((ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC+((-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC).*2.0)./lC+(((IlAxz.*1.0./lA.^2.*(dx-LA.*dq1.*sin(alphaAA).*sin(q1)))./2.0-(IlAxx.*1.0./lA.^2.*(dz+LA.*dq1.*cos(q1)))./2.0).*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA).*2.0)./lA+(((IlBxz.*1.0./lB.^2.*(dx-LB.*dq2.*sin(alphaAB).*sin(q2)))./2.0-(IlBxx.*1.0./lB.^2.*(dz+LB.*dq2.*cos(q2)))./2.0).*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB).*2.0)./lB+(((IlCxz.*1.0./lC.^2.*(dx-LC.*dq3.*sin(alphaAC).*sin(q3)))./2.0-(IlCxx.*1.0./lC.^2.*(dz+LC.*dq3.*cos(q3)))./2.0).*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC).*2.0)./lC-(((IlAyz.*1.0./lA.^2.*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./2.0-(IlAxy.*1.0./lA.^2.*(z+LA.*sin(q1)))./2.0).*(((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA-((dz+LA.*dq1.*cos(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA).*2.0)./lA-(((IlByz.*1.0./lB.^2.*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./2.0-(IlBxy.*1.0./lB.^2.*(z+LB.*sin(q2)))./2.0).*(((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB-((dz+LB.*dq2.*cos(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB).*2.0)./lB-(((IlCyz.*1.0./lC.^2.*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./2.0-(IlCxy.*1.0./lC.^2.*(z+LC.*sin(q3)))./2.0).*(((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC-((dz+LC.*dq3.*cos(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC).*2.0)./lC-((((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*((IlAzz.*1.0./lA.^2.*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./2.0-(IlAxz.*1.0./lA.^2.*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./2.0).*2.0)./lA-((((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*((IlBzz.*1.0./lB.^2.*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./2.0-(IlBxz.*1.0./lB.^2.*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./2.0).*2.0)./lB-((((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*((IlCzz.*1.0./lC.^2.*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./2.0-(IlCxz.*1.0./lC.^2.*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./2.0).*2.0)./lC+((((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*((IlAyz.*1.0./lA.^2.*(dx-LA.*dq1.*sin(alphaAA).*sin(q1)))./2.0-(IlAxy.*1.0./lA.^2.*(dz+LA.*dq1.*cos(q1)))./2.0).*2.0)./lA+((((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*((IlByz.*1.0./lB.^2.*(dx-LB.*dq2.*sin(alphaAB).*sin(q2)))./2.0-(IlBxy.*1.0./lB.^2.*(dz+LB.*dq2.*cos(q2)))./2.0).*2.0)./lB+((((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*((IlCyz.*1.0./lC.^2.*(dx-LC.*dq3.*sin(alphaAC).*sin(q3)))./2.0-(IlCxy.*1.0./lC.^2.*(dz+LC.*dq3.*cos(q3)))./2.0).*2.0)./lC-1.0./lA.^2.*((IlAxz.*(((dz+LA.*dq1.*cos(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA))./(lA.*2.0)-(IlAzz.*(-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA+((ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA+((-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA))./(lA.*2.0)+(IlAyz.*(((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA-((dz+LA.*dq1.*cos(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA))./(lA.*2.0)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)).*2.0-1.0./lB.^2.*((IlBxz.*(((dz+LB.*dq2.*cos(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB))./(lB.*2.0)-(IlBzz.*(-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB+((ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB+((-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB))./(lB.*2.0)+(IlByz.*(((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB-((dz+LB.*dq2.*cos(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB))./(lB.*2.0)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)).*2.0-1.0./lC.^2.*((IlCxz.*(((dz+LC.*dq3.*cos(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC))./(lC.*2.0)-(IlCzz.*(-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC+((ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC+((-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC))./(lC.*2.0)+(IlCyz.*(((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC-((dz+LC.*dq3.*cos(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC))./(lC.*2.0)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)).*2.0+1.0./lA.^2.*(z+LA.*sin(q1)).*((IlAxx.*(((dz+LA.*dq1.*cos(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA))./(lA.*2.0)-(IlAxz.*(-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA+((ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA+((-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA))./(lA.*2.0)+(IlAxy.*(((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA-((dz+LA.*dq1.*cos(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA))./(lA.*2.0)).*2.0+1.0./lB.^2.*(z+LB.*sin(q2)).*((IlBxx.*(((dz+LB.*dq2.*cos(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB))./(lB.*2.0)-(IlBxz.*(-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB+((ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB+((-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB))./(lB.*2.0)+(IlBxy.*(((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB-((dz+LB.*dq2.*cos(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB))./(lB.*2.0)).*2.0+1.0./lC.^2.*(z+LC.*sin(q3)).*((IlCxx.*(((dz+LC.*dq3.*cos(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC))./(lC.*2.0)-(IlCxz.*(-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC+((ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC+((-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC))./(lC.*2.0)+(IlCxy.*(((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC-((dz+LC.*dq3.*cos(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC))./(lC.*2.0)).*2.0+(((IlAxz.*1.0./lA.^2.*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./2.0-(IlAxx.*1.0./lA.^2.*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./2.0).*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA).*2.0)./lA+(((IlBxz.*1.0./lB.^2.*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./2.0-(IlBxx.*1.0./lB.^2.*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./2.0).*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB).*2.0)./lB;ddz.*mTCP-g.*mTCP-g.*mlA-g.*mlB-g.*mlC+(mlA.*(ddz./2.0-(LA.*ddq1.*sin(pi./2.0).*cos(q1))./2.0+(LA.*dq1.^2.*sin(pi./2.0).*sin(q1))./2.0))./2.0+mlA.*(ddz./4.0-(LA.*ddq1.*sin(pi./2.0).*cos(q1))./4.0+(LA.*dq1.^2.*sin(pi./2.0).*sin(q1))./4.0)+(mlB.*(ddz./2.0-(LB.*ddq2.*sin(pi./2.0).*cos(q2))./2.0+(LB.*dq2.^2.*sin(pi./2.0).*sin(q2))./2.0))./2.0+mlB.*(ddz./4.0-(LB.*ddq2.*sin(pi./2.0).*cos(q2))./4.0+(LB.*dq2.^2.*sin(pi./2.0).*sin(q2))./4.0)+(mlC.*(ddz./2.0-(LC.*ddq3.*sin(pi./2.0).*cos(q3))./2.0+(LC.*dq3.^2.*sin(pi./2.0).*sin(q3))./2.0))./2.0+mlC.*(ddz./4.0-(LC.*ddq3.*sin(pi./2.0).*cos(q3))./4.0+(LC.*dq3.^2.*sin(pi./2.0).*sin(q3))./4.0)+1.0./lA.^2.*(dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*((IlAxz.*(((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*(-1.0./2.0))./lA+(IlAxx.*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA))./(lA.*2.0)+(IlAxy.*(((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA))./(lA.*2.0)).*2.0+1.0./lB.^2.*(dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*((IlBxz.*(((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*(-1.0./2.0))./lB+(IlBxx.*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB))./(lB.*2.0)+(IlBxy.*(((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB))./(lB.*2.0)).*2.0+1.0./lC.^2.*(dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*((IlCxz.*(((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*(-1.0./2.0))./lC+(IlCxx.*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC))./(lC.*2.0)+(IlCxy.*(((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC))./(lC.*2.0)).*2.0-1.0./lA.^2.*(dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*((IlAyz.*(((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*(-1.0./2.0))./lA+(IlAxy.*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA))./(lA.*2.0)+(IlAyy.*(((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA))./(lA.*2.0)).*2.0-1.0./lB.^2.*(dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*((IlByz.*(((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*(-1.0./2.0))./lB+(IlBxy.*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB))./(lB.*2.0)+(IlByy.*(((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB))./(lB.*2.0)).*2.0-1.0./lC.^2.*(dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*((IlCyz.*(((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*(-1.0./2.0))./lC+(IlCxy.*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC))./(lC.*2.0)+(IlCyy.*(((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC))./(lC.*2.0)).*2.0+(((IlAyz.*1.0./lA.^2.*(dx-LA.*dq1.*sin(alphaAA).*sin(q1)))./2.0-(IlAxz.*1.0./lA.^2.*(dy-LA.*dq1.*cos(alphaAA).*sin(q1)))./2.0).*(((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*2.0)./lA+(((IlByz.*1.0./lB.^2.*(dx-LB.*dq2.*sin(alphaAB).*sin(q2)))./2.0-(IlBxz.*1.0./lB.^2.*(dy-LB.*dq2.*cos(alphaAB).*sin(q2)))./2.0).*(((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*2.0)./lB+(((IlCyz.*1.0./lC.^2.*(dx-LC.*dq3.*sin(alphaAC).*sin(q3)))./2.0-(IlCxz.*1.0./lC.^2.*(dy-LC.*dq3.*cos(alphaAC).*sin(q3)))./2.0).*(((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*2.0)./lC-1.0./lA.^2.*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*((IlAyz.*(((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*(-1.0./2.0))./lA+(IlAxy.*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA))./(lA.*2.0)+(IlAyy.*(((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA))./(lA.*2.0)).*2.0-1.0./lB.^2.*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*((IlByz.*(((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*(-1.0./2.0))./lB+(IlBxy.*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB))./(lB.*2.0)+(IlByy.*(((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB))./(lB.*2.0)).*2.0-1.0./lC.^2.*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*((IlCyz.*(((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*(-1.0./2.0))./lC+(IlCxy.*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC))./(lC.*2.0)+(IlCyy.*(((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC))./(lC.*2.0)).*2.0+(((IlAxz.*1.0./lA.^2.*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./2.0-(IlAyz.*1.0./lA.^2.*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./2.0).*(-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA+((ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA+((-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA).*2.0)./lA+(((IlBxz.*1.0./lB.^2.*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./2.0-(IlByz.*1.0./lB.^2.*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./2.0).*(-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB+((ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB+((-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB).*2.0)./lB+(((IlCxz.*1.0./lC.^2.*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./2.0-(IlCyz.*1.0./lC.^2.*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./2.0).*(-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC+((ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC+((-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC).*2.0)./lC+1.0./lA.^2.*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*((IlAxz.*(((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*(-1.0./2.0))./lA+(IlAxx.*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA))./(lA.*2.0)+(IlAxy.*(((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA))./(lA.*2.0)).*2.0+1.0./lB.^2.*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*((IlBxz.*(((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*(-1.0./2.0))./lB+(IlBxx.*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB))./(lB.*2.0)+(IlBxy.*(((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB))./(lB.*2.0)).*2.0+1.0./lC.^2.*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*((IlCxz.*(((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*(-1.0./2.0))./lC+(IlCxx.*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC))./(lC.*2.0)+(IlCxy.*(((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC))./(lC.*2.0)).*2.0-(((IlAxy.*1.0./lA.^2.*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./2.0-(IlAyy.*1.0./lA.^2.*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./2.0).*(((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA-((dz+LA.*dq1.*cos(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA).*2.0)./lA-(((IlBxy.*1.0./lB.^2.*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./2.0-(IlByy.*1.0./lB.^2.*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./2.0).*(((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB-((dz+LB.*dq2.*cos(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB).*2.0)./lB-(((IlCxy.*1.0./lC.^2.*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./2.0-(IlCyy.*1.0./lC.^2.*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./2.0).*(((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC-((dz+LC.*dq3.*cos(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC).*2.0)./lC-(((IlAxy.*1.0./lA.^2.*(dx-LA.*dq1.*sin(alphaAA).*sin(q1)))./2.0-(IlAxx.*1.0./lA.^2.*(dy-LA.*dq1.*cos(alphaAA).*sin(q1)))./2.0).*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA).*2.0)./lA-(((IlBxy.*1.0./lB.^2.*(dx-LB.*dq2.*sin(alphaAB).*sin(q2)))./2.0-(IlBxx.*1.0./lB.^2.*(dy-LB.*dq2.*cos(alphaAB).*sin(q2)))./2.0).*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB).*2.0)./lB-(((IlCxy.*1.0./lC.^2.*(dx-LC.*dq3.*sin(alphaAC).*sin(q3)))./2.0-(IlCxx.*1.0./lC.^2.*(dy-LC.*dq3.*cos(alphaAC).*sin(q3)))./2.0).*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC).*2.0)./lC+((((dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*((IlAyz.*1.0./lA.^2.*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./2.0-(IlAxz.*1.0./lA.^2.*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./2.0).*2.0)./lA+((((dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*((IlByz.*1.0./lB.^2.*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./2.0-(IlBxz.*1.0./lB.^2.*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./2.0).*2.0)./lB+((((dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*((IlCyz.*1.0./lC.^2.*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./2.0-(IlCxz.*1.0./lC.^2.*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./2.0).*2.0)./lC-((((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*((IlAyy.*1.0./lA.^2.*(dx-LA.*dq1.*sin(alphaAA).*sin(q1)))./2.0-(IlAxy.*1.0./lA.^2.*(dy-LA.*dq1.*cos(alphaAA).*sin(q1)))./2.0).*2.0)./lA-((((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*((IlByy.*1.0./lB.^2.*(dx-LB.*dq2.*sin(alphaAB).*sin(q2)))./2.0-(IlBxy.*1.0./lB.^2.*(dy-LB.*dq2.*cos(alphaAB).*sin(q2)))./2.0).*2.0)./lB-((((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*((IlCyy.*1.0./lC.^2.*(dx-LC.*dq3.*sin(alphaAC).*sin(q3)))./2.0-(IlCxy.*1.0./lC.^2.*(dy-LC.*dq3.*cos(alphaAC).*sin(q3)))./2.0).*2.0)./lC-1.0./lA.^2.*((IlAxx.*(((dz+LA.*dq1.*cos(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA))./(lA.*2.0)-(IlAxz.*(-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA+((ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA+((-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA))./(lA.*2.0)+(IlAxy.*(((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA-((dz+LA.*dq1.*cos(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA))./(lA.*2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)).*2.0-1.0./lB.^2.*((IlBxx.*(((dz+LB.*dq2.*cos(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB))./(lB.*2.0)-(IlBxz.*(-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB+((ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB+((-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB))./(lB.*2.0)+(IlBxy.*(((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB-((dz+LB.*dq2.*cos(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB))./(lB.*2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)).*2.0-1.0./lC.^2.*((IlCxx.*(((dz+LC.*dq3.*cos(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC))./(lC.*2.0)-(IlCxz.*(-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC+((ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC+((-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC))./(lC.*2.0)+(IlCxy.*(((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC-((dz+LC.*dq3.*cos(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC))./(lC.*2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)).*2.0-(((IlAxy.*1.0./lA.^2.*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./2.0-(IlAxx.*1.0./lA.^2.*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./2.0).*(((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA-((z+LA.*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA).*2.0)./lA-(((IlBxy.*1.0./lB.^2.*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./2.0-(IlBxx.*1.0./lB.^2.*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./2.0).*(((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB-((z+LB.*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB).*2.0)./lB-(((IlCxy.*1.0./lC.^2.*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./2.0-(IlCxx.*1.0./lC.^2.*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./2.0).*(((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC-((z+LC.*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC).*2.0)./lC-((((z+LA.*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA-((dz+LA.*dq1.*sin(pi./2.0).*cos(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA).*((IlAyy.*1.0./lA.^2.*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./2.0-(IlAxy.*1.0./lA.^2.*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./2.0).*2.0)./lA-((((z+LB.*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB-((dz+LB.*dq2.*sin(pi./2.0).*cos(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB).*((IlByy.*1.0./lB.^2.*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./2.0-(IlBxy.*1.0./lB.^2.*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./2.0).*2.0)./lB-((((z+LC.*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC-((dz+LC.*dq3.*sin(pi./2.0).*cos(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC).*((IlCyy.*1.0./lC.^2.*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./2.0-(IlCxy.*1.0./lC.^2.*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./2.0).*2.0)./lC+1.0./lA.^2.*((IlAxy.*(((dz+LA.*dq1.*cos(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA))./(lA.*2.0)-(IlAyz.*(-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA+((ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA+((-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA))./(lA.*2.0)+(IlAyy.*(((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./lA-((dz+LA.*dq1.*cos(q1)).*(dx-LA.*dq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA+((dx-LA.*dq1.*sin(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(-ddx+LA.*dq1.^2.*cos(q1).*cos(alphaAA-pi./2.0)+LA.*ddq1.*sin(q1).*cos(alphaAA-pi./2.0)))./lA))./(lA.*2.0)).*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)).*2.0+1.0./lB.^2.*((IlBxy.*(((dz+LB.*dq2.*cos(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB))./(lB.*2.0)-(IlByz.*(-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB+((ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB+((-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB))./(lB.*2.0)+(IlByy.*(((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./lB-((dz+LB.*dq2.*cos(q2)).*(dx-LB.*dq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB+((dx-LB.*dq2.*sin(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(-ddx+LB.*dq2.^2.*cos(q2).*cos(alphaAB-pi./2.0)+LB.*ddq2.*sin(q2).*cos(alphaAB-pi./2.0)))./lB))./(lB.*2.0)).*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)).*2.0+1.0./lC.^2.*((IlCxy.*(((dz+LC.*dq3.*cos(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC))./(lC.*2.0)-(IlCyz.*(-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC+((ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC+((-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC))./(lC.*2.0)+(IlCyy.*(((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./lC-((dz+LC.*dq3.*cos(q3)).*(dx-LC.*dq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC+((dx-LC.*dq3.*sin(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(-ddx+LC.*dq3.^2.*cos(q3).*cos(alphaAC-pi./2.0)+LC.*ddq3.*sin(q3).*cos(alphaAC-pi./2.0)))./lC))./(lC.*2.0)).*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)).*2.0-(((IlAxx.*1.0./lA.^2.*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./2.0-(IlAxy.*1.0./lA.^2.*(x+aA.*sin(alphaAA)-eA.*sin(alphaAA)+LA.*cos(q1).*sin(alphaAA)))./2.0).*(((dz+LA.*dq1.*cos(q1)).*(dy+LA.*dq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((dy-LA.*dq1.*cos(alphaAA).*sin(q1)).*(dz+LA.*dq1.*sin(pi./2.0).*cos(q1)))./lA+((z+LA.*sin(q1)).*(ddy+LA.*dq1.^2.*sin(alphaAA-pi./2.0).*cos(q1)+LA.*ddq1.*sin(alphaAA-pi./2.0).*sin(q1)))./lA-((ddz+LA.*ddq1.*sin(pi./2.0).*cos(q1)-LA.*dq1.^2.*sin(pi./2.0).*sin(q1)).*(y+aA.*cos(alphaAA)-eA.*cos(alphaAA)+LA.*cos(alphaAA).*cos(q1)))./lA).*2.0)./lA-(((IlBxx.*1.0./lB.^2.*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./2.0-(IlBxy.*1.0./lB.^2.*(x+aB.*sin(alphaAB)-eB.*sin(alphaAB)+LB.*cos(q2).*sin(alphaAB)))./2.0).*(((dz+LB.*dq2.*cos(q2)).*(dy+LB.*dq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((dy-LB.*dq2.*cos(alphaAB).*sin(q2)).*(dz+LB.*dq2.*sin(pi./2.0).*cos(q2)))./lB+((z+LB.*sin(q2)).*(ddy+LB.*dq2.^2.*sin(alphaAB-pi./2.0).*cos(q2)+LB.*ddq2.*sin(alphaAB-pi./2.0).*sin(q2)))./lB-((ddz+LB.*ddq2.*sin(pi./2.0).*cos(q2)-LB.*dq2.^2.*sin(pi./2.0).*sin(q2)).*(y+aB.*cos(alphaAB)-eB.*cos(alphaAB)+LB.*cos(alphaAB).*cos(q2)))./lB).*2.0)./lB-(((IlCxx.*1.0./lC.^2.*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./2.0-(IlCxy.*1.0./lC.^2.*(x+aC.*sin(alphaAC)-eC.*sin(alphaAC)+LC.*cos(q3).*sin(alphaAC)))./2.0).*(((dz+LC.*dq3.*cos(q3)).*(dy+LC.*dq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((dy-LC.*dq3.*cos(alphaAC).*sin(q3)).*(dz+LC.*dq3.*sin(pi./2.0).*cos(q3)))./lC+((z+LC.*sin(q3)).*(ddy+LC.*dq3.^2.*sin(alphaAC-pi./2.0).*cos(q3)+LC.*ddq3.*sin(alphaAC-pi./2.0).*sin(q3)))./lC-((ddz+LC.*ddq3.*sin(pi./2.0).*cos(q3)-LC.*dq3.^2.*sin(pi./2.0).*sin(q3)).*(y+aC.*cos(alphaAC)-eC.*cos(alphaAC)+LC.*cos(alphaAC).*cos(q3)))./lC).*2.0)./lC];
        Mddthetaf = @(IL_1xx,IL_1xy,IL_1yy,IL_2xx,IL_2xy,IL_2yy,IL_3xx,IL_3xy,IL_3yy,Il_1xx,Il_1xy,Il_1xz,Il_1yy,Il_1yz,Il_1zz,Il_2xx,Il_2xy,Il_2xz,Il_2yy,Il_2yz,Il_2zz,Il_3xx,Il_3xy,Il_3xz,Il_3yy,Il_3yz,Il_3zz,L1,L2,L3,a1,a2,a3,alpha11,alpha12,alpha13,ddq1,ddq2,ddq3,ddx,ddy,ddz,dtheta1,dtheta2,dtheta3,dx,dy,dz,e1,e2,e3,g,l1,l2,l3,mL1,mL2,mL3,mlA,ml2,ml3,pi,theta1,theta2,theta3,x,y,z)[ddq1.*(IL_1xx.*sin(alpha11-pi./2.0).^2+IL_1yy.*cos(alpha11-pi./2.0).^2+IL_1xy.*sin(alpha11-pi./2.0).*cos(alpha11-pi./2.0).*2.0)+((((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xz.*((L1.*cos(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*cos(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*(-1.0./2.0))./l1+(Il_1yz.*((L1.*cos(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+(L1.*sin(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1-(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1zz.*((L1.*cos(alpha11).*sin(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-(L1.*sin(alpha11).*sin(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1+((((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xz.*((L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dz+L1.*dtheta1.*cos(theta1)))./l1-(L1.*sin(pi./2.0).*cos(theta1).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)+(Il_1yz.*((L1.*sin(pi./2.0).*cos(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dz+L1.*dtheta1.*cos(theta1)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1-(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)-(Il_1zz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1-(((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xz.*(((dz+L1.*dtheta1.*cos(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((z+L1.*sin(theta1)).*(ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)+L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((ddz+L1.*ddq1.*sin(pi./2.0).*cos(theta1)-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)-(Il_1zz.*(-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+((ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)+L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+((-ddx+L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)+L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)+(Il_1yz.*(((ddz+L1.*ddq1.*sin(pi./2.0).*cos(theta1)-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1-((dz+L1.*dtheta1.*cos(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((z+L1.*sin(theta1)).*(-ddx+L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)+L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1))./(l1.*2.0)).*2.0)./l1+(((Il_1xz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xx.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1xy.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*cos(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*cos(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1-(((Il_1xz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xx.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1xy.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dz+L1.*dtheta1.*cos(theta1)))./l1-(L1.*sin(pi./2.0).*cos(theta1).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1-(((Il_1xz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1xx.*((L1.*sin(pi./2.0).*cos(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(z+L1.*sin(theta1)))./l1))./(l1.*2.0)-(Il_1xy.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*sin(pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*(((dz+L1.*dtheta1.*cos(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((z+L1.*sin(theta1)).*(ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)+L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((ddz+L1.*ddq1.*sin(pi./2.0).*cos(theta1)-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1-(((Il_1yz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xy.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1yy.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*cos(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+(L1.*sin(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1-(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*2.0)./l1-(((Il_1yz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xy.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1yy.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*sin(pi./2.0).*cos(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dz+L1.*dtheta1.*cos(theta1)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1-(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*2.0)./l1-((((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1).*((Il_1xx.*((L1.*cos(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*cos(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*(-1.0./2.0))./l1+(Il_1xy.*((L1.*cos(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+(L1.*sin(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1-(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1xz.*((L1.*cos(alpha11).*sin(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-(L1.*sin(alpha11).*sin(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1-((((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1).*((Il_1xx.*((L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dz+L1.*dtheta1.*cos(theta1)))./l1-(L1.*sin(pi./2.0).*cos(theta1).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)+(Il_1xy.*((L1.*sin(pi./2.0).*cos(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dz+L1.*dtheta1.*cos(theta1)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1-(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)-(Il_1xz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1-(((L1.*sin(pi./2.0).*cos(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(z+L1.*sin(theta1)))./l1).*((Il_1xx.*(((dz+L1.*dtheta1.*cos(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((z+L1.*sin(theta1)).*(ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)+L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((ddz+L1.*ddq1.*sin(pi./2.0).*cos(theta1)-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)-(Il_1xz.*(-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+((ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)+L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+((-ddx+L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)+L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)+(Il_1xy.*(((ddz+L1.*ddq1.*sin(pi./2.0).*cos(theta1)-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1-((dz+L1.*dtheta1.*cos(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((z+L1.*sin(theta1)).*(-ddx+L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)+L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1))./(l1.*2.0)).*2.0)./l1-(((Il_1zz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xz.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1yz.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*cos(alpha11).*sin(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-(L1.*sin(alpha11).*sin(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1+(((Il_1zz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xz.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1yz.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1-((((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xy.*((L1.*cos(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*cos(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*(-1.0./2.0))./l1+(Il_1yy.*((L1.*cos(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+(L1.*sin(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1-(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1yz.*((L1.*cos(alpha11).*sin(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-(L1.*sin(alpha11).*sin(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1-((((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xy.*((L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dz+L1.*dtheta1.*cos(theta1)))./l1-(L1.*sin(pi./2.0).*cos(theta1).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)+(Il_1yy.*((L1.*sin(pi./2.0).*cos(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dz+L1.*dtheta1.*cos(theta1)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1-(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)-(Il_1yz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1+(((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*sin(pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xy.*(((dz+L1.*dtheta1.*cos(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((z+L1.*sin(theta1)).*(ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)+L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((ddz+L1.*ddq1.*sin(pi./2.0).*cos(theta1)-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)-(Il_1yz.*(-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+((ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)+L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+((-ddx+L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)+L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)+(Il_1yy.*(((ddz+L1.*ddq1.*sin(pi./2.0).*cos(theta1)-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1-((dz+L1.*dtheta1.*cos(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((z+L1.*sin(theta1)).*(-ddx+L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)+L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1))./(l1.*2.0)).*2.0)./l1+(((Il_1zz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1xz.*((L1.*sin(pi./2.0).*cos(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(z+L1.*sin(theta1)))./l1))./(l1.*2.0)-(Il_1yz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*sin(pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*(-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+((ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)+L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+((-ddx+L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)+L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1-(((Il_1yz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1xy.*((L1.*sin(pi./2.0).*cos(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(z+L1.*sin(theta1)))./l1))./(l1.*2.0)-(Il_1yy.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*sin(pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*(((ddz+L1.*ddq1.*sin(pi./2.0).*cos(theta1)-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1-((dz+L1.*dtheta1.*cos(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((z+L1.*sin(theta1)).*(-ddx+L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)+L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1).*2.0)./l1+(L1.*g.*mL1.*cos(theta1))./2.0+L1.*g.*mlA.*cos(theta1)-(L1.*mlA.*sin(pi./2.0).*cos(theta1).*(ddz./2.0-(L1.*ddq1.*sin(pi./2.0).*cos(theta1))./2.0+(L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1))./2.0))./2.0-L1.*mlA.*sin(pi./2.0).*cos(theta1).*(ddz./4.0-(L1.*ddq1.*sin(pi./2.0).*cos(theta1))./4.0+(L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1))./4.0)+(L1.*mlA.*sin(theta1).*cos(alpha11-pi./2.0).*(ddx./2.0+(L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0))./2.0+(L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0))./2.0))./2.0+L1.*mlA.*sin(theta1).*cos(alpha11-pi./2.0).*(ddx./4.0+(L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0))./4.0+(L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0))./4.0)+(L1.*mlA.*sin(alpha11-pi./2.0).*sin(theta1).*(ddy.*(-1.0./2.0)+(L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1))./2.0+(L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1))./2.0))./2.0+L1.*mlA.*sin(alpha11-pi./2.0).*sin(theta1).*(ddy.*(-1.0./4.0)+(L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1))./4.0+(L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1))./4.0);ddq2.*(IL_2xx.*sin(alpha12-pi./2.0).^2+IL_2yy.*cos(alpha12-pi./2.0).^2+IL_2xy.*sin(alpha12-pi./2.0).*cos(alpha12-pi./2.0).*2.0)+((((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xz.*((L2.*cos(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*cos(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*(-1.0./2.0))./l2+(Il_2yz.*((L2.*cos(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+(L2.*sin(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2-(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2zz.*((L2.*cos(alpha12).*sin(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-(L2.*sin(alpha12).*sin(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2+((((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xz.*((L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dz+L2.*dtheta2.*cos(theta2)))./l2-(L2.*sin(pi./2.0).*cos(theta2).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)+(Il_2yz.*((L2.*sin(pi./2.0).*cos(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dz+L2.*dtheta2.*cos(theta2)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2-(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)-(Il_2zz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2-(((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xz.*(((dz+L2.*dtheta2.*cos(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((z+L2.*sin(theta2)).*(ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)+L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((ddz+L2.*ddq2.*sin(pi./2.0).*cos(theta2)-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)-(Il_2zz.*(-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+((ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)+L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+((-ddx+L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)+L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)+(Il_2yz.*(((ddz+L2.*ddq2.*sin(pi./2.0).*cos(theta2)-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2-((dz+L2.*dtheta2.*cos(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((z+L2.*sin(theta2)).*(-ddx+L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)+L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2))./(l2.*2.0)).*2.0)./l2+(((Il_2xz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xx.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2xy.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*cos(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*cos(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2-(((Il_2xz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xx.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2xy.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dz+L2.*dtheta2.*cos(theta2)))./l2-(L2.*sin(pi./2.0).*cos(theta2).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2-(((Il_2xz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2xx.*((L2.*sin(pi./2.0).*cos(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(z+L2.*sin(theta2)))./l2))./(l2.*2.0)-(Il_2xy.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*sin(pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*(((dz+L2.*dtheta2.*cos(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((z+L2.*sin(theta2)).*(ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)+L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((ddz+L2.*ddq2.*sin(pi./2.0).*cos(theta2)-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2-(((Il_2yz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xy.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2yy.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*cos(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+(L2.*sin(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2-(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*2.0)./l2-(((Il_2yz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xy.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2yy.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*sin(pi./2.0).*cos(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dz+L2.*dtheta2.*cos(theta2)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2-(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*2.0)./l2-((((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2).*((Il_2xx.*((L2.*cos(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*cos(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*(-1.0./2.0))./l2+(Il_2xy.*((L2.*cos(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+(L2.*sin(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2-(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2xz.*((L2.*cos(alpha12).*sin(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-(L2.*sin(alpha12).*sin(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2-((((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2).*((Il_2xx.*((L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dz+L2.*dtheta2.*cos(theta2)))./l2-(L2.*sin(pi./2.0).*cos(theta2).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)+(Il_2xy.*((L2.*sin(pi./2.0).*cos(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dz+L2.*dtheta2.*cos(theta2)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2-(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)-(Il_2xz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2-(((L2.*sin(pi./2.0).*cos(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(z+L2.*sin(theta2)))./l2).*((Il_2xx.*(((dz+L2.*dtheta2.*cos(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((z+L2.*sin(theta2)).*(ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)+L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((ddz+L2.*ddq2.*sin(pi./2.0).*cos(theta2)-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)-(Il_2xz.*(-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+((ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)+L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+((-ddx+L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)+L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)+(Il_2xy.*(((ddz+L2.*ddq2.*sin(pi./2.0).*cos(theta2)-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2-((dz+L2.*dtheta2.*cos(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((z+L2.*sin(theta2)).*(-ddx+L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)+L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2))./(l2.*2.0)).*2.0)./l2-(((Il_2zz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xz.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2yz.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*cos(alpha12).*sin(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-(L2.*sin(alpha12).*sin(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2+(((Il_2zz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xz.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2yz.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2-((((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xy.*((L2.*cos(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*cos(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*(-1.0./2.0))./l2+(Il_2yy.*((L2.*cos(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+(L2.*sin(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2-(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2yz.*((L2.*cos(alpha12).*sin(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-(L2.*sin(alpha12).*sin(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2-((((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xy.*((L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dz+L2.*dtheta2.*cos(theta2)))./l2-(L2.*sin(pi./2.0).*cos(theta2).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)+(Il_2yy.*((L2.*sin(pi./2.0).*cos(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dz+L2.*dtheta2.*cos(theta2)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2-(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)-(Il_2yz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2+(((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*sin(pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xy.*(((dz+L2.*dtheta2.*cos(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((z+L2.*sin(theta2)).*(ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)+L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((ddz+L2.*ddq2.*sin(pi./2.0).*cos(theta2)-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)-(Il_2yz.*(-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+((ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)+L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+((-ddx+L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)+L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)+(Il_2yy.*(((ddz+L2.*ddq2.*sin(pi./2.0).*cos(theta2)-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2-((dz+L2.*dtheta2.*cos(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((z+L2.*sin(theta2)).*(-ddx+L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)+L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2))./(l2.*2.0)).*2.0)./l2+(((Il_2zz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2xz.*((L2.*sin(pi./2.0).*cos(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(z+L2.*sin(theta2)))./l2))./(l2.*2.0)-(Il_2yz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*sin(pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*(-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+((ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)+L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+((-ddx+L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)+L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2-(((Il_2yz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2xy.*((L2.*sin(pi./2.0).*cos(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(z+L2.*sin(theta2)))./l2))./(l2.*2.0)-(Il_2yy.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*sin(pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*(((ddz+L2.*ddq2.*sin(pi./2.0).*cos(theta2)-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2-((dz+L2.*dtheta2.*cos(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((z+L2.*sin(theta2)).*(-ddx+L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)+L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2).*2.0)./l2+(L2.*g.*mL2.*cos(theta2))./2.0+L2.*g.*ml2.*cos(theta2)-(L2.*ml2.*sin(pi./2.0).*cos(theta2).*(ddz./2.0-(L2.*ddq2.*sin(pi./2.0).*cos(theta2))./2.0+(L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2))./2.0))./2.0-L2.*ml2.*sin(pi./2.0).*cos(theta2).*(ddz./4.0-(L2.*ddq2.*sin(pi./2.0).*cos(theta2))./4.0+(L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2))./4.0)+(L2.*ml2.*sin(theta2).*cos(alpha12-pi./2.0).*(ddx./2.0+(L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0))./2.0+(L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0))./2.0))./2.0+L2.*ml2.*sin(theta2).*cos(alpha12-pi./2.0).*(ddx./4.0+(L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0))./4.0+(L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0))./4.0)+(L2.*ml2.*sin(alpha12-pi./2.0).*sin(theta2).*(ddy.*(-1.0./2.0)+(L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2))./2.0+(L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2))./2.0))./2.0+L2.*ml2.*sin(alpha12-pi./2.0).*sin(theta2).*(ddy.*(-1.0./4.0)+(L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2))./4.0+(L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2))./4.0);ddq3.*(IL_3xx.*sin(alpha13-pi./2.0).^2+IL_3yy.*cos(alpha13-pi./2.0).^2+IL_3xy.*sin(alpha13-pi./2.0).*cos(alpha13-pi./2.0).*2.0)+((((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xz.*((L3.*cos(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*cos(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*(-1.0./2.0))./l3+(Il_3yz.*((L3.*cos(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+(L3.*sin(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3-(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3zz.*((L3.*cos(alpha13).*sin(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-(L3.*sin(alpha13).*sin(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3+((((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xz.*((L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dz+L3.*dtheta3.*cos(theta3)))./l3-(L3.*sin(pi./2.0).*cos(theta3).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)+(Il_3yz.*((L3.*sin(pi./2.0).*cos(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dz+L3.*dtheta3.*cos(theta3)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3-(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)-(Il_3zz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3-(((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xz.*(((dz+L3.*dtheta3.*cos(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((z+L3.*sin(theta3)).*(ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)+L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((ddz+L3.*ddq3.*sin(pi./2.0).*cos(theta3)-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)-(Il_3zz.*(-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+((ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)+L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+((-ddx+L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)+L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)+(Il_3yz.*(((ddz+L3.*ddq3.*sin(pi./2.0).*cos(theta3)-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3-((dz+L3.*dtheta3.*cos(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((z+L3.*sin(theta3)).*(-ddx+L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)+L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3))./(l3.*2.0)).*2.0)./l3+(((Il_3xz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xx.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3xy.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*cos(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*cos(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3-(((Il_3xz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xx.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3xy.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dz+L3.*dtheta3.*cos(theta3)))./l3-(L3.*sin(pi./2.0).*cos(theta3).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3-(((Il_3xz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3xx.*((L3.*sin(pi./2.0).*cos(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(z+L3.*sin(theta3)))./l3))./(l3.*2.0)-(Il_3xy.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*sin(pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*(((dz+L3.*dtheta3.*cos(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((z+L3.*sin(theta3)).*(ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)+L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((ddz+L3.*ddq3.*sin(pi./2.0).*cos(theta3)-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3-(((Il_3yz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xy.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3yy.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*cos(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+(L3.*sin(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3-(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*2.0)./l3-(((Il_3yz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xy.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3yy.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*sin(pi./2.0).*cos(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dz+L3.*dtheta3.*cos(theta3)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3-(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*2.0)./l3-((((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3).*((Il_3xx.*((L3.*cos(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*cos(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*(-1.0./2.0))./l3+(Il_3xy.*((L3.*cos(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+(L3.*sin(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3-(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3xz.*((L3.*cos(alpha13).*sin(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-(L3.*sin(alpha13).*sin(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3-((((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3).*((Il_3xx.*((L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dz+L3.*dtheta3.*cos(theta3)))./l3-(L3.*sin(pi./2.0).*cos(theta3).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)+(Il_3xy.*((L3.*sin(pi./2.0).*cos(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dz+L3.*dtheta3.*cos(theta3)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3-(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)-(Il_3xz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3-(((L3.*sin(pi./2.0).*cos(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(z+L3.*sin(theta3)))./l3).*((Il_3xx.*(((dz+L3.*dtheta3.*cos(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((z+L3.*sin(theta3)).*(ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)+L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((ddz+L3.*ddq3.*sin(pi./2.0).*cos(theta3)-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)-(Il_3xz.*(-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+((ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)+L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+((-ddx+L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)+L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)+(Il_3xy.*(((ddz+L3.*ddq3.*sin(pi./2.0).*cos(theta3)-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3-((dz+L3.*dtheta3.*cos(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((z+L3.*sin(theta3)).*(-ddx+L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)+L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3))./(l3.*2.0)).*2.0)./l3-(((Il_3zz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xz.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3yz.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*cos(alpha13).*sin(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-(L3.*sin(alpha13).*sin(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3+(((Il_3zz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xz.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3yz.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3-((((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xy.*((L3.*cos(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*cos(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*(-1.0./2.0))./l3+(Il_3yy.*((L3.*cos(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+(L3.*sin(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3-(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3yz.*((L3.*cos(alpha13).*sin(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-(L3.*sin(alpha13).*sin(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3-((((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xy.*((L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dz+L3.*dtheta3.*cos(theta3)))./l3-(L3.*sin(pi./2.0).*cos(theta3).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)+(Il_3yy.*((L3.*sin(pi./2.0).*cos(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dz+L3.*dtheta3.*cos(theta3)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3-(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)-(Il_3yz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3+(((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*sin(pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xy.*(((dz+L3.*dtheta3.*cos(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((z+L3.*sin(theta3)).*(ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)+L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((ddz+L3.*ddq3.*sin(pi./2.0).*cos(theta3)-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)-(Il_3yz.*(-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+((ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)+L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+((-ddx+L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)+L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)+(Il_3yy.*(((ddz+L3.*ddq3.*sin(pi./2.0).*cos(theta3)-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3-((dz+L3.*dtheta3.*cos(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((z+L3.*sin(theta3)).*(-ddx+L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)+L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3))./(l3.*2.0)).*2.0)./l3+(((Il_3zz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3xz.*((L3.*sin(pi./2.0).*cos(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(z+L3.*sin(theta3)))./l3))./(l3.*2.0)-(Il_3yz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*sin(pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*(-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+((ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)+L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+((-ddx+L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)+L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3-(((Il_3yz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3xy.*((L3.*sin(pi./2.0).*cos(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(z+L3.*sin(theta3)))./l3))./(l3.*2.0)-(Il_3yy.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*sin(pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*(((ddz+L3.*ddq3.*sin(pi./2.0).*cos(theta3)-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3-((dz+L3.*dtheta3.*cos(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((z+L3.*sin(theta3)).*(-ddx+L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)+L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3).*2.0)./l3+(L3.*g.*mL3.*cos(theta3))./2.0+L3.*g.*ml3.*cos(theta3)-(L3.*ml3.*sin(pi./2.0).*cos(theta3).*(ddz./2.0-(L3.*ddq3.*sin(pi./2.0).*cos(theta3))./2.0+(L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3))./2.0))./2.0-L3.*ml3.*sin(pi./2.0).*cos(theta3).*(ddz./4.0-(L3.*ddq3.*sin(pi./2.0).*cos(theta3))./4.0+(L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3))./4.0)+(L3.*ml3.*sin(theta3).*cos(alpha13-pi./2.0).*(ddx./2.0+(L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0))./2.0+(L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0))./2.0))./2.0+L3.*ml3.*sin(theta3).*cos(alpha13-pi./2.0).*(ddx./4.0+(L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0))./4.0+(L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0))./4.0)+(L3.*ml3.*sin(alpha13-pi./2.0).*sin(theta3).*(ddy.*(-1.0./2.0)+(L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3))./2.0+(L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3))./2.0))./2.0+L3.*ml3.*sin(alpha13-pi./2.0).*sin(theta3).*(ddy.*(-1.0./4.0)+(L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3))./4.0+(L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3))./4.0)];
        dGammasdthetaf = @(L1,L2,L3,a1,a2,a3,alpha11,alpha12,alpha13,e1,e2,e3,theta1,theta2,theta3,x,y,z)reshape([L1.*cos(theta1).*(z+L1.*sin(theta1)).*2.0-L1.*cos(alpha11).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)).*2.0-L1.*sin(alpha11).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)).*2.0,0.0,0.0,0.0,L2.*cos(theta2).*(z+L2.*sin(theta2)).*2.0-L2.*cos(alpha12).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)).*2.0-L2.*sin(alpha12).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)).*2.0,0.0,0.0,0.0,L3.*cos(theta3).*(z+L3.*sin(theta3)).*2.0-L3.*cos(alpha13).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)).*2.0-L3.*sin(alpha13).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)).*2.0],[3,3]);
        tauf = @(IL_1xx,IL_1xy,IL_1yy,IL_2xx,IL_2xy,IL_2yy,IL_3xx,IL_3xy,IL_3yy,Il_1xx,Il_1xy,Il_1xz,Il_1yy,Il_1yz,Il_1zz,Il_2xx,Il_2xy,Il_2xz,Il_2yy,Il_2yz,Il_2zz,Il_3xx,Il_3xy,Il_3xz,Il_3yy,Il_3yz,Il_3zz,L1,L2,L3,a1,a2,a3,alpha11,alpha12,alpha13,ddq1,ddq2,ddq3,ddx,ddy,ddz,dtheta1,dtheta2,dtheta3,dx,dy,dz,e1,e2,e3,g,l1,l2,l3,lambda1,lambda2,lambda3,mL1,mL2,mL3,mlA,ml2,ml3,pi,theta1,theta2,theta3,x,y,z)[lambda1.*(L1.*cos(theta1).*(z+L1.*sin(theta1)).*-2.0+L1.*cos(alpha11).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)).*2.0+L1.*sin(alpha11).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)).*2.0)+ddq1.*(IL_1xx.*sin(alpha11-pi./2.0).^2+IL_1yy.*cos(alpha11-pi./2.0).^2+IL_1xy.*sin(alpha11-pi./2.0).*cos(alpha11-pi./2.0).*2.0)+((((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xz.*((L1.*cos(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*cos(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*(-1.0./2.0))./l1+(Il_1yz.*((L1.*cos(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+(L1.*sin(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1-(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1zz.*((L1.*cos(alpha11).*sin(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-(L1.*sin(alpha11).*sin(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1+((((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xz.*((L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dz+L1.*dtheta1.*cos(theta1)))./l1-(L1.*sin(pi./2.0).*cos(theta1).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)+(Il_1yz.*((L1.*sin(pi./2.0).*cos(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dz+L1.*dtheta1.*cos(theta1)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1-(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)-(Il_1zz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1-(((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xz.*(((dz+L1.*dtheta1.*cos(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((z+L1.*sin(theta1)).*(ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)+L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((ddz+L1.*ddq1.*sin(pi./2.0).*cos(theta1)-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)-(Il_1zz.*(-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+((ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)+L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+((-ddx+L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)+L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)+(Il_1yz.*(((ddz+L1.*ddq1.*sin(pi./2.0).*cos(theta1)-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1-((dz+L1.*dtheta1.*cos(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((z+L1.*sin(theta1)).*(-ddx+L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)+L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1))./(l1.*2.0)).*2.0)./l1+(((Il_1xz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xx.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1xy.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*cos(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*cos(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1-(((Il_1xz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xx.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1xy.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dz+L1.*dtheta1.*cos(theta1)))./l1-(L1.*sin(pi./2.0).*cos(theta1).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1-(((Il_1xz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1xx.*((L1.*sin(pi./2.0).*cos(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(z+L1.*sin(theta1)))./l1))./(l1.*2.0)-(Il_1xy.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*sin(pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*(((dz+L1.*dtheta1.*cos(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((z+L1.*sin(theta1)).*(ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)+L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((ddz+L1.*ddq1.*sin(pi./2.0).*cos(theta1)-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1-(((Il_1yz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xy.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1yy.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*cos(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+(L1.*sin(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1-(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*2.0)./l1-(((Il_1yz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xy.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1yy.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*sin(pi./2.0).*cos(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dz+L1.*dtheta1.*cos(theta1)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1-(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*2.0)./l1-((((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1).*((Il_1xx.*((L1.*cos(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*cos(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*(-1.0./2.0))./l1+(Il_1xy.*((L1.*cos(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+(L1.*sin(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1-(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1xz.*((L1.*cos(alpha11).*sin(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-(L1.*sin(alpha11).*sin(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1-((((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1).*((Il_1xx.*((L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dz+L1.*dtheta1.*cos(theta1)))./l1-(L1.*sin(pi./2.0).*cos(theta1).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)+(Il_1xy.*((L1.*sin(pi./2.0).*cos(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dz+L1.*dtheta1.*cos(theta1)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1-(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)-(Il_1xz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1-(((L1.*sin(pi./2.0).*cos(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(z+L1.*sin(theta1)))./l1).*((Il_1xx.*(((dz+L1.*dtheta1.*cos(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((z+L1.*sin(theta1)).*(ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)+L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((ddz+L1.*ddq1.*sin(pi./2.0).*cos(theta1)-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)-(Il_1xz.*(-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+((ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)+L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+((-ddx+L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)+L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)+(Il_1xy.*(((ddz+L1.*ddq1.*sin(pi./2.0).*cos(theta1)-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1-((dz+L1.*dtheta1.*cos(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((z+L1.*sin(theta1)).*(-ddx+L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)+L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1))./(l1.*2.0)).*2.0)./l1-(((Il_1zz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xz.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1yz.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*cos(alpha11).*sin(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-(L1.*sin(alpha11).*sin(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1+(((Il_1zz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xz.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1yz.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1-((((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xy.*((L1.*cos(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*cos(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*(-1.0./2.0))./l1+(Il_1yy.*((L1.*cos(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+(L1.*sin(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1-(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1yz.*((L1.*cos(alpha11).*sin(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-(L1.*sin(alpha11).*sin(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1-((((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xy.*((L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dz+L1.*dtheta1.*cos(theta1)))./l1-(L1.*sin(pi./2.0).*cos(theta1).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)+(Il_1yy.*((L1.*sin(pi./2.0).*cos(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dz+L1.*dtheta1.*cos(theta1)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1-(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)-(Il_1yz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1+(((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*sin(pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xy.*(((dz+L1.*dtheta1.*cos(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((z+L1.*sin(theta1)).*(ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)+L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((ddz+L1.*ddq1.*sin(pi./2.0).*cos(theta1)-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)-(Il_1yz.*(-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+((ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)+L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+((-ddx+L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)+L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)+(Il_1yy.*(((ddz+L1.*ddq1.*sin(pi./2.0).*cos(theta1)-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1-((dz+L1.*dtheta1.*cos(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((z+L1.*sin(theta1)).*(-ddx+L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)+L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1))./(l1.*2.0)).*2.0)./l1+(((Il_1zz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1xz.*((L1.*sin(pi./2.0).*cos(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(z+L1.*sin(theta1)))./l1))./(l1.*2.0)-(Il_1yz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*sin(pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*(-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+((ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)+L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+((-ddx+L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)+L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1-(((Il_1yz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1xy.*((L1.*sin(pi./2.0).*cos(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(z+L1.*sin(theta1)))./l1))./(l1.*2.0)-(Il_1yy.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*sin(pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*(((ddz+L1.*ddq1.*sin(pi./2.0).*cos(theta1)-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1-((dz+L1.*dtheta1.*cos(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((z+L1.*sin(theta1)).*(-ddx+L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)+L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1).*2.0)./l1+(L1.*g.*mL1.*cos(theta1))./2.0+L1.*g.*mlA.*cos(theta1)-(L1.*mlA.*sin(pi./2.0).*cos(theta1).*(ddz./2.0-(L1.*ddq1.*sin(pi./2.0).*cos(theta1))./2.0+(L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1))./2.0))./2.0-L1.*mlA.*sin(pi./2.0).*cos(theta1).*(ddz./4.0-(L1.*ddq1.*sin(pi./2.0).*cos(theta1))./4.0+(L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1))./4.0)+(L1.*mlA.*sin(theta1).*cos(alpha11-pi./2.0).*(ddx./2.0+(L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0))./2.0+(L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0))./2.0))./2.0+L1.*mlA.*sin(theta1).*cos(alpha11-pi./2.0).*(ddx./4.0+(L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0))./4.0+(L1.*ddq1.*sin(theta1).*cos(alpha11-pi./2.0))./4.0)+(L1.*mlA.*sin(alpha11-pi./2.0).*sin(theta1).*(ddy.*(-1.0./2.0)+(L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1))./2.0+(L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1))./2.0))./2.0+L1.*mlA.*sin(alpha11-pi./2.0).*sin(theta1).*(ddy.*(-1.0./4.0)+(L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1))./4.0+(L1.*ddq1.*sin(alpha11-pi./2.0).*sin(theta1))./4.0);lambda2.*(L2.*cos(theta2).*(z+L2.*sin(theta2)).*-2.0+L2.*cos(alpha12).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)).*2.0+L2.*sin(alpha12).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)).*2.0)+ddq2.*(IL_2xx.*sin(alpha12-pi./2.0).^2+IL_2yy.*cos(alpha12-pi./2.0).^2+IL_2xy.*sin(alpha12-pi./2.0).*cos(alpha12-pi./2.0).*2.0)+((((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xz.*((L2.*cos(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*cos(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*(-1.0./2.0))./l2+(Il_2yz.*((L2.*cos(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+(L2.*sin(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2-(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2zz.*((L2.*cos(alpha12).*sin(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-(L2.*sin(alpha12).*sin(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2+((((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xz.*((L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dz+L2.*dtheta2.*cos(theta2)))./l2-(L2.*sin(pi./2.0).*cos(theta2).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)+(Il_2yz.*((L2.*sin(pi./2.0).*cos(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dz+L2.*dtheta2.*cos(theta2)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2-(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)-(Il_2zz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2-(((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xz.*(((dz+L2.*dtheta2.*cos(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((z+L2.*sin(theta2)).*(ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)+L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((ddz+L2.*ddq2.*sin(pi./2.0).*cos(theta2)-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)-(Il_2zz.*(-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+((ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)+L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+((-ddx+L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)+L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)+(Il_2yz.*(((ddz+L2.*ddq2.*sin(pi./2.0).*cos(theta2)-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2-((dz+L2.*dtheta2.*cos(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((z+L2.*sin(theta2)).*(-ddx+L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)+L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2))./(l2.*2.0)).*2.0)./l2+(((Il_2xz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xx.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2xy.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*cos(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*cos(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2-(((Il_2xz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xx.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2xy.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dz+L2.*dtheta2.*cos(theta2)))./l2-(L2.*sin(pi./2.0).*cos(theta2).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2-(((Il_2xz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2xx.*((L2.*sin(pi./2.0).*cos(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(z+L2.*sin(theta2)))./l2))./(l2.*2.0)-(Il_2xy.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*sin(pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*(((dz+L2.*dtheta2.*cos(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((z+L2.*sin(theta2)).*(ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)+L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((ddz+L2.*ddq2.*sin(pi./2.0).*cos(theta2)-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2-(((Il_2yz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xy.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2yy.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*cos(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+(L2.*sin(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2-(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*2.0)./l2-(((Il_2yz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xy.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2yy.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*sin(pi./2.0).*cos(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dz+L2.*dtheta2.*cos(theta2)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2-(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*2.0)./l2-((((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2).*((Il_2xx.*((L2.*cos(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*cos(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*(-1.0./2.0))./l2+(Il_2xy.*((L2.*cos(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+(L2.*sin(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2-(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2xz.*((L2.*cos(alpha12).*sin(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-(L2.*sin(alpha12).*sin(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2-((((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2).*((Il_2xx.*((L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dz+L2.*dtheta2.*cos(theta2)))./l2-(L2.*sin(pi./2.0).*cos(theta2).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)+(Il_2xy.*((L2.*sin(pi./2.0).*cos(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dz+L2.*dtheta2.*cos(theta2)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2-(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)-(Il_2xz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2-(((L2.*sin(pi./2.0).*cos(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(z+L2.*sin(theta2)))./l2).*((Il_2xx.*(((dz+L2.*dtheta2.*cos(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((z+L2.*sin(theta2)).*(ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)+L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((ddz+L2.*ddq2.*sin(pi./2.0).*cos(theta2)-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)-(Il_2xz.*(-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+((ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)+L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+((-ddx+L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)+L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)+(Il_2xy.*(((ddz+L2.*ddq2.*sin(pi./2.0).*cos(theta2)-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2-((dz+L2.*dtheta2.*cos(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((z+L2.*sin(theta2)).*(-ddx+L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)+L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2))./(l2.*2.0)).*2.0)./l2-(((Il_2zz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xz.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2yz.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*cos(alpha12).*sin(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-(L2.*sin(alpha12).*sin(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2+(((Il_2zz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xz.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2yz.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2-((((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xy.*((L2.*cos(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*cos(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*(-1.0./2.0))./l2+(Il_2yy.*((L2.*cos(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+(L2.*sin(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2-(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2yz.*((L2.*cos(alpha12).*sin(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-(L2.*sin(alpha12).*sin(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2-((((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xy.*((L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dz+L2.*dtheta2.*cos(theta2)))./l2-(L2.*sin(pi./2.0).*cos(theta2).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)+(Il_2yy.*((L2.*sin(pi./2.0).*cos(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dz+L2.*dtheta2.*cos(theta2)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2-(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)-(Il_2yz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2+(((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*sin(pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xy.*(((dz+L2.*dtheta2.*cos(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((z+L2.*sin(theta2)).*(ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)+L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((ddz+L2.*ddq2.*sin(pi./2.0).*cos(theta2)-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)-(Il_2yz.*(-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+((ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)+L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+((-ddx+L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)+L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)+(Il_2yy.*(((ddz+L2.*ddq2.*sin(pi./2.0).*cos(theta2)-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2-((dz+L2.*dtheta2.*cos(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((z+L2.*sin(theta2)).*(-ddx+L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)+L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2))./(l2.*2.0)).*2.0)./l2+(((Il_2zz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2xz.*((L2.*sin(pi./2.0).*cos(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(z+L2.*sin(theta2)))./l2))./(l2.*2.0)-(Il_2yz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*sin(pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*(-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+((ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)+L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+((-ddx+L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)+L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2-(((Il_2yz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2xy.*((L2.*sin(pi./2.0).*cos(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(z+L2.*sin(theta2)))./l2))./(l2.*2.0)-(Il_2yy.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*sin(pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*(((ddz+L2.*ddq2.*sin(pi./2.0).*cos(theta2)-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2-((dz+L2.*dtheta2.*cos(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((z+L2.*sin(theta2)).*(-ddx+L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)+L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2).*2.0)./l2+(L2.*g.*mL2.*cos(theta2))./2.0+L2.*g.*ml2.*cos(theta2)-(L2.*ml2.*sin(pi./2.0).*cos(theta2).*(ddz./2.0-(L2.*ddq2.*sin(pi./2.0).*cos(theta2))./2.0+(L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2))./2.0))./2.0-L2.*ml2.*sin(pi./2.0).*cos(theta2).*(ddz./4.0-(L2.*ddq2.*sin(pi./2.0).*cos(theta2))./4.0+(L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2))./4.0)+(L2.*ml2.*sin(theta2).*cos(alpha12-pi./2.0).*(ddx./2.0+(L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0))./2.0+(L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0))./2.0))./2.0+L2.*ml2.*sin(theta2).*cos(alpha12-pi./2.0).*(ddx./4.0+(L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0))./4.0+(L2.*ddq2.*sin(theta2).*cos(alpha12-pi./2.0))./4.0)+(L2.*ml2.*sin(alpha12-pi./2.0).*sin(theta2).*(ddy.*(-1.0./2.0)+(L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2))./2.0+(L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2))./2.0))./2.0+L2.*ml2.*sin(alpha12-pi./2.0).*sin(theta2).*(ddy.*(-1.0./4.0)+(L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2))./4.0+(L2.*ddq2.*sin(alpha12-pi./2.0).*sin(theta2))./4.0);lambda3.*(L3.*cos(theta3).*(z+L3.*sin(theta3)).*-2.0+L3.*cos(alpha13).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)).*2.0+L3.*sin(alpha13).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)).*2.0)+ddq3.*(IL_3xx.*sin(alpha13-pi./2.0).^2+IL_3yy.*cos(alpha13-pi./2.0).^2+IL_3xy.*sin(alpha13-pi./2.0).*cos(alpha13-pi./2.0).*2.0)+((((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xz.*((L3.*cos(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*cos(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*(-1.0./2.0))./l3+(Il_3yz.*((L3.*cos(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+(L3.*sin(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3-(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3zz.*((L3.*cos(alpha13).*sin(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-(L3.*sin(alpha13).*sin(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3+((((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xz.*((L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dz+L3.*dtheta3.*cos(theta3)))./l3-(L3.*sin(pi./2.0).*cos(theta3).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)+(Il_3yz.*((L3.*sin(pi./2.0).*cos(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dz+L3.*dtheta3.*cos(theta3)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3-(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)-(Il_3zz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3-(((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xz.*(((dz+L3.*dtheta3.*cos(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((z+L3.*sin(theta3)).*(ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)+L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((ddz+L3.*ddq3.*sin(pi./2.0).*cos(theta3)-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)-(Il_3zz.*(-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+((ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)+L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+((-ddx+L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)+L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)+(Il_3yz.*(((ddz+L3.*ddq3.*sin(pi./2.0).*cos(theta3)-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3-((dz+L3.*dtheta3.*cos(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((z+L3.*sin(theta3)).*(-ddx+L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)+L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3))./(l3.*2.0)).*2.0)./l3+(((Il_3xz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xx.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3xy.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*cos(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*cos(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3-(((Il_3xz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xx.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3xy.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dz+L3.*dtheta3.*cos(theta3)))./l3-(L3.*sin(pi./2.0).*cos(theta3).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3-(((Il_3xz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3xx.*((L3.*sin(pi./2.0).*cos(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(z+L3.*sin(theta3)))./l3))./(l3.*2.0)-(Il_3xy.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*sin(pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*(((dz+L3.*dtheta3.*cos(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((z+L3.*sin(theta3)).*(ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)+L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((ddz+L3.*ddq3.*sin(pi./2.0).*cos(theta3)-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3-(((Il_3yz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xy.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3yy.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*cos(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+(L3.*sin(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3-(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*2.0)./l3-(((Il_3yz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xy.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3yy.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*sin(pi./2.0).*cos(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dz+L3.*dtheta3.*cos(theta3)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3-(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*2.0)./l3-((((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3).*((Il_3xx.*((L3.*cos(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*cos(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*(-1.0./2.0))./l3+(Il_3xy.*((L3.*cos(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+(L3.*sin(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3-(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3xz.*((L3.*cos(alpha13).*sin(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-(L3.*sin(alpha13).*sin(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3-((((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3).*((Il_3xx.*((L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dz+L3.*dtheta3.*cos(theta3)))./l3-(L3.*sin(pi./2.0).*cos(theta3).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)+(Il_3xy.*((L3.*sin(pi./2.0).*cos(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dz+L3.*dtheta3.*cos(theta3)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3-(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)-(Il_3xz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3-(((L3.*sin(pi./2.0).*cos(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(z+L3.*sin(theta3)))./l3).*((Il_3xx.*(((dz+L3.*dtheta3.*cos(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((z+L3.*sin(theta3)).*(ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)+L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((ddz+L3.*ddq3.*sin(pi./2.0).*cos(theta3)-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)-(Il_3xz.*(-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+((ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)+L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+((-ddx+L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)+L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)+(Il_3xy.*(((ddz+L3.*ddq3.*sin(pi./2.0).*cos(theta3)-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3-((dz+L3.*dtheta3.*cos(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((z+L3.*sin(theta3)).*(-ddx+L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)+L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3))./(l3.*2.0)).*2.0)./l3-(((Il_3zz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xz.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3yz.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*cos(alpha13).*sin(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-(L3.*sin(alpha13).*sin(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3+(((Il_3zz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xz.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3yz.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3-((((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xy.*((L3.*cos(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*cos(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*(-1.0./2.0))./l3+(Il_3yy.*((L3.*cos(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+(L3.*sin(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3-(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3yz.*((L3.*cos(alpha13).*sin(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-(L3.*sin(alpha13).*sin(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3-((((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xy.*((L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dz+L3.*dtheta3.*cos(theta3)))./l3-(L3.*sin(pi./2.0).*cos(theta3).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)+(Il_3yy.*((L3.*sin(pi./2.0).*cos(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dz+L3.*dtheta3.*cos(theta3)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3-(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)-(Il_3yz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3+(((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*sin(pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xy.*(((dz+L3.*dtheta3.*cos(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((z+L3.*sin(theta3)).*(ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)+L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((ddz+L3.*ddq3.*sin(pi./2.0).*cos(theta3)-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)-(Il_3yz.*(-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+((ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)+L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+((-ddx+L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)+L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)+(Il_3yy.*(((ddz+L3.*ddq3.*sin(pi./2.0).*cos(theta3)-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3-((dz+L3.*dtheta3.*cos(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((z+L3.*sin(theta3)).*(-ddx+L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)+L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3))./(l3.*2.0)).*2.0)./l3+(((Il_3zz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3xz.*((L3.*sin(pi./2.0).*cos(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(z+L3.*sin(theta3)))./l3))./(l3.*2.0)-(Il_3yz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*sin(pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*(-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+((ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)+L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+((-ddx+L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)+L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3-(((Il_3yz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3xy.*((L3.*sin(pi./2.0).*cos(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(z+L3.*sin(theta3)))./l3))./(l3.*2.0)-(Il_3yy.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*sin(pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*(((ddz+L3.*ddq3.*sin(pi./2.0).*cos(theta3)-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3-((dz+L3.*dtheta3.*cos(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((z+L3.*sin(theta3)).*(-ddx+L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)+L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3).*2.0)./l3+(L3.*g.*mL3.*cos(theta3))./2.0+L3.*g.*ml3.*cos(theta3)-(L3.*ml3.*sin(pi./2.0).*cos(theta3).*(ddz./2.0-(L3.*ddq3.*sin(pi./2.0).*cos(theta3))./2.0+(L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3))./2.0))./2.0-L3.*ml3.*sin(pi./2.0).*cos(theta3).*(ddz./4.0-(L3.*ddq3.*sin(pi./2.0).*cos(theta3))./4.0+(L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3))./4.0)+(L3.*ml3.*sin(theta3).*cos(alpha13-pi./2.0).*(ddx./2.0+(L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0))./2.0+(L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0))./2.0))./2.0+L3.*ml3.*sin(theta3).*cos(alpha13-pi./2.0).*(ddx./4.0+(L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0))./4.0+(L3.*ddq3.*sin(theta3).*cos(alpha13-pi./2.0))./4.0)+(L3.*ml3.*sin(alpha13-pi./2.0).*sin(theta3).*(ddy.*(-1.0./2.0)+(L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3))./2.0+(L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3))./2.0))./2.0+L3.*ml3.*sin(alpha13-pi./2.0).*sin(theta3).*(ddy.*(-1.0./4.0)+(L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3))./4.0+(L3.*ddq3.*sin(alpha13-pi./2.0).*sin(theta3))./4.0)];
        Mf = @(IL_1xx,IL_1xy,IL_1yy,IL_2xx,IL_2xy,IL_2yy,IL_3xx,IL_3xy,IL_3yy,Il_1xx,Il_1xy,Il_1xz,Il_1yy,Il_1yz,Il_1zz,Il_2xx,Il_2xy,Il_2xz,Il_2yy,Il_2yz,Il_2zz,Il_3xx,Il_3xy,Il_3xz,Il_3yy,Il_3yz,Il_3zz,L1,L2,L3,a1,a2,a3,alpha11,alpha12,alpha13,e1,e2,e3,l1,l2,l3,mlA,ml2,ml3,pi,theta1,theta2,theta3,x,y,z)reshape([IL_1xx.*sin(alpha11-pi./2.0).^2+IL_1yy.*cos(alpha11-pi./2.0).^2+(((L1.*sin(pi./2.0).*cos(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(z+L1.*sin(theta1)))./l1).*((Il_1xz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1xx.*((L1.*sin(pi./2.0).*cos(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(z+L1.*sin(theta1)))./l1))./(l1.*2.0)-(Il_1xy.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*sin(pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*4.0)./l1-(((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*sin(pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1yz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1xy.*((L1.*sin(pi./2.0).*cos(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(z+L1.*sin(theta1)))./l1))./(l1.*2.0)-(Il_1yy.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*sin(pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*4.0)./l1+IL_1xy.*sin(alpha11-pi./2.0).*cos(alpha11-pi./2.0).*2.0+(((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1zz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1xz.*((L1.*sin(pi./2.0).*cos(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(z+L1.*sin(theta1)))./l1))./(l1.*2.0)-(Il_1yz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*sin(pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*4.0)./l1+(L1.^2.*mlA.*sin(theta1).^2.*cos(alpha11-pi./2.0).^2)./2.0+(L1.^2.*mlA.*sin(alpha11-pi./2.0).^2.*sin(theta1).^2)./2.0+(L1.^2.*mlA.*sin(pi./2.0).^2.*cos(theta1).^2)./2.0,0.0,0.0,0.0,IL_2xx.*sin(alpha12-pi./2.0).^2+IL_2yy.*cos(alpha12-pi./2.0).^2+(((L2.*sin(pi./2.0).*cos(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(z+L2.*sin(theta2)))./l2).*((Il_2xz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2xx.*((L2.*sin(pi./2.0).*cos(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(z+L2.*sin(theta2)))./l2))./(l2.*2.0)-(Il_2xy.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*sin(pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*4.0)./l2-(((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*sin(pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2yz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2xy.*((L2.*sin(pi./2.0).*cos(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(z+L2.*sin(theta2)))./l2))./(l2.*2.0)-(Il_2yy.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*sin(pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*4.0)./l2+IL_2xy.*sin(alpha12-pi./2.0).*cos(alpha12-pi./2.0).*2.0+(((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2zz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2xz.*((L2.*sin(pi./2.0).*cos(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(z+L2.*sin(theta2)))./l2))./(l2.*2.0)-(Il_2yz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*sin(pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*4.0)./l2+(L2.^2.*ml2.*sin(theta2).^2.*cos(alpha12-pi./2.0).^2)./2.0+(L2.^2.*ml2.*sin(alpha12-pi./2.0).^2.*sin(theta2).^2)./2.0+(L2.^2.*ml2.*sin(pi./2.0).^2.*cos(theta2).^2)./2.0,0.0,0.0,0.0,IL_3xx.*sin(alpha13-pi./2.0).^2+IL_3yy.*cos(alpha13-pi./2.0).^2+(((L3.*sin(pi./2.0).*cos(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(z+L3.*sin(theta3)))./l3).*((Il_3xz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3xx.*((L3.*sin(pi./2.0).*cos(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(z+L3.*sin(theta3)))./l3))./(l3.*2.0)-(Il_3xy.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*sin(pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*4.0)./l3-(((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*sin(pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3yz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3xy.*((L3.*sin(pi./2.0).*cos(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(z+L3.*sin(theta3)))./l3))./(l3.*2.0)-(Il_3yy.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*sin(pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*4.0)./l3+IL_3xy.*sin(alpha13-pi./2.0).*cos(alpha13-pi./2.0).*2.0+(((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3zz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3xz.*((L3.*sin(pi./2.0).*cos(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(z+L3.*sin(theta3)))./l3))./(l3.*2.0)-(Il_3yz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*sin(pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*4.0)./l3+(L3.^2.*ml3.*sin(theta3).^2.*cos(alpha13-pi./2.0).^2)./2.0+(L3.^2.*ml3.*sin(alpha13-pi./2.0).^2.*sin(theta3).^2)./2.0+(L3.^2.*ml3.*sin(pi./2.0).^2.*cos(theta3).^2)./2.0],[3,3]);
        ff = @(Il_1xx,Il_1xy,Il_1xz,Il_1yy,Il_1yz,Il_1zz,Il_2xx,Il_2xy,Il_2xz,Il_2yy,Il_2yz,Il_2zz,Il_3xx,Il_3xy,Il_3xz,Il_3yy,Il_3yz,Il_3zz,L1,L2,L3,a1,a2,a3,alpha11,alpha12,alpha13,ddx,ddy,ddz,dtheta1,dtheta2,dtheta3,dx,dy,dz,e1,e2,e3,g,l1,l2,l3,lambda1,lambda2,lambda3,mL1,mL2,mL3,mlA,ml2,ml3,pi,theta1,theta2,theta3,x,y,z)[-lambda1.*(L1.*cos(theta1).*(z+L1.*sin(theta1)).*-2.0+L1.*cos(alpha11).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)).*2.0+L1.*sin(alpha11).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)).*2.0)-((((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xz.*((L1.*cos(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*cos(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*(-1.0./2.0))./l1+(Il_1yz.*((L1.*cos(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+(L1.*sin(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1-(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1zz.*((L1.*cos(alpha11).*sin(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-(L1.*sin(alpha11).*sin(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1-((((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xz.*((L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dz+L1.*dtheta1.*cos(theta1)))./l1-(L1.*sin(pi./2.0).*cos(theta1).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)+(Il_1yz.*((L1.*sin(pi./2.0).*cos(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dz+L1.*dtheta1.*cos(theta1)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1-(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)-(Il_1zz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1+(((Il_1xz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1xx.*((L1.*sin(pi./2.0).*cos(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(z+L1.*sin(theta1)))./l1))./(l1.*2.0)-(Il_1xy.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*sin(pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*(((dz+L1.*dtheta1.*cos(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)).*(z+L1.*sin(theta1)))./l1-((ddz-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1-(((Il_1xz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xx.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1xy.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*cos(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*cos(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1+(((Il_1xz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xx.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1xy.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dz+L1.*dtheta1.*cos(theta1)))./l1-(L1.*sin(pi./2.0).*cos(theta1).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1+(((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xz.*(((dz+L1.*dtheta1.*cos(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)).*(z+L1.*sin(theta1)))./l1-((ddz-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)-(Il_1yz.*(((dz+L1.*dtheta1.*cos(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((z+L1.*sin(theta1)).*(ddx-L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)))./l1-((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1-((ddz-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1zz.*(((ddx-L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*2.0)./l1-(((Il_1yz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1xy.*((L1.*sin(pi./2.0).*cos(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(z+L1.*sin(theta1)))./l1))./(l1.*2.0)-(Il_1yy.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*sin(pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*(((dz+L1.*dtheta1.*cos(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((z+L1.*sin(theta1)).*(ddx-L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)))./l1-((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1-((ddz-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*2.0)./l1+(((Il_1zz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1xz.*((L1.*sin(pi./2.0).*cos(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(z+L1.*sin(theta1)))./l1))./(l1.*2.0)-(Il_1yz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*sin(pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*(((ddx-L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*2.0)./l1+(((Il_1yz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xy.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1yy.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*cos(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+(L1.*sin(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1-(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*2.0)./l1+(((Il_1yz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xy.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1yy.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*sin(pi./2.0).*cos(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dz+L1.*dtheta1.*cos(theta1)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1-(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*2.0)./l1+((((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1).*((Il_1xx.*((L1.*cos(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*cos(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*(-1.0./2.0))./l1+(Il_1xy.*((L1.*cos(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+(L1.*sin(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1-(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1xz.*((L1.*cos(alpha11).*sin(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-(L1.*sin(alpha11).*sin(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1+((((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1).*((Il_1xx.*((L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dz+L1.*dtheta1.*cos(theta1)))./l1-(L1.*sin(pi./2.0).*cos(theta1).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)+(Il_1xy.*((L1.*sin(pi./2.0).*cos(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dz+L1.*dtheta1.*cos(theta1)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1-(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)-(Il_1xz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1+(((Il_1zz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xz.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1yz.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*cos(alpha11).*sin(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-(L1.*sin(alpha11).*sin(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1-(((Il_1zz.*(((dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*(-1.0./2.0))./l1+(Il_1xz.*(((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-((z+L1.*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1))./(l1.*2.0)+(Il_1yz.*(((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*2.0)./l1+((((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xy.*((L1.*cos(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*cos(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1).*(-1.0./2.0))./l1+(Il_1yy.*((L1.*cos(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+(L1.*sin(alpha11).*sin(theta1).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1-(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1yz.*((L1.*cos(alpha11).*sin(theta1).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-(L1.*sin(alpha11).*sin(theta1).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1+((((z+L1.*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xy.*((L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dz+L1.*dtheta1.*cos(theta1)))./l1-(L1.*sin(pi./2.0).*cos(theta1).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(z+L1.*sin(theta1)))./l1+(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)+(Il_1yy.*((L1.*sin(pi./2.0).*cos(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dz+L1.*dtheta1.*cos(theta1)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1-(L1.*dtheta1.*sin(pi./2.0).*sin(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)-(Il_1yz.*((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)))./l1+(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)))./l1+(L1.*dtheta1.*sin(alpha11-pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1+(L1.*dtheta1.*cos(theta1).*cos(alpha11-pi./2.0).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)).*2.0)./l1+(((L1.*sin(pi./2.0).*cos(theta1).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1-(L1.*sin(alpha11-pi./2.0).*sin(theta1).*(z+L1.*sin(theta1)))./l1).*((Il_1xx.*(((dz+L1.*dtheta1.*cos(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)).*(z+L1.*sin(theta1)))./l1-((ddz-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)-(Il_1xy.*(((dz+L1.*dtheta1.*cos(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((z+L1.*sin(theta1)).*(ddx-L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)))./l1-((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1-((ddz-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1xz.*(((ddx-L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*2.0)./l1-(((L1.*sin(theta1).*cos(alpha11-pi./2.0).*(z+L1.*sin(theta1)))./l1+(L1.*sin(pi./2.0).*cos(theta1).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1).*((Il_1xy.*(((dz+L1.*dtheta1.*cos(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1+((ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)).*(z+L1.*sin(theta1)))./l1-((ddz-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1))./(l1.*2.0)-(Il_1yy.*(((dz+L1.*dtheta1.*cos(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1+((z+L1.*sin(theta1)).*(ddx-L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)))./l1-((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dz+L1.*dtheta1.*sin(pi./2.0).*cos(theta1)))./l1-((ddz-L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)+(Il_1yz.*(((ddx-L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0)).*(y+a1.*cos(alpha11)-e1.*cos(alpha11)+L1.*cos(alpha11).*cos(theta1)))./l1+((dy-L1.*dtheta1.*cos(alpha11).*sin(theta1)).*(dx-L1.*dtheta1.*sin(theta1).*cos(alpha11-pi./2.0)))./l1-((dx-L1.*dtheta1.*sin(alpha11).*sin(theta1)).*(dy+L1.*dtheta1.*sin(alpha11-pi./2.0).*sin(theta1)))./l1-((ddy+L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1)).*(x+a1.*sin(alpha11)-e1.*sin(alpha11)+L1.*cos(theta1).*sin(alpha11)))./l1))./(l1.*2.0)).*2.0)./l1-(L1.*g.*mL1.*cos(theta1))./2.0-L1.*g.*mlA.*cos(theta1)-(L1.*mlA.*sin(theta1).*cos(alpha11-pi./2.0).*(ddx./2.0+(L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0))./2.0))./2.0-L1.*mlA.*sin(theta1).*cos(alpha11-pi./2.0).*(ddx./4.0+(L1.*dtheta1.^2.*cos(theta1).*cos(alpha11-pi./2.0))./4.0)+(L1.*mlA.*sin(alpha11-pi./2.0).*sin(theta1).*(ddy./2.0-(L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1))./2.0))./2.0+L1.*mlA.*sin(alpha11-pi./2.0).*sin(theta1).*(ddy./4.0-(L1.*dtheta1.^2.*sin(alpha11-pi./2.0).*cos(theta1))./4.0)+(L1.*mlA.*sin(pi./2.0).*cos(theta1).*(ddz./2.0+(L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1))./2.0))./2.0+L1.*mlA.*sin(pi./2.0).*cos(theta1).*(ddz./4.0+(L1.*dtheta1.^2.*sin(pi./2.0).*sin(theta1))./4.0);-lambda2.*(L2.*cos(theta2).*(z+L2.*sin(theta2)).*-2.0+L2.*cos(alpha12).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)).*2.0+L2.*sin(alpha12).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)).*2.0)-((((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xz.*((L2.*cos(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*cos(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*(-1.0./2.0))./l2+(Il_2yz.*((L2.*cos(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+(L2.*sin(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2-(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2zz.*((L2.*cos(alpha12).*sin(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-(L2.*sin(alpha12).*sin(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2-((((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xz.*((L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dz+L2.*dtheta2.*cos(theta2)))./l2-(L2.*sin(pi./2.0).*cos(theta2).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)+(Il_2yz.*((L2.*sin(pi./2.0).*cos(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dz+L2.*dtheta2.*cos(theta2)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2-(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)-(Il_2zz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2+(((Il_2xz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2xx.*((L2.*sin(pi./2.0).*cos(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(z+L2.*sin(theta2)))./l2))./(l2.*2.0)-(Il_2xy.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*sin(pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*(((dz+L2.*dtheta2.*cos(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)).*(z+L2.*sin(theta2)))./l2-((ddz-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2-(((Il_2xz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xx.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2xy.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*cos(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*cos(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2+(((Il_2xz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xx.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2xy.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dz+L2.*dtheta2.*cos(theta2)))./l2-(L2.*sin(pi./2.0).*cos(theta2).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2+(((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xz.*(((dz+L2.*dtheta2.*cos(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)).*(z+L2.*sin(theta2)))./l2-((ddz-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)-(Il_2yz.*(((dz+L2.*dtheta2.*cos(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((z+L2.*sin(theta2)).*(ddx-L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)))./l2-((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2-((ddz-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2zz.*(((ddx-L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*2.0)./l2-(((Il_2yz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2xy.*((L2.*sin(pi./2.0).*cos(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(z+L2.*sin(theta2)))./l2))./(l2.*2.0)-(Il_2yy.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*sin(pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*(((dz+L2.*dtheta2.*cos(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((z+L2.*sin(theta2)).*(ddx-L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)))./l2-((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2-((ddz-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*2.0)./l2+(((Il_2zz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2xz.*((L2.*sin(pi./2.0).*cos(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(z+L2.*sin(theta2)))./l2))./(l2.*2.0)-(Il_2yz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*sin(pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*(((ddx-L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*2.0)./l2+(((Il_2yz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xy.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2yy.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*cos(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+(L2.*sin(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2-(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*2.0)./l2+(((Il_2yz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xy.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2yy.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*sin(pi./2.0).*cos(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dz+L2.*dtheta2.*cos(theta2)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2-(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*2.0)./l2+((((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2).*((Il_2xx.*((L2.*cos(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*cos(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*(-1.0./2.0))./l2+(Il_2xy.*((L2.*cos(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+(L2.*sin(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2-(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2xz.*((L2.*cos(alpha12).*sin(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-(L2.*sin(alpha12).*sin(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2+((((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2).*((Il_2xx.*((L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dz+L2.*dtheta2.*cos(theta2)))./l2-(L2.*sin(pi./2.0).*cos(theta2).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)+(Il_2xy.*((L2.*sin(pi./2.0).*cos(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dz+L2.*dtheta2.*cos(theta2)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2-(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)-(Il_2xz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2+(((Il_2zz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xz.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2yz.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*cos(alpha12).*sin(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-(L2.*sin(alpha12).*sin(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2-(((Il_2zz.*(((dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*(-1.0./2.0))./l2+(Il_2xz.*(((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-((z+L2.*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2))./(l2.*2.0)+(Il_2yz.*(((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*2.0)./l2+((((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xy.*((L2.*cos(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*cos(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2).*(-1.0./2.0))./l2+(Il_2yy.*((L2.*cos(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+(L2.*sin(alpha12).*sin(theta2).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2-(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2yz.*((L2.*cos(alpha12).*sin(theta2).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-(L2.*sin(alpha12).*sin(theta2).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2+((((z+L2.*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xy.*((L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dz+L2.*dtheta2.*cos(theta2)))./l2-(L2.*sin(pi./2.0).*cos(theta2).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(z+L2.*sin(theta2)))./l2+(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)+(Il_2yy.*((L2.*sin(pi./2.0).*cos(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dz+L2.*dtheta2.*cos(theta2)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2-(L2.*dtheta2.*sin(pi./2.0).*sin(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)-(Il_2yz.*((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)))./l2+(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)))./l2+(L2.*dtheta2.*sin(alpha12-pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2+(L2.*dtheta2.*cos(theta2).*cos(alpha12-pi./2.0).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)).*2.0)./l2+(((L2.*sin(pi./2.0).*cos(theta2).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2-(L2.*sin(alpha12-pi./2.0).*sin(theta2).*(z+L2.*sin(theta2)))./l2).*((Il_2xx.*(((dz+L2.*dtheta2.*cos(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)).*(z+L2.*sin(theta2)))./l2-((ddz-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)-(Il_2xy.*(((dz+L2.*dtheta2.*cos(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((z+L2.*sin(theta2)).*(ddx-L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)))./l2-((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2-((ddz-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2xz.*(((ddx-L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*2.0)./l2-(((L2.*sin(theta2).*cos(alpha12-pi./2.0).*(z+L2.*sin(theta2)))./l2+(L2.*sin(pi./2.0).*cos(theta2).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2).*((Il_2xy.*(((dz+L2.*dtheta2.*cos(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2+((ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)).*(z+L2.*sin(theta2)))./l2-((ddz-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2))./(l2.*2.0)-(Il_2yy.*(((dz+L2.*dtheta2.*cos(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2+((z+L2.*sin(theta2)).*(ddx-L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)))./l2-((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dz+L2.*dtheta2.*sin(pi./2.0).*cos(theta2)))./l2-((ddz-L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)+(Il_2yz.*(((ddx-L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0)).*(y+a2.*cos(alpha12)-e2.*cos(alpha12)+L2.*cos(alpha12).*cos(theta2)))./l2+((dy-L2.*dtheta2.*cos(alpha12).*sin(theta2)).*(dx-L2.*dtheta2.*sin(theta2).*cos(alpha12-pi./2.0)))./l2-((dx-L2.*dtheta2.*sin(alpha12).*sin(theta2)).*(dy+L2.*dtheta2.*sin(alpha12-pi./2.0).*sin(theta2)))./l2-((ddy+L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2)).*(x+a2.*sin(alpha12)-e2.*sin(alpha12)+L2.*cos(theta2).*sin(alpha12)))./l2))./(l2.*2.0)).*2.0)./l2-(L2.*g.*mL2.*cos(theta2))./2.0-L2.*g.*ml2.*cos(theta2)-(L2.*ml2.*sin(theta2).*cos(alpha12-pi./2.0).*(ddx./2.0+(L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0))./2.0))./2.0-L2.*ml2.*sin(theta2).*cos(alpha12-pi./2.0).*(ddx./4.0+(L2.*dtheta2.^2.*cos(theta2).*cos(alpha12-pi./2.0))./4.0)+(L2.*ml2.*sin(alpha12-pi./2.0).*sin(theta2).*(ddy./2.0-(L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2))./2.0))./2.0+L2.*ml2.*sin(alpha12-pi./2.0).*sin(theta2).*(ddy./4.0-(L2.*dtheta2.^2.*sin(alpha12-pi./2.0).*cos(theta2))./4.0)+(L2.*ml2.*sin(pi./2.0).*cos(theta2).*(ddz./2.0+(L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2))./2.0))./2.0+L2.*ml2.*sin(pi./2.0).*cos(theta2).*(ddz./4.0+(L2.*dtheta2.^2.*sin(pi./2.0).*sin(theta2))./4.0);-lambda3.*(L3.*cos(theta3).*(z+L3.*sin(theta3)).*-2.0+L3.*cos(alpha13).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)).*2.0+L3.*sin(alpha13).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)).*2.0)-((((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xz.*((L3.*cos(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*cos(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*(-1.0./2.0))./l3+(Il_3yz.*((L3.*cos(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+(L3.*sin(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3-(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3zz.*((L3.*cos(alpha13).*sin(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-(L3.*sin(alpha13).*sin(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3-((((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xz.*((L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dz+L3.*dtheta3.*cos(theta3)))./l3-(L3.*sin(pi./2.0).*cos(theta3).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)+(Il_3yz.*((L3.*sin(pi./2.0).*cos(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dz+L3.*dtheta3.*cos(theta3)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3-(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)-(Il_3zz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3+(((Il_3xz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3xx.*((L3.*sin(pi./2.0).*cos(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(z+L3.*sin(theta3)))./l3))./(l3.*2.0)-(Il_3xy.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*sin(pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*(((dz+L3.*dtheta3.*cos(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)).*(z+L3.*sin(theta3)))./l3-((ddz-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3-(((Il_3xz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xx.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3xy.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*cos(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*cos(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3+(((Il_3xz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xx.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3xy.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dz+L3.*dtheta3.*cos(theta3)))./l3-(L3.*sin(pi./2.0).*cos(theta3).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3+(((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xz.*(((dz+L3.*dtheta3.*cos(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)).*(z+L3.*sin(theta3)))./l3-((ddz-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)-(Il_3yz.*(((dz+L3.*dtheta3.*cos(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((z+L3.*sin(theta3)).*(ddx-L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)))./l3-((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3-((ddz-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3zz.*(((ddx-L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*2.0)./l3-(((Il_3yz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3xy.*((L3.*sin(pi./2.0).*cos(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(z+L3.*sin(theta3)))./l3))./(l3.*2.0)-(Il_3yy.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*sin(pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*(((dz+L3.*dtheta3.*cos(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((z+L3.*sin(theta3)).*(ddx-L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)))./l3-((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3-((ddz-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*2.0)./l3+(((Il_3zz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3xz.*((L3.*sin(pi./2.0).*cos(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(z+L3.*sin(theta3)))./l3))./(l3.*2.0)-(Il_3yz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*sin(pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*(((ddx-L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*2.0)./l3+(((Il_3yz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xy.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3yy.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*cos(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+(L3.*sin(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3-(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*2.0)./l3+(((Il_3yz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xy.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3yy.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*sin(pi./2.0).*cos(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dz+L3.*dtheta3.*cos(theta3)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3-(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*2.0)./l3+((((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3).*((Il_3xx.*((L3.*cos(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*cos(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*(-1.0./2.0))./l3+(Il_3xy.*((L3.*cos(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+(L3.*sin(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3-(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3xz.*((L3.*cos(alpha13).*sin(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-(L3.*sin(alpha13).*sin(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3+((((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3).*((Il_3xx.*((L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dz+L3.*dtheta3.*cos(theta3)))./l3-(L3.*sin(pi./2.0).*cos(theta3).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)+(Il_3xy.*((L3.*sin(pi./2.0).*cos(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dz+L3.*dtheta3.*cos(theta3)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3-(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)-(Il_3xz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3+(((Il_3zz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xz.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3yz.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*cos(alpha13).*sin(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-(L3.*sin(alpha13).*sin(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3-(((Il_3zz.*(((dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*(-1.0./2.0))./l3+(Il_3xz.*(((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-((z+L3.*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3))./(l3.*2.0)+(Il_3yz.*(((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*2.0)./l3+((((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xy.*((L3.*cos(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*cos(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3).*(-1.0./2.0))./l3+(Il_3yy.*((L3.*cos(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+(L3.*sin(alpha13).*sin(theta3).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3-(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3yz.*((L3.*cos(alpha13).*sin(theta3).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-(L3.*sin(alpha13).*sin(theta3).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3+((((z+L3.*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xy.*((L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dz+L3.*dtheta3.*cos(theta3)))./l3-(L3.*sin(pi./2.0).*cos(theta3).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(z+L3.*sin(theta3)))./l3+(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)+(Il_3yy.*((L3.*sin(pi./2.0).*cos(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dz+L3.*dtheta3.*cos(theta3)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3-(L3.*dtheta3.*sin(pi./2.0).*sin(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)-(Il_3yz.*((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)))./l3+(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)))./l3+(L3.*dtheta3.*sin(alpha13-pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3+(L3.*dtheta3.*cos(theta3).*cos(alpha13-pi./2.0).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)).*2.0)./l3+(((L3.*sin(pi./2.0).*cos(theta3).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3-(L3.*sin(alpha13-pi./2.0).*sin(theta3).*(z+L3.*sin(theta3)))./l3).*((Il_3xx.*(((dz+L3.*dtheta3.*cos(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)).*(z+L3.*sin(theta3)))./l3-((ddz-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)-(Il_3xy.*(((dz+L3.*dtheta3.*cos(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((z+L3.*sin(theta3)).*(ddx-L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)))./l3-((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3-((ddz-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3xz.*(((ddx-L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*2.0)./l3-(((L3.*sin(theta3).*cos(alpha13-pi./2.0).*(z+L3.*sin(theta3)))./l3+(L3.*sin(pi./2.0).*cos(theta3).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3).*((Il_3xy.*(((dz+L3.*dtheta3.*cos(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3+((ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)).*(z+L3.*sin(theta3)))./l3-((ddz-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3))./(l3.*2.0)-(Il_3yy.*(((dz+L3.*dtheta3.*cos(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3+((z+L3.*sin(theta3)).*(ddx-L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)))./l3-((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dz+L3.*dtheta3.*sin(pi./2.0).*cos(theta3)))./l3-((ddz-L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)+(Il_3yz.*(((ddx-L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0)).*(y+a3.*cos(alpha13)-e3.*cos(alpha13)+L3.*cos(alpha13).*cos(theta3)))./l3+((dy-L3.*dtheta3.*cos(alpha13).*sin(theta3)).*(dx-L3.*dtheta3.*sin(theta3).*cos(alpha13-pi./2.0)))./l3-((dx-L3.*dtheta3.*sin(alpha13).*sin(theta3)).*(dy+L3.*dtheta3.*sin(alpha13-pi./2.0).*sin(theta3)))./l3-((ddy+L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3)).*(x+a3.*sin(alpha13)-e3.*sin(alpha13)+L3.*cos(theta3).*sin(alpha13)))./l3))./(l3.*2.0)).*2.0)./l3-(L3.*g.*mL3.*cos(theta3))./2.0-L3.*g.*ml3.*cos(theta3)-(L3.*ml3.*sin(theta3).*cos(alpha13-pi./2.0).*(ddx./2.0+(L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0))./2.0))./2.0-L3.*ml3.*sin(theta3).*cos(alpha13-pi./2.0).*(ddx./4.0+(L3.*dtheta3.^2.*cos(theta3).*cos(alpha13-pi./2.0))./4.0)+(L3.*ml3.*sin(alpha13-pi./2.0).*sin(theta3).*(ddy./2.0-(L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3))./2.0))./2.0+L3.*ml3.*sin(alpha13-pi./2.0).*sin(theta3).*(ddy./4.0-(L3.*dtheta3.^2.*sin(alpha13-pi./2.0).*cos(theta3))./4.0)+(L3.*ml3.*sin(pi./2.0).*cos(theta3).*(ddz./2.0+(L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3))./2.0))./2.0+L3.*ml3.*sin(pi./2.0).*cos(theta3).*(ddz./4.0+(L3.*dtheta3.^2.*sin(pi./2.0).*sin(theta3))./4.0)];

        % for the 3d grid depicitons
        bit
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        function obj = GTPR(a1, a2, a3, L1, L2, L3, l1, l2, l3, e1, e2, e3, alpha12, alpha13, varargin)
        % varargin input can be multiple ones as well.
        % Input should be in the following order, if not empty:
        % a1, a2, a3, L1, L2, L3, l1, l2, l3, e1, e2, e3, alpha12, alpha13
        
            if contains(varargin, 'zeroTolerance')
               idx = find(contains(varargin, 'zeroTolerance'));
               obj.zeroTolerance = varargin{idx+1};
            else
                obj.zeroTolerance = 1e-7;
            end
        
            obj.a1 = a1;
            obj.a2 = a2;
            obj.a3 = a3;

            obj.e1 = e1;
            obj.e2 = e2;
            obj.e3 = e3;

            obj.L1 = L1;
            obj.L2 = L2;
            obj.L3 = L3;

            obj.l1 = l1;
            obj.l2 = l2;
            obj.l3 = l3;

            obj.alpha12 = alpha12;
            obj.alpha13 = alpha13;
            
            obj.getBaseTrianglePts(a1, a2, a3, alpha12, alpha13);
            
            %Double BitMap (0-64 bits are present in a double)
            obj.bit(1).Name = "W";
            obj.bit(1).Position = 1;
            obj.bit(1).Color = [1 0 1];
            obj.bit(1).Cursor = 's';
            obj.bit(1).Opacity = 0.5;

            obj.bit(2).Name = "A1";
            obj.bit(2).Position = 2;
            obj.bit(2).Color = [1 0 0];
            obj.bit(2).Cursor = 'o';
            obj.bit(2).Opacity = 0.5;

            obj.bit(3).Name = "A2";
            obj.bit(3).Position = 3;
            obj.bit(3).Color = [0 1 0];
            obj.bit(3).Cursor = 'o';
            obj.bit(3).Opacity = 0.5;

            obj.bit(4).Name = "A3";
            obj.bit(4).Position = 4;
            obj.bit(4).Color = [0 0 1];
            obj.bit(4).Cursor = 'o';
            obj.bit(4).Opacity = 0.5;

            obj.bit(5).Name = "B1";
            obj.bit(5).Position = 5;
            obj.bit(5).Color = [1 0 0];
            obj.bit(5).Cursor = 'x';
            obj.bit(5).Opacity = 0.5;

            obj.bit(6).Name = "B2";
            obj.bit(6).Position = 6;
            obj.bit(6).Color = [0 1 0];
            obj.bit(6).Cursor = 'x';
            obj.bit(6).Opacity = 0.5;

            obj.bit(7).Name = "B3";
            obj.bit(7).Position = 7;
            obj.bit(7).Color = [0 0 1];
            obj.bit(7).Cursor = 'x';
            obj.bit(7).Opacity = 0.5;

            obj.bit(8).Name = "C1";
            obj.bit(8).Position = 8;
            obj.bit(8).Color = [1 0 0];
            obj.bit(8).Cursor = '+';
            obj.bit(8).Opacity = 0.5;

            obj.bit(9).Name = "C2";
            obj.bit(9).Position = 9;
            obj.bit(9).Color = [0 1 0];
            obj.bit(9).Cursor = '+';
            obj.bit(9).Opacity = 0.5;

            obj.bit(10).Name = "C3";
            obj.bit(10).Position = 10;
            obj.bit(10).Color = [0 0 1];
            obj.bit(10).Cursor = '+';
            obj.bit(10).Opacity = 0.5;

            obj.bit(11).Name = "O";
            obj.bit(11).Position = 11;
            obj.bit(11).Color = [1 1 1];
            obj.bit(11).Cursor = 'o';
            obj.bit(11).Opacity = 0.5;

            obj.bit(12).Name = "A1A2";
            obj.bit(12).Position = 12;
            obj.bit(12).Color = [1 1 1];
            obj.bit(12).lineStyle = '-';
            obj.bit(12).Opacity = 0.5;

            obj.bit(13).Name = "A2A3";
            obj.bit(13).Position = 13;
            obj.bit(13).Color = [1 1 1];
            obj.bit(13).lineStyle = '-';
            obj.bit(13).Opacity = 0.5;

            obj.bit(14).Name = "A1A3";
            obj.bit(14).Position = 14;
            obj.bit(14).Color = [1 1 1];
            obj.bit(14).lineStyle = '-';
            obj.bit(14).Opacity = 0.5;

            obj.bit(15).Name = "A1B1";
            obj.bit(15).Position = 15;
            obj.bit(15).Color = [1 0 0];
            obj.bit(15).lineStyle = '-';
            obj.bit(15).Opacity = 0.2;

            obj.bit(16).Name = "A2B2";
            obj.bit(16).Position = 16;
            obj.bit(16).Color = [0 1 0];
            obj.bit(16).lineStyle = '-';
            obj.bit(16).Opacity = 0.2;

            obj.bit(17).Name = "A3B3";
            obj.bit(17).Position = 17;
            obj.bit(17).Color = [0 0 1];
            obj.bit(17).lineStyle = '-';
            obj.bit(17).Opacity = 0.2;

            obj.bit(18).Name = "B1C1";
            obj.bit(18).Position = 18;
            obj.bit(18).Color = [1 0 0];
            obj.bit(18).lineStyle = '-';
            obj.bit(18).Opacity = 0.2;

            obj.bit(19).Name = "B2C2";
            obj.bit(19).Position = 19;
            obj.bit(19).Color = [0 1 0];
            obj.bit(19).lineStyle = '-';
            obj.bit(19).Opacity = 0.2;

            obj.bit(20).Name = "B3C3";
            obj.bit(20).Position = 20;
            obj.bit(20).Color = [0 0 1];
            obj.bit(20).lineStyle = '-';
            obj.bit(20).Opacity = 0.2;

            obj.bit(21).Name = "C1C2";
            obj.bit(21).Position = 21;
            obj.bit(21).Color = [0 0 0];
            obj.bit(21).lineStyle = '-';
            obj.bit(21).Opacity = 0.2;

            obj.bit(22).Name = "C2C3";
            obj.bit(22).Position = 22;
            obj.bit(22).Color = [0 0 0];
            obj.bit(22).lineStyle = '-';
            obj.bit(22).Opacity = 0.2;

            obj.bit(23).Name = "C1C3";
            obj.bit(23).Position = 23;
            obj.bit(23).Color = [0 0 0];
            obj.bit(23).lineStyle = '-';
            obj.bit(23).Opacity = 0.2;

            obj.bit(24).Name = "B1L";
            obj.bit(24).Position = 24;
            obj.bit(24).Color = [1 0 0];
            obj.bit(24).Cursor = 'x';
            obj.bit(24).Opacity = 0.5;

            obj.bit(25).Name = "B1R";
            obj.bit(25).Position = 25;
            obj.bit(25).Color = [1 0 0];
            obj.bit(25).Cursor = 'x';
            obj.bit(25).Opacity = 0.5;

            obj.bit(26).Name = "B2L";
            obj.bit(26).Position = 26;
            obj.bit(26).Color = [0 1 0];
            obj.bit(26).Cursor = 'x';
            obj.bit(26).Opacity = 0.5;

            obj.bit(27).Name = "B2R";
            obj.bit(27).Position = 27;
            obj.bit(27).Color = [0 1 0];
            obj.bit(27).Cursor = 'x';
            obj.bit(27).Opacity = 0.5;

            obj.bit(28).Name = "B3L";
            obj.bit(28).Position = 28;
            obj.bit(28).Color = [0 0 1];
            obj.bit(28).Cursor = 'x';
            obj.bit(28).Opacity = 0.5;

            obj.bit(29).Name = "B3R";
            obj.bit(29).Position = 29;
            obj.bit(29).Color = [0 0 1];
            obj.bit(29).Cursor = 'x';
            obj.bit(29).Opacity = 0.5;

            obj.bit(30).Name = "C1L";
            obj.bit(30).Position = 30;
            obj.bit(30).Color = [1 0 0];
            obj.bit(30).Cursor = '+';
            obj.bit(30).Opacity = 0.5;

            obj.bit(31).Name = "C1R";
            obj.bit(31).Position = 31;
            obj.bit(31).Color = [1 0 0];
            obj.bit(31).Cursor = '+';
            obj.bit(31).Opacity = 0.5;

            obj.bit(32).Name = "C2L";
            obj.bit(32).Position = 32;
            obj.bit(32).Color = [0 1 0];
            obj.bit(32).Cursor = '+';
            obj.bit(32).Opacity = 0.5;

            obj.bit(33).Name = "C2R";
            obj.bit(33).Position = 33;
            obj.bit(33).Color = [0 1 0];
            obj.bit(33).Cursor = '+';
            obj.bit(33).Opacity = 0.5;

            obj.bit(34).Name = "C3L";
            obj.bit(34).Position = 34;
            obj.bit(34).Color = [0 0 1];
            obj.bit(34).Cursor = '+';
            obj.bit(34).Opacity = 0.5;

            obj.bit(35).Name = "C3R";
            obj.bit(35).Position = 35;
            obj.bit(35).Color = [0 0 1];
            obj.bit(35).Cursor = '+';
            obj.bit(35).Opacity = 0.5;

            obj.bit(36).Name = "B1LC1L";
            obj.bit(36).Position = 36;
            obj.bit(36).Color = [1 0 0];
            obj.bit(36).lineStyle = '-';
            obj.bit(36).Opacity = 0.2;

            obj.bit(37).Name = "B1RC1R";
            obj.bit(37).Position = 37;
            obj.bit(37).Color = [1 0 0];
            obj.bit(37).lineStyle = '-';
            obj.bit(37).Opacity = 0.2;
            
            obj.bit(38).Name = "B2LC2L";
            obj.bit(38).Position = 38;
            obj.bit(38).Color = [0 1 0];
            obj.bit(38).lineStyle = '-';
            obj.bit(38).Opacity = 0.2;
            
            obj.bit(39).Name = "B2RC2R";
            obj.bit(39).Position = 39;
            obj.bit(39).Color = [0 1 0];
            obj.bit(39).lineStyle = '-';
            obj.bit(39).Opacity = 0.2;
            
            obj.bit(40).Name = "B3LC3L";
            obj.bit(40).Position = 40;
            obj.bit(40).Color = [0 0 1];
            obj.bit(40).lineStyle = '-';
            obj.bit(40).Opacity = 0.2;

            obj.bit(41).Name = "B3RC3R";
            obj.bit(41).Position = 41;
            obj.bit(41).Color = [0 0 1];
            obj.bit(41).lineStyle = '-';
            obj.bit(41).Opacity = 0.2;
            
            obj.bit(42).Name = "B1LB1R";
            obj.bit(42).Position = 42;
            obj.bit(42).Color = [1 0 0];
            obj.bit(42).lineStyle = '-';
            obj.bit(42).Opacity = 0.2;
            
            obj.bit(43).Name = "B2LB2R";
            obj.bit(43).Position = 43;
            obj.bit(43).Color = [0 1 0];
            obj.bit(43).lineStyle = '-';
            obj.bit(43).Opacity = 0.2;
            
            obj.bit(44).Name = "B3LB3R";
            obj.bit(44).Position = 44;
            obj.bit(44).Color = [0 0 1];
            obj.bit(44).lineStyle = '-';
            obj.bit(44).Opacity = 0.2;
            
            obj.bit(45).Name = "C1LC1R";
            obj.bit(45).Position = 45;
            obj.bit(45).Color = [1 0 0];
            obj.bit(45).lineStyle = '-';
            obj.bit(45).Opacity = 0.2;
            
            obj.bit(46).Name = "C2LC2R";
            obj.bit(46).Position = 46;
            obj.bit(46).Color = [0 1 0];
            obj.bit(46).lineStyle = '-';
            obj.bit(46).Opacity = 0.2;
            
            obj.bit(47).Name = "C3LC3R";
            obj.bit(47).Position = 47;
            obj.bit(47).Color = [0 0 1];
            obj.bit(47).lineStyle = '-';
            obj.bit(47).Opacity = 0.2;
            
            obj.bit(48).Name = "WorkTriangle";
            obj.bit(48).Position = 48;
            obj.bit(48).Color = [1 1 0];
            obj.bit(48).Opacity = 0.2;

            obj.bit(49).Name = "BaseTriangle";
            obj.bit(49).Position = 49;
            obj.bit(49).Color = [1 1 1];
            obj.bit(49).Opacity = 0.2;
        end
        
        function W = directGeometry(obj, theta, varargin)
           % Forward geometry solution for the GTPR
           % Input: (theta)
           
%             obj.theta = theta;
            
            if size(theta,2) ~= 3
                theta = theta';
            end
            
            obj.theta = theta;
            
		    obj.B1 = [   0   ;
                    (-obj.a1 - obj.L1 * cos(theta(1)))'  ;
                    (-obj.L1*sin(theta(1)))'];
		   
		    obj.B2 =    [   -obj.a2*sin(obj.alpha12) - obj.L2*sin(obj.alpha12)*cos(theta(2));
                            -obj.a2*cos(obj.alpha12) - obj.L2*cos(obj.alpha12)*cos(theta(2));
                        -obj.L2*sin(theta(2))];
		   
			obj.B3 =    [   -obj.a3*sin(obj.alpha13) - obj.L3*sin(obj.alpha13)*cos(theta(3));
                            -obj.a3*cos(obj.alpha13) - obj.L3*cos(obj.alpha13)*cos(theta(3));
                            -obj.L3*sin(theta(3))];
							
			
							
           obj.B1v =    obj.B1 + [0; obj.e1; 0];

           obj.B2v =    obj.B2 + [obj.e2*sin(obj.alpha12);
                                  obj.e2*cos(obj.alpha12);
                                  0];
                    
           obj.B3v =    obj.B3 + [obj.e3*sin(obj.alpha13);
                                  obj.e3*cos(obj.alpha13);
                                  0];
           
           [~, obj.W] = obj.threeSpheresIntersection(obj.B1v,obj.B2v,obj.B3v,obj.l1,obj.l2,obj.l3);     
           
           if ~isreal(obj.W)
               W = [];
               obj.W = [];
               return;
           end

           if size(obj.W,1) ~= 3
               obj.W = obj.W';
           end

           obj.getWorkTrianglePts(obj.W, obj.e1, obj.e2, obj.e3, obj.alpha12, obj.alpha13);

           W = obj.W;
        end
        
        function theta = inverseGeometry(obj, W)
            
            

            if size(W,1) ~= 3
                W = W';
            end
            
            obj.W = W;
            
            obj.getWorkTrianglePts(obj.W, obj.e1, obj.e2, obj.e3, obj.alpha12, obj.alpha13); 
            
            nUnit1 = cross([0;0;0]-obj.A1, obj.A1-[0;0;obj.L1]);
            nUnit1 = nUnit1 / norm(nUnit1);
            B1 = obj.sphereCircleIntersection(obj.A1, obj.L1, nUnit1, obj.C1, obj.l1); %#ok<PROPLC>
            
            nUnit2 = cross([0;0;0]-obj.A2, obj.A2-[0;0;obj.L2]);
            nUnit2 = nUnit2 / norm(nUnit2);
            B2 = obj.sphereCircleIntersection(obj.A2, obj.L2, nUnit2, obj.C2, obj.l2); %#ok<PROPLC>
            
            nUnit3 = cross([0;0;0]-obj.A3, obj.A3-[0;0;obj.L3]);
            nUnit3 = nUnit3 / norm(nUnit3);
            B3 = obj.sphereCircleIntersection(obj.A3, obj.L3, nUnit3, obj.C3, obj.l3); %#ok<PROPLC>
            
            if ~isempty(B1.p0) && ~isempty(B2.p0) && ~isempty(B3.p0)
                obj.B1 = B1.p0; %#ok<PROPLC>
                obj.B2 = B2.p0; %#ok<PROPLC>
                obj.B3 = B3.p0; %#ok<PROPLC>
            else
                theta = [complex(0,0);complex(0,0);complex(0,0)];
                return
            end
            OA1 = [0;0;0]-obj.A1;
            B1A1 = obj.B1-obj.A1;
            theta(1) = pi - acos(dot(B1A1,OA1)/(norm(OA1)*norm(B1A1)));
            
            if ~isreal(theta(1)) && abs(imag(theta(1))) < obj.zeroTolerance
                theta(1) = real(theta(1));
            end
            
            OA2 = [0;0;0]-obj.A2;
            B2A2 = obj.B2-obj.A2;
            theta(2) = pi - acos(dot(B2A2,OA2)/(norm(OA2)*norm(B2A2)));
            
            if ~isreal(theta(2)) && abs(imag(theta(2))) < obj.zeroTolerance
                theta(2) = real(theta(2));
            end
            
            OA3 = [0;0;0]-obj.A3;
            B3A3 = obj.B3-obj.A3;
            theta(3) = pi - acos(dot(B3A3,OA3)/(norm(OA3)*norm(B3A3)));

            if ~isreal(theta(3)) && abs(imag(theta(3))) < obj.zeroTolerance
                theta(3) = real(theta(3));
            end
            
            
            if B1A1(3) > 0
                theta(1) = -theta(1);
            end
            
            if B2A2(3) > 0
                theta(2) = -theta(2);
            end
            
            if B3A3(3) > 0
                theta(3) = -theta(3);
            end
            
            theta = [theta(1);theta(2);theta(3)];
            obj.theta = theta;
        end
        
        function S = calculateScrews(obj, varargin)
			% Getting the Screws of the system:
			% Arm A:
			%obj.S01_1 = zeros(6,1);
            
            if size(obj.A1,1) ~= 3
                obj.A1 = obj.A1';
            end
            
            if size(obj.A2,1) ~= 3
                obj.A2 = obj.A2';
            end
            
            if size(obj.A3,1) ~= 3
                obj.A3 = obj.A3';
            end
            
            
            if size(obj.B1,1) ~= 3
                obj.B1 = obj.B1';
            end
            
            if size(obj.B2,1) ~= 3
                obj.B2 = obj.B2';
            end
            
            if size(obj.B3,1) ~= 3
                obj.B3 = obj.B3';
            end
            
            
            if size(obj.C1,1) ~= 3
                obj.C1 = obj.C1';
            end
            
            if size(obj.C2,1) ~= 3
                obj.C2 = obj.C2';
            end
            
            if size(obj.C3,1) ~= 3
                obj.C3 = obj.C3';
            end
            

            % Arm System 1:
            obj.S12_1 = zeros(6,1);
            obj.S23_1 = zeros(6,1);
            obj.S34_1 = zeros(6,1);
            obj.S45_1 = zeros(6,1);
            obj.S56_1 = zeros(6,1);
            
            obj.S12_1(1:3) = obj.rotz(pi/2) * ((obj.A1)/norm(obj.A1));
            obj.S23_1(1:3) = obj.S12_1(1:3);
            obj.S34_1(1:3) = cross(obj.S23_1(1:3),(obj.C1-obj.B1)); 
            obj.S34_1(1:3) = obj.S34_1(1:3)/norm(obj.S34_1(1:3));
            obj.S45_1(1:3) = obj.S23_1(1:3);
            obj.S56_1(1:3) = obj.S34_1(1:3);

            obj.S12_1(4:6) = cross(obj.A1, obj.S12_1(1:3));
            obj.S23_1(4:6) = cross(obj.B1, obj.S23_1(1:3));
            obj.S34_1(4:6) = cross(obj.B1, obj.S34_1(1:3));
            obj.S45_1(4:6) = cross(obj.C1, obj.S45_1(1:3));
            obj.S56_1(4:6) = cross(obj.C1, obj.S56_1(1:3));

            obj.S01_1 = zeros(6,1);
            
            % Arm System 2:
            obj.S12_2 = zeros(6,1);
            obj.S23_2 = zeros(6,1);
            obj.S34_2 = zeros(6,1);
            obj.S45_2 = zeros(6,1);
            obj.S56_2 = zeros(6,1);

            obj.S12_2(1:3) = obj.rotz(pi/2) * ((obj.A2)/norm(obj.A2));
            obj.S23_2(1:3) = obj.S12_2(1:3);
            obj.S34_2(1:3) = cross(obj.S23_2(1:3), (obj.C2-obj.B2)); 
            obj.S34_2(1:3) = obj.S34_2(1:3)/norm(obj.S34_2(1:3));
            obj.S45_2(1:3) = obj.S23_2(1:3);
            obj.S56_2(1:3) = obj.S34_2(1:3);

            obj.S12_2(4:6) = cross(obj.A2, obj.S12_2(1:3));
            obj.S23_2(4:6) = cross(obj.B2, obj.S23_2(1:3));
            obj.S34_2(4:6) = cross(obj.B2, obj.S34_2(1:3));
            obj.S45_2(4:6) = cross(obj.C2, obj.S45_2(1:3));
            obj.S56_2(4:6) = cross(obj.C2, obj.S56_2(1:3));

            obj.S01_2 = zeros(6,1);

            % Arm System 3:
            obj.S12_3 = zeros(6,1);
            obj.S23_3 = zeros(6,1);
            obj.S34_3 = zeros(6,1);
            obj.S45_3 = zeros(6,1);
            obj.S56_3 = zeros(6,1);

            obj.S12_3(1:3) = obj.rotz(pi/2) * ((obj.A3)/norm(obj.A3));
            obj.S23_3(1:3) = obj.S12_3(1:3);
            obj.S34_3(1:3) = cross(obj.S23_3(1:3), (obj.C3-obj.B3)); 
            obj.S34_3(1:3) = obj.S34_3(1:3)/norm(obj.S34_3(1:3));
            obj.S45_3(1:3) = obj.S23_3(1:3);
            obj.S56_3(1:3) = obj.S34_3(1:3);

            obj.S12_3(4:6) = cross(obj.A3, obj.S12_3(1:3));
            obj.S23_3(4:6) = cross(obj.B3, obj.S23_3(1:3));
            obj.S34_3(4:6) = cross(obj.B3, obj.S34_3(1:3));
            obj.S45_3(4:6) = cross(obj.C3, obj.S45_3(1:3));
            obj.S56_3(4:6) = cross(obj.C3, obj.S56_3(1:3));

            obj.S01_3 = zeros(6,1);


            % Reciprocal Screws (lower arms, passive joints):
            % Arm system 1
            obj.Sl1 = zeros(6,1);
            obj.Sl1(1:3) = (obj.C1-obj.B1) / norm(obj.C1-obj.B1);
            obj.Sl1(4:6) = cross(obj.C1, obj.Sl1(1:3));

            % Arm System 2
            obj.Sl2 = zeros(6,1);
            obj.Sl2(1:3) = (obj.C2-obj.B2) / norm(obj.C2-obj.B2);
            obj.Sl2(4:6) = cross(obj.B2, obj.Sl2(1:3));

            % Arm System 3
            obj.Sl3 = zeros(6,1);
            obj.Sl3(1:3) = (obj.C3-obj.B3) / norm(obj.C3-obj.B3);
            obj.Sl3(4:6) = cross(obj.C3, obj.Sl3(1:3));
            
            obj.J_A = [obj.Sl1(1:3), obj.Sl2(1:3), obj.Sl3(1:3)]';
            BB1 = obj.kleinForm(obj.S12_1,obj.Sl1);
            BB2 = obj.kleinForm(obj.S12_2,obj.Sl2);
            BB3 = obj.kleinForm(obj.S12_3,obj.Sl3);
            obj.J_B = diag([BB1, BB2, BB3]);

            % Screws of upper arms:
            obj.SL1 = zeros(6,1);
            obj.SL1(1:3) = (obj.B1-obj.A1)/(norm(obj.B1-obj.A1));
            obj.SL1(4:6) = cross(obj.A1,obj.SL1(1:3));

            obj.SL2 = zeros(6,1);
            obj.SL2(1:3) = (obj.B2-obj.A2)/(norm(obj.B2-obj.A2));
            obj.SL2(4:6) = cross(obj.A2,obj.SL2(1:3));

            obj.SL3 = zeros(6,1);
            obj.SL3(1:3) = (obj.B3-obj.A3)/(norm(obj.B3-obj.A3));
            obj.SL3(4:6) = cross(obj.A3,obj.SL3(1:3));
			
			S.S01_1 = obj.S01_1;
			S.S12_1 = obj.S12_1;
			S.S23_1 = obj.S23_1;
			S.S34_1 = obj.S34_1;
			S.S45_1 = obj.S45_1;
			S.S56_1 = obj.S56_1;
			S.SL1 = obj.SL1;
			
			S.S01_2 = obj.S01_2;
			S.S12_2 = obj.S12_2;
			S.S23_2 = obj.S23_2;
			S.S34_2 = obj.S34_2;
			S.S45_2 = obj.S45_2;
			S.S56_2 = obj.S56_2;
			S.SL2 = obj.SL2;

			S.S01_3 = obj.S01_3;
			S.S12_3 = obj.S12_3;
			S.S23_3 = obj.S23_3;
			S.S34_3 = obj.S34_3;
			S.S45_3 = obj.S45_3;
			S.S56_3 = obj.S56_3;
			S.SL3 = obj.SL3;

			S.J_A = obj.J_A;
			S.J_B = obj.J_B;
%             a = inv(S.J_A)*S.J_B;
        end
        
        
        
        function J = Jacobian(obj, varargin)
            % Returns the Jacobian of the current configuration
            % The default convention is used in a default setup:
            % S = [w;v];
            
			
			
			%if any(strcmp(varargin, 'O')) % To change the viewpoint of the system from the origin
			%	idx = find(strcmp(varargin,'O'));
			%	O = varargin{idx+1};
			%end
            
			obj.calculateScrews(varargin);
			
			J = obj.J_A\obj.J_B;

%             J(J<1e-6) = 0;
            
            obj.J = J;
            
			obj.Jinv = obj.J_B\obj.J_A;

%             obj.Jinv(obj.Jinv<1e-6) = 0;
        end
        
        function setDynamicParameters(obj, mL1, mL2, mL3, ml1, ml2, ml3, mW, IL0_1, IL0_2, IL0_3, Il0_1, Il0_2, Il0_3, g)
        
            
            
            obj.mL1 = mL1;
            obj.mL2 = mL2;
            obj.mL3 = mL3;
            obj.ml1 = ml1;
            obj.ml2 = ml2;
            obj.ml3 = ml3;
            obj.mW = mW;
            obj.IL0_1 = IL0_1;
            obj.IL0_2 = IL0_2;
            obj.IL0_3 = IL0_3;
            obj.Il0_1 = Il0_1;
            obj.Il0_2 = Il0_2;
            obj.Il0_3 = Il0_3;
            obj.g = g;
            
        end
        
        function [tau, ddtheta] = inverseDynamics(obj, v,a,varargin)
            if any(strcmp(varargin, 'fext'))
               fextidx = find(strcmp(varargin, 'fext'));
               Fext = varargin{fextidx+1};
            else
                if any(strcmp(varargin, 'Fext'))
                   fextidx = find(strcmp(varargin, "Fext"));
                   Fext = varargin{fextidx+1};
                else
                    Fext = [0;0;0];
                end
            end
            
            if any(strcmp(varargin, 'theta'))
               thetaidx = find(varargin, 'theta');
               obj.theta = varargin{thetaidx+1};
               obj.directGeometry(obj.theta);
            end
            obj.a = a;
            obj.dtheta = obj.inverseKinematics(v);
            obj.ddtheta = obj.inverseKinematics_Acceleration(a,v);
            
            % make a function for calculating w34x
            % Upper arms velocity:
            obj.Jacobian_arms_acceleration(v);

            
            obj.dynamicModelParameters(Fext);
            
            obj.tau = obj.Mddtheta - obj.dGammasdtheta*obj.lambda;
            
            tau = obj.tau;
                                                        
            ddtheta = obj.ddtheta;
        end
        
        function [ddtheta, a] = directDynamics(obj, v, tau, varargin)
            if any(strcmp(varargin, 'fext'))
               fextidx = find(varargin, 'fext');
               fext = varargin{fextidx+1};
            else
                fext = [0;0;0];
            end
            
            if any(strcmp(varargin, 'theta'))
               qidx = find(varargin, 'theta');
               obj.theta = varargin{qidx+1};
               obj.directGeometry(obj.theta);
            end
            
            obj.dtheta = obj.inverseKinematics(v);
            
            obj.calculateArmsInertiaMatrix(); % upper and lower arms inertia matrix calculator
            error("check! not finished!");
            % make a function for calculating w34x
            % Upper arms velocity:
            [~, ~, ~, ~, phiA, phiB, phiC] = obj.Jacobian_arms_acceleration(v);

            obj.w23A = phiA(3);
            obj.w23B = phiB(3);
            obj.w23C = phiC(3);

            obj.w34A = phiA(4);
            obj.w34B = phiB(4);
            obj.w34C = phiC(4);
           
            obj.dynamicModelParameters();
            
            ddtheta = obj.M\(tau - obj.C*obj.dtheta - obj.taug + obj.J'*fext);
            
            a = obj.directKinematics_Acceleration(obj.dtheta, ddtheta);
            
        end
            
            
        
        function calculateArmsInertiaMatrix(obj)
            % upper and lower arms inertia matrix calculator
            
            [obj.IL_1, obj.IL_2, obj.IL_3] = obj.upperArmsInertiaMatrix(obj.mL1, obj.mL2, obj.mL3, obj.IL0_1, obj.IL0_2, obj.IL0_3);
            obj.arms_rotationalMatrices(); % with rotational matrices
            [obj.Il_1, obj.Il_2, obj.Il_3] = obj.lowerArmsInertiaMatrix(obj.ml1, obj.ml2, obj.ml3, obj.Il0_1, obj.Il0_2, obj.Il0_3);
        
        end
        
        function dynamicModelParameters(obj, Fext)
            % based on: http://www.scholarpedia.org/article/Robot_dynamics
            % https://www.sciencedirect.com/science/article/pii/S089812211630147X
            % The base equation is the following: (M can also be as H):
            % tau = H(theta)*ddtheta + C(theta,dtheta)*dtheta + taug(theta) + J(theta)'*fext
            % tau: moment needed by the motors
            % H(theta): Hessian, mass matrix (M(theta))
            % C(theta,dtheta): Coriolis forces matrix
            % taug(theta): moments due to gravitational effects
            % J(theta)'*fext: external forces acting on W
            % J(theta)': transpose of Jacobian of the system
            
            
            obj.calculateArmsInertiaMatrix(); % upper and lower arms inertia matrix calculator
            
            obj.dGammasdxyz = obj.dGammasdxyzf(obj.L1,obj.L2,obj.L3,obj.a1,obj.a2,obj.a3,0,obj.alpha12,obj.alpha13,obj.e1,obj.e2,obj.e3,obj.theta(1),obj.theta(2),obj.theta(3),obj.W(1),obj.W(2),obj.W(3));
            obj.Mxyz = obj.Mxyzf(obj.Il0_1(1,1),obj.Il0_1(1,2),obj.Il0_1(1,3),obj.Il0_1(2,2),obj.Il0_1(2,3),obj.Il0_1(3,3),obj.Il0_2(1,1),obj.Il0_2(1,2),obj.Il0_2(1,3),obj.Il0_2(2,2),obj.Il0_2(2,3),obj.Il0_2(3,3),obj.Il0_3(1,1),obj.Il0_3(1,2),obj.Il0_3(1,3),obj.Il0_3(2,2),obj.Il0_3(2,3),obj.Il0_3(3,3),obj.L1,obj.L2,obj.L3,obj.a1,obj.a2,obj.a3,0,obj.alpha12,obj.alpha13,obj.ddtheta(1),obj.ddtheta(2),obj.ddtheta(3),obj.a(1),obj.a(2),obj.a(3),obj.dtheta(1),obj.dtheta(2),obj.dtheta(3),obj.v(1),obj.v(2),obj.v(3),obj.e1,obj.e2,obj.e3,obj.g,obj.l1,obj.l2,obj.l3,obj.mW,obj.ml1,obj.ml2,obj.ml3,pi,obj.theta(1),obj.theta(2),obj.theta(3),obj.W(1),obj.W(2),obj.W(3));
            obj.Mddtheta = obj.Mddthetaf(obj.IL0_1(1,1),obj.IL0_1(1,2),obj.IL0_1(2,2),obj.IL0_2(1,1),obj.IL0_2(1,2),obj.IL0_2(2,2),obj.IL0_3(1,1),obj.IL0_3(1,2),obj.IL0_3(2,2),obj.Il0_1(1,1),obj.Il0_1(1,2),obj.Il0_1(1,3),obj.Il0_1(2,2),obj.Il0_1(2,3),obj.Il0_1(3,3),obj.Il0_2(1,1),obj.Il0_2(1,2),obj.Il0_2(1,3),obj.Il0_2(2,2),obj.Il0_2(2,3),obj.Il0_2(3,3),obj.Il0_3(1,1),obj.Il0_3(1,2),obj.Il0_3(1,3),obj.Il0_3(2,2),obj.Il0_3(2,3),obj.Il0_3(3,3),obj.L1,obj.L2,obj.L3,obj.a1,obj.a2,obj.a3,0,obj.alpha12,obj.alpha13,obj.ddtheta(1),obj.ddtheta(2),obj.ddtheta(3),obj.a(1),obj.a(2),obj.a(3),obj.dtheta(1),obj.dtheta(2),obj.dtheta(3),obj.v(1),obj.v(2),obj.v(3),obj.e1,obj.e2,obj.e3,obj.g,obj.l1,obj.l2,obj.l3,obj.mL1,obj.mL2,obj.mL3,obj.ml1,obj.ml2,obj.ml3,pi,obj.theta(1),obj.theta(2),obj.theta(3),obj.W(1),obj.W(2),obj.W(3));
            obj.dGammasdtheta = obj.dGammasdthetaf(obj.L1,obj.L2,obj.L3,obj.a1,obj.a2,obj.a3,0,obj.alpha12,obj.alpha13,obj.e1,obj.e2,obj.e3,obj.theta(1),obj.theta(2),obj.theta(3),obj.W(1),obj.W(2),obj.W(3));

     

            obj.lambda = obj.dGammasdxyz\(obj.Mxyz - obj.J'*Fext);
            
        end
        
        
        function draw_VelocityEllipsoid(obj, varargin)
            % varargin can be:
            % 'Jacobian' => Input jacobian
            % 'center' => 3x1 vector giving the ellipsoid center point,
            %                default is 0
            % 'ratio' => a number, which magnifies or shrinks the plotted
            %                ellipsoid
            % 'gridSize => default is 20, number of points to calculate for 
            %                the ellipse
           
            if any(strcmp(varargin, 'id'))
                idx = find(strcmp(varargin, 'id'));
                id = varargin{idx+1};
            else
                id = 1; 
            end
            
            if any(strcmp(varargin, 'Jacobian'))
                idx = find(strcmp(varargin, 'Jacobian'));
                Jacobian = varargin{idx+1};
                obj.getVelocityfs('Jacobian', Jacobian);
            elseif obj.emptyJacobian == true    
                obj.getVelocityEigens();
            end
            
            if any(strcmp(varargin, 'center'))
                idx = find(strcmp(varargin, 'center'));
                center = varargin{idx+1};
            else
                center = zeros(3,1); 
            end
            
            if any(strcmp(varargin, 'ratio'))
                idx = find(strcmp(varargin, 'ratio'));
                ratio = varargin{idx+1};
            else
                ratio = 1; 
            end
            
            if any(strcmp(varargin, 'gridSize'))
                idx = find(strcmp(varargin, 'gridSize'));
                N = varargin{idx+1};
            else
                N = 20; 
            end
            
            if any(strcmp(varargin, 'FaceColor'))
                idx = find(strcmp(varargin, 'FaceColor'));
                faceColor = varargin{idx+1};
            else
                faceColor = 'b'; 
            end
            
            if any(strcmp(varargin, 'EdgeColor'))
                idx = find(strcmp(varargin, 'EdgeColor'));
                edgeColor = varargin{idx+1};
            else
                edgeColor = 'w'; 
            end
            
            if any(strcmp(varargin, 'FaceAlpha'))
                idx = find(strcmp(varargin, 'FaceAlpha'));
                faceAlpha = varargin{idx+1};
            else
                faceAlpha = 0.5; 
            end
              
            if any(strcmp(varargin, 'rotate'))
                idx = find(strcmp(varargin, 'rotate'));
                rotate = varargin{idx+1};
                switch rotate
                    case 'on'
                        rotate = true;
                    otherwise
                        rotate = false;
                end
            else
                rotate = true; 
            end
            
            a = 1\sqrt(obj.eigValues(1,1,id)); %#ok<PROPLC>
            b = 1\sqrt(obj.eigValues(2,2,id));
            c = 1\sqrt(obj.eigValues(3,3,id));
            
            if imag(a) ~= 0%#ok<PROPLC>
                a = real(a);%#ok<PROPLC>
            end
            
            if imag(b) ~= 0
                b = real(b);
            end
            
            if imag(c) ~= 0
                c = real(c);
            end
            
            if size(center) ~= size(obj.Velocity_EllipsePts.center)
                center = center';
            end
            
            if ~isempty(obj.Velocity_EllipsePts.center)
                if obj.Velocity_EllipsePts.a ~= a || obj.Velocity_EllipsePts.b ~= b || obj.Velocity_EllipsePts.c ~= c || all(obj.Velocity_EllipsePts.center ~= center) || obj.Velocity_EllipsePts.ratio ~= ratio%#ok<PROPLC>

                    [X, Y, Z] = ellipsoid(0,0,0,a,b,c,N);%#ok<PROPLC>

                    XX = zeros(N+1,N+1);
                    YY = XX;
                    ZZ = XX;

                    for k = 1:length(X)
                        for j = 1:length(X)
                            point = [X(k,j) Y(k,j) Z(k,j)]';

                            if rotate == true
                                if size(obj.eigVecs,3) ~= 1
                                    P = obj.eigVecs(:,:,idx) * point;
                                else
                                    P = obj.eigVecs * point;
                                end
                            else
                                P = point;
                            end

                            XX(k,j) = P(1)*ratio+center(1);
                            YY(k,j) = P(2)*ratio+center(2);
                            ZZ(k,j) = P(3)*ratio+center(3);
                        end
                    end

                    obj.Velocity_EllipsePts = struct('X', XX, 'Y', YY, 'Z', ZZ, 'a', a, 'b', b, 'c', c, 'center', center, 'ratio', ratio);%#ok<PROPLC>
                    mesh(XX,YY,ZZ, 'FaceColor', faceColor, 'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha);
                else
                    mesh(obj.Velocity_EllipsePts.X,obj.Velocity_EllipsePts.Y,obj.Velocity_EllipsePts.Z, 'FaceColor', faceColor, 'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha);
                end
            else
                [X, Y, Z] = ellipsoid(0,0,0,a,b,c,N);%#ok<PROPLC>

                XX = zeros(N+1,N+1);
                YY = XX;
                ZZ = XX;

                for k = 1:length(X)
                    for j = 1:length(X)
                        point = [X(k,j) Y(k,j) Z(k,j)]';

                        if rotate == true
                            if size(obj.eigVecs,3) ~= 1
                                P = obj.eigVecs(:,:,idx) * point;
                            else
                                P = obj.eigVecs * point;
                            end
                        else
                            P = point;
                        end

                        XX(k,j) = P(1)*ratio+center(1);
                        YY(k,j) = P(2)*ratio+center(2);
                        ZZ(k,j) = P(3)*ratio+center(3);
                    end
                end

                obj.Velocity_EllipsePts = struct('X', XX, 'Y', YY, 'Z', ZZ, 'a', a, 'b', b, 'c', c, 'center', center, 'ratio', ratio);%#ok<PROPLC>
                mesh(XX,YY,ZZ, 'FaceColor', faceColor, 'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha);
            end
            
            
           
            hidden off;
            
        end
        
        function [eigVectors, eigValues, rotation] = getVelocityEigens(obj, varargin)
            
            if any(strcmp(varargin, 'Jacobian'))
                idx = find(strcmp(varargin, 'Jacobian'));
                obj.J = varargin{idx+1};
            end 
            
%             if size(obj.J,3) == 1
                M = obj.J * obj.J';
%             else
%                 jinv = permute(obj.J, [2,1,3]);
%                 M = mtimesx(obj.J, jinv);
%             end
            
            [eigVecs, eigValues, rotation] = svd(M);    %#ok<PROPLC> 
            
            obj.eigVecs = eigVecs; %#ok<PROPLC>
            obj.eigValues = eigValues;
            
            eigVectors = obj.eigVecs;
        end
        
        function [ws_cube, mapData] = mapWsUntil_wci(obj, resolution, wci_des, varargin)
    
            mappedThetas = zeros(100, 3);

            ws_cube_init = obj.getCubeAboutRobot(resolution);

        %     ws_cube = gtpr.createCubeIdxMap('ws_cube', ws_cube);

            thetas = [];

            wci = 0;
            ctr = 1;
            tic;
            while (wci < wci_des) || isinf(wci) || isnan(wci)

                theta = rand(1,3)*pi; %#ok<PROPLC>

                if ismember(theta, mappedThetas) %#ok<PROPLC>
                    continue;
                end

                mappedThetas(ctr,:) =  theta; %#ok<PROPLC>

                if ctr == 1
                    ws_cube = ws_cube_init;
                end

                theta = rand(100,3)*pi; %#ok<PROPLC>
                thetas = [thetas; theta]; %#ok<PROPLC>
                for i = 1:size(theta,1) %#ok<PROPLC>

                    try
                        W(ctr,:) = obj.directGeometry(theta(i,:)); %#ok<PROPLC>

                    catch ME

                        if strcmpi(ME.message, "Impossible Configuration present!")
                            continue;
                        else
                            disp("a");
                        end
                    end


                    mapData.Jacobian{ctr} = obj.Jacobian();

                    [mapData.mu1{ctr}, mapData.mu2{ctr}, mapData.mu3{ctr}] = obj.getManipulabilityMeasures();

                    [mapData.eigVectors{ctr}, mapData.eigValues{ctr}] = obj.getVelocityEigens();

                    ws_cube = obj.mapBitPointInSpace(ws_cube, W(ctr,:), obj.bit(1).Position, obj.bit(1).Name, 1);  %#ok<PROPLC>
                    ws_cube.W.points = W; %#ok<PROPLC>
                    ctr = ctr + 1;

                end

                ws_cube = obj.getCubeAboutRobot(resolution, 'thetas', thetas);
                ws_cube.W.points = W; %#ok<PROPLC>

                for i = 1:size(W,1) %#ok<PROPLC>
                    ws_cube = obj.mapBitPointInSpace(ws_cube, W(i,:), obj.bit(1).Position, obj.bit(1).Name, 1); %#ok<PROPLC>
                end

                wci = obj.workspaceCoverageIndex(ws_cube);%, "plotBoundary");

                if isinf(wci) || isnan(wci)
                    break;
                end

            end
            findingTime = toc;

            mapData.wci = wci;
            mapData.thetas = nonzeros(mappedThetas);
            mapData.runTime = findingTime;
            mapData.avgTime = findingTime/ctr;
            mapData.dataPtsMapped = ctr;

            ws_cube.W.theta = mapData.thetas;

        end
        
        function [mu1, mu2, mu3, eigenVectors, eigenValues, rotation] = getManipulabilityMeasures(obj, varargin)
            % Returns the 3 basic Manipulability Measures:
            % - mu1 => ratio of longest and shortest axes
            % - mu2 => condition number of M (J*J')
            % - mu3 => proportional volume of ellipsoid
           
            
            [eigenVectors, eigenValues, rotation] = obj.getVelocityEigens();
           
            mu1 = [];
            mu2 = [];
            mu3 = [];
            for i = 1:size(obj.eigValues,3)
                eigValues = nonzeros(obj.eigValues(:,:,i)); %#ok<PROPLC>
                mu1(i) = sqrt(max(eigValues))/sqrt(min(eigValues)); %#ok<PROPLC>
                mu2(i) = max(eigValues)/min(eigValues); %#ok<PROPLC>
                mu3(i) = sqrt(det(obj.J*obj.J')); %#ok<PROPLC>
            end

            eigenVectors = obj.eigVecs;
            eigenValues = obj.eigValues;
        end


        
        function Pts = realGeometry(obj)
        
            Pts.B1L = obj.B1 + norm(obj.A1)/10 * obj.rotz(pi/2)*[sin(0); cos(0); 0];
            Pts.B1R = obj.B1 - norm(obj.A1)/10 * obj.rotz(pi/2)*[sin(0); cos(0); 0];
            Pts.B2L = obj.B2 + norm(obj.A2)/10 * obj.rotz(pi/2)*[sin(obj.alpha12); cos(obj.alpha12); 0];
            Pts.B2R = obj.B2 - norm(obj.A2)/10 * obj.rotz(pi/2)*[sin(obj.alpha12); cos(obj.alpha12); 0];
            Pts.B3L = obj.B3 + norm(obj.A3)/10 * obj.rotz(pi/2)*[sin(obj.alpha13); cos(obj.alpha13); 0];
            Pts.B3R = obj.B3 - norm(obj.A3)/10 * obj.rotz(pi/2)*[sin(obj.alpha13); cos(obj.alpha13); 0];

            Pts.C1L = obj.C1 + norm(obj.A1)/10 * obj.rotz(pi/2)*[sin(0); cos(0); 0];
            Pts.C1R = obj.C1 - norm(obj.A1)/10 * obj.rotz(pi/2)*[sin(0); cos(0); 0];
            Pts.C2L = obj.C2 + norm(obj.A2)/10 * obj.rotz(pi/2)*[sin(obj.alpha12); cos(obj.alpha12); 0];
            Pts.C2R = obj.C2 - norm(obj.A2)/10 * obj.rotz(pi/2)*[sin(obj.alpha12); cos(obj.alpha12); 0];
            Pts.C3L = obj.C3 + norm(obj.A3)/10 * obj.rotz(pi/2)*[sin(obj.alpha13); cos(obj.alpha13); 0];
            Pts.C3R = obj.C3 - norm(obj.A3)/10 * obj.rotz(pi/2)*[sin(obj.alpha13); cos(obj.alpha13); 0];
        
            obj.B1L = Pts.B1L;
            obj.B2L = Pts.B2L;
            obj.B3L = Pts.B3L;

            obj.C1L = Pts.C1L;
            obj.C2L = Pts.C2L;
            obj.C3L = Pts.C3L;

            obj.B1R = Pts.B1R;
            obj.B2R = Pts.B2R;
            obj.B3R = Pts.B3R;

            obj.C1R = Pts.C1R;
            obj.C2R = Pts.C2R;
            obj.C3R = Pts.C3R;
            
        end
        
        function plotHandles = plotMechanism(obj, varargin)
            
            
            hold on;
            O = plot3(0,0,0,'ko');
            
            A1A2 = plot3([obj.A1(1) obj.A2(1)], [obj.A1(2)  obj.A2(2)], [obj.A1(3)  obj.A2(3)], 'k');
			A1A3 = plot3([obj.A1(1) obj.A3(1)], [obj.A1(2)  obj.A3(2)], [obj.A1(3)  obj.A3(3)], 'k');
			A2A3 = plot3([obj.A2(1) obj.A3(1)], [obj.A2(2)  obj.A3(2)], [obj.A2(3)  obj.A3(3)], 'k');
            
			C1C2 = plot3([obj.C1(1) obj.C2(1)], [obj.C1(2)  obj.C2(2)], [obj.C1(3)  obj.C2(3)], 'k');
			C2C3 = plot3([obj.C2(1) obj.C3(1)], [obj.C2(2)  obj.C3(2)], [obj.C2(3)  obj.C3(3)], 'k');
			C1C3 = plot3([obj.C1(1) obj.C3(1)], [obj.C1(2)  obj.C3(2)], [obj.C1(3)  obj.C3(3)], 'k');
            
            plotHandles.O = O;
            plotHandles.A1A2 = A1A2;
            plotHandles.A1A3 = A1A3;
            plotHandles.A2A3 = A2A3;
            
            plotHandles.C1C2 = C1C2;
            plotHandles.C1C3 = C1C3;
            plotHandles.C2C3 = C2C3;
            
            
            if any(strcmp('real', varargin))
                
                obj.realGeometry();
                
                plotHandles.B1LB1R = plot3([obj.B1L(1) obj.B1R(1)], [obj.B1L(2) obj.B1R(2)], [obj.B1L(3) obj.B1R(3)], 'k');
                plotHandles.B2LB2R = plot3([obj.B2L(1) obj.B2R(1)], [obj.B2L(2) obj.B2R(2)], [obj.B2L(3) obj.B2R(3)], 'k');
                plotHandles.B3LB3R = plot3([obj.B3L(1) obj.B3R(1)], [obj.B3L(2) obj.B3R(2)], [obj.B3L(3) obj.B3R(3)], 'k');

                plotHandles.C1LC1R = plot3([obj.C1L(1) obj.C1R(1)], [obj.C1L(2) obj.C1R(2)], [obj.C1L(3) obj.C1R(3)], 'k');
                plotHandles.C2LC2R = plot3([obj.C2L(1) obj.C2R(1)], [obj.C2L(2) obj.C2R(2)], [obj.C2L(3) obj.C2R(3)], 'k');
                plotHandles.C3LC3R = plot3([obj.C3L(1) obj.C3R(1)], [obj.C3L(2) obj.C3R(2)], [obj.C3L(3) obj.C3R(3)], 'k');

                plotHandles.B1LC1L = plot3([obj.B1L(1) obj.C1L(1)], [obj.B1L(2) obj.C1L(2)], [obj.B1L(3) obj.C1L(3)], 'k');
                plotHandles.B2LC2R = plot3([obj.B2L(1) obj.C2L(1)], [obj.B2L(2) obj.C2L(2)], [obj.B2L(3) obj.C2L(3)], 'k');
                plotHandles.B3LC3R = plot3([obj.B3L(1) obj.C3L(1)], [obj.B3L(2) obj.C3L(2)], [obj.B3L(3) obj.C3L(3)], 'k');

                plotHandles.B1RC1R = plot3([obj.B1R(1) obj.C1R(1)], [obj.B1R(2) obj.C1R(2)], [obj.B1R(3) obj.C1R(3)], 'k');
                plotHandles.B2RC2R = plot3([obj.B2R(1) obj.C2R(1)], [obj.B2R(2) obj.C2R(2)], [obj.B2R(3) obj.C2R(3)], 'k');
                plotHandles.B3RC3R = plot3([obj.B3R(1) obj.C3R(1)], [obj.B3R(2) obj.C3R(2)], [obj.B3R(3) obj.C3R(3)], 'k');
                
%                 C = [178, 178, 178]/255;
%                     
%                 plotHandles.BaseTriangle = fill3([obj.A1(1), obj.A2(1), obj.A3(1)], [obj.A1(2), obj.A2(2), obj.A3(2)], [obj.A1(3), obj.A2(3), obj.A3(3)], C, 'o');
%        
%                 C = [217, 217, 217]/255;
%                 plotHandles.workTriangle = fill3([obj.C1(1), obj.C2(1), obj.C3(1)], [obj.C1(2), obj.C2(2), obj.C3(2)], [obj.C1(3), obj.C2(3), obj.C3(3)], C);
                    
                if any(strcmpi('LineWidth', varargin))
                    idx = find(strcmpi('LineWidth', varargin));
                    pt = varargin{idx+1};

                    plotHandles.A1A2.LineWidth = pt;
                    plotHandles.A1A3.LineWidth = pt;
                    plotHandles.A2A3.LineWidth = pt;	
                    plotHandles.C1C2.LineWidth = pt;
                    plotHandles.C2C3.LineWidth = pt;
                    plotHandles.C1C3.LineWidth = pt;

                    plotHandles.B1LB1R.LineWidth = pt;
                    plotHandles.B2LB2R.LineWidth = pt;
                    plotHandles.B3LB3R.LineWidth = pt;
                    
                    plotHandles.C1LC1R.LineWidth = pt;
                    plotHandles.C2LC2R.LineWidth = pt;
                    plotHandles.C3LC3R.LineWidth = pt;
                    
                    plotHandles.B1LC1L.LineWidth = pt;
                    plotHandles.B2LC2R.LineWidth = pt;
                    plotHandles.B3LC3R.LineWidth = pt;
                    
                    plotHandles.B1RC1R.LineWidth = pt;
                    plotHandles.B2RC2R.LineWidth = pt;
                    plotHandles.B3RC3R.LineWidth = pt;
                    
                    plotHandles.BaseTriangle.LineWidth = pt;
                    plotHandles.workTriangle.LineWidth = pt;

                end

            else
                
                B1C1 = plot3([obj.B1(1) obj.C1(1)], [obj.B1(2) obj.C1(2)], [obj.B1(3) obj.C1(3)], 'k');
                B2C2 = plot3([obj.B2(1) obj.C2(1)], [obj.B2(2) obj.C2(2)], [obj.B2(3) obj.C2(3)], 'k');
                B3C3 = plot3([obj.B3(1) obj.C3(1)], [obj.B3(2) obj.C3(2)], [obj.B3(3) obj.C3(3)], 'k');
                
                plotHandles.B1C1 = B1C1;
                plotHandles.B2C2 = B2C2;
                plotHandles.B3C3 = B3C3;
                
            end
            
            if any(strcmp('annotate', varargin))
               A1_text = text(obj.A1(1), obj.A1(2), obj.A1(3), '$\bf{A_1}$', 'Interpreter', 'latex', 'Color', 'r', 'FontSize', 20);
               A2_text = text(obj.A2(1), obj.A2(2), obj.A2(3), '$\bf{A_2}$', 'Interpreter', 'latex', 'Color', 'r', 'FontSize', 20);
               A3_text = text(obj.A3(1), obj.A3(2), obj.A3(3), '$\bf{A_3}$', 'Interpreter', 'latex', 'Color', 'r', 'FontSize', 20);
               
               B1_text = text(obj.B1(1), obj.B1(2), obj.B1(3), '$\bf{B_1}$', 'Interpreter', 'latex', 'Color', 'g', 'FontSize', 20);
               B2_text = text(obj.B2(1), obj.B2(2), obj.B2(3), '$\bf{B_2}$', 'Interpreter', 'latex', 'Color', 'g', 'FontSize', 20);
               B3_text = text(obj.B3(1), obj.B3(2), obj.B3(3), '$\bf{B_3}$', 'Interpreter', 'latex', 'Color', 'g', 'FontSize', 20);
               
               C1_text = text(obj.C1(1), obj.C1(2), obj.C1(3), '$\bf{C_1}$', 'Interpreter', 'latex', 'Color', 'b', 'FontSize', 20);
               C2_text = text(obj.C2(1), obj.C2(2), obj.C2(3), '$\bf{C_2}$', 'Interpreter', 'latex', 'Color', 'b', 'FontSize', 20);
               C3_text = text(obj.C3(1), obj.C3(2), obj.C3(3), '$\bf{C_3}$', 'Interpreter', 'latex', 'Color', 'b', 'FontSize', 20);
               
               plotHandles.A1_text = A1_text;
               plotHandles.A2_text = A2_text;
               plotHandles.A3_text = A3_text;
               
               plotHandles.B1_text = B1_text;
               plotHandles.B2_text = B2_text;
               plotHandles.B3_text = B3_text;
               
               plotHandles.C1_text = C1_text;
               plotHandles.C2_text = C2_text;
               plotHandles.C3_text = C3_text;
               
            end
            
            A1B1 = plot3([obj.A1(1) obj.B1(1)], [obj.A1(2) obj.B1(2)], [obj.A1(3) obj.B1(3)], 'k');
            A2B2 = plot3([obj.A2(1) obj.B2(1)], [obj.A2(2) obj.B2(2)], [obj.A2(3) obj.B2(3)], 'k');
            A3B3 = plot3([obj.A3(1) obj.B3(1)], [obj.A3(2) obj.B3(2)], [obj.A3(3) obj.B3(3)], 'k');
            
            plotHandles.A1B1 = A1B1;
            plotHandles.A2B2 = A2B2;
            plotHandles.A3B3 = A3B3;
            
            if any(strcmpi('LineWidth', varargin))
                idx = find(strcmpi('LineWidth', varargin));
                pt = varargin{idx+1};
                plotHandles.A1B1.LineWidth = pt;
                plotHandles.A2B2.LineWidth = pt;
                plotHandles.A3B3.LineWidth = pt;
            end

            plotHandles.W = plot3(obj.W(1), obj.W(2), obj.W(3), 'ms'); 
            
            hold off;
        end
        
        function CoM = CoMs(obj, varargin)
            % Calculate Center of masses of linkages, and W
            
            CoM_L1 = ((obj.B1-obj.A1)/2)+obj.A1;
            CoM_L2 = ((obj.B2-obj.A2)/2)+obj.A2;
            CoM_L3 = ((obj.B3-obj.A3)/2)+obj.A3;
            
            if obj.e3 == obj.e2 && obj.e3 == obj.e1
                CoM_W = obj.W;
            else
                CoM_W = (1/3)*(obj.A3 + obj.B3 + obj.C3);
            end
            
            if strcmp('real', varargin)
                obj.B1L = obj.B1 + norm(obj.A1)/10 * obj.rotz(pi/2)*[sin(0); cos(0); 0];
                obj.B1R = obj.B1 - norm(obj.A1)/10 * obj.rotz(pi/2)*[sin(0); cos(0); 0];
                obj.B2L = obj.B2 + norm(obj.A2)/10 * obj.rotz(pi/2)*[sin(obj.alpha12); cos(obj.alpha12); 0];
                obj.B2R = obj.B2 - norm(obj.A2)/10 * obj.rotz(pi/2)*[sin(obj.alpha12); cos(obj.alpha12); 0];
                obj.B3L = obj.B3 + norm(obj.A3)/10 * obj.rotz(pi/2)*[sin(obj.alpha13); cos(obj.alpha13); 0];
                obj.B3R = obj.B3 - norm(obj.A3)/10 * obj.rotz(pi/2)*[sin(obj.alpha13); cos(obj.alpha13); 0];

                obj.C1L = obj.C1 + norm(obj.A1)/10 * obj.rotz(pi/2)*[sin(0); cos(0); 0];
                obj.C1R = obj.C1 - norm(obj.A1)/10 * obj.rotz(pi/2)*[sin(0); cos(0); 0];
                obj.C2L = obj.C2 + norm(obj.A2)/10 * obj.rotz(pi/2)*[sin(obj.alpha12); cos(obj.alpha12); 0];
                obj.C2R = obj.C2 - norm(obj.A2)/10 * obj.rotz(pi/2)*[sin(obj.alpha12); cos(obj.alpha12); 0];
                obj.C3L = obj.C3 + norm(obj.A3)/10 * obj.rotz(pi/2)*[sin(obj.alpha13); cos(obj.alpha13); 0];
                obj.C3R = obj.C3 - norm(obj.A3)/10 * obj.rotz(pi/2)*[sin(obj.alpha13); cos(obj.alpha13); 0];

                % CoMs of parallel mechanism:

                CoM_B1C1L = ((obj.C1L-obj.B1L)/2)+obj.B1L;
                CoM_B2C2L = ((obj.C2L-obj.B2L)/2)+obj.B2L;
                CoM_B3C3L = ((obj.C3L-obj.B3L)/2)+obj.B3L;
                
                CoM_B1C1R = ((obj.C1R-obj.B1R)/2)+obj.B1R;
                CoM_B2C2R = ((obj.C2R-obj.B2R)/2)+obj.B2R;
                CoM_B3C3R = ((obj.C3R-obj.B3R)/2)+obj.B3R;
                
                CoM = struct("W", CoM_W, ...
                              "B1", obj.B1, ...
                              "B2", obj.B2, ...
                              "B3", obj.B3, ...
                              "C1", obj.C1, ...
                              "C2", obj.C2, ...
                              "C3", obj.C3, ...
                              "B1C1L", CoM_B1C1L, ...
                              "B2C2L", CoM_B2C2L, ...
                              "B3C3L", CoM_B3C3L, ...
                              "B1C1R", CoM_B1C1R, ...
                              "B2C2R", CoM_B2C2R, ...
                              "B3C3R", CoM_B3C3R, ...
                              "parallelPoints", struct("B1L", obj.B1L, "B1R", obj.B1R, ...
                                                       "B2L", obj.B2L, "B2R", obj.B2R, ...
                                                       "B3L", obj.B3L, "B3R", obj.B3R, ...
                                                       "C1L", obj.C1L, "C1R", obj.C1R, ...
                                                       "C2L", obj.C2L, "C2R", obj.C2R, ...
                                                       "C3L", obj.C3L, "C3R", obj.C3R));
            else
                CoM_l1 = ((obj.C1-obj.B1)/2)+obj.B1;
                CoM_l2 = ((obj.C2-obj.B2)/2)+obj.B2;
                CoM_l3 = ((obj.C3-obj.B3)/2)+obj.B3;
                
                CoM = struct("W", CoM_W, ...
                              "L1", CoM_L1, ...
                              "L2", CoM_L2, ...
                              "L3", CoM_L3, ...
                              "l1", CoM_l1, ...
                              "l2", CoM_l2, ...
                              "l3", CoM_l3);
            end
            
            obj.CoM = CoM;
            
        end
        
        function [J_arm1, J_arm2, J_arm3, J_C, phi1, phi2, phi3] = Jacobian_arms_acceleration(obj,V)
            
            % Withing the phi-s lay the angular velocities of each
            % respective arms, for the given screws:
            % phi = [w01, w12, w23, w34, w45, w56]';
            
            if size(V,1) ~= 6
                V = [0;0;0;V];
            end
            
            % Jacobian by arms:            
            J_arm1 = [ones(6,1), obj.S12_1, obj.S23_1, obj.S34_1, obj.S45_1, obj.S56_1];
            J_arm2 = [ones(6,1), obj.S12_2, obj.S23_2, obj.S34_2, obj.S45_2, obj.S56_2];
            J_arm3 = [ones(6,1), obj.S12_3, obj.S23_3, obj.S34_3, obj.S45_3, obj.S56_3];

            phi1 = J_arm1\V;
            phi2 = J_arm2\V;
            phi3 = J_arm3\V;

            phi1(phi1 <= obj.zeroTolerance) = 0;
            phi2(phi2 <= obj.zeroTolerance) = 0;
            phi3(phi3 <= obj.zeroTolerance) = 0;
            
            Lie1 = obj.lieProduct(phi1(1)*obj.S01_1, phi1(2)*obj.S12_1 + phi1(3)*obj.S23_1 + phi1(4)*obj.S34_1 + phi1(5)*obj.S45_1 + phi1(6)*obj.S56_1) + ... 
                   obj.lieProduct(phi1(2)*obj.S12_1, phi1(3)*obj.S23_1 + phi1(4)*obj.S34_1 + phi1(5)*obj.S45_1 + phi1(6)*obj.S56_1) + ...
                   obj.lieProduct(phi1(3)*obj.S23_1, phi1(4)*obj.S34_1 + phi1(5)*obj.S45_1 + phi1(6)*obj.S56_1) + ...
                   obj.lieProduct(phi1(4)*obj.S34_1, phi1(5)*obj.S45_1 + phi1(6)*obj.S56_1) + ...
                   obj.lieProduct(phi1(5)*obj.S45_1, phi1(6)*obj.S56_1);

            Lie2 = obj.lieProduct(phi2(1)*obj.S01_2, phi2(2)*obj.S12_2 + phi2(3)*obj.S23_2 + phi2(4)*obj.S34_2 + phi2(5)*obj.S45_2 + phi2(6)*obj.S56_2) + ... 
                   obj.lieProduct(phi2(2)*obj.S12_2, phi2(3)*obj.S23_2 + phi2(4)*obj.S34_2 + phi2(5)*obj.S45_2 + phi2(6)*obj.S56_2) + ...
                   obj.lieProduct(phi2(3)*obj.S23_2, phi2(4)*obj.S34_2 + phi2(5)*obj.S45_2 + phi2(6)*obj.S56_2) + ...
                   obj.lieProduct(phi2(4)*obj.S34_2, phi2(5)*obj.S45_2 + phi2(6)*obj.S56_2) + ...
                   obj.lieProduct(phi2(5)*obj.S45_2, phi2(6)*obj.S56_2);

            Lie3 = obj.lieProduct(phi3(1)*obj.S01_3, phi3(2)*obj.S12_3 + phi3(3)*obj.S23_3 + phi3(4)*obj.S34_3 + phi3(5)*obj.S45_3 + phi3(6)*obj.S56_3) + ... 
                   obj.lieProduct(phi3(2)*obj.S12_3, phi3(3)*obj.S23_3 + phi3(4)*obj.S34_3 + phi3(5)*obj.S45_3 + phi3(6)*obj.S56_3) + ...
                   obj.lieProduct(phi3(3)*obj.S23_3, phi3(4)*obj.S34_3 + phi3(5)*obj.S45_3 + phi3(6)*obj.S56_3) + ...
                   obj.lieProduct(phi3(4)*obj.S34_3, phi3(5)*obj.S45_3 + phi3(6)*obj.S56_3) + ...
                   obj.lieProduct(phi3(5)*obj.S45_3, phi3(6)*obj.S56_3);

            obj.J_C = [obj.kleinForm(Lie1, obj.Sl1); obj.kleinForm(Lie2, obj.Sl2); obj.kleinForm(Lie3, obj.Sl3)];
            
            J_C = obj.J_C;
            
            obj.J_arm1 = J_arm1;
            obj.J_arm2 = J_arm2;
            obj.J_arm3 = J_arm3;
            obj.phi1 = phi1;
            obj.phi2 = phi2;
            obj.phi3 = phi3;
        end
        
        function a = directKinematics_Acceleration(obj, ddtheta, V, varargin)
        
            obj.Jacobian();
            obj.armAngles();
            
            % Lie screw of the system, according to: chapter 14.4.2 of
            % Kinematic Analysis of Parallel Manipulators by Algebraic
            % Screw Theory, by: Jaime Gallardo-Alvarado
            
            % Jacobian by arms:
            obj.Jacobian_arms_acceleration(V);
            
            a = obj.J_A\(obj.J_B*ddtheta + obj.J_C);
            
            a(a <= obj.zeroTolerance) = 0;
        end
        
        function ddtheta = inverseKinematics_Acceleration(obj, a, V)
            obj.Jacobian();
            obj.armAngles();
            
            % Lie screw of the system, according to: chapter 14.4.2 of
            % Kinematic Analysis of Parallel Manipulators by Algebraic
            % Screw Theory, by: Jaime Gallardo-Alvarado
            
            % Jacobian by arms:
            obj.Jacobian_arms_acceleration(V);
                ddtheta = obj.J_B\(obj.J_A*a - obj.J_C);
        end
        
        function V = directKinematics(obj, dtheta, varargin)
			obj.Jacobian();
            V = obj.J * dtheta;
            
            V(abs(V)<=obj.zeroTolerance) = 0;
        end 
        
        function dtheta = inverseKinematics(obj, V, varargin)
            obj.v = V;
            obj.Jacobian();
            obj.dtheta = obj.J\V;
            dtheta = obj.dtheta;
        end
        
        function [IL_1, IL_2, IL_3] = upperArmsInertiaMatrix(obj,mL1, mL2, mL3, IL_1, IL_2, IL_3, varargin)
            % Calculates the inertia matrix of each arm mainly based off
            % of: https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec26.pdf
            % https://www.boschrexroth.com/en/xc/products/product-groups/assembly-technology/basic-mechanic-elements/strut-profiles/strut-profiles-slot-6-modular-dimensions-20/20x20
            
            % I is given as on the Bosch Rexroth site, so Izz = 0, that is
            % the main axis at begginning
            
%             if any(strcmp(varargin, 'I')) % If inertia matrix is given at input
%                idx = find(strcmp(varargin, 'I'));
%                obj.I = varargin{idx+1};
%            end

            % THE INERTIA MATRIX IS GIVEN AT THE RODS ORIGIN, SO WE CAN USE
            % S12 (qi) OMEGAS DIRECTLY
            rA = [0;0;0];
            rB = [0;0;0];
            rC = [0;0;0];

            % Arm A
            TA = obj.rotx(90 + obj.theta(1));
%             rA = obj.A1 + (obj.A2-obj.A1); 
            IL_1 = (TA*IL_1*TA') + mL1 * [rA(2)^2 + rA(3)^2, rA(1)*rA(2), rA(1)*rA(3);
                                      rA(2)*rA(1), rA(1)^2 + rA(3)^2, rA(2)*rA(3);
                                      rA(3)*rA(1), rA(3)*rA(2), rA(1)^2 + rA(2)^2];
                                  
            % Arm B
            TB = obj.rotz(-obj.alpha12)*obj.rotx(90+obj.theta(2));
%             rB = obj.B1 + (obj.B2-obj.B1);
            IL_2 = (TB*IL_2*TB') + mL2 * [rB(2)^2 + rB(3)^2, rB(1)*rB(2), rB(1)*rB(3);
                                      rB(2)*rB(1), rB(1)^2 + rB(3)^2, rB(2)*rB(3);
                                      rB(3)*rB(1), rB(3)*rB(2), rB(1)^2 + rB(2)^2];
                                  
            % Arm C
            TC = obj.rotz(-obj.alpha13)*obj.rotx(90+obj.theta(3));
%             rC = obj.C1 + (obj.C2-obj.C1);
            IL_3 = (TC*IL_3*TC') + mL3 * [rC(2)^2 + rC(3)^2, rC(1)*rC(2), rC(1)*rC(3);
                                      rC(2)*rC(1), rC(1)^2 + rC(3)^2, rC(2)*rC(3);
                                      rC(3)*rC(1), rC(3)*rC(2), rC(1)^2 + rC(2)^2];
        end
        
        function [Il_1, Il_2, Il_3] = lowerArmsInertiaMatrix(obj,ml1, ml2, ml3, Il_1, Il_2, Il_3, varargin)
            % Calculates the inertia matrix of each arm mainly based off
            % of: https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec26.pdf
            % https://www.boschrexroth.com/en/xc/products/product-groups/assembly-technology/basic-mechanic-elements/strut-profiles/strut-profiles-slot-6-modular-dimensions-20/20x20
            
            % I is given as on the Bosch Rexroth site, so Izz = 0, that is
            % the main axis at begginning
            
            % THE INERTIA MATRIX IS GIVEN AT THE RODS ORIGIN, SO WE CAN USE
            % S23 S34 OMEGAS DIRECTLY
            rA = [0;0;0];
            rB = [0;0;0];
            rC = [0;0;0];
            
%             if any(strcmp(varargin, 'I')) % If inertia matrix is given at input
%                idx = find(strcmp(varargin, 'I'));
%                obj.I = varargin{idx+1};
%            end
            % Arm A
            TA = obj.Rl1; % Rotational matrix, previously calculated in: obj.arms_rotationalMatrices()
%             rA = obj.A2 + ((obj.A3-obj.A2) / 2); 
            Il_1 = (TA*Il_1*TA') + ml1 * [rA(2)^2 + rA(3)^2, rA(1)*rA(2), rA(1)*rA(3);
                                      rA(2)*rA(1), rA(1)^2 + rA(3)^2, rA(2)*rA(3);
                                      rA(3)*rA(1), rA(3)*rA(2), rA(1)^2 + rA(2)^2];
                                  
            % Arm B
            TB = obj.Rl2;
%             rB = obj.B2 + ((obj.B3-obj.B2)/2);
            Il_2 = (TB*Il_2*TB') + ml2 * [rB(2)^2 + rB(3)^2, rB(1)*rB(2), rB(1)*rB(3);
                                      rB(2)*rB(1), rB(1)^2 + rB(3)^2, rB(2)*rB(3);
                                      rB(3)*rB(1), rB(3)*rB(2), rB(1)^2 + rB(2)^2];
                                  
            % Arm C
            TC = obj.Rl3;
%             rC = obj.C2 + (obj.C3-obj.C2);
            Il_3 = (TC*Il_3*TC') + ml3 * [rC(2)^2 + rC(3)^2, rC(1)*rC(2), rC(1)*rC(3);
                                      rC(2)*rC(1), rC(1)^2 + rC(3)^2, rC(2)*rC(3);
                                      rC(3)*rC(1), rC(3)*rC(2), rC(1)^2 + rC(2)^2];
        end
        
        function out = lieProduct(~,S1,S2)
           out = zeros(6,1);
           out(1:3) = cross(S1(1:3), S2(1:3));
           out(4:6) = cross(S1(1:3),S2(4:6)) - cross(S2(1:3), S1(4:6));
        end
        
        function plot_rotFrames(obj)
            
            quiver3(0,0,0, 1,0,0, 'r', 'AutoScale', 'off');
            quiver3(0,0,0, 0,1,0, 'g', 'AutoScale', 'off');
            quiver3(0,0,0, 0,0,1, 'b', 'AutoScale', 'off');
            
            quiver3(obj.A1(1), obj.A1(2), obj.A1(3), obj.RL1(1,1), obj.RL1(2,1), obj.RL1(3,1), 'r', 'AutoScale', 'off');
            quiver3(obj.A1(1), obj.A1(2), obj.A1(3), obj.RL1(1,2), obj.RL1(2,2), obj.RL1(3,2), 'g', 'AutoScale', 'off');
            quiver3(obj.A1(1), obj.A1(2), obj.A1(3), obj.RL1(1,3), obj.RL1(2,3), obj.RL1(3,3), 'b', 'AutoScale', 'off');
            
            quiver3(obj.B1(1), obj.B1(2), obj.B1(3), obj.RL2(1,1), obj.RL2(2,1), obj.RL2(3,1), 'r', 'AutoScale', 'off');
            quiver3(obj.B1(1), obj.B1(2), obj.B1(3), obj.RL2(1,2), obj.RL2(2,2), obj.RL2(3,2), 'g', 'AutoScale', 'off');
            quiver3(obj.B1(1), obj.B1(2), obj.B1(3), obj.RL2(1,3), obj.RL2(2,3), obj.RL2(3,3), 'b', 'AutoScale', 'off');
            
            quiver3(obj.C1(1), obj.C1(2), obj.C1(3), obj.RL3(1,1), obj.RL3(2,1), obj.RL3(3,1), 'r', 'AutoScale', 'off');
            quiver3(obj.C1(1), obj.C1(2), obj.C1(3), obj.RL3(1,2), obj.RL3(2,2), obj.RL3(3,2), 'g', 'AutoScale', 'off');
            quiver3(obj.C1(1), obj.C1(2), obj.C1(3), obj.RL3(1,3), obj.RL3(2,3), obj.RL3(3,3), 'b', 'AutoScale', 'off');
            
            quiver3(obj.A2(1), obj.A2(2), obj.A2(3), obj.Rl1(1,1), obj.Rl1(2,1), obj.Rl1(3,1), 'r', 'AutoScale', 'off');
            quiver3(obj.A2(1), obj.A2(2), obj.A2(3), obj.Rl1(1,2), obj.Rl1(2,2), obj.Rl1(3,2), 'g', 'AutoScale', 'off');
            quiver3(obj.A2(1), obj.A2(2), obj.A2(3), obj.Rl1(1,3), obj.Rl1(2,3), obj.Rl1(3,3), 'b', 'AutoScale', 'off');
            
            quiver3(obj.B2(1), obj.B2(2), obj.B2(3), obj.Rl2(1,1), obj.Rl2(2,1), obj.Rl2(3,1), 'r', 'AutoScale', 'off');
            quiver3(obj.B2(1), obj.B2(2), obj.B2(3), obj.Rl2(1,2), obj.Rl2(2,2), obj.Rl2(3,2), 'g', 'AutoScale', 'off');
            quiver3(obj.B2(1), obj.B2(2), obj.B2(3), obj.Rl2(1,3), obj.Rl2(2,3), obj.Rl2(3,3), 'b', 'AutoScale', 'off');
            
            quiver3(obj.C2(1), obj.C2(2), obj.C2(3), obj.Rl3(1,1), obj.Rl3(2,1), obj.Rl3(3,1), 'r', 'AutoScale', 'off');
            quiver3(obj.C2(1), obj.C2(2), obj.C2(3), obj.Rl3(1,2), obj.Rl3(2,2), obj.Rl3(3,2), 'g', 'AutoScale', 'off');
            quiver3(obj.C2(1), obj.C2(2), obj.C2(3), obj.Rl3(1,3), obj.Rl3(2,3), obj.Rl3(3,3), 'b', 'AutoScale', 'off');
            
            quiver3(obj.W(1), obj.W(2), obj.W(3), 1, 0, 0, 'r', 'AutoScale', 'off');
            quiver3(obj.W(1), obj.W(2), obj.W(3), 0, 1, 0, 'g', 'AutoScale', 'off');
            quiver3(obj.W(1), obj.W(2), obj.W(3), 0, 0, 1, 'b', 'AutoScale', 'off');
        end
        
        function [RLA, RLB, RLC, RlA, RlB, RlC] = arms_rotationalMatrices(obj)
            % Returns rotational matrices, where the Z axis is aligned with
            % the direction of the arms.
            % from: https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
            % Visualize it by using:
            % quiver3(gtpr.A2(1), gtpr.A2(2), gtpr.A2(3), gtpr.RlA(1,1), gtpr.RlA(2,1), gtpr.RlA(3,1), 'r');
            % quiver3(gtpr.A2(1), gtpr.A2(2), gtpr.A2(3), gtpr.RlA(1,2), gtpr.RlA(2,2), gtpr.RlA(3,2), 'g');
            % quiver3(gtpr.A2(1), gtpr.A2(2), gtpr.A2(3), gtpr.RlA(1,3), gtpr.RlA(2,3), gtpr.RlA(3,3), 'b');
    
            
            RLA = obj.rotx(obj.theta(1)+90);
            RLB = obj.rotz(-obj.alpha12)*obj.rotx(obj.theta(2)+90);
            RLC = obj.rotz(-obj.alpha13)*obj.rotx(obj.theta(3)+90);
            
            % lower limb A
            rlAv = cross((obj.A2-obj.A1)/norm(obj.A2-obj.A1), (obj.A3-obj.A2)/norm(obj.A3-obj.A2));
            rlAvx = [0, -rlAv(3), rlAv(2);
                     rlAv(3), 0, -rlAv(1);
                     -rlAv(2), rlAv(1), 0];

            rlAc = dot((obj.A2-obj.A1)/norm(obj.A2-obj.A1), (obj.A3-obj.A2)/norm(obj.A3-obj.A2));
            rlAs = norm(rlAv);

            RlA = eye(3) + rlAvx + rlAvx^2*((1-rlAc)/rlAs^2);

            RlA = RlA*RLA;

            % lower limb B
            rlBv = cross((obj.B2-obj.B1)/norm(obj.B2-obj.B1), (obj.B3-obj.B2)/norm(obj.B3-obj.B2));
            rlBvx = [0, -rlBv(3), rlBv(2);
                     rlBv(3), 0, -rlBv(1);
                     -rlBv(2), rlBv(1), 0];

            rlBc = dot((obj.B2-obj.B1)/norm(obj.B2-obj.B1), (obj.B3-obj.B2)/norm(obj.B3-obj.B2));
            rlBs = norm(rlBv);

            RlB = eye(3) + rlBvx + rlBvx^2*((1-rlBc)/rlBs^2);

            RlB = RlB*RLB;

            % lower limb C
            rlCv = cross((obj.C2-obj.C1)/norm(obj.C2-obj.C1), (obj.C3-obj.C2)/norm(obj.C3-obj.C2));
            rlCvx = [0, -rlCv(3), rlCv(2);
                     rlCv(3), 0, -rlCv(1);
                     -rlCv(2), rlCv(1), 0];

            rlCc = dot((obj.C2-obj.C1)/norm(obj.C2-obj.C1), (obj.C3-obj.C2)/norm(obj.C3-obj.C2));
            rlCs = norm(rlCv);

            RlC = eye(3) + rlCvx + rlCvx^2*((1-rlCc)/rlCs^2);
            RlC = RlC*RLC;
        
            obj.RL1 = RLA;
            obj.RL2 = RLB;
            obj.RL3 = RLC;
            
            obj.Rl1 = RlA;
            obj.Rl2 = RlB;
            obj.Rl3 = RlC;
            
        end
        
        function [phi, theta, d, e_perp, m_perp, p, XiL, SlXf] = armAngles(obj, varargin)
            % Angles of arms on the passive joints, from the Screw
            % relations
            %
            % Outputs:
            % w23A = obj.phi.A 
            % w45A = -obj.phi.A  
            % w34A = obj.theta.A
            % w56A = -obj.theta.A
            %   -   p:  Intersection point of common perpendicular. Can have 2 set of 3
            %           element vectors, for 2 separate points, or 1 set of 3 element vector
            %           for a common intersection point, which result in both screws lieing in
            %           the same plane
            %   -   d:  Distance between lines, at common perpendicular points
            %   -   e_perp: cross product of e1,e2 orthonormal of directions
            %   -   m_perp: orthonormal of moments m1,m2
            %   -   alpha:  angle enclosed by the two screws
            
            % sideways angles
            if isempty(obj.S12_1)
                obj.calculateScrews();
            end
            
            [p_S12ASlA, d_S12ASlA, e_perp_S12ASlA, m_perp_S12ASlA, phi_S12ASlA] = obj.SpatialRelationShipOf2Screws(obj.S12_1,obj.Sl1);
%             phi_S12ASlA = phi_S12ASlA;%+90;
            
            [p_S12BSlB, d_S12BSlB, e_perp_S12BSlB, m_perp_S12BSlB, phi_S12BSlB] = obj.SpatialRelationShipOf2Screws(obj.S12_2,obj.Sl2);
%             phi_S12BSlB = phi_S12BSlB;%+90;

            [p_S12CSlC, d_S12CSlC, e_perp_S12CSlC, m_perp_S12CSlC, phi_S12CSlC] = obj.SpatialRelationShipOf2Screws(obj.S12_3,obj.Sl3);
%             phi_S12CSlC = phi_S12CSlC;%+90;
            
            if ~isreal(phi_S12ASlA) && abs(imag(phi_S12ASlA)) < obj.zeroTolerance
                phi_S12ASlA = real(phi_S12ASlA);
            end
            
            if ~isreal(phi_S12BSlB) && abs(imag(phi_S12BSlB)) < obj.zeroTolerance
                phi_S12BSlB = real(phi_S12BSlB);
            end
            
            if ~isreal(phi_S12CSlC) && abs(imag(phi_S12CSlC)) < obj.zeroTolerance
                phi_S12CSlC = real(phi_S12CSlC);
            end
            
            % To calculate the angle between the upper and lower arms, in
            % the plane of the upper arms movement, the bottom screw is
            % offset, so that they do not have a common intersection point
            A2L = obj.A2 + norm(obj.A1)/10 * obj.rotz(pi/2)*[sin(0); cos(0); 0];
            B2L = obj.B2 + norm(obj.B1)/10 * obj.rotz(pi/2)*[sin(obj.alpha12); cos(obj.alpha12); 0]; %#ok<PROPLC>
            C2L = obj.C2 + norm(obj.C1)/10 * obj.rotz(pi/2)*[sin(obj.alpha13); cos(obj.alpha13); 0]; %#ok<PROPLC>
            
            A3L = obj.A3 + norm(obj.A1)/10 * obj.rotz(pi/2)*[sin(0); cos(0); 0];
            B3L = obj.B3 + norm(obj.B1)/10 * obj.rotz(pi/2)*[sin(obj.alpha12); cos(obj.alpha12); 0]; %#ok<PROPLC>
            C3L = obj.C3 + norm(obj.C1)/10 * obj.rotz(pi/2)*[sin(obj.alpha13); cos(obj.alpha13); 0]; %#ok<PROPLC>
            
            % Fictional screw construction:
            SlAf = zeros(6,1);
            SlAf(1:3) = A3L-A2L;% / norm(A3L-A2L);
            SlAf(4:6) = cross(A3L,SlAf(1:3));
            
            SlBf = zeros(6,1);
            SlBf(1:3) = B3L-B2L; %#ok<PROPLC>% / norm(B3L-B2L);
            SlBf(4:6) = cross(B3L,SlBf(1:3)); %#ok<PROPLC>
            
            SlCf = zeros(6,1);
            SlCf(1:3) = C3L-C2L; %#ok<PROPLC> % / norm(C3L-C2L);
            SlCf(4:6) = cross(C3L,SlCf(1:3)); %#ok<PROPLC>

            [p_SLASlAf, d_SLASlAf, e_perp_SLASlAf, m_perp_SLASlAf, theta_SLASlAf] = obj.SpatialRelationShipOf2Screws(SlAf,obj.SL1);
            [p_SLBSlBf, d_SLBSlBf, e_perp_SLBSlBf, m_perp_SLBSlBf, theta_SLBSlBf] = obj.SpatialRelationShipOf2Screws(SlBf,obj.SL2);
            [p_SLCSlCf, d_SLCSlCf, e_perp_SLCSlCf, m_perp_SLCSlCf, theta_SLCSlCf] = obj.SpatialRelationShipOf2Screws(SlCf,obj.SL3);
            
            
            %
            obj.phi = struct("A", phi_S12ASlA, "B", phi_S12BSlB, "C", phi_S12CSlC);
            obj.beta = struct("A", theta_SLASlAf, "B", theta_SLBSlBf, "C", theta_SLCSlCf);
            
            obj.p = struct("S12A_SlA", p_S12ASlA, "S12B_SlB", p_S12BSlB, "S12C_SlC", p_S12CSlC,...
                       "SLA_SlAf", p_SLASlAf, "SLB_SlBf", p_SLBSlBf, "SLC_SlCf", p_SLCSlCf);
            
            obj.d = struct("S12A_SlA", d_S12ASlA, "S12B_SlB", d_S12BSlB, "S12C_SlC", d_S12CSlC,...
                       "SLA_SlAf", d_SLASlAf, "SLB_SlBf", d_SLBSlBf, "SLC_SlCf", d_SLCSlCf);
                          
            obj.e_perp = struct("S12A_SlA", e_perp_S12ASlA, "S12B_SlB", e_perp_S12BSlB, "S12C_SlC", e_perp_S12CSlC,...
                            "SLA_SlAf", e_perp_SLASlAf, "SLB_SlBf", e_perp_SLBSlBf, "SLC_SlCf", e_perp_SLCSlCf);
                          
            obj.m_perp = struct("S12A_SlA", m_perp_S12ASlA, "S12B_SlB", m_perp_S12BSlB, "S12C_SlC", m_perp_S12CSlC,...
                            "SLA_SlAf", m_perp_SLASlAf, "SLB_SlBf", m_perp_SLBSlBf, "SLC_SlCf", m_perp_SLCSlCf);
            
            obj.XiL = struct("A2L", A2L, "A3L", A3L, "B2L", B2L, "B3L", B3L, "C2L", C2L, "C3L", C3L); %#ok<PROPLC>
            obj.SlXf = struct("SlAf", SlAf, "SlBf", SlBf, "SlCf", SlCf);
            
            phi = obj.phi; theta = obj.theta; d = obj.d; e_perp = obj.e_perp; m_perp = obj.m_perp; p = obj.p; XiL = obj.XiL; SlXf = obj.SlXf;
        end
        
        function figureHandle = plotJointSpace(obj, varargin)
           % Plots the joint space of the given obj.
           % If no 'theta' varargin argument is given, it works with whatever
           % is inside the object (obj.theta).
           
           if any(strcmp(varargin, 'theta'))
               idx = find(strcmp(varargin, 'theta'));
               obj.theta = varargin{idx+1};
           end
           
%            props = obj.extractPlotArguments(varargin);
           if any(strcmp(varargin, 'figureNo'))
               idx = find(strcmp(varargin, 'figureNo'));
               figureHandle = figure(varargin{idx+1});
           else
               figureHandle = figure();
           end
           
           
           hold on;
           
           plotProps = obj.extractPlotArguments(varargin);
           
           cmd = [];
           
           if ~isempty(plotProps.Color)
               cmd = strcat(cmd, "'Color', plotProps.Color,");
           end
           
           if ~isempty(plotProps.LineStyle)
               cmd = strcat(cmd, "'LineStyle', plotProps.LineStyle,");
           end
           
           if ~isempty(plotProps.LineWidth)
               cmd = strcat(cmd, "'LineWidth', plotProps.LineWidth,");
           end
           
           if ~isempty(plotProps.Marker)
               cmd = strcat(cmd, "'Marker', plotProps.Marker,");
           end
           
           if ~isempty(plotProps.MarkerIndices)
               cmd = strcat(cmd, "'MarkerIndices', plotProps.MarkerIndices,");
           end
           
           if ~isempty(plotProps.MarkerEdgeColor)
               cmd = strcat(cmd, "'MarkerEdgeColor', plotProps.MarkerEdgeColor,");
           end
           
           if ~isempty(plotProps.MarkerFaceColor)
               cmd = strcat(cmd, "'MarkerFaceColor', plotProps.MarkerFaceColor,");
           end
           
           if ~isempty(plotProps.MarkerSize)
               cmd = strcat(cmd, "'MarkerSize', plotProps.MarkerSize,");
           end
           
           if ~isempty(plotProps.DatetimeTickFormat)
               cmd = strcat(cmd, "'DatetimeTickFormat', plotProps.DatetimeTickFormat,");
           end
           
           if ~isempty(plotProps.DurationTickFormat)
               cmd = strcat(cmd, "'DurationTickFormat', plotProps.DurationTickFormat,");
           end
           formatCmd = cmd;
           
           for i = 1:size(obj.theta,1)
               if isempty(formatCmd)
                   cmd = strcat("plot3(obj.theta(i,1), obj.theta(i,2), obj.theta(i,3));");
               else
                cmd = strcat("plot3(obj.theta(i,1), obj.theta(i,2), obj.theta(i,3), ", strip(formatCmd, 'right', ','), ");");
               end
               eval(cmd);
           end

           
%            plot3(obj.theta(1,:), obj.theta(2,:), obj.theta(3,:), varargin);
           hold off;
           
        end
        
        
        
        
    end
    
    methods
        
        
        
        
        function baseTrianglePts = getBaseTrianglePts(obj,a1, a2, a3, alpha12, alpha13)
            obj.A1 = [0; -a1; 0]; 
            obj.A2 = [-a2*sin(alpha12); -a2*cos(alpha12); 0]; 
            obj.A3 = [-a3*sin(alpha13); -a3*cos(alpha13); 0]; 
            
            baseTrianglePts = struct("A1", obj.A1, "A2", obj.A2, "A3", obj.A3); 
        end
        
        function workTrianglePts = getWorkTrianglePts(obj, W, e1, e2, e3, alpha12, alpha13)
            obj.C1 = (W - [0; e1; 0]);
            obj.C2 = (W - [e2*sin(alpha12); e2*cos(alpha12); 0]);
            obj.C3 = (W - [e3*sin(alpha13); e3*cos(alpha13); 0]);
            
            if size(obj.C1,1) ~= 3
                obj.C1 = obj.C1';
            end
            
            if size(obj.C2,1) ~= 3
                obj.C2 = obj.C2';
            end
            
            if size(obj.C3,1) ~= 3
                obj.C3 = obj.C3';
            end
            
            workTrianglePts = struct("C1", obj.C1, "C2", obj.C2, "C3", obj.C3);
        end
        
        function [UpperIntersectionPoint, LowerIntersectionPoint] = threeSpheresIntersection(obj, c1, c2, c3, r1, r2, r3)
            % Three-speheres intersection algorithm from:
            % The Delta Parallel Robot: Kinematics Solutions
            % Robert L. Williams II, Ph.D., williar4@ohio.edu
            % Mechanical Engineering, Ohio University, October 2016
            % Appendix A-B.
            % This function uses both algorithms, and chooses the one that is solvable.

            % UpperIntersectionPoint = nan;
            % LowerIntersectionPoint = interx(c1,c2,c3,r1,r2,r3,0);
            
            % Precision for the comparison
%             e = 1e-8;

            if size(c1,2) ~= 3
                c1 = c1';
            end

            if size(c2,2) ~= 3
                c2 = c2';
            end
            
            if size(c3,2) ~= 3
                c3 = c3';
            end
            
    
            if c1(3) == c2(3) && c2(3) == c3(3)
                [LowerIntersectionPoint, UpperIntersectionPoint] = obj.SimplifiedThreeSpheresIntersectionAlgorithm(c1, c2, c3, r1, r2, r3);
            else
                [UpperIntersectionPoint, LowerIntersectionPoint] = obj.threeSpheresIntersectionAlgorithm(c1, c2, c3, r1, r2, r3);
                
                if ~any(isreal(LowerIntersectionPoint)) || any(abs(LowerIntersectionPoint) > 100) || any(isnan(UpperIntersectionPoint)) || any(isnan(LowerIntersectionPoint))
                    clear UpperIntersectionPoint LowerIntersectionPoint;
                    [UpperIntersectionPoint, LowerIntersectionPoint] = obj.threeSpheresIntersectionAlgorithm(c1, c3, c2, r1, r3, r2);
                end
                
                if ~any(isreal(LowerIntersectionPoint)) || any(abs(LowerIntersectionPoint) > 100) || any(isnan(UpperIntersectionPoint)) || any(isnan(LowerIntersectionPoint))
                    clear UpperIntersectionPoint LowerIntersectionPoint;
                    [UpperIntersectionPoint, LowerIntersectionPoint] = obj.threeSpheresIntersectionAlgorithm(c2, c1, c3, r2, r1, r3);
                end
                
                if ~any(isreal(LowerIntersectionPoint)) || any(abs(LowerIntersectionPoint) > 100) || any(isnan(UpperIntersectionPoint)) || any(isnan(LowerIntersectionPoint))
                    clear UpperIntersectionPoint LowerIntersectionPoint;
                    [UpperIntersectionPoint, LowerIntersectionPoint] = obj.threeSpheresIntersectionAlgorithm(c2, c3, c1, r2, r3, r1);
                end
                
                if ~any(isreal(LowerIntersectionPoint)) || any(abs(LowerIntersectionPoint) > 100) || any(isnan(UpperIntersectionPoint)) || any(isnan(LowerIntersectionPoint))
                    clear UpperIntersectionPoint LowerIntersectionPoint;
                    [UpperIntersectionPoint, LowerIntersectionPoint] = obj.threeSpheresIntersectionAlgorithm(c3, c1, c2, r3, r1, r2);
                end
                
                if ~any(isreal(LowerIntersectionPoint)) || any(abs(LowerIntersectionPoint) > 100) || any(isnan(UpperIntersectionPoint)) || any(isnan(LowerIntersectionPoint))
                    clear UpperIntersectionPoint LowerIntersectionPoint;
                    [UpperIntersectionPoint, LowerIntersectionPoint] = obj.threeSpheresIntersectionAlgorithm(c3, c2, c1, r3, r2, r1);
                end
                
                if ~any(isreal(LowerIntersectionPoint)) || any(abs(LowerIntersectionPoint) > 100) || any(isnan(UpperIntersectionPoint)) || any(isnan(LowerIntersectionPoint))
                    error("Impossible Configuration present!");
                end
                
            end
            
        end
        
        function [UpperIntersectionPoint, LowerIntersectionPoint] = SimplifiedThreeSpheresIntersectionAlgorithm(obj, c1, c2, c3, r1, r2, r3)
            % Simplified Three-speheres intersection algorithm from:
            % The Delta Parallel Robot: Kinematics Solutions
            % Robert L. Williams II, Ph.D., williar4@ohio.edu
            % Mechanical Engineering, Ohio University, October 2016
            % Appendix B.

            x1 = c1(1);
            y1 = c1(2);
            z1 = c1(3);

            x2 = c2(1);
            y2 = c2(2);
            z2 = c2(3);

            x3 = c3(1);
            y3 = c3(2);
            z3 = c3(3);

            if (all(z1 ~= z2) && all(z2 ~= z3) && all(z3 ~= z1)) && obj.errors == true
                error("Sphere centers are not in the same height!");
            end

            a = 2*(x3-x1); %#ok<PROPLC>
            b = 2*(y3-y1);
            c = power(r1,2)-power(r3,2)-power(x1,2)-power(y1,2)+power(x3,2)+power(y3,2);
            d = 2*(x3-x2); %#ok<PROPLC>
            e = 2*(y3-y2);
            f = power(r2,2)-power(r3,2)-power(x2,2)-power(y2,2)+power(x3,2)+power(y3,2);

            aebd = a.*e-b.*d; %#ok<PROPLC>

            if any(aebd == 0) && obj.errors == true
                error("First singularity condition fulfilled! See the appendix B. of The Delta Parallel Robot: Kinematics Solutions Robert L. Williams II, Ph.D., williar4@ohio.edu, Mechanical Engineering, Ohio University, October 2016 article for explanation!");
            end

            x = (c.*e-b.*f)./aebd;
            y = (a.*f-c.*d)./aebd; %#ok<PROPLC>

            A = 1;
            B = -2*z1;
            C = power(z1,2)-power(r1,2)+power(x-x1,2)+power(y-y1,2);

            z = roots([A,B,C]);
            
            [warnMsg, ~] = lastwarn;
            if strcmp(warnMsg,'Imaginary solution for sqrt(power(b,2) - 4 * a * c).')
                warning('Not all 3 spheres intersect! Robot assembly is not fulfilled!')
            end
            
            UpperIntersectionPoint = [x, y, z(2,:)]';
            LowerIntersectionPoint = [x, y, z(1,:)]';

        end
        
        
        
        function [UpperIntersectionPoint, LowerIntersectionPoint] = threeSpheresIntersectionAlgorithm(~, c1, c2, c3, r1, r2, r3)
            % Three-speheres intersection algorithm from:
            % The Delta Parallel Robot: Kinematics Solutions
            % Robert L. Williams II, Ph.D., williar4@ohio.edu
            % Mechanical Engineering, Ohio University, October 2016
            % Appendix A.

            x1 = c1(1);
            y1 = c1(2);
            z1 = c1(3);

            x2 = c2(1);
            y2 = c2(2);
            z2 = c2(3);

            x3 = c3(1);
            y3 = c3(2);
            z3 = c3(3);


            a11 = 2*(x3-x1);
            a12 = 2*(y3-y1);
            a13 = 2*(z3-z1);

            a21 = 2*(x3-x2);
            a22 = 2*(y3-y2);
            a23 = 2*(z3-z2);

            b1 = power(r1,2)-power(r3,2)-power(x1,2)-power(y1,2)-power(z1,2)+power(x3,2)+power(y3,2)+power(z3,2);
            b2 = power(r2,2)-power(r3,2)-power(x2,2)-power(y2,2)-power(z2,2)+power(x3,2)+power(y3,2)+power(z3,2);

            a1 = (a11./a13)-(a21./a23); %#ok<PROPLC>
            a2 = (a12./a13)-(a22./a23); %#ok<PROPLC>
            a3 = (b2./a23)-(b1./a13); %#ok<PROPLC>
            a4 = -(a2./a1); %#ok<PROPLC>
            a5 = -(a3./a1); %#ok<PROPLC>
            a6 = (-a21.*a4-a22)./a23;
            a7 = (b2-a21.*a5)./a23;


            a = power(a4,2) + 1 + power(a6,2); %#ok<PROPLC>

            % Singularity conditions:
            if (any(a13 == 0) || any(a23 == 0) || any(a1 == 0) || any(a == 0)) %#ok<PROPLC>
                warning('Singularity condition present! Sphere centers are in same Z height!');
            end

            b = 2.*a4.*(a5-x1)-2.*y1+2.*a6.*(a7-z1);
            c = a5.*(a5-2.*x1)+a7.*(a7-2.*z1)+power(x1,2)+power(y1,2)+power(z1,2)-power(r1,2);

            
            try
                y = roots([a,b,c]); %#ok<PROPLC>
            catch ME
                if contains(ME.message, "Input to ROOTS must not contain NaN or Inf.")
                    UpperIntersectionPoint = 1i;
                    LowerIntersectionPoint = 1i;
                    return;
                end
            end
            if size(y,2) ~= 2
                y = y';
            end
            
            [warnMsg, ~] = lastwarn;

            if strcmp(warnMsg,'Imaginary solution for sqrt(power(b,2) - 4 * a * c).') 
                warning('Not all 3 spheres intersect!')
            end
            
%             x = nan(length(a),2);
%             z = nan(length(a),2);

            x(:,1) = a4.*y(:,1)+a5;
            x(:,2) = a4.*y(:,2)+a5;

            z(:,1) = a6.*y(:,1)+a7;
            z(:,2) = a6.*y(:,2)+a7;

            [rows, ~] = find(z(:,2) > z(:,1));
            
            temp = z(rows,2);
            z(rows,2) = z(rows,1);
            z(rows,1) = temp;
            
            temp = x(rows,2);
            x(rows,2) = x(rows,1);
            x(rows,1) = temp;
            
            temp = y(rows,2);
            y(rows,2) = y(rows,1);
            y(rows,1) = temp;
            
            UpperIntersectionPoint = [x(:,1), y(:,1), z(:,1)];
            LowerIntersectionPoint = [x(:,2), y(:,2), z(:,2)];
            
        end
        
        function props = extractPlotArguments(~,inputs)
           
            Color = [];
            LineStyle = [];
            LineWidth = [];
            Marker = [];
            MarkerIndices = [];
            MarkerEdgeColor = [];
            MarkerFaceColor = [];
            MarkerSize = [];
            DatetimeTickFormat = [];
            DurationTickFormat = [];
            
            propName = 'Color';
            if any(strcmp(inputs, propName))
                idx = find(strcmp(inputs, propName));
                Color = inputs{idx+1};
            end
            
            propName = 'LineStyle';
            if any(strcmp(inputs, propName))
                idx = find(strcmp(inputs, propName));
                LineStyle = inputs{idx+1};
            end
            
            propName = 'LineWidth';
            if any(strcmp(inputs, propName))
                idx = find(strcmp(inputs, propName));
                LineWidth = inputs{idx+1};
            end
            
            propName = 'Marker';
            if any(strcmp(inputs, propName))
                idx = find(strcmp(inputs, propName));
                Marker = inputs{idx+1};
            end
            
            propName = 'MarkerIndices';
            if any(strcmp(inputs, propName))
                idx = find(strcmp(inputs, propName));
                MarkerIndices = inputs{idx+1};
            end
            
            propName = 'MarkerEdgeColor';
            if any(strcmp(inputs, propName))
                idx = find(strcmp(inputs, propName));
                MarkerEdgeColor = inputs{idx+1};
            end
            
            propName = 'MarkerFaceColor';
            if any(strcmp(inputs, propName))
                idx = find(strcmp(inputs, propName));
                MarkerFaceColor = inputs{idx+1};
            end
            
            propName = 'MarkerSize';
            if any(strcmp(inputs, propName))
                idx = find(strcmp(inputs, propName));
                MarkerSize = inputs{idx+1};
            end
            
            propName = 'DatetimeTickFormat';
            if any(strcmp(inputs, propName))
                idx = find(strcmp(inputs, propName));
                DatetimeTickFormat = inputs{idx+1};
            end
            
            propName = 'DurationTickFormat';
            if any(strcmp(inputs, propName))
                idx = find(strcmp(inputs, propName));
                DurationTickFormat = inputs{idx+1};
            end
            
            props = struct( 'Color'             ,   Color           ,...
                            'LineStyle'         ,   LineStyle       ,...
                            'LineWidth'         ,   LineWidth       ,...
                            'Marker'            ,   Marker          ,...
                            'MarkerIndices'     ,   MarkerIndices   ,...
                            'MarkerEdgeColor'   ,   MarkerEdgeColor ,...
                            'MarkerFaceColor'   ,   MarkerFaceColor ,...
                            'MarkerSize'        ,   MarkerSize      ,...
                            'DatetimeTickFormat',   DatetimeTickFormat ,...
                            'DurationTickFormat',   DurationTickFormat);
        end
        
        function [p, d, e_perp, m_perp, alpha] = SpatialRelationShipOf2Screws(obj,S1,S2)
        % Based off of:
        % http://web.cs.iastate.edu/~cs577/handouts/plucker-coordinates.pdf
        % Chapter 5
        %
        % Inputs:
        %   - S1,S2: Screws, with normalized directions 2 set of 6 element vectors
        % Outputs:
        %   -   p:  Intersection point of common perpendicular. Can have 2 set of 3
        %           element vectors, for 2 separate points, or 1 set of 3 element vector
        %           for a common intersection point, which result in both screws lieing in
        %           the same plane
        %   -   d:  Distance between lines, at common perpendicular points
        %   -   e_perp: cross product of e1,e2 orthonormal of directions
        %   -   m_perp: orthonormal of moments m1,m2
        %   -   alpha:  angle enclosed by the two screws


            if size(S1,1) ~= 6
                S1 = S1';
            end

            if size(S2,1) ~= 6
                S2 = S2';
            end
            l1 = S1(1:3); %#ok<PROPLC>
            m1 = S1(4:6);
            l2 = S2(1:3); %#ok<PROPLC>
            m2 = S2(4:6);

            p = [];
            d = [];
            e_perp = [];
            m_perp = [];
            alpha = [];

            % 1) Check if e1xe2 = 0 holds:
            I = norm(cross(l1, l2)); %#ok<PROPLC>

            if I == 0 %they are parallel with distance d (eq 10)
                d = norm(cross(l1, (m1-m2/1)))/norm(l1)^2; %#ok<PROPLC> % s = 1
            else
                % 2) Check reciprocal stuff
                reciprocal = obj.kleinForm([l1;m1],[l2;m2]); %#ok<PROPLC>
                if reciprocal == 0 %2a) the 2 lines intersect at point p* given in eq 29)
                    p = (dot(m1,l2)*eye(3) + l1*m2' - l2*m1')*(cross(l1,l2)/(norm(cross(l1,l2)))^2); %#ok<PROPLC>
                else
                    % 2b) the 2 lines are not coplanar...
                    d = norm(reciprocal)/norm(cross(l1,m1)); %#ok<PROPLC>
                    e_perp = cross(l1,m1); %#ok<PROPLC>
                    m_perp = cross(m1,m1) - cross(m2,l1) + ((dot(S1,S2)*(dot(l1,m1)))/(norm(cross(l1,m1))^2))*cross(l1,m1); %#ok<PROPLC>

                    p(1:3,1) = (cross(-m1, cross(l2, cross(l1,l2))) + dot(m2,cross(l1, l2))*l1)/(norm(cross(l1,l2))^2); %#ok<PROPLC>
                    p(1:3,2) = (cross(m2,cross(l1,cross(l1, l2))) - dot(m1, cross(l1, l2))*l2)/(norm(cross(l1,l2))^2); %#ok<PROPLC>
                end
            end


            if size(p,2) > 1
                m21 = m2-cross(p(1:3,2)-p(1:3,1),l2); %#ok<PROPLC>

                alpha = -asin((dot(l1,m21)/norm(m21))); %#ok<PROPLC>

            end


        end

        function S = kleinForm(~,S1, S2)

            S = dot(S1(1:3),S2(4:6)) + dot(S1(4:6),S2(1:3));

        end

        function E = rotx(~, theta)
            
            c = cos(theta);
            s = sin(theta);

            E = [ 1 0  0
                  0 c -s;
                  0 s  c];
        end

        function  E = rotz(~, theta )

            % rz  3x3 coordinate rotation (Z-axis)
            % rz(theta)  calculates the 3x3 rotational coordinate transform matrix from
            % A to B coordinates, where coordinate frame B is rotated by an angle theta
            % (radians) relative to frame A about their common Z axis.

            c = cos(theta);
            s = sin(theta);

            E = [ c  -s  0;
                  s   c  0;
                  0   0  1 ];

        end
        
        function data = sphereCircleIntersection(~, cc, rc, nUnit, cs, rs)
            % https://gamedev.stackexchange.com/questions/75756/sphere-sphere-intersection-and-circle-sphere-intersection
            if size(cc) ~= [3,1] %#ok<BDSCA>
                cc = cc';
            end

            if size(cs) ~= [3,1] %#ok<BDSCA>
                cs = cs';
            end

            if size(nUnit) ~= [3,1] %#ok<BDSCA>
                nUnit = nUnit';
            end

            d = dot(nUnit, cc-cs); %#ok<PROPLC>

            if abs(d) > rs %#ok<PROPLC>
                data.p0 = [];
                data.p1 = [];
                data.p2 = [];
                return % no intersection
            end

            cp = cs + d*nUnit; %#ok<PROPLC>

            rp = sqrt(rs^2 - d^2); %#ok<PROPLC>

            % now the circle circle intersection can be introduced
            dc = norm(cc - cp);

            h = 0.5 + ((rc^2-rp^2)/(2*dc^2));

            ci = cc + h*(cc-cp);

            ri = sqrt(rc^2 - h^2*dc^2);

            % ni = (cc - cp)/d;

            t = cross(cp-cc,nUnit)/norm(cross(cp-cc,nUnit));

            p0 = cc - t*ri + (cc-ci); % Added part for compensation and ci has been modified to cc in both cases
            p1 = cc + t*ri + (cc-ci);

            data.ci = ci;
            data.p0 = p0;
            data.p1 = p1;

        end
        
        function [ws_cube, PtsMapped] = mapWS(obj, ws_cube, theta, varargin)
            
            resolution = ws_cube.resolution/2;
            
            % Map Origin
            [ws_cube, ws_O] = obj.mapBitPointInSpace(ws_cube, [0 0 0], obj.bit(11).Position, obj.bit(11).Name, 1);
            plotcube(resolution*[1 1 1], [0 0 0], obj.bit(11).Opacity, obj.bit(11).Color);
            hold on;
            
            obj.directGeometry(theta);

            if any(strcmpi(varargin, "depiction"))
                idx = find(strcmpi(varargin, "depiction"));
                depiction = varargin{idx+1};
            else
                depiction = "";
            end

            if any(strcmpi(varargin, "plotMechanism"))
                obj.plotMechanism(depiction);
            end
            

            % Map the lines to the points:
            [ws_cube, ws_TCPValuesMapped] = obj.mapBitPointInSpace(ws_cube, obj.W, obj.bit(1).Position, obj.bit(1).Name, 1);

            [ws_cube, ws_A1PtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.A1, obj.bit(2).Position, obj.bit(2).Name, 1);
            [ws_cube, ws_A2PtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.A2, obj.bit(3).Position, obj.bit(3).Name, 1);
            [ws_cube, ws_A3PtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.A3, obj.bit(4).Position, obj.bit(4).Name, 1);

            [ws_cube, ws_C1PtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.C1, obj.bit(8).Position, obj.bit(8).Name, 1);
            [ws_cube, ws_C2PtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.C2, obj.bit(9).Position, obj.bit(9).Name, 1);
            [ws_cube, ws_C3PtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.C3, obj.bit(10).Position, obj.bit(10).Name, 1);
            
            [ws_cube, ws_B1PtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.B1, obj.bit(5).Position, obj.bit(5).Name, 1);
            [ws_cube, ws_B2PtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.B2, obj.bit(6).Position, obj.bit(6).Name, 1);
            [ws_cube, ws_B3PtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.B3, obj.bit(7).Position, obj.bit(7).Name, 1);
            
            if any(strcmpi(varargin, "plotCubes"))
                plotcube(resolution*[1 1 1], obj.A1', obj.bit(2).Opacity, obj.bit(2).Color);
                plotcube(resolution*[1 1 1], obj.A2', obj.bit(3).Opacity, obj.bit(3).Color);
                plotcube(resolution*[1 1 1], obj.A3', obj.bit(4).Opacity, obj.bit(4).Color);

                plotcube(resolution*[1 1 1], obj.W', obj.bit(1).Opacity, obj.bit(1).Color);
                
                obj.plotBresenhamLine([ws_A1PtsMapped.x.idx, ws_A1PtsMapped.y.idx, ws_A1PtsMapped.z.idx], ...
                [ws_B1PtsMapped.x.idx, ws_B1PtsMapped.y.idx, ws_B1PtsMapped.z.idx], ...
                ws_cube, resolution, obj.bit(15).Opacity, obj.bit(15).Color);

                obj.plotBresenhamLine([ws_A2PtsMapped.x.idx, ws_A2PtsMapped.y.idx, ws_A2PtsMapped.z.idx], ...
                [ws_B2PtsMapped.x.idx, ws_B2PtsMapped.y.idx, ws_B2PtsMapped.z.idx], ...
                ws_cube, resolution, obj.bit(16).Opacity, obj.bit(16).Color);

                obj.plotBresenhamLine([ws_A3PtsMapped.x.idx, ws_A3PtsMapped.y.idx, ws_A3PtsMapped.z.idx], ...
                [ws_B3PtsMapped.x.idx, ws_B3PtsMapped.y.idx, ws_B3PtsMapped.z.idx], ...
                ws_cube, resolution, obj.bit(17).Opacity, obj.bit(17).Color);
            end
            
            PtsMapped.O = ws_O;

            PtsMapped.A1 = ws_A1PtsMapped;
            PtsMapped.A2 = ws_A2PtsMapped;
            PtsMapped.A3 = ws_A3PtsMapped;


            PtsMapped.C1 = ws_C1PtsMapped;
            PtsMapped.C2 = ws_C2PtsMapped;
            PtsMapped.C3 = ws_C3PtsMapped;
            
            PtsMapped.W = ws_TCPValuesMapped;
            
            obj.realGeometry();
            
            if any(strcmpi(depiction, 'real')) || any(strcmpi(varargin, "plotCubes"))
                [ws_cube, ws_B1LPtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.B1L, obj.bit(24).Position, obj.bit(24).Name, 1);
                [ws_cube, ws_B1RPtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.B1R, obj.bit(25).Position, obj.bit(25).Name, 1);

                [ws_cube, ws_B2LPtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.B2L, obj.bit(26).Position, obj.bit(24).Name, 1);
                [ws_cube, ws_B2RPtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.B2R, obj.bit(27).Position, obj.bit(25).Name, 1);

                [ws_cube, ws_B3LPtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.B3L, obj.bit(28).Position, obj.bit(24).Name, 1);
                [ws_cube, ws_B3RPtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.B3R, obj.bit(29).Position, obj.bit(25).Name, 1);

                [ws_cube, ws_C1LPtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.C1L, obj.bit(30).Position, obj.bit(30).Name, 1);
                [ws_cube, ws_C1RPtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.C1R, obj.bit(31).Position, obj.bit(31).Name, 1);

                [ws_cube, ws_C2LPtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.C2L, obj.bit(32).Position, obj.bit(32).Name, 1);
                [ws_cube, ws_C2RPtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.C2R, obj.bit(33).Position, obj.bit(33).Name, 1);

                [ws_cube, ws_C3LPtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.C3L, obj.bit(34).Position, obj.bit(34).Name, 1);
                [ws_cube, ws_C3RPtsMapped] = obj.mapBitPointInSpace(ws_cube, obj.C3R, obj.bit(35).Position, obj.bit(35).Name, 1);

                
                PtsMapped.B1L = ws_B1LPtsMapped;
                PtsMapped.B1R = ws_B1RPtsMapped;

                PtsMapped.B2L = ws_B2LPtsMapped;
                PtsMapped.B2R = ws_B2RPtsMapped;

                PtsMapped.B3L = ws_B3LPtsMapped;
                PtsMapped.B3R = ws_B3RPtsMapped;

                PtsMapped.C1L = ws_C1LPtsMapped;
                PtsMapped.C1R = ws_C1RPtsMapped;

                PtsMapped.C2L = ws_C2LPtsMapped;
                PtsMapped.C2R = ws_C2RPtsMapped;

                PtsMapped.C3L = ws_C3LPtsMapped;
                PtsMapped.C3R = ws_C3RPtsMapped;
                
                obj.plotBresenhamLine([ws_A1PtsMapped.x.idx, ws_A1PtsMapped.y.idx, ws_A1PtsMapped.z.idx], ...
                [ws_O.x.idx, ws_O.y.idx, ws_O.z.idx], ...
                ws_cube, resolution, obj.bit(49).Opacity, obj.bit(49).Color);

                obj.plotBresenhamLine([ws_A2PtsMapped.x.idx, ws_A2PtsMapped.y.idx, ws_A2PtsMapped.z.idx], ...
                [ws_O.x.idx, ws_O.y.idx, ws_O.z.idx], ...
                ws_cube, resolution, obj.bit(49).Opacity, obj.bit(49).Color);

                obj.plotBresenhamLine([ws_A3PtsMapped.x.idx, ws_A3PtsMapped.y.idx, ws_A3PtsMapped.z.idx], ...
                [ws_O.x.idx, ws_O.y.idx, ws_O.z.idx], ...
                ws_cube, resolution, obj.bit(49).Opacity, obj.bit(49).Color);

                obj.plotBresenhamLine([ws_A1PtsMapped.x.idx, ws_A1PtsMapped.y.idx, ws_A1PtsMapped.z.idx], ...
                [ws_A2PtsMapped.x.idx, ws_A2PtsMapped.y.idx, ws_A2PtsMapped.z.idx], ...
                ws_cube, resolution, obj.bit(49).Opacity, obj.bit(49).Color);

                obj.plotBresenhamLine([ws_A2PtsMapped.x.idx, ws_A2PtsMapped.y.idx, ws_A2PtsMapped.z.idx], ...
                [ws_A3PtsMapped.x.idx, ws_A3PtsMapped.y.idx, ws_A3PtsMapped.z.idx], ...
                ws_cube, resolution, obj.bit(49).Opacity, obj.bit(49).Color);

                obj.plotBresenhamLine([ws_A1PtsMapped.x.idx, ws_A1PtsMapped.y.idx, ws_A1PtsMapped.z.idx], ...
                [ws_A3PtsMapped.x.idx, ws_A3PtsMapped.y.idx, ws_A3PtsMapped.z.idx], ...
                ws_cube, resolution, obj.bit(49).Opacity, obj.bit(49).Color);

                obj.plotBresenhamLine([ws_C1PtsMapped.x.idx, ws_C1PtsMapped.y.idx, ws_C1PtsMapped.z.idx], ...
                [ws_TCPValuesMapped.x.idx, ws_TCPValuesMapped.y.idx, ws_TCPValuesMapped.z.idx], ...
                ws_cube, resolution, obj.bit(48).Opacity, obj.bit(48).Color);

                obj.plotBresenhamLine([ws_C2PtsMapped.x.idx, ws_C2PtsMapped.y.idx, ws_C2PtsMapped.z.idx], ...
                [ws_TCPValuesMapped.x.idx, ws_TCPValuesMapped.y.idx, ws_TCPValuesMapped.z.idx], ...
                ws_cube, resolution, obj.bit(48).Opacity, obj.bit(48).Color);

                obj.plotBresenhamLine([ws_C3PtsMapped.x.idx, ws_C3PtsMapped.y.idx, ws_C3PtsMapped.z.idx], ...
                [ws_TCPValuesMapped.x.idx, ws_TCPValuesMapped.y.idx, ws_TCPValuesMapped.z.idx], ...
                ws_cube, resolution, obj.bit(48).Opacity, obj.bit(48).Color);

                obj.plotBresenhamLine([ws_C1PtsMapped.x.idx, ws_C1PtsMapped.y.idx, ws_C1PtsMapped.z.idx], ...
                [ws_C2PtsMapped.x.idx, ws_C2PtsMapped.y.idx, ws_C2PtsMapped.z.idx], ...
                ws_cube, resolution, obj.bit(48).Opacity, obj.bit(48).Color);

                obj.plotBresenhamLine([ws_C2PtsMapped.x.idx, ws_C2PtsMapped.y.idx, ws_C2PtsMapped.z.idx], ...
                [ws_C3PtsMapped.x.idx, ws_C3PtsMapped.y.idx, ws_C3PtsMapped.z.idx], ...
                ws_cube, resolution, obj.bit(48).Opacity, obj.bit(48).Color);

                obj.plotBresenhamLine([ws_C1PtsMapped.x.idx, ws_C1PtsMapped.y.idx, ws_C1PtsMapped.z.idx], ...
                [ws_C3PtsMapped.x.idx, ws_C3PtsMapped.y.idx, ws_C3PtsMapped.z.idx], ...
                ws_cube, resolution, obj.bit(48).Opacity, obj.bit(48).Color);
                
                if any(strcmpi(varargin, "plotCubes"))
                    plotcube(resolution*[1 1 1], obj.B1L', obj.bit(24).Opacity, obj.bit(24).Color);
                    plotcube(resolution*[1 1 1], obj.B1R', obj.bit(25).Opacity, obj.bit(25).Color);
                    plotcube(resolution*[1 1 1], obj.B2L', obj.bit(26).Opacity, obj.bit(26).Color);
                    plotcube(resolution*[1 1 1], obj.B2R', obj.bit(27).Opacity, obj.bit(27).Color);
                    plotcube(resolution*[1 1 1], obj.B3L', obj.bit(28).Opacity, obj.bit(28).Color);
                    plotcube(resolution*[1 1 1], obj.B3R', obj.bit(29).Opacity, obj.bit(29).Color);

                    plotcube(resolution*[1 1 1], obj.C1L', obj.bit(30).Opacity, obj.bit(30).Color);
                    plotcube(resolution*[1 1 1], obj.C1R', obj.bit(31).Opacity, obj.bit(31).Color);
                    plotcube(resolution*[1 1 1], obj.C2L', obj.bit(32).Opacity, obj.bit(32).Color);
                    plotcube(resolution*[1 1 1], obj.C2R', obj.bit(33).Opacity, obj.bit(33).Color);
                    plotcube(resolution*[1 1 1], obj.C3L', obj.bit(34).Opacity, obj.bit(34).Color);
                    plotcube(resolution*[1 1 1], obj.C3R', obj.bit(35).Opacity, obj.bit(35).Color);
                    
                    obj.plotBresenhamLine([ws_B1LPtsMapped.x.idx, ws_B1LPtsMapped.y.idx, ws_B1LPtsMapped.z.idx], ...
                    [ws_C1LPtsMapped.x.idx, ws_C1LPtsMapped.y.idx, ws_C1LPtsMapped.z.idx], ...
                    ws_cube, resolution, obj.bit(36).Opacity, obj.bit(36).Color);
                
                    obj.plotBresenhamLine([ws_B1RPtsMapped.x.idx, ws_B1RPtsMapped.y.idx, ws_B1RPtsMapped.z.idx], ...
                    [ws_C1RPtsMapped.x.idx, ws_C1RPtsMapped.y.idx, ws_C1RPtsMapped.z.idx], ...
                    ws_cube, resolution, obj.bit(37).Opacity, obj.bit(37).Color);
                
                    obj.plotBresenhamLine([ws_B2LPtsMapped.x.idx, ws_B2LPtsMapped.y.idx, ws_B2LPtsMapped.z.idx], ...
                    [ws_C2LPtsMapped.x.idx, ws_C2LPtsMapped.y.idx, ws_C2LPtsMapped.z.idx], ...
                    ws_cube, resolution, obj.bit(38).Opacity, obj.bit(38).Color);
                
                    obj.plotBresenhamLine([ws_B2RPtsMapped.x.idx, ws_B2RPtsMapped.y.idx, ws_B2RPtsMapped.z.idx], ...
                    [ws_C2RPtsMapped.x.idx, ws_C2RPtsMapped.y.idx, ws_C2RPtsMapped.z.idx], ...
                    ws_cube, resolution, obj.bit(39).Opacity, obj.bit(39).Color);
                
                    obj.plotBresenhamLine([ws_B3LPtsMapped.x.idx, ws_B3LPtsMapped.y.idx, ws_B3LPtsMapped.z.idx], ...
                    [ws_C3LPtsMapped.x.idx, ws_C3LPtsMapped.y.idx, ws_C3LPtsMapped.z.idx], ...
                    ws_cube, resolution, obj.bit(40).Opacity, obj.bit(40).Color);
                
                    obj.plotBresenhamLine([ws_B3RPtsMapped.x.idx, ws_B3RPtsMapped.y.idx, ws_B3RPtsMapped.z.idx], ...
                    [ws_C3RPtsMapped.x.idx, ws_C3RPtsMapped.y.idx, ws_C3RPtsMapped.z.idx], ...
                    ws_cube, resolution, obj.bit(41).Opacity, obj.bit(41).Color);
                
                    obj.plotBresenhamLine([ws_B1LPtsMapped.x.idx, ws_B1LPtsMapped.y.idx, ws_B1LPtsMapped.z.idx], ...
                    [ws_B1RPtsMapped.x.idx, ws_B1RPtsMapped.y.idx, ws_B1RPtsMapped.z.idx], ...
                    ws_cube, resolution, obj.bit(42).Opacity, obj.bit(42).Color);
                
                    obj.plotBresenhamLine([ws_B2LPtsMapped.x.idx, ws_B2LPtsMapped.y.idx, ws_B2LPtsMapped.z.idx], ...
                    [ws_B2RPtsMapped.x.idx, ws_B2RPtsMapped.y.idx, ws_B2RPtsMapped.z.idx], ...
                    ws_cube, resolution, obj.bit(43).Opacity, obj.bit(43).Color);
                
                    obj.plotBresenhamLine([ws_B3LPtsMapped.x.idx, ws_B3LPtsMapped.y.idx, ws_B3LPtsMapped.z.idx], ...
                    [ws_B3RPtsMapped.x.idx, ws_B3RPtsMapped.y.idx, ws_B3RPtsMapped.z.idx], ...
                    ws_cube, resolution, obj.bit(44).Opacity, obj.bit(44).Color);
                
                    obj.plotBresenhamLine([ws_C1LPtsMapped.x.idx, ws_C1LPtsMapped.y.idx, ws_C1LPtsMapped.z.idx], ...
                    [ws_C1RPtsMapped.x.idx, ws_C1RPtsMapped.y.idx, ws_C1RPtsMapped.z.idx], ...
                    ws_cube, resolution, obj.bit(45).Opacity, obj.bit(45).Color);
                
                    obj.plotBresenhamLine([ws_C2LPtsMapped.x.idx, ws_C2LPtsMapped.y.idx, ws_C2LPtsMapped.z.idx], ...
                    [ws_C2RPtsMapped.x.idx, ws_C2RPtsMapped.y.idx, ws_C2RPtsMapped.z.idx], ...
                    ws_cube, resolution, obj.bit(46).Opacity, obj.bit(46).Color);
                
                    obj.plotBresenhamLine([ws_C3LPtsMapped.x.idx, ws_C3LPtsMapped.y.idx, ws_C3LPtsMapped.z.idx], ...
                    [ws_C3RPtsMapped.x.idx, ws_C3RPtsMapped.y.idx, ws_C3RPtsMapped.z.idx], ...
                    ws_cube, resolution, obj.bit(47).Opacity, obj.bit(47).Color);
                end
            else
                
                PtsMapped.B1 = ws_B1PtsMapped;
                PtsMapped.B2 = ws_B2PtsMapped;
                PtsMapped.B3 = ws_B3PtsMapped;
                
                if any(strcmpi(varargin, "plotCubes"))
                    plotcube(resolution*[1 1 1], obj.B1', obj.bit(5).Opacity, obj.bit(5).Color);
                    plotcube(resolution*[1 1 1], obj.B2', obj.bit(6).Opacity, obj.bit(6).Color);
                    plotcube(resolution*[1 1 1], obj.B3', obj.bit(7).Opacity, obj.bit(7).Color);

                    plotcube(resolution*[1 1 1], obj.C1', obj.bit(8).Opacity, obj.bit(8).Color);
                    plotcube(resolution*[1 1 1], obj.C2', obj.bit(9).Opacity, obj.bit(9).Color);
                    plotcube(resolution*[1 1 1], obj.C3', obj.bit(10).Opacity, obj.bit(10).Color);
                
                    obj.plotBresenhamLine([ws_B1PtsMapped.x.idx, ws_B1PtsMapped.y.idx, ws_B1PtsMapped.z.idx], ...
                    [ws_C1PtsMapped.x.idx, ws_C1PtsMapped.y.idx, ws_C1PtsMapped.z.idx], ...
                    ws_cube, resolution, obj.bit(18).Opacity, obj.bit(18).Color);

                    obj.plotBresenhamLine([ws_B2PtsMapped.x.idx, ws_B2PtsMapped.y.idx, ws_B2PtsMapped.z.idx], ...
                    [ws_C2PtsMapped.x.idx, ws_C2PtsMapped.y.idx, ws_C2PtsMapped.z.idx], ...
                    ws_cube, resolution, obj.bit(19).Opacity, obj.bit(19).Color);

                    obj.plotBresenhamLine([ws_B3PtsMapped.x.idx, ws_B3PtsMapped.y.idx, ws_B3PtsMapped.z.idx], ...
                    [ws_C3PtsMapped.x.idx, ws_C3PtsMapped.y.idx, ws_C3PtsMapped.z.idx], ...
                    ws_cube, resolution, obj.bit(20).Opacity, obj.bit(20).Color);
                end
            end
            
            
        end

        function [ws_cube, PtsMapped, anglesMapped, coordsMapped] = MapWorkspace(obj, resolution, thetas)

            ws_cube = obj.getCubeAboutRobot(resolution, 'thetas', thetas);
            
            i = 1;
            while i <= size(thetas, 1)
                try
                    
                    TCP(i,:) = obj.directGeometry(thetas(i,:));
                    
                    obj.Jacobian();
    
                    [ws_cube.manipulability.mu1(i), ws_cube.manipulability.mu2(i), ws_cube.manipulability.mu3(i)] = obj.getManipulabilityMeasures();


                    % find the TCP where it belongs to:
                        [~, closestIndex_x_TCP] = min(abs(TCP(i,1) - ws_cube.grid.xAxis.points));
                        [~, closestIndex_y_TCP] = min(abs(TCP(i,2) - ws_cube.grid.yAxis.points));
                        [~, closestIndex_z_TCP] = min(abs(TCP(i,3) - ws_cube.grid.zAxis.points));
                        % mark the cube appropriately
                        ws_cube.grid.idxMap(closestIndex_x_TCP, closestIndex_y_TCP, closestIndex_z_TCP) = bitset(ws_cube.grid.idxMap(closestIndex_x_TCP, closestIndex_y_TCP, closestIndex_z_TCP), obj.bit(1).Position, 1);
                        ws_cube.grid.mappedIndeces(i,:) = [closestIndex_x_TCP, closestIndex_y_TCP, closestIndex_z_TCP];
                
                    if i == 1
                        % find point A1
                        [~, closestIndex_x_A1] = min(abs(obj.A1(1) - ws_cube.grid.xAxis.points));
                        [~, closestIndex_y_A1] = min(abs(obj.A1(2) - ws_cube.grid.yAxis.points));
                        [~, closestIndex_z_A1] = min(abs(obj.A1(3) - ws_cube.grid.zAxis.points));
                        % mark the cube appropriately
                        ws_cube.grid.idxMap(closestIndex_x_A1, closestIndex_y_A1, closestIndex_z_A1) = bitset(ws_cube.grid.idxMap(closestIndex_x_A1, closestIndex_y_A1, closestIndex_z_A1), obj.bit(2).Position, 1);


                        % find point A2
                        [~, closestIndex_x_A2] = min(abs(obj.A2(1) - ws_cube.grid.xAxis.points));
                        [~, closestIndex_y_A2] = min(abs(obj.A2(2) - ws_cube.grid.yAxis.points));
                        [~, closestIndex_z_A2] = min(abs(obj.A2(3) - ws_cube.grid.zAxis.points));
                        % mark the cube appropriately
                        ws_cube.grid.idxMap(closestIndex_x_A2, closestIndex_y_A2, closestIndex_z_A2) = bitset(ws_cube.grid.idxMap(closestIndex_x_A2, closestIndex_y_A2, closestIndex_z_A2), obj.bit(3).Position, 1);


                        % find point A3
                        [~, closestIndex_x_A3] = min(abs(obj.A3(1) - ws_cube.grid.xAxis.points));
                        [~, closestIndex_y_A3] = min(abs(obj.A3(2) - ws_cube.grid.yAxis.points));
                        [~, closestIndex_z_A3] = min(abs(obj.A3(3) - ws_cube.grid.zAxis.points));
                        % mark the cube appropriately
                        ws_cube.grid.idxMap(closestIndex_x_A3, closestIndex_y_A3, closestIndex_z_A3) = bitset(ws_cube.grid.idxMap(closestIndex_x_A3, closestIndex_y_A3, closestIndex_z_A3), obj.bit(4).Position, 1);


                        % find point o
                        [~, closestIndex_x_O] = min(abs(obj.O(1) - ws_cube.grid.xAxis.points));
                        [~, closestIndex_y_O] = min(abs(obj.O(2) - ws_cube.grid.yAxis.points));
                        [~, closestIndex_z_O] = min(abs(obj.O(3) - ws_cube.grid.zAxis.points));
                        % mark the cube appropriately
                        ws_cube.grid.idxMap(closestIndex_x_O, closestIndex_y_O, closestIndex_z_O) = bitset(ws_cube.grid.idxMap(closestIndex_x_O, closestIndex_y_O, closestIndex_z_O), obj.bit(11).Position, 1);


%                         precision = ceil(abs(log10(resolution)));
                        % find lines A1A2
                        pts_A1 = obj.mapLine([closestIndex_x_A1; closestIndex_y_A1; closestIndex_z_A1], [closestIndex_x_A2; closestIndex_y_A2; closestIndex_z_A2]);
                        for j = 1:length(pts_A1.x)
                            ws_cube.grid.idxMap(pts_A1.x(j), pts_A1.y(j), pts_A1.z(j)) = bitset(ws_cube.grid.idxMap(pts_A1.x(j), pts_A1.y(j), pts_A1.z(j)), obj.bit(12).Position, 1);
                        end

                        % find lines A2A3
                        pts_A2A3 = obj.mapLine([closestIndex_x_A2; closestIndex_y_A2; closestIndex_z_A2], [closestIndex_x_A3; closestIndex_y_A3; closestIndex_z_A3]);
                        for j = 1:size(pts_A2A3,1)
                            ws_cube.grid.idxMap(pts_A2A3.x(j), pts_A2A3.y(j), pts_A2A3.z(j)) = bitset(ws_cube.grid.idxMap(pts_A2A3.x(j), pts_A2A3.y(j), pts_A2A3.z(j)), obj.bit(13).Position, 1);
                        end

                        % find lines A1A3
                        pts_A1A3 = obj.mapLine([closestIndex_x_A1; closestIndex_y_A1; closestIndex_z_A1], [closestIndex_x_A3; closestIndex_y_A3; closestIndex_z_A3]);
                        for j = 1:size(pts_A1A3,1)
                            ws_cube.grid.idxMap(pts_A1A3.x(j), pts_A1A3.y(j), pts_A1A3.z(j)) = bitset(ws_cube.grid.idxMap(pts_A1A3.x(j), pts_A1A3.y(j), pts_A1A3.z(j)), obj.bit(14).Position, 1);
                        end
                    end

                    % find point B1
                    [~, closestIndex_x_B1] = min(abs(obj.B1(1) - ws_cube.grid.xAxis.points));
                    [~, closestIndex_y_B1] = min(abs(obj.B1(2) - ws_cube.grid.yAxis.points));
                    [~, closestIndex_z_B1] = min(abs(obj.B1(3) - ws_cube.grid.zAxis.points));
                    % mark the cube appropriately
                    ws_cube.grid.idxMap(closestIndex_x_B1, closestIndex_y_B1, closestIndex_z_B1) = bitset(ws_cube.grid.idxMap(closestIndex_x_B1, closestIndex_y_B1, closestIndex_z_B1), obj.bit(5).Position, 1);


                    % find point B2
                    [~, closestIndex_x_B2] = min(abs(obj.B2(1) - ws_cube.grid.xAxis.points));
                    [~, closestIndex_y_B2] = min(abs(obj.B2(2) - ws_cube.grid.yAxis.points));
                    [~, closestIndex_z_B2] = min(abs(obj.B2(3) - ws_cube.grid.zAxis.points));
                    % mark the cube appropriately
                    ws_cube.grid.idxMap(closestIndex_x_B2, closestIndex_y_B2, closestIndex_z_B2) = bitset(ws_cube.grid.idxMap(closestIndex_x_B2, closestIndex_y_B2, closestIndex_z_B2), obj.bit(6).Position, 1);


                    % find point B3
                    [~, closestIndex_x_B3] = min(abs(obj.B3(1) - ws_cube.grid.xAxis.points));
                    [~, closestIndex_y_B3] = min(abs(obj.B3(2) - ws_cube.grid.yAxis.points));
                    [~, closestIndex_z_B3] = min(abs(obj.B3(3) - ws_cube.grid.zAxis.points));
                    % mark the cube appropriately
                    ws_cube.grid.idxMap(closestIndex_x_B3, closestIndex_y_B3, closestIndex_z_B3) = bitset(ws_cube.grid.idxMap(closestIndex_x_B3, closestIndex_y_B3, closestIndex_z_B3), obj.bit(7).Position, 1);


                    % find point C1
                    [~, closestIndex_x_C1] = min(abs(obj.C1(1) - ws_cube.grid.xAxis.points));
                    [~, closestIndex_y_C1] = min(abs(obj.C1(2) - ws_cube.grid.yAxis.points));
                    [~, closestIndex_z_C1] = min(abs(obj.C1(3) - ws_cube.grid.zAxis.points));
                    % mark the cube appropriately
                    ws_cube.grid.idxMap(closestIndex_x_C1, closestIndex_y_C1, closestIndex_z_C1) = bitset(ws_cube.grid.idxMap(closestIndex_x_C1, closestIndex_y_C1, closestIndex_z_C1), obj.bit(8).Position, 1);

                    
                    % find point C2
                    [~, closestIndex_x_C2] = min(abs(obj.C2(1) - ws_cube.grid.xAxis.points));
                    [~, closestIndex_y_C2] = min(abs(obj.C2(2) - ws_cube.grid.yAxis.points));
                    [~, closestIndex_z_C2] = min(abs(obj.C2(3) - ws_cube.grid.zAxis.points));
                    % mark the cube appropriately
                    ws_cube.grid.idxMap(closestIndex_x_C2, closestIndex_y_C2, closestIndex_z_C2) = bitset(ws_cube.grid.idxMap(closestIndex_x_C2, closestIndex_y_C2, closestIndex_z_C2), obj.bit(9).Position, 1);


                    % find point C3
                    [~, closestIndex_x_C3] = min(abs(obj.C3(1) - ws_cube.grid.xAxis.points));
                    [~, closestIndex_y_C3] = min(abs(obj.C3(2) - ws_cube.grid.yAxis.points));
                    [~, closestIndex_z_C3] = min(abs(obj.C3(3) - ws_cube.grid.zAxis.points));
                    % mark the cube appropriately
                    ws_cube.grid.idxMap(closestIndex_x_C3, closestIndex_y_C3, closestIndex_z_C3) = bitset(ws_cube.grid.idxMap(closestIndex_x_C3, closestIndex_y_C3, closestIndex_z_C3), obj.bit(10).Position, 1);


                    % find lines A1B1
                    pts_A1B1 = obj.mapLine([closestIndex_x_A1; closestIndex_y_A1; closestIndex_z_A1], [closestIndex_x_B1; closestIndex_y_B1; closestIndex_z_B1]);
                    for j = 1:size(pts_A1B1,1)
                        ws_cube.grid.idxMap(pts_A1B1.x(j), pts_A1B1.y(j), pts_A1B1.z(j)) = bitset(ws_cube.grid.idxMap(pts_A1B1.x(j), pts_A1B1.y(j), pts_A1B1.z(j)), obj.bit(15).Position, 1);
                    end

                    % find lines A2B2
                    pts_A2B2 = obj.mapLine([closestIndex_x_A2; closestIndex_y_A2; closestIndex_z_A2], [closestIndex_x_B2; closestIndex_y_B2; closestIndex_z_B2]);
                    for j = 1:size(pts_A2B2,1)
                        ws_cube.grid.idxMap(pts_A2B2.x(j), pts_A2B2.y(j), pts_A2B2.z(j)) = bitset(ws_cube.grid.idxMap(pts_A2B2.x(j), pts_A2B2.y(j), pts_A2B2.z(j)), obj.bit(16).Position, 1);
                    end

                    % find lines A3B3
                    pts_A3B3 = obj.mapLine([closestIndex_x_A3; closestIndex_y_A3; closestIndex_z_A3], [closestIndex_x_B3; closestIndex_y_B3; closestIndex_z_B3]);
                    for j = 1:size(pts_A3B3,1)
                        ws_cube.grid.idxMap(pts_A3B3.x(j), pts_A3B3.y(j), pts_A3B3.z(j)) = bitset(ws_cube.grid.idxMap(pts_A3B3.x(j), pts_A3B3.y(j), pts_A3B3.z(j)), obj.bit(17).Position, 1);
                    end

                    % find lines B1C1
                    pts_B1C1 = obj.mapLine([closestIndex_x_B1; closestIndex_y_B1; closestIndex_z_B1], [closestIndex_x_C1; closestIndex_y_C1; closestIndex_z_C1]);
                    for j = 1:size(pts_B1C1,1)
                        ws_cube.grid.idxMap(pts_B1C1.x(j), pts_B1C1.y(j), pts_B1C1.z(j)) = bitset(ws_cube.grid.idxMap(pts_B1C1.x(j), pts_B1C1.y(j), pts_B1C1.z(j)), obj.bit(18).Position, 1);
                    end

                    % find lines B2C2
                    pts_B2C2 = obj.mapLine([closestIndex_x_B2; closestIndex_y_B2; closestIndex_z_B2], [closestIndex_x_C2; closestIndex_y_C2; closestIndex_z_C2]);
                    for j = 1:size(pts_B2C2,1)
                        ws_cube.grid.idxMap(pts_B2C2.x(j), pts_B2C2.y(j), pts_B2C2.z(j)) = bitset(ws_cube.grid.idxMap(pts_B2C2.x(j), pts_B2C2.y(j), pts_B2C2.z(j)), obj.bit(19).Position, 1);
                    end

                    % find lines B3C3
                    pts_B3C3 = obj.mapLine([closestIndex_x_B3; closestIndex_y_B3; closestIndex_z_B3], [closestIndex_x_C3; closestIndex_y_C3; closestIndex_z_C3]);
                    for j = 1:size(pts_B3C3,1)
                        ws_cube.grid.idxMap(pts_B3C3.x(j), pts_B3C3.y(j), pts_B3C3.z(j)) = bitset(ws_cube.grid.idxMap(pts_B3C3.x(j), pts_B3C3.y(j), pts_B3C3.z(j)), obj.bit(20).Position, 1);
                    end

                    % find lines C1C2
                    pts_C1C2 = obj.mapLine([closestIndex_x_C1; closestIndex_y_C1; closestIndex_z_C1], [closestIndex_x_C2; closestIndex_y_C2; closestIndex_z_C2]);
                    for j = 1:size(pts_C1C2,1)
                        ws_cube.grid.idxMap(pts_C1C2.x(j), pts_C1C2.y(j), pts_C1C2.z(j)) = bitset(ws_cube.grid.idxMap(pts_C1C2.x(j), pts_C1C2.y(j), pts_C1C2.z(j)), obj.bit(21).Position, 1);
                    end

                    % find lines C2C3
                    pts_C2C3 = obj.mapLine([closestIndex_x_C2; closestIndex_y_C2; closestIndex_z_C2], [closestIndex_x_C3; closestIndex_y_C3; closestIndex_z_C3]);
                    for j = 1:size(pts_C2C3,1)
                        ws_cube.grid.idxMap(pts_C2C3.x(j), pts_C2C3.y(j), pts_C2C3.z(j)) = bitset(ws_cube.grid.idxMap(pts_C2C3.x(j), pts_C2C3.y(j), pts_C2C3.z(j)), obj.bit(22).Position, 1);
                    end

                    % find lines C1C3
                    pts_C1C3 = obj.mapLine([closestIndex_x_C1; closestIndex_y_C1; closestIndex_z_C1], [closestIndex_x_C3; closestIndex_y_C3; closestIndex_z_C3]);
                    for j = 1:size(pts_C1C3,1)
                        ws_cube.grid.idxMap(pts_C1C3.x(j), pts_C1C3.y(j), pts_C1C3.z(j)) = bitset(ws_cube.grid.idxMap(pts_C1C3.x(j), pts_C1C3.y(j), pts_C1C3.z(j)), obj.bit(23).Position, 1);
                    end

%                     ws_cubePlot(i) = gtpr.ws_cubePlotByPiece([closestIndex_x, closestIndex_y, closestIndex_z], ws_resolution, ws_cube);
                
                catch ME
                    if strcmpi(ME.message, 'Impossible Configuration present!') || ...
                            (strcmpi(ME.message, 'Unable to perform assignment because the size of the left side is 1-by-3 and the size of the right side is 0-by-0.') || strcmpi(ME.message, 'Unable to perform assignment because the indices on the left side are not compatible with the size of the right side.') && ME.stack.line == 2766)
                        thetas(i,:) = rand(1,3);
                        continue;
                    else
                        error("error!!!")
                    end
                end

                pts(:,i) = [closestIndex_x_TCP; closestIndex_y_TCP; closestIndex_z_TCP];
                i = i + 1;

                    
            end


            PtsMapped = pts;
            coordsMapped = TCP';
            anglesMapped = thetas;

        end



%         function [ws_cube, PtsMapped] = MapWorkspace(obj, resolution, varargin)
%            
%             switch length(varargin)
%                 case 0
%                     arg1 = [];
%                     arg2 = [];
%                     arg3 = [];
%                 case 1
%                     arg1 = varargin{1};
%                     arg2 = [];
%                     arg3 = [];
%                 case 2
%                     arg1 = varargin{1};
%                     arg2 = varargin{2};
%                     arg3 = [];
%                 case 3
%                     arg1 = varargin{1};
%                     arg2 = varargin{2};
%                     arg3 = varargin{3};
%                 otherwise
%                     error("Case not implemented yet, implement it.");
%             end
%             
%             if any(strcmpi("thetas", varargin))
%                 idx = find(strcmpi("thetas", varargin));
%                 thetas = varargin{idx+1};
%             elseif any(strcmpi("NoOfPts", varargin))
%                 idx = find(strcmpi("NoOfPts", varargin));
%                 thetas = rand(varargin{idx+1},3)*pi-pi/2;
%             else
%                 thetas = rand(1000,3)*pi-pi/2;
%             end
%             
%             ws_cube = obj.getCubeAboutRobot(resolution, arg1, arg2, arg3);
%             
%             
%             thetas = ws_cube.theta;
%             
%             
%             
%             for i = size(thetas,1):-1:1
%                 obj.directGeometry(thetas(i,:));
%                 [ws_cube, PtsMapped(i)] = obj.plotWS(ws_cube, thetas(i,:), resolution, arg1, arg2, arg3);
%                 
%                 ws_cube.W.points(i,:) = obj.W;
%                 ws_cube.W.theta(i,:) = thetas(i,:);
%                 
%                 [ws_cube.W.manipulability.mu1(i), ws_cube.W.manipulability.mu2(i), ws_cube.W.manipulability.mu3(i)] = obj.getManipulabilityMeasures();
%                 [ws_cube.W.manipulability.eigVectors{i}, ws_cube.W.manipulability.eigValues{i}] = obj.getVelocityEigens();
%                 ws_cube.W.manipulability.Jacobian{i} = obj.Jacobian();
%             end
%             ws_cube.resolution = resolution;
%             ws_cube.cubicResolution = resolution^3;
%         end
        


        function plotWSCube(obj, ws_cube, varargin)
            resolution = ws_cube.resolution;
%             resolution = resolution/2;
            hold on;
            grid on;
           
            if any(strcmpi(varargin, "skipIds"))
                idx = find(strcmpi(varargin, 'skipIds'));
                skipIdxs = varargin{idx+1};
            else
                skipIdxs = [];
            end
           
            idxMap = ws_cube.grid.idxMap;
           
            for i = 1:length(obj.bit)
               
                if any(skipIdxs == i)
                    continue;
                end
               
                currBitMap = bitget(idxMap, obj.bit(i).Position);
               
                idxs = find(currBitMap==1);
               
                [idx_x, idx_y, idx_z] = ind2sub(size(currBitMap), idxs);
               
                for j = 1:length(idx_x)
                    
%                    plotcube(resolution*[1 1 1], obj.B1', obj.bit(5).Opacity, obj.bit(5).Color);myPlotCube(~, pts, edgeResolution, varargin)
                    obj.myPlotCube(...
                        [ws_cube.grid.xAxis.points(idx_x(j)), ws_cube.grid.yAxis.points(idx_y(j)), ws_cube.grid.zAxis.points(idx_z(j))]', ...
                        resolution*[1 1 1]/2, 'edgePlot');%, ...
%                         obj.bit(i).Opacity, ...
%                         obj.bit(i).Color);    
                end
            end
        end

        function [wci, data] = workspaceCoverageIndex(obj, ws_cube, varargin)
    
            [ws_boundary, ws_volume] = boundary(ws_cube.W.points(:,1), ws_cube.W.points(:,2), ws_cube.W.points(:,3));

            currentWsVol = length(ws_cube.W.points)*ws_cube.cubicResolution;
            
            wci = (currentWsVol/ws_volume)*100;

        %     ws_boundary_coords = [ws_cube.W.points(ws_boundary(:,1),1), ws_cube.W.points(ws_boundary(:,2),2), ws_cube.W.points(ws_boundary(:,3),3)];

        %     [ws_convHull, ws_convHull_vol] = convhull(ws_cube.W.points(:,1), ws_cube.W.points(:,2), ws_cube.W.points(:,3), 'Simplify', true);

            hold on;
            if any(strcmpi("plotBoundary", varargin))
                trisurf(ws_boundary, ws_cube.W.points(:,1), ws_cube.W.points(:,2), ws_cube.W.points(:,3), 'FaceColor', obj.bit(1).Color, 'FaceAlpha', obj.bit(1).Opacity);
            end

        %     if any(strcmpi("convhull", varargin))
        %         trisurf(ws_convHull, ws_cube.W.points(:,1), ws_cube.W.points(:,2), ws_cube.W.points(:,3), 'FaceColor', 'green', 'FaceAlpha', gtpr.bit(1).Opacity/2);   
        %     end
            data.ws_boundaryIdx = ws_boundary;
            data.full_ws_volume = ws_volume;
            data.current_ws_volume = currentWsVol;
        end


        
        function plotBresenhamLine(obj, Pt1_idxs, Pt2_idxs, cube, resolution, Opacity, Color)
            idxs = obj.mapLine(Pt1_idxs, Pt2_idxs);

            coords = [ cube.grid.xAxis.points(idxs.x);...
                       cube.grid.yAxis.points(idxs.y);...
                       cube.grid.zAxis.points(idxs.z)]';

            for j = 1:size(coords,1)
               plotcube(resolution*[1 1 1], coords(j,:), Opacity, Color);
            end
        end

        function idxs = mapLine(~, P1, P2)
            [idxs.x, idxs.y, idxs.z] = bresenham_line3d(P1,P2);
        end
        
        function cube = createCubeIdxMap(obj, varargin)
            
            if any(strcmpi(varargin, 'ws_cube'))
                idx = find(strcmpi(varargin, 'ws_cube'));
                cube = varargin{idx+1};
            else
                cube = zeros(10,10,10);
            end
            
            cube.grid.idxMap = zeros(length(cube.grid.xAxis.points), length(cube.grid.yAxis.points), length(cube.grid.zAxis.points));
            
        end

        function cube = CubeAboutPoints(obj, pts, edgeResolution, varargin)
            %pts should be the size of n x 3
            %varargin should include the plot style related stuff

            if any(strcmpi(varargin, 'scale'))
                idx = find(strcmpi(varargin, 'scale'));
                scale = varargin{idx+1};
            else
                scale = 1;
            end
            
            [min_x, max_x] = bounds(pts(:,1));
            [min_y, max_y] = bounds(pts(:,2));
            [min_z, max_z] = bounds(pts(:,3));

            min_x = min_x * scale;
            max_x = max_x * scale;
            min_y = min_y * scale;
            max_y = max_y * scale;
            min_z = min_z * scale;
            max_z = max_z * scale;

            pt(1,:) = [min_x, min_y, min_z];
            pt(2,:) = [min_x, min_y, max_z];
            pt(3,:) = [min_x, max_y, min_z];
            pt(4,:) = [min_x, max_y, max_z];
            pt(5,:) = [max_x, min_y, min_z];
            pt(6,:) = [max_x, min_y, max_z];
            pt(7,:) = [max_x, max_y, min_z];
            pt(8,:) = [max_x, max_y, max_z];

            switch length(varargin)
                case 0
                    arg1 = [];
                    arg2 = [];
                    arg3 = [];
                    arg4 = [];
                case 1
                    arg1 = varargin{1};
                    arg2 = [];
                    arg3 = [];
                    arg4 = [];
                case 2
                    arg1 = varargin{1};
                    arg2 = varargin{2};
                    arg3 = [];
                    arg4 = [];
                case 3
                    arg1 = varargin{1};
                    arg2 = varargin{2};
                    arg3 = varargin{3};
                    arg4 = [];
                case 4
                    arg1 = varargin{1};
                    arg2 = varargin{2};
                    arg3 = varargin{3};
                    arg4 = varargin{4};
                    arg5 = [];
                case 5
                    arg1 = varargin{1};
                    arg2 = varargin{2};
                    arg3 = varargin{3};
                    arg4 = varargin{4};
                    arg5 = varargin{5};
                otherwise
                    error("Edit this part, to allow more argumets!");
            end
            
            if any(strcmpi(varargin, 'ws_cube'))
                idx = find(strcmpi(varargin, 'ws_cube'));
                cube = varargin{idx+1};
                
            else
                cube = obj.myPlotCube(pts, edgeResolution, arg1, arg2, arg3, arg4, arg5);
            end
            cube.pts = pt;

            cube.xAxis.min = min_x;
            cube.xAxis.max = max_x;

            cube.yAxis.min = min_y;
            cube.yAxis.max = max_y;

            cube.zAxis.min = min_z;
            cube.zAxis.max = max_z;
            
            if any(strcmpi(varargin, 'ws_cube'))
                idx = find(strcmpi(varargin, 'ws_cube'));
                ws_cube = varargin{idx+1};
                if exist('cube.grid.idxMap', 'var') && exist('ws_cube.grid.idxMap', 'var')
                    cube.grid.idxMap = cube.grid.idxMap | ws_cube.grid.idxMap;
                end
                
            end
        end

        function [cubeplotH, points] = ws_cubePlotByPiece(obj, idx, res, ws_cube, varargin)
%             pts1 = pt + [res/2, res/2, res/2];
%             pts2 = pt + [res/2, res/2, -res/2];
%             pts3 = pt + [res/2, -res/2, res/2];
%             pts4 = pt + [res/2, -res/2, -res/2];
%             pts5 = pt + [-res/2, res/2, res/2];
%             pts6 = pt + [-res/2, res/2, -res/2];
%             pts7 = pt + [-res/2, -res/2, res/2];
%             pts8 = pt + [-res/2, -res/2, -res/2];

%             pt = [ws_cube.grid.xAxis(idx(1)), ws_cube.grid.yAxis(idx(2)), ws_cube.grid.zAxis(idx(3))];
            idx_x_p1 = idx(1) + 1;
            idx_y_p1 = idx(2) + 1;
            idx_z_p1 = idx(3) + 1;

            idx_x_m1 = idx(1) - 1;
            idx_y_m1 = idx(2) - 1;
            idx_z_m1 = idx(3) - 1;
            
            if idx_x_p1 > length(ws_cube.grid.xAxis.points)
                idx_x_p1 = length(ws_cube.grid.xAxis.points);
            else
                if idx_x_m1 < 1
                    idx_x_m1 = 1;
                end
            end

            if idx_y_p1 > length(ws_cube.grid.yAxis.points)
                idx_y_p1 = length(ws_cube.grid.yAxis.points);
            else
                if idx_y_m1 < 1
                    idx_y_m1 = 1;
                end
            end


            if idx_z_p1 > length(ws_cube.grid.zAxis.points)
                idx_z_p1 = length(ws_cube.grid.zAxis.points);
            else
                if idx_z_m1 < 1
                    idx_z_m1 = 1;
                end
            end

            pts1 = [ws_cube.grid.xAxis.points(idx_x_p1)-res/2, ws_cube.grid.yAxis.points(idx_y_p1)-res/2, ws_cube.grid.zAxis.points(idx_z_p1)-res/2];
            pts2 = [ws_cube.grid.xAxis.points(idx_x_p1)-res/2, ws_cube.grid.yAxis.points(idx_y_p1)-res/2, ws_cube.grid.zAxis.points(idx_z_m1)+res/2];
            pts3 = [ws_cube.grid.xAxis.points(idx_x_p1)-res/2, ws_cube.grid.yAxis.points(idx_y_m1)+res/2, ws_cube.grid.zAxis.points(idx_z_p1)-res/2];
            pts4 = [ws_cube.grid.xAxis.points(idx_x_p1)-res/2, ws_cube.grid.yAxis.points(idx_y_m1)+res/2, ws_cube.grid.zAxis.points(idx_z_m1)+res/2];
            pts5 = [ws_cube.grid.xAxis.points(idx_x_m1)+res/2, ws_cube.grid.yAxis.points(idx_y_p1)-res/2, ws_cube.grid.zAxis.points(idx_z_p1)-res/2];
            pts6 = [ws_cube.grid.xAxis.points(idx_x_m1)+res/2, ws_cube.grid.yAxis.points(idx_y_p1)-res/2, ws_cube.grid.zAxis.points(idx_z_m1)+res/2];
            pts7 = [ws_cube.grid.xAxis.points(idx_x_m1)+res/2, ws_cube.grid.yAxis.points(idx_y_m1)+res/2, ws_cube.grid.zAxis.points(idx_z_p1)-res/2];
            pts8 = [ws_cube.grid.xAxis.points(idx_x_m1)+res/2, ws_cube.grid.yAxis.points(idx_y_m1)+res/2, ws_cube.grid.zAxis.points(idx_z_m1)+res/2];

%             pts1 = pt + [res/4, res/4, res/4];
%             pts2 = pt + [res/4, res/4, -res/4];
%             pts3 = pt + [res/4, -res/4, res/4];
%             pts4 = pt + [res/4, -res/4, -res/4];
%             pts5 = pt + [-res/4, res/4, res/4];
%             pts6 = pt + [-res/4, res/4, -res/4];
%             pts7 = pt + [-res/4, -res/4, res/4];
%             pts8 = pt + [-res/4, -res/4, -res/4];

            pts = [pts1; pts2; pts3; pts4; pts5; pts6; pts7; pts8];

            if ~any(strcmp(varargin, 'noPlot'))
                switch length(varargin)
                    case 0
                        cubeplotH = obj.myPlotCube(pts, res, 'edgePlot');
                    case 1
                        cubeplotH = obj.myPlotCube(pts, res, 'edgePlot', varargin{1});
                    case 2
                        cubeplotH = obj.myPlotCube(pts, res, 'edgePlot', varargin{1}, varargin{2});
                    case 3
                    cubeplotH = obj.myPlotCube(pts, res, 'edgePlot', varargin{1}, varargin{2}, varargin{3});
                    otherwise
                        error("Incorrect amount of arguments given.");
                end
            else
                cubeplotH = double.empty();
            end
            points = pts;
        end

        function cube = myPlotCube(~, pts, edgeResolution, varargin)
            
            edgeColor = 'r';

            pt = 0.5;
            
            if any(strcmpi('edgeColor', varargin))
                idx = find(strcmpi('edgeColor', varargin));
                edgeColor = varargin{idx+1};
            end

            if all(size(edgeResolution) == [1 1]) 
                edgeResolution = [1 1 1] * edgeResolution;
            end

            if all(size(pts) == [3 1]) || all(size(pts) == [1 3]) 
                
                pt1 =   pts + [edgeResolution(1); edgeResolution(2); edgeResolution(3)];
                pt2 =   pts + [edgeResolution(1); edgeResolution(2); -edgeResolution(3)];
                pt3 =   pts + [edgeResolution(1); -edgeResolution(2); edgeResolution(3)];
                pt4 =   pts + [edgeResolution(1); -edgeResolution(2); -edgeResolution(3)];
    
                pt5 =   pts + [-edgeResolution(1); edgeResolution(2); edgeResolution(3)];
                pt6 =   pts + [-edgeResolution(1); edgeResolution(2); -edgeResolution(3)];
                pt7 =   pts + [-edgeResolution(1); -edgeResolution(2); edgeResolution(3)];
                pt8 =   pts + [-edgeResolution(1); -edgeResolution(2); -edgeResolution(3)];
    
                pts = [pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8]';
            end

            if any(strcmpi('linewidth', varargin))
                idx = find(strcmpi('linewidth', varargin));
                pt = varargin{idx+1};
            end

            cube_pts = gobjects(size(pts,1),1);
            
           if any(strcmpi(varargin, 'PointsPlot'))
                for i = 1:size(pts,1)
                    cube_pts(i) = plot3(pts(i,1), pts(i,2), pts(i,3), 'Marker', 'x', 'Color', 'r');
                end
                cube.corners = cube_pts;
            end
            
            if any(strcmpi(varargin, 'edgePlot'))
                
                min_x = min(pts(:,1));
                max_x = max(pts(:,1));
                
                min_y = min(pts(:,2));
                max_y = max(pts(:,2));
                
                min_z = min(pts(:,3));
                max_z = max(pts(:,3));
                
                edge1 = plot3([min_x min_x], [min_y max_y], [min_z min_z], edgeColor, 'LineWidth', pt);
                edge2 = plot3([min_x max_x], [min_y min_y], [min_z min_z], edgeColor, 'LineWidth', pt);
                edge3 = plot3([max_x max_x], [min_y max_y], [min_z min_z], edgeColor, 'LineWidth', pt);
                edge4 = plot3([min_x max_x], [max_y max_y], [min_z min_z], edgeColor, 'LineWidth', pt);

                edge5 = plot3([min_x min_x], [min_y max_y], [max_z max_z], edgeColor, 'LineWidth', pt);
                edge6 = plot3([min_x max_x], [min_y min_y], [max_z max_z], edgeColor, 'LineWidth', pt);
                edge7 = plot3([max_x max_x], [min_y max_y], [max_z max_z], edgeColor, 'LineWidth', pt);
                edge8 = plot3([min_x max_x], [max_y max_y], [max_z max_z], edgeColor, 'LineWidth', pt);

                edge9 = plot3([min_x min_x], [min_y min_y], [min_z max_z], edgeColor, 'LineWidth', pt);
                edge10 = plot3([max_x max_x], [min_y min_y], [min_z max_z], edgeColor, 'LineWidth', pt);
                edge11 = plot3([min_x min_x], [max_y max_y], [min_z max_z], edgeColor, 'LineWidth', pt);
                edge12 = plot3([max_x max_x], [max_y max_y], [min_z max_z], edgeColor, 'LineWidth', pt);
                
                cube.edges.edge1 = edge1;
                cube.edges.edge2 = edge2;
                cube.edges.edge3 = edge3;
                cube.edges.edge4 = edge4;
                cube.edges.edge5 = edge5;
                cube.edges.edge6 = edge6;
                cube.edges.edge7 = edge7;
                cube.edges.edge8 = edge8;
                cube.edges.edge9 = edge9;
                cube.edges.edge10 = edge10;
                cube.edges.edge11 = edge11;
                cube.edges.edge12 = edge12;
            end
            
            if any(strcmpi(varargin, 'scale'))
                idx = find(strcmpi(varargin{:}, 'scale'));
                scale = varargin{idx+1};
            else
                scale = 1.1;
            end

            [min_x, max_x] = bounds(pts(:,1));
            [min_y, max_y] = bounds(pts(:,2));
            [min_z, max_z] = bounds(pts(:,3));

            min_x = (min_x - edgeResolution) * scale;
            max_x = (max_x + edgeResolution) * scale;
            min_y = (min_y - edgeResolution) * scale;
            max_y = (max_y + edgeResolution) * scale;
            min_z = (min_z - edgeResolution) * scale;
            max_z = (max_z + edgeResolution) * scale;

            
    
            cube.center = [min_x + ((max_x - min_x)/2); ...
                           min_y + ((max_y - min_y)/2); ...
                           min_z + ((max_z - min_z)/2)];
        end

        function cube = prepareGridFromCube(~, cube, edgeResolution, varargin)

            cube.grid.xAxis.points =    cube.xAxis.min:edgeResolution:cube.xAxis.max;
            cube.grid.yAxis.points = cube.yAxis.min:edgeResolution:cube.yAxis.max;
            cube.grid.zAxis.points = cube.zAxis.min:edgeResolution:cube.zAxis.max;

            cube.grid.xAxis.idx_length = length(cube.grid.xAxis.points);
            cube.grid.yAxis.idx_length = length(cube.grid.yAxis.points);
            cube.grid.zAxis.idx_length = length(cube.grid.zAxis.points);
    
          
            cube.grid.xAxis.length_m = range(cube.grid.xAxis.points);
            cube.grid.yAxis.length_m = range(cube.grid.yAxis.points);
            cube.grid.zAxis.length_m = range(cube.grid.zAxis.points);
        end

        function [cube, gridValuesMapped] = mapBitPointInSpace(~, cube, Coordinate_point, bitPosition, bitName, bitValue)

            [x_minValue, x_closestIndex] = min(abs(cube.grid.xAxis.points - Coordinate_point(1)));
            [y_minValue, y_closestIndex] = min(abs(cube.grid.yAxis.points - Coordinate_point(2)));
            [z_minValue, z_closestIndex] = min(abs(cube.grid.zAxis.points - Coordinate_point(3)));

        %     [x_minValue, x_closestIndex] = min(cube.grid.xAxis.points - Coordinate_point(1));
        %     [y_minValue, y_closestIndex] = min(cube.grid.yAxis.points - Coordinate_point(2));
        %     [z_minValue, z_closestIndex] = min(cube.grid.zAxis.points - Coordinate_point(3));

            cube.grid.idxMap(x_closestIndex, y_closestIndex, z_closestIndex) = bitset(cube.grid.idxMap(x_closestIndex, y_closestIndex, z_closestIndex), bitPosition, bitValue, 'uint64');

            gridValuesMapped.x.coords = cube.grid.xAxis.points(x_closestIndex);
            gridValuesMapped.y.coords = cube.grid.yAxis.points(y_closestIndex);
            gridValuesMapped.z.coords = cube.grid.zAxis.points(z_closestIndex);

            gridValuesMapped.x.idx = x_closestIndex;
            gridValuesMapped.y.idx = y_closestIndex;
            gridValuesMapped.z.idx = z_closestIndex;

            gridValuesMapped.Name = bitName;

        end
        
        function ws_cube = getCubeAboutRobot(obj, resolution, varargin)
            
            if any(strcmpi(varargin, "thetas"))
                idx = find(strcmpi(varargin, "thetas"));
                
                theta_cuboid = varargin{idx+1};
            else
                theta_cuboid = [-pi/4,  -pi/4,  -pi/4;
                            pi/2,   pi/2,   pi/2;
                            pi/2,   0,      0;
                            0,      pi/2,   0;
                            0,      0,      pi/2];
            end
            
            pts_W = zeros(size(theta_cuboid,1), 3);
            
            pts_A1 = zeros(size(theta_cuboid,1), 3);
            pts_A2 = zeros(size(theta_cuboid,1), 3);
            pts_A3 = zeros(size(theta_cuboid,1), 3);
            
            pts_B1 = zeros(size(theta_cuboid,1), 3);
            pts_B2 = zeros(size(theta_cuboid,1), 3);
            pts_B3 = zeros(size(theta_cuboid,1), 3);
            
            pts_C1 = zeros(size(theta_cuboid,1), 3);
            pts_C2 = zeros(size(theta_cuboid,1), 3);
            pts_C3 = zeros(size(theta_cuboid,1), 3);
            
            wrongIdx = [];
            
            for i = 1:size(theta_cuboid,1)
                try
                    pts_W(i,:) = obj.directGeometry(theta_cuboid(i,:));

%                     if pts_W(i,3) > 0
%                         thetas(i,:) = thetas(i-1,:);
%                         pts_W(i,:) = pts_W(i-1,:);
%                         continue
%                     end

                    thetas(i,:) = theta_cuboid(i,:);
                catch ME
                    if contains(ME.message, "Impossible Configuration present!")
                        wrongIdx = [wrongIdx, i];
                        continue
                    end
                end
                
                pts_A1(i,:) = obj.A1;
                pts_A2(i,:) = obj.A2;
                pts_A3(i,:) = obj.A3;

                pts_B1(i,:) = obj.B1;
                pts_B2(i,:) = obj.B2;
                pts_B3(i,:) = obj.B3;

                pts_C1(i,:) = obj.C1;
                pts_C2(i,:) = obj.C2;
                pts_C3(i,:) = obj.C3; 
            end

            pts = [ pts_W; ...
                    pts_A1; pts_A2; pts_A3; ...
                    pts_B1; pts_B2; pts_B3; ...
                    pts_C1; pts_C2; pts_C3];
            
            if any(strcmpi(varargin, "ws_cube"))
                idx = find(strcmpi(varargin, "ws_cube"));
                ws_cube = varargin{idx+1};
                
                pts = [pts; ws_cube.pts]; 
            end
            
            
            switch length(varargin)
                case 0
                    arg1 = [];
                    arg2 = [];
                    arg3 = [];
                    arg4 = [];
                    arg5 = [];
                case 1
                    arg1 = varargin{1};
                    arg2 = [];
                    arg3 = [];
                    arg4 = [];
                    arg5 = [];
                case 2
                    arg1 = varargin{1};
                    arg2 = varargin{2};
                    arg3 = [];
                    arg4 = [];
                    arg5 = [];
                case 3
                    arg1 = varargin{1};
                    arg2 = varargin{2};
                    arg3 = varargin{3};
                    arg4 = [];
                    arg5 = [];
                case 4
                    arg1 = varargin{1};
                    arg2 = varargin{2};
                    arg3 = varargin{3};
                    arg4 = varargin{4};
                    arg5 = [];
                case 5
                    arg1 = varargin{1};
                    arg2 = varargin{2};
                    arg3 = varargin{3};
                    arg4 = varargin{4};
                    arg5 = varargin{5};
                otherwise
                    error("Edit this part, to allow more argumets!");
            end
               
            ws_cube = obj.CubeAboutPoints(pts, resolution, arg1, arg2, arg3, arg4, arg5);
            
            ws_cube = obj.prepareGridFromCube(ws_cube, resolution, 'ws_cube');
            
            ws_cube.wrongIdx = wrongIdx;
%             ws_cube.theta = thetas;
            
            ws_cube.resolution = resolution;
            ws_cube.cubicResolution = resolution^3;
%             if any(strcmpi(varargin, 'tall'))
%                 ws_cube.grid.idxMap = tall(zeros(length(ws_cube.grid.xAxis.points), length(ws_cube.grid.yAxis.points), length(ws_cube.grid.zAxis.points)));
%             else
                ws_cube.grid.idxMap = zeros(length(ws_cube.grid.xAxis.points), length(ws_cube.grid.yAxis.points), length(ws_cube.grid.zAxis.points));
%             end
        end
        
        function [t, tau, ddq, dq, q, CirclePath, vt, a, J11, J12, J13, J21, J22, J23, J31, J32, J33] = CircularMotion_Dynamic(obj, NoOfCircles, q0_center, radi, normalToPlane, phi0, omega0, omega_max, beta, t_step)
            % GTPR Circle Motion of TCP - Open-Loop
            % NoOfCircles: Number of circles to write
            % center: Center of the circle
            % radi: Radius of the circle
            % normalToPlane: Normal to the plane in which to write the circle
            % omega: the angular velocity to go on the circle
            
            normalToPlane = normalToPlane/norm(normalToPlane);
            
            % Path Parameters: (one Circle):
            CircleCenter = obj.directGeometry(q0_center);

            pathLength_arc = 2*pi*NoOfCircles;
            
            t_ramp = omega_max/beta;
            
            phi_ramp = omega0*t_ramp + 1/2*beta*t_ramp^2;
            
            pathConstOmega_arc = pathLength_arc - 2*phi_ramp; % For acceleration and deceleration;
            
            phi_total = 2*phi_ramp + pathConstOmega_arc;
            
            t_const = pathConstOmega_arc/omega_max;
            
            t_end = 2*t_ramp + t_const;
            
            t = 0:t_step:t_end;
            
            tau = zeros(3,length(t));
            ddq = zeros(3,length(t));
            dq = zeros(3,length(t));
            q = zeros(3,length(t));
            
            phi = zeros(1,length(t)); %#ok<PROPLC>
            omega = zeros(3,length(t));
            CirclePath = zeros(3,length(t));
            beta0 = beta;
            
            phi(1) = phi0; %#ok<PROPLC>
            
            vt = zeros(3,length(t));
            at = zeros(3,length(t));
            acp = zeros(3,length(t));
            a = zeros(3,length(t));
            
            J11 = zeros(1,length(t));
            J12 = zeros(1,length(t));
            J13 = zeros(1,length(t));
            J21 = zeros(1,length(t));
            J22 = zeros(1,length(t));
            J23 = zeros(1,length(t));
            J31 = zeros(1,length(t));
            J32 = zeros(1,length(t));
            J33 = zeros(1,length(t));
            
            R = vrrotvec2mat(vrrotvec(normalToPlane,[0;0;1]));
            
            CirclePath(:,1) = CircleCenter + (R*(radi * [cos(phi(1)); sin(phi(1)); 0])); %#ok<PROPLC>
            
            q(:,1) = obj.inverseGeometry(CirclePath(:,1));
            
            
            for i = 1:length(t)
                omega(i+1) = omega(i) + beta*t_step;
                phi(i+1) = phi(i) + omega(i+1)*t_step + 1/2*beta*t_step^2; %#ok<PROPLC>
                
                if omega(i) > omega_max
                    omega(i) = omega_max;
                    beta = 0;
                elseif phi(i) > phi_total-phi_ramp %#ok<PROPLC>
                    beta = -beta0;
                else
                    % do nothing
                end
                
%                 obj.directGeometry(q(:,i));
                CirclePath(:,i) = CircleCenter + (R*(radi * [cos(phi(i)); sin(phi(i)); 0])); %#ok<PROPLC>
                
                
                if phi(i) > phi_total - phi_ramp && phi(i) > phi_ramp %#ok<PROPLC>
                    beta = -beta0;
                end
                
%                 scalar velocity, accelerations:
                    % centripetal Direction:
%                 dir_cp = (CircleCenter - obj.TCP) / norm(CircleCenter - obj.TCP);
                dir_cp = (CircleCenter - CirclePath(:,i)) / norm(CircleCenter - CirclePath(:,i));
                dir_tg = cross(dir_cp, normalToPlane) / norm(cross(dir_cp, normalToPlane));
               
%                 v_tan = omega(i)*radi;
                a_tan = beta*radi;
                a_centp = omega(i)^2*radi;
                
%                 vt(:,i) = dir_tg*v_tan;
                at(:,i) = dir_tg*a_tan;
                acp(:,i) = dir_cp*a_centp;
                
                a(:,i) = at(:,i) + acp(:,i);
                v_tan = omega(i)*radi;
                
                vt(:,i) = v_tan*dir_tg;
%                 obj.directGeometry(q(:,i));
                
                q(:,i) = obj.inverseGeometry(CirclePath(:,i));
%                 dq(:,i) = obj.inverseKinematics(vt(:,i));
%                 vt(:,i) = obj.directKinematics(dq(:,i));
                
                
                
                [tau(:,i), ddq(:,i)] = obj.inverseDynamics(vt(:,i), a(:,i));
                
                dq(:,i+1) = dq(:,i) + ddq(:,i)*t_step;
%                 q(:,i+1) = q(:,i) + dq(:,i)*t_step;
                
                J = obj.Jacobian(); %#ok<PROPLC>
                
                J11(i) = J(1,1); %#ok<PROPLC>
                J12(i) = J(1,2); %#ok<PROPLC>
                J13(i) = J(1,3); %#ok<PROPLC>
                J21(i) = J(2,1); %#ok<PROPLC>
                J22(i) = J(2,2); %#ok<PROPLC>
                J23(i) = J(2,3); %#ok<PROPLC>
                J31(i) = J(3,1); %#ok<PROPLC>
                J32(i) = J(3,2); %#ok<PROPLC>
                J33(i) = J(3,3); %#ok<PROPLC>
%                 figure(1);
%                 clf;
%                 hold on;
%                 grid on;
%                 axis equal;
%                 xlabel("X");
%                 ylabel("Y");
%                 zlabel("Z");
%                 title(t(i));
%                 
%                 scale = 0.05;
                
%                 plot3(CirclePath(:,1), CirclePath(:,2), CirclePath(:,3), 'k');

%                 plot3(CircleCenter(1), CircleCenter(2), CircleCenter(3), 'ro');
%                 plot3(CirclePath(1,1:i), CirclePath(2,1:i), CirclePath(3,1:i), 'kx');
%                 plot3([CircleCenter(1) CirclePath(1,i)], [CircleCenter(2) CirclePath(2,i)], [CircleCenter(3) CirclePath(3,i)], 'r');
%                 plot3(obj.W(1), obj.W(2), obj.W(3), 'ms');
% %                 quiver3(CirclePath(1,i), CirclePath(2,i), CirclePath(3,i), dir_cp(1), dir_cp(2), dir_cp(3), 'g');
% %                 quiver3(CirclePath(1,i), CirclePath(2,i), CirclePath(3,i), dir_tg(1), dir_tg(2), dir_tg(3), 'm');
%                 
%                 quiver3(CirclePath(1,i), CirclePath(2,i), CirclePath(3,i), scale*vt(1,i), scale*vt(2,i), scale*vt(3,i), 'm:');
%                 quiver3(CirclePath(1,i), CirclePath(2,i), CirclePath(3,i), scale*at(1,i), scale*at(2,i), scale*at(3,i), 'm');
%                 quiver3(CirclePath(1,i), CirclePath(2,i), CirclePath(3,i), scale*acp(1,i), scale*acp(2,i), scale*acp(3,i), 'g');
%                 quiver3(CirclePath(1,i), CirclePath(2,i), CirclePath(3,i), scale*a(1,i), scale*a(2,i), scale*a(3,i), 'c');
%                 view(30,30);
% 
%                 figure(2);
%                 clf;
%                 hold on;
%                 grid on;
%                 plot(1:i, q(1,1:i), 'r');
%                 plot(1:i, q(2,1:i), 'g');
%                 plot(1:i, q(3,1:i), 'b');
%                 
%                 pause(1e-5);
            end
            
            
            
        end
        
    end
end

