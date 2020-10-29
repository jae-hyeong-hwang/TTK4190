% Project in TTK4190 Guidance and Control of Vehicles 
%
% Author:           Edvard Grødem
% Study program:    Mttk

clear;
clc;

curretEnabled = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h  = 0.1;    % sampling time [s]
Ns = 10000;  % no. of samples

psi_ref = [ones(1, Ns-5000)*10*pi/180, -ones(1, 5001)*20*pi/180];%10 * pi/180;  % desired yaw angle (rad)
U_d = 7;                % desired cruise speed (m/s)
               
% ship parameters 
m = 17.0677e6;          % mass (kg)
Iz = 2.1732e10;         % yaw moment of inertia about CO (kg m^3)
xg = -3.7;              % CG x-ccordinate (m)
L = 161;                % length (m)
B = 21.8;               % beam (m)
T = 8.9;                % draft (m)
KT = 0.7;               % propeller coefficient (-)
Dia = 3.3;              % propeller diameter (m)
rho = 1025;             % density of water (kg/m^3)
visc = 1e-6;            % kinematic viscousity at 20 degrees (m/s^2)
eps = 0.001;            % a small number added to ensure that the denominator of Cf is well defined at u=0
k = 0.1;                % form factor giving a viscous correction
t_thr = 0.05;           % thrust deduction number

% rudder limitations
delta_max  = 40 * pi/180;        % max rudder angle      (rad)
Ddelta_max = 5  * pi/180;        % max rudder derivative (rad/s)

% added mass matrix about CO
Xudot = -8.9830e5;
Yvdot = -5.1996e6;
Yrdot =  9.3677e5;
Nvdot =  Yrdot;
Nrdot = -2.4283e10;
MA = -[ Xudot 0    0 
        0 Yvdot Yrdot
        0 Nvdot Nrdot ];

% rigid-body mass matrix
MRB = [ m 0    0 
        0 m    m*xg
        0 m*xg Iz ];
    
Minv = inv(MRB + MA); % Added mass is included to give the total inertia

% ocean current in NED
Vc = 1;                                         % current speed (m/s)
betaVc = deg2rad(45);                           % current direction (rad)

% wind expressed in NED
Vw = 10;                   % wind speed (m/s)
betaVw = deg2rad(135);     % wind direction (rad)
rho_a = 1.247;             % air density at 10 deg celsius
cy = 0.95;                 % wind coefficient in sway
cn = 0.15;                 % wind coefficient in yaw
A_Lw = 10 * L;             % projected lateral area

% linear damping matrix (only valid for zero speed)
T1 = 20; T2 = 20; T6 = 10;

Xu = -(m - Xudot) / T1;
Yv = -(m - Yvdot) / T2;
Nr = -(Iz - Nrdot)/ T6;
D = diag([-Xu -Yv -Nr]);         % zero speed linear damping

% rudder coefficients (Section 9.5)
b = 2;
AR = 8;
CB = 0.8;

lambda = b^2 / AR;
tR = 0.45 - 0.28*CB;
CN = 6.13*lambda / (lambda + 2.25);
aH = 0.75;
xH = -0.4 * L;
xR = -0.5 * L;

X_delta2 = 0.5 * (1 - tR) * rho * AR * CN;
Y_delta = 0.25 * (1 + aH) * rho * AR * CN; 
N_delta = 0.25 * (xR + aH*xH) * rho * AR * CN;   

% input matrix
Bu = @(u_r,delta) [ (1-t_thr)  -u_r^2 * X_delta2 * delta
                        0      -u_r^2 * Y_delta
                        0      -u_r^2 * N_delta            ];
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
% Heading Controller
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rudder control law
wb = 0.06;
zeta = 1;
wn = 1 / sqrt( 1 - 2*zeta^2 + sqrt( 4*zeta^4 - 4*zeta^2 + 2) ) * wb;

% linearized sway-yaw model (see (7.15)-(7.19) in Fossen (2021)) used
% for controller design. The code below should be modified.

C_A_lin = [ 0           0               0;
            0           0            -Xudot*U_d;
            0   (Xudot-Yvdot)*U_d    -Yrdot*U_d];
C_RB_lin =[ 0           0               0;
            0           0             m*U_d;
            0           0           m*xg*U_d];
N_lin =C_RB_lin(2:3, 2:3)+  C_A_lin(2:3, 2:3) + D(2:3, 2:3);
M_lin = (MA(2:3, 2:3) + MRB(2:3, 2:3));
M_lin_inv = inv(M_lin); 
b_lin = [-2*U_d*Y_delta -2*U_d*N_delta]';
[btf, atf] = ss2tf(-M_lin_inv*N_lin, M_lin_inv*b_lin,  [0, 1], [0]);
lamda1 = (-atf(2)+sqrt(atf(2)^2-4*atf(3)))/2;
lamda2 = (-atf(2)-sqrt(atf(2)^2-4*atf(3)))/2;
T_1 = -1/lamda1;
T_2 = -1/lamda2;
T_3 = btf(2)/btf(3);
K_sys = btf(3)/(lamda1*lamda2);
T_sys = T_1 + T_2- T_3;

% tfLin = tf(btf, atf);
% tfLinNom = tf(K_sys, [T_sys, 1]);
% figure();
% step(tfLin);
% hold on;
% step(tfLinNom);

wb = 0.06;
chi = 1;
wref = 0.03;
wn = 1/sqrt(1-2*chi^2+sqrt(4*chi^4-4*chi^2+2))*wb;
mpid = T_sys;
dpid = 1;
kpid = 0;
K_p = mpid*wn^2-kpid;
K_d = 2*chi*wn*mpid-dpid;
K_i = wn*K_p/10;

[~, ~, Cr1, Dr1] = tf2ss([wref^3], [1, (2*chi+1)*wref, (2*chi+1)*wref^2, wref^3]); % yaw referance
[~, ~, Cr2, Dr2] = tf2ss([wref^3, 0], [1, (2*chi+1)*wref, (2*chi+1)*wref^2, wref^3]); % yaw rate referance
[Ar, Br, Cr3, Dr3] = tf2ss([wref^3, 0, 0], [1, (2*chi+1)*wref, (2*chi+1)*wref^2, wref^3]); % yaw rate rate feed forward

Xref = zeros(3, 1);
e_i = 0;
% initial states
eta = [0 0 0]';
nu  = [0.1 0 0]';
delta = 0;
n = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simdata = zeros(Ns+1,15);                % table of simulation data

for i=1:Ns+1

    t = (i-1) * h;                      % time (s)
    R = Rzyx(0,0,eta(3));               % R_n_b
    
    % current 
    if curretEnabled
        nu_n_c = [-cos(betaVc) -sin(betaVc) 0].'*Vc;    % current speed decomposed
        nu_r = nu - R.'*nu_n_c;                             
    else
        nu_r = nu;
    end
    

    % wind 
    if t > 200
        gamma_rw = eta(3) -betaVw -pi;
        C_Y = cy*sin(gamma_rw);
        C_N = cn*sin(2*gamma_rw);
        Ywind = 0.5 * rho_a * Vw^2 * C_Y * A_Lw;    % Wind moment in sway
        Nwind = 0.5 * rho_a * Vw^2 * C_N * A_Lw*L;  % Wind moment in yaw
    else
        Ywind = 0;
        Nwind = 0;
    end
    tau_env = [0 Ywind Nwind]';
    
    % state-dependent time-varying matrices
    CRB = m * nu(3) * [ 0 -1 -xg 
                        1  0  0 
                        xg 0  0  ];
                    
    % coriolis due to added mass
    CA = [  0   0   Yvdot * nu_r(2) + Yrdot * nu_r(3)
            0   0   -Xudot * nu_r(1) 
          -Yvdot * nu_r(2) - Yrdot * nu_r(3)    Xudot * nu_r(1)   0];
    N = CRB + CA + D;
    
    % nonlinear surge damping
    Rn = L/visc * abs(nu_r(1));
    Cf = 0.075 / ( (log(Rn) - 2)^2 + eps);
    Xns = -0.5 * rho * (B*L) * (1 + k) * Cf * abs(nu_r(1)) * nu_r(1);
    
    % cross-flow drag
    Ycf = 0;
    Ncf = 0;
    dx = L/10;
    Cd_2D = Hoerner(B,T);
    for xL = -L/2:dx:L/2
        vr = nu_r(2);
        r = nu_r(3);
        Ucf = abs(vr + xL * r) * (vr + xL * r);
        Ycf = Ycf - 0.5 * rho * T * Cd_2D * Ucf * dx;
        Ncf = Ncf - 0.5 * rho * T * Cd_2D * xL * Ucf * dx;
    end
    d = -[Xns Ycf Ncf]';
    
    % reference models
    Xref = Xref + (Ar*Xref + Br*psi_ref(i))*h;
    
    psi_d = Cr1*Xref+ Dr1*psi_ref(i);
    r_d = Cr2*Xref+ Dr2*psi_ref(i);
    u_d = U_d; 
    
    % thrust 
    thr = rho * Dia^4 * KT * abs(n) * n;    % thrust command (N)
        
    % control law
    psi = eta(3);
    e_psi = (psi-psi_d);
    r = nu(3);
    feedForward = 1/K_sys*(T_sys*(Cr3*Xref+ Dr3*psi_ref(i))+r_d);
    delta_c = feedForward - 1/K_sys*(K_p*e_psi + K_d*r + K_i*e_i);  
    %delta_c = 0.1;              % rudder angle command (rad)
    
    % ship dynamics
    u = [ thr delta ]';
    tau = Bu(nu_r(1),delta) * u;
    nu_dot = Minv * (tau_env + tau - N * nu_r - d); 
    eta_dot = R * nu;    
    toIntegrate = true;
    % Rudder saturation and dynamics (Sections 9.5.2)
    if abs(delta_c) >= delta_max
        delta_c = sign(delta_c)*delta_max;
        toIntegrate = false;
    end
    
    delta_dot = delta_c - delta;
    if abs(delta_dot) >= Ddelta_max
        delta_dot = sign(delta_dot)*Ddelta_max;
        toIntegrate = false;
    end
    if toIntegrate
       e_i = e_i +  e_psi*h;
    end

    % propeller dynamics
    Im = 100000; Tm = 10; Km = 0.6;         % propulsion parameters
    n_c = 10;                               % propeller speed (rps)
    n_dot = (1/10) * (n_c - n);             % should be changed in Part 3
    
    % store simulation data in a table (for testing)
    beta = asin(nu_r(2)/norm(nu_r));
    simdata(i,:) = [t n_c delta_c n delta eta' nu' u_d psi_d r_d, beta];       
     
    % Euler integration
    eta = euler2(eta_dot,eta,h);
    nu  = euler2(nu_dot,nu,h);
    delta = euler2(delta_dot,delta,h);   
    n  = euler2(n_dot,n,h);    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t       = simdata(:,1);                 % s
n_c     = 60 * simdata(:,2);            % rpm
delta_c = (180/pi) * simdata(:,3);      % deg
n       = 60 * simdata(:,4);            % rpm
delta   = (180/pi) * simdata(:,5);      % deg
x       = simdata(:,6);                 % m
y       = simdata(:,7);                 % m
psi     = (180/pi) * simdata(:,8);      % deg
u       = simdata(:,9);                 % m/s
v       = simdata(:,10);                % m/s
r       = (180/pi) * simdata(:,11);     % deg/s
u_d     = simdata(:,12);                % m/s
psi_d   = (180/pi) * simdata(:,13);     % deg
r_d     = (180/pi) * simdata(:,14);     % deg/s
beta    = (180/pi) * simdata(:,15);     % deg
beta_c  = (180/pi) * atan2(v, u);       % deg


figure(1)
figure(gcf)
subplot(311)
plot(y,x,'linewidth',2); axis('equal')
title('North-East positions'); xlabel('x (m)'); ylabel('y (m)'); 
subplot(312)
plot(t,psi,t,psi_d,'linewidth',2);
title('Actual and desired yaw angles (deg)'); xlabel('time (s)');
subplot(313)
plot(t,r,t,r_d,'linewidth',2);
title('Actual and desired yaw rates (deg/s)'); xlabel('time (s)');

figure(2)
figure(gcf)
subplot(311)
plot(t,u,t,u_d,'linewidth',2);
title('Actual and desired surge velocities (m/s)'); xlabel('time (s)');
subplot(312)
plot(t,n,t,n_c,'linewidth',2);
title('Actual and commanded propeller speed (rpm)'); xlabel('time (s)');
subplot(313)
plot(t,delta,t,delta_c,'linewidth',2);
title('Actual and commanded rudder angles (deg)'); xlabel('time (s)');

figure(3) 
figure(gcf)
subplot(211)
plot(t,u,'linewidth',2);
title('Actual surge velocity (m/s)'); xlabel('time (s)');
subplot(212)
plot(t,v,'linewidth',2);
title('Actual sway velocity (m/s)'); xlabel('time (s)');

figure(4)
figure(gcf)
subplot(211)
plot(t,beta,'linewidth',2);
title('Actual sideslip angle \beta (deg) '); xlabel('time (s)');
subplot(212)
plot(t,beta_c,'linewidth',2);
title('Actual crab angle \beta_c (deg)'); xlabel('time (s)');
