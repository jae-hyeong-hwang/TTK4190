%%

close all; 
A = [   -0.322  0.052   0.028   -1.12   0.002;
        0       0       1       -0.001  0;
        -10.6   0       -2.87   0.46    -0.65;
        6.87    0       -0.04   -0.32   -0.02;
        0       0       0       0       -7.5];
B = [0 0 0 0 7.5].';
C = [eye(4), zeros(4,1)];
D = [0 0 0 0].';
a_p1 = -A(3,3);
a_p2 = A(3,5);



V_g = 580;
V_a = V_g;
g = 9.81;

d_max = deg2rad(30);
e_max = deg2rad(15);

zeta = 0.70;
wp = sqrt(abs(a_p2)*d_max/e_max);

k_dp = (2*zeta*wp - a_p1)/a_p2;
k_pp = d_max/e_max*sign(a_p2);

k_ip = 0:0.005:1;
phiController = tf([k_ip], [1, 0]) + k_pp + tf([k_dp, 0],[1]);
phiSystem = tf([a_p2],[1, a_p1]);
phiFeedbackLoop = phiController*phiSystem/(1+phiController*phiSystem);
figure();
p = pole(phiFeedbackLoop);
plot3(real(p), imag(p), k_ip, '*');
xlabel('Real');
ylabel('Imag');
zlabel('K_ip Gain');
grid on;
title('phi Controller loop');

%This gives k_pi = 0.04;
k_ip = -0.04;
% wp = 1.14 this means that wc = 0.114
%%
n = 20;
wc = 1/n*wp;

zeta = 0.5;
k_pc = 2*zeta*wc*V_g/g;
k_ic = 0.3*wc^2*V_g/g;

Simulation_Time = 1000;
sim('planeModel',Simulation_Time)
figure();
plot(out.inputX)
hold on;
plot(out.outputX)
title('Simple system');

sim('planeModelFullModel',Simulation_Time)
figure();
plot(out.inputX)
hold on;
plot(maout.outputX)
title('Full state space system');



