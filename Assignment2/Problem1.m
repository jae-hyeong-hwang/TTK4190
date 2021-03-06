
A = [   -0.322  0.052   0.028   -1.12   0.002;
        0       0       1       -0.001  0;
        -10.6   0       -2.87   0.46    -0.65;
        6.87    0       -0.04   -0.32   -0.02;
        0       0       0       0       -7.5];
    
B = [0 0 0 0 7.5].';

C = [eye(4), zeros(4,1)];

%Add comment 
%% Task c
A1 = [  A(1,1), A(1,4);
        A(4,1), A(4,4)];

B = [];
C = eye(2)
D = []
sys = ss(A1, B, C, D)
damp(sys)


%% Task d
w = -A(4,1)*A(3,4)/A(3,1) + A(4,4)
Ts = 1/(2*pi*-w)

%% Task e
Tr = 1/(2*pi*2.87)