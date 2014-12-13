% test hinfdesign
% requires Yalmip and SeDuMi

disp 'Running an example'

a = 10;
b = 5;

A = [-a 0;1 -1];
B1 = [0 0; 0 1];
B2 = [1;0];
C1 = [0,0;b,1-b];
D11 = [0 0;-1, b];
D12 = [1;0];
D21 = [1 0;0 b];
C2 = C1;

Pss.A = A;
Pss.B = [B1 B2];
Pss.C = [C1;C2];
Pss.D = [D11,D12;D21,[0;0]];

Pcontrol = hinfsyn6241(Pss, 1, 1);

