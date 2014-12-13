function subOptimalController = hinfdesign(openLoopPlant, nMeasurements, nControls)

eps = 5E-4;

A = openLoopPlant.a;
B = openLoopPlant.b;
C = openLoopPlant.c;
D = openLoopPlant.d;
% 
nNoises = size(B,2)-nControls;
nCosts = size(C,1)-nMeasurements;
dim = size(A,1);

% Plant
A = A;
B1 = B(1:dim,1:(nNoises));
B2 = B(1:dim,(nNoises+1):end);
C1 = C(1:nCosts,1:dim);
C2 = C((nCosts+1):end,1:dim);
D12 = D(1:nCosts,(nNoises+1):end);
D21 = D((nCosts+1):end,1:nNoises);
D11 = D(1:nCosts,1:nNoises);
D22 = D((nCosts+1):end,(nNoises+1):end);

% Define SDP variables with Yalmip
h = sdpvar(1,1);
g = sdpvar(1,1);
P = sdpvar(dim,dim);
Q = sdpvar(dim,dim);

% Set up PD constraints for the SDP
Xp = [-g*C1'*C1-P*A-A'*P, -g*C1'*D11-P*B1;...
    -g*D11'*C1-B1'*P, h*eye(size(D11'*D11))-g*D11'*D11];

Lz = null([C2 D21]);

Xq = [-g*B1*B1'-Q*A'-A*Q, -g*B1*D11'-Q*C1';...
      -g*D11*B1'-C1*Q, h*eye(size(D11*D11'))-g*D11*D11'];

Lf = null([B2' D12']);

PIQ = [P, eye(dim);...
       eye(dim), Q];
   
h1g = [h, 1;...
       1, g];

% These are the four positive definite constraints
warning('off','YALMIP:strict')
PDs = [Lz'*Xp*Lz>0, Lf'*Xq*Lf>0, PIQ>0, h1g>0];
%    
% We must actively select SeDuMi as the solve Yalmip should use 
OPTIONS = sdpsettings('solver', 'sedumi','verbose',0);

% % SeDuMi via Yalmip to solve SDP, minimizing for h
solvesdp(PDs,h,OPTIONS);
warning('on','YALMIP:strict')

Pd = double(P);
Qd = double(Q);
hd = double(h);
gd = 1/hd;

disp(['hinfsyn6241 finds h as: ' num2str(hd)]);

% Solve for symmetric H
iQd = Qd\eye(size(Qd));
Hd = [Pd, iQd-Pd;...
    iQd-Pd, Pd-iQd];

% Intermediate parameters to simplify the derivation
R = [zeros(nNoises,dim), eye(nNoises), zeros(nNoises,dim)];
S = [C1, D11, zeros(size(C1,1),dim)];
T = [zeros(size(D12,1),dim), D12];
U = [eye(dim), zeros(dim,nNoises), zeros(dim);...
     zeros(dim), zeros(dim,nNoises), eye(dim)];
V = [A, B1, zeros(dim);...
     zeros(dim, dim+dim+size(B1,2))];
W = [zeros(size(B2,1),dim), B2;...
     eye(dim), zeros(dim,size(B2,2))];

% Collected quadratic form
cA = hd*R'*R-gd*S'*S-U'*Hd*V-V'*Hd*U; % This is ensured to be symmetric
cB = -gd*S'*T-U'*Hd*W;
cC = -gd*T'*T;

% connection with Parrotts lemma
C_pl = [zeros(dim), zeros(dim,nNoises), eye(dim);...
        C2, D21, zeros(size(C2,1),dim)];

% Regularized parameters
ceA = cA - eps*eye(size(cA,1));
ceC = cC - eps*eye(size(cC,1));

% Solve for optimal controller after having obtained H as composite storage function
iceA = ceA\eye(size(ceA));
Fi = (ceC - cB'*iceA*cB);
iFi = Fi\eye(size(Fi));
Se = (C_pl*iceA*C_pl');
iSe = Se\eye(size(Se));
lZ = iFi*cB'*iceA*C_pl'*iSe;
offset = lZ*C_pl*iceA*cB;
F = (eye(size(offset))-offset)\lZ;

% Return the controller
subOptimalController.A = F(1:dim,1:dim);
subOptimalController.B = F(1:dim,(dim+1):end);
subOptimalController.C = F((dim+1):end,1:dim);
subOptimalController.D = F((dim+1):end,(dim+1):end);

end
