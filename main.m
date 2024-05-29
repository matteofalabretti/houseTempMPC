% Definizioni costanti
T_ext = 278;
k_ext = 9;

C = [6300;4600;4200];
k = [16;18;19];
tau = [580;520;540];

u = [0; 0 ; 0 ; 100 ; 100 ; 100 ];

% ODE del sistema
dxdt = @(t,x) tempCasa(t, x, k, C, tau, T_ext, k_ext, u);


x0 = [280 285 290 0 0 0; 250 255 260 25 50 100; 285 290 295 25 50 100]';

% for i= 1:3
%     [tt, xx] = ode45(dxdt, [0 5000], x0(:, i));
%     figure
%     plot(tt, xx)
% end

[tt, xx] = ode45(dxdt, [0 5000], x0(:, 1));
figure
plot(tt, xx);

% Linearizzazione
A_lin = [-eye(3)*k_ext eye(3); zeros(3,3) -eye(3)./tau];
B_lin = [0 0 0 1/tau(1) 1/tau(2) 1/tau(3)]';
C_lin = eye(6);
D_lin = zeros(6,1);

sys_lineare = ss(A_lin, B_lin, C_lin, D_lin);

% Discretizziamo
Ts = 60;
sys_discretizzato = c2d(sys_lineare, Ts);
pzmap(sys_discretizzato)

%PUSH