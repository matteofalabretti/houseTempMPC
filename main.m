%%

clear;
clc;
close all
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
% set(0,'DefaultFigureWindowStyle', 'normal')
set(0,'defaulttextInterpreter','latex')
rng('default');

%% Definizioni costanti
T_ext = 278;
k_ext = 9;

C = [6300;4600;4200];
k = [16;18;19];
tau = [580;520;540];


%forse meglio mettere un vettore da 3 (abbiamo 3 ingressi)
u = [0; 0 ; 0 ; 100 ; 100 ; 100 ];

%% ODE del sistema
dxdt = @(t,x) tempCasa(t, x, k, C, tau, T_ext, k_ext, u);


x0 = [280 285 290 0 0 0; 250 255 260 25 50 100; 285 290 295 25 50 100]';

% for i= 1:3
%     [tt, xx] = ode45(dxdt, [0 5000], x0(:, i));
%     figure
%     plot(tt, xx)
% end

[tt, xx] = ode45(dxdt, [0 5000], x0(:, 1));
figure
hold on
subplot(2,1,1)
plot(tt, xx(: , 1:3));
title("Temperature delle stanze")
legend(["T1" "T2" "T3"])

subplot(2, 1,2)
plot(tt , xx(: , 4:6));
title("Potenza termica dei termosifoni")
legend(["Q1" "Q2" "Q3"])
hold off

%% Linearizzazione

A_11 = [ [-45 -18 18] ./C(1); [18 -47 20] ./C(2) ; [18 -20 -47] ./ C(3)];
A_12 = eye(3) ./ C;
A_21 = zeros(3,3);
A_22 = -eye(3) ./tau;

A_lin = [A_11 , A_12 ; A_21 , A_22];
B_lin = [0 0 0 1/tau(1) 1/tau(2) 1/tau(3)]';
C_lin = eye(6);
D_lin = zeros(6,1);

sys_lineare = ss(A_lin, B_lin, C_lin, D_lin);

% verifica stabilità
disp("Autovalori di della matrice A linearizzata:")
disp(eig(A_lin));


%% Discretizziamo
Ts = 60;
sys_discretizzato = c2d(sys_lineare, Ts);
figure
pzmap(sys_discretizzato)

% verifica stabilità
disp("Modulo degli autovalori di della matrice A discretizzata:")
disp(abs(eig(sys_discretizzato.A)));


%% analisi della raggiungibilità

Mr_lineare = ctrb(sys_lineare);
Mr_discretizzato = ctrb(sys_discretizzato);


disp("Rango matrice di raggiungibilità sitema linearizzato:")
disp(rank(Mr_lineare))
disp("Dimensioni:")
disp(width(Mr_lineare) + " x " + height(Mr_lineare))

disp(" ")
disp("Rango matrice di raggiungibilità sitema discretizzato:")
disp(rank(Mr_discretizzato))
disp("Dimensioni:")
disp(width(Mr_discretizzato) + " x " + height(Mr_discretizzato))
