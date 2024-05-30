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

% Definizione Obbiettivi di Controllo
x_ref = [289 289 289 100 100 100];
u = [0; 0 ; 0 ; 100 ; 100 ; 100 ];


%% ODE del sistema
dxdt = @(t,x) tempCasa(t, x, k, C, tau, T_ext, k_ext, u);

% simulazione del comportamento del sistema al variare degli stati
x0 = [280 285 290 0 0 0; 250 255 260 25 50 100; 285 290 295 25 50 100]';

for i= 1:3
    [tt, xx] = ode45(dxdt, [0 5000], x0(:, i));
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

end

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

% Verifica della Stabilità del sistema lineare
disp("Autovalori di della matrice A linearizzata:")
disp(eig(A_lin));


%% Discretizziamo
Ts = 60;
sys_discretizzato = c2d(sys_lineare, Ts);
figure
pzmap(sys_discretizzato)

% Verifica della Stabilità del sistema discretizzato
disp("Modulo degli autovalori di della matrice A discretizzata:")
disp(abs(eig(sys_discretizzato.A)));


%% analisi della raggiungibilità

% computazione delle matici di raggiungibilità
Mr_lineare = ctrb(sys_lineare);
Mr_discretizzato = ctrb(sys_discretizzato);

disp("Rango matrice di raggiungibilità sitema linearizzato: " + rank(Mr_lineare))
disp("Dimensioni: " + width(Mr_lineare) + " x " + height(Mr_lineare))

disp("------------------------------------------------------------")
disp("Rango matrice di raggiungibilità sitema discretizzato: " + rank(Mr_discretizzato))
disp("Dimensioni: " + width(Mr_discretizzato) + " x " + height(Mr_discretizzato))


%% Definizione dei vincoli su stato e su ingressi
%U_vinc = [0 150]; % [W]
%X_vinc = [282.5 300]; % [K]

Hx = [eye(6);-eye(6)];
hx = [300*ones(3,1); 150*ones(3,1); 282.5*ones(3,1); zeros(3,1)];
Hu = [eye(6); -eye(6)];
hu = [150*ones(3,1); zeros(3,1)];

%% definizione delle matrici del costo quadratico
Q = eye(6);
R = 1;

%% Verifica dell'esistenza del Controllable Invariant Set
[G, g]= CIS(sys_discretizzato.A, sys_discretizzato.B, x_ref, u, Hx, hx, Hu, hu, Q, R);

%% Verifica della fattibilità del n-step controllable invariant set
Np = 10;
[Np_steps_H, Np_steps_h] = controllable_set(Hx, hx, Hu, hu, G, g, A, B, Np);