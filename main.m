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

c = [6300;4600;4200];
k = [16;18;19];
tau = [580;520;540];

% Definizione Obbiettivi di Controllo
x_ref = [289; 289; 289; 100; 100; 100];
% u a 3 ingressi
u = [100 ; 100 ; 100 ];


%% ODE del sistema
dxdt = @(t,x) tempCasa(t, x, k, c, tau, T_ext, k_ext, u);

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

A_11 = [ [-45 -18 18] ./c(1); [18 -47 20] ./c(2) ; [18 -20 -47] ./ c(3)];
A_12 = eye(3) ./ c;
A_21 = zeros(3,3);
A_22 = -eye(3) ./tau;

A = [A_11 , A_12 ; A_21 , A_22];
B = [zeros(3,3); eye(3) ./ tau];
C = eye(6);
D = zeros(6,3);

sys_lineare = ss(A, B, C, D);

% Verifica della Stabilità del sistema lineare
disp("Autovalori di della matrice A linearizzata:")
disp(eig(A));


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
Hu = [eye(3); -eye(3)]; % vettore colonna così poi quando viene 
% moltiplicato per K (vettore riga) si forma la matrice
hu = [150*ones(3,1); zeros(3,1)]; %vettore dei vincoli

%% definizione delle matrici del costo quadratico
Q = eye(6);
R = 1;

%% Verifica dell'esistenza del Controllable Invariant Set

[G, g]= CIS(sys_discretizzato.A, sys_discretizzato.B, x_ref, u, Hx, hx, Hu, hu, Q, R);

cis = Polyhedron(G , g);

figure(3)
cis.plot();
