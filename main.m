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
x_ref = [289 289 289 100 100 100]';
u = [100 100  100 ]';


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

syms T1 T2 T3 Q1 Q2 Q3 Q1r Q2r Q3r real
 
k12 = k(1) + (4)/(1+exp(-0.5*abs(T1 - T2)));
k13 = k(1) + (4)/(1+exp(-0.5*abs(T1 - T3)));
k23 = k(2) + (4)/(1+exp(-0.5*abs(T2 - T3)));

T1_dot = (Q1 - k12*(T1 - T2) - k13*(T1-T3) - k_ext*(T1 - T_ext))/C(1);
T2_dot = (Q2 + k12*(T1 - T2) - k23*(T2-T3) - k_ext*(T2 - T_ext))/C(2);
T3_dot = (Q3 + k13*(T1 - T3) + k23*(T2-T3) - k_ext*(T3 - T_ext))/C(3);

Q1_dot = (Q1r - Q1)/tau(1);   
Q2_dot = (Q2r - Q2)/tau(2);
Q3_dot = (Q3r - Q3)/tau(3);

F = [T1_dot; T2_dot; T3_dot; Q1_dot; Q2_dot; Q3_dot];
Stati = [T1 T2 T3 Q1 Q2 Q3];
Ingressi = [Q1r Q2r Q3r];

%jacobiane
A(T1, T2, T3, Q1, Q2, Q3) = jacobian(F , Stati);
B(Q1r, Q2r, Q3r) = jacobian(F , Ingressi);

%Valutiamo nell'equilibrio
A_lin = double(A(x_ref(1) , x_ref(2) , x_ref(3) ,x_ref(4) ,x_ref(5) ,x_ref(6)));
B_lin = double(B(u(1) , u(2) , u(3)));

C_lin = eye(6);
D_lin = zeros(6,3);

sys_lineare = ss(A_lin, B_lin, C_lin, D_lin);
% 
% % Verifica della Stabilità del sistema lineare
disp("Autovalori di della matrice A linearizzata:")
disp(eig(A_lin));

disp("Stati con ingresso costante a: " + u(1));
disp(-sys_lineare.A^-1 * sys_lineare.B * zeros(3,1));

%% Simulazione Sistema

tt_sim = 0:1:4999;
uu_sim = ones(3,5000) * 100;
xx_0sim = [284 285 284 0 10 0]';


lsim(sys_lineare , uu_sim , tt_sim , xx_0sim)

%% Discretizziamo
Ts = 60;
sys_discretizzato = c2d(sys_lineare, Ts);
figure
pzmap(sys_discretizzato)

% Verifica della Stabilità del sistema discretizzato
disp("Modulo degli autovalori di della matrice A discretizzata:")
disp(abs(eig(sys_discretizzato.A)));

%punto di equilibrio
x_eq = (eye(width(sys_discretizzato.A)) - sys_discretizzato.A)^-1 * sys_discretizzato.B * [100 100 100]'; 

disp("Punto di equilibiro con ingresso pari a 100:")
disp(x_eq);

%% simulazione sistema discretizzato

tt_sim = 0:Ts:4999;
uu_sim = ones(3,length(tt_sim)) * 100;
xx_0sim = [284 285 284 0 10 0]';

lsim(sys_discretizzato , uu_sim , tt_sim , xx_0sim)

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

Hx = [eye(6); -eye(6)];
hx = [11*ones(3,1); 50*ones(3,1); 6.5*ones(3,1); 100*ones(3,1)];
Hu = [eye(3); -eye(3)];
hu = [50*ones(3,1); 100*ones(3,1)];

%% definizione delle matrici del costo quadratico
Q = eye(6);
R = 1e-10;

%% Verifica dell'esistenza del Controllable Invariant Set
[G, g]= CIS(sys_discretizzato.A, sys_discretizzato.B, zeros(6,1), zeros(3,1), Hx, hx, Hu, hu, Q, R);

CIS_G = Polyhedron(G, g);
CIS_G = minHRep(CIS_G);
disp("Il Control invariant set è un insieme vuoto? " + boolean(CIS_G.isEmptySet));

CIS_G_T = projection(CIS_G , 1:3);
CIS_G_Q = projection(CIS_G , 4:6);

figure
CIS_G_T.plot();
figure
CIS_G_Q.plot();

%% Verifica della fattibilità del n-step controllable invariant set
Np = 10;
[Np_steps_H, Np_steps_h] = controllable_set(Hx, hx, Hu, hu, G, g, sys_discretizzato.A, sys_discretizzato.B, Np);