% In questo script valutiamo MPC con vincolo di disuguaglianza, quindi con un control invariant set (CIS)

% * Se Q > 1e2 probabilmente il calcolo del N Steps controllable set si
% blocca quindi in quel caso saltare la parte e impostare manualmente i
% passi facendo: Np = #passi * 

clear;
clc;
close all

%% Impostazioni dell script
%Impostiamo il tempo di campionamento
Ts = 60; % [secondi]

%% Richiamiamo lo script di inizzializzazione
inizializzazione

%% Definizione delle matrici del costo quadratico
Q = 1.e3*eye(6);
R = 1e1*eye(3);
% S come soluzione di Riccati
[~ , S] = dlqr(sys_discretizzato.A , sys_discretizzato.B , Q , R); 

%% Verifica dell'esistenza del Controllable Invariant Set
[G, g]= CIS(sys_discretizzato.A, sys_discretizzato.B, zeros(6,1), zeros(3,1), Hx, hx, Hu, hu, Q, R);


%% Plot del CIS
CIS_G = Polyhedron(G, g);
CIS_G = minHRep(CIS_G);

CIS_G_T = projection(CIS_G , 1:3);
CIS_G_Q = projection(CIS_G , 4:6);

figure
subplot(1 , 2 , 1)
CIS_G_T.plot();
title("Proiezione del CIS delle temperature nelle stanze")
limitiTemp = [X_vinc_lin(7) X_vinc_lin(1)];
xlim(limitiTemp)
ylim(limitiTemp)
zlim(limitiTemp)
xlabel("T1 $[ ^{\circ}C]$" , Interpreter="latex")
ylabel("T2 $[ ^{\circ}C]$" , Interpreter="latex")
zlabel("T3 $[ ^{\circ}C]$" , Interpreter="latex")

subplot(1 , 2 , 2);
CIS_G_Q.plot();
title("Proiezione del CIS della potenza termica dei termosifoni")
limitiQ = [X_vinc_lin(10) X_vinc_lin(4)];
xlim(limitiQ)
ylim(limitiQ)
zlim(limitiQ)
xlabel("Q1 $[W]$" , Interpreter="latex")
ylabel("Q2 $[W]$" , Interpreter="latex")
zlabel("Q3 $[W]$" , Interpreter="latex")

%% N-step controllable set

[Np_steps_H, Np_steps_h , Np] = controllable_set(Hx, hx, Hu, hu, G, g, sys_discretizzato.A, sys_discretizzato.B, x0_centrato);

disp(" ")
disp("Passi minimi per entrare nel CIS: " + Np);
%% Verifica fattibilita dal punto di partenza
trasp = 0.3; %impostiamo la trasparenza delle figure

Np_step = Polyhedron(Np_steps_H , Np_steps_h);
Np_step = Np_step.minHRep();

Np_steps_T = projection(Np_step , 1:3);
Np_steps_Q = projection(Np_step , 4:6);


figure
subplot(1 , 2 ,1)
Np_steps_T.plot();
title("Temperature nelle stanze")
limitiTemp = [X_vinc_lin(7) X_vinc_lin(1)];
xlim(limitiTemp)
ylim(limitiTemp)
zlim(limitiTemp)
trasparenzaFigura(trasp)
hold on
plot3(x0_centrato(1) ,x0_centrato(2), x0_centrato(3) , "." , MarkerSize=50)
xlabel("T1 $[ ^{\circ}C]$" , Interpreter="latex")
ylabel("T2 $[ ^{\circ}C]$" , Interpreter="latex")
zlabel("T3 $[ ^{\circ}C]$" , Interpreter="latex")

% legend(["n-steps" , "Punto di partenza"])

subplot(1 , 2 ,2)
Np_steps_Q.plot();
title("Potenza termica dei termosifoni")
limitiQ = [X_vinc_lin(10) X_vinc_lin(4)];
xlim(limitiQ)
ylim(limitiQ)
zlim(limitiQ)
trasparenzaFigura(trasp)
hold on
plot3(x0_centrato(4) ,x0_centrato(5), x0_centrato(6) , "." , MarkerSize=50)
xlabel("Q1 $[W]$" , Interpreter="latex")
ylabel("Q2 $[W]$" , Interpreter="latex")
zlabel("Q3 $[W]$" , Interpreter="latex")

% legend(["n-steps" , "Punto di partenza"])

%% Simulazione a tempo continuo con il controllo

if ~exist("Np" , "var")
    Np = 10;
end

htt=[];
hxx = [];
u_online = [];
x_ini = [284 285 284 0 10 0]';
flag = zeros(1 , Np);
n_sim = 100;

[A_cal , A_cal_n , B_cal , B_cal_n,  Q_cal , R_cal] = Calligrafica(sys_discretizzato.A , sys_discretizzato.B , Q , R , S , Np);


for i = 1:n_sim

    if i == 1
        x_run = x_ini-x_ref(1:6);
    else
        x_run = hxx(: , end)-x_ref(1:6);
    end

    [controlAction , flag(i)]= MPC(x_run , A_cal , A_cal_n , B_cal , B_cal_n,  Q_cal , R_cal , Np, G,g, X_vinc_lin, U_vinc_lin);
    tempo = linspace(Ts*(i-1), Ts*i , Ts);
    controlAction = controlAction(1:3) + [100; 100; 100];
    u_online = [u_online,repmat(controlAction , 1 , Ts)];
    dxdt = @(t,x) tempCasa(t, x, k, C, tau, T_ext, k_ext, controlAction);
    [tt, xx] = ode45(dxdt , tempo , x_run+x_ref);
    htt = [htt,tt'];
    hxx = [hxx,xx'];

end

%% Plot

[plot_T  , plot_Q , plot_U] = plotSimulazione(htt , hxx , u_online , x_ref);

