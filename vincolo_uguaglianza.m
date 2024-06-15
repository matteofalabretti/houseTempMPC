% In questo script valutiamo MPC con vincolo di uguaglianza, con R e Q uguali
% A tentativi abbiamo trovato che Np deve essere di almeno 20 passi con:
% Ts = 60s e Q e R identità

clear;
clc;
close all

%% Impostazioni dell script
%Impostiamo il tempo di campionamento
Ts = 60; % [secondi]
% Definizione delle matrici del costo quadratico
Q = 1.e1*eye(6);
R = 1e1*eye(3);

%% Richiamiamo lo script di inizzializzazione
inizializzazione

%% Creaiamo i vincoli iniziali
% poniamo: 0 <= x <= 0
% quindi: x = 0
G = Hx;
g = [zeros(12,1)];

setIniziale = Polyhedron(G , g);
disp("Il set iniziale è vuoto? " + setIniziale.isEmptySet)

%% LA PROSSIMA SEZIONE FA BLOCCARE MATLAB
%% Non Far andare!!! dopo 3 step si blocca matlab
% [Np_steps_H, Np_steps_h , Np] = controllable_set(Hx, hx, Hu, hu, G, g, sys_discretizzato.A, sys_discretizzato.B, x0_centrato);

%% Verifica fattibilita dal punto di partenza
% trasp = 0.3; %impostiamo la trasparenza delle figure
% 
% Np_step = Polyhedron(Np_steps_H , Np_steps_h);
% Np_step = Np_step.minHRep();
% 
% Np_steps_T = projection(Np_step , 1:3);
% Np_steps_Q = projection(Np_step , 4:6);
% 
% 
% figure
% Np_steps_T.plot();
% title("Proiezione del dominio di attrazione delle temperature nelle stanze")
% limitiTemp = [X_vinc_lin(7) X_vinc_lin(1)];
% xlim(limitiTemp)
% ylim(limitiTemp)
% zlim(limitiTemp)
% trasparenzaFigura(trasp)
% hold on
% plot3(x0_centrato(1) ,x0_centrato(2), x0_centrato(3) , "." , MarkerSize=50)
% 
% legend(["n-steps" , "Punto di partenza"])
% 
% figure;
% Np_steps_Q.plot();
% title("Proiezione del dominio di attrazione della potenza termica dei termosifoni")
% limitiQ = [X_vinc_lin(10) X_vinc_lin(4)];
% xlim(limitiQ)
% ylim(limitiQ)
% zlim(limitiQ)
% trasparenzaFigura(trasp)
% hold on
% plot3(x0_centrato(4) ,x0_centrato(5), x0_centrato(6) , "." , MarkerSize=50)
% 
% legend(["n-steps" , "Punto di partenza"])

%% simulazione a tempo continuo con il controllo e vincolo terminale
Np = 20;

htt=[];
hxx = [];
u_online = [];
x_ini = [284 285 284 0 10 0]';

n_sim = 100;

for i = 1:n_sim

    if i == 1
        x_run = x_ini-x_ref(1:6);
    else
        x_run = hxx(: , end)-x_ref(1:6);
    end

    controlAction = MPC_Uguaglianza(x_run, sys_discretizzato, Q, R, Np, G,g, X_vinc_lin, U_vinc_lin);
    tempo = linspace(Ts*(i-1), Ts*i , Ts);
    controlAction = controlAction(1:3) + [100; 100; 100];
    u_online = [u_online,repmat(controlAction , 1 , Ts)];
    dxdt = @(t,x) tempCasa(t, x, k, C, tau, T_ext, k_ext, controlAction);
    [tt, xx] = ode45(dxdt , tempo , x_run+x_ref);
    htt = [htt,tt'];
    hxx = [hxx,xx'];

end

%% Plot

plotSimulazione(htt , hxx , u_online , x_ref);
