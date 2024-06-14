% In questo script valutiamo MPC con vincolo di uguaglianza, con R e Q uguali

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
[Np_steps_H, Np_steps_h , Np] = controllable_set(Hx, hx, Hu, hu, G, g, sys_discretizzato.A, sys_discretizzato.B, x0_centrato);

%% Verifica fattibilita dal punto di partenza
trasp = 0.3; %impostiamo la trasparenza delle figure

Np_step = Polyhedron(Np_steps_H , Np_steps_h);
Np_step = Np_step.minHRep();

Np_steps_T = projection(Np_step , 1:3);
Np_steps_Q = projection(Np_step , 4:6);


figure
Np_steps_T.plot();
title("Proiezione del dominio di attrazione delle temperature nelle stanze")
limitiTemp = [X_vinc_lin(7) X_vinc_lin(1)];
xlim(limitiTemp)
ylim(limitiTemp)
zlim(limitiTemp)
trasparenzaFigura(trasp)
hold on
plot3(x0_centrato(1) ,x0_centrato(2), x0_centrato(3) , "." , MarkerSize=50)

legend(["n-steps" , "Punto di partenza"])

figure;
Np_steps_Q.plot();
title("Proiezione del dominio di attrazione della potenza termica dei termosifoni")
limitiQ = [X_vinc_lin(10) X_vinc_lin(4)];
xlim(limitiQ)
ylim(limitiQ)
zlim(limitiQ)
trasparenzaFigura(trasp)
hold on
plot3(x0_centrato(4) ,x0_centrato(5), x0_centrato(6) , "." , MarkerSize=50)

legend(["n-steps" , "Punto di partenza"])

%% simulazione a tempo continuo con il controllo e vincolo terminale

htt=[];
hxx = [];
u_online = [];
x_ini = [284 285 284 0 10 0]';
n_sim = 100;

for i = 1:n_sim


    if i == 1
        x_run = x_ini-x_ref(1:6);
    else
        x_run = xx(height(xx), 1:6)'-x_ref(1:6);
    end

    controlAction = MPC_Uguaglianza(x_run, sys_discretizzato, Q, R, Np, G,g, X_vinc_lin, U_vinc_lin);
    u_online = [u_online;controlAction(1:3)'];
    [tt, xx] = ode45(@(t,x) tempCasa(t, x, k, C, tau, T_ext, k_ext, [100; 100; 100] + controlAction(1:3)), [60*(i-1) 60*i], x_run+x_ref);
    htt = [htt;tt];
    hxx = [hxx;xx];

end

%% Plot
tempo = (1:n_sim) * Ts/60; %[min]

figure

sgtitle("Evoluzioni degli stati")

subplot(2 , 1, 1)
plot( tempo, storia_x(1:3 , :) + x_ref(1:3))
yline(x_ref(1))
legend(["T1" , "T2" , "T3" ,"Obbiettivo"])
ylabel("Temperatura $[^{\circ}C]$" , Interpreter="latex");
xlabel("Tempo $[min]$" , Interpreter="latex");
title("Temperatura")

subplot(2 , 1, 2)
plot(tempo, storia_x(4:end , :) + x_ref(4:end))
yline(x_ref(4))
ylim([0 , 150])
legend(["Q1"  "Q2"  "Q3" "Obbiettivo"])
ylabel("Potenza $[W]$" , Interpreter="latex");
xlabel("Tempo $[min]$" , Interpreter="latex");
title("Potenza dei termosifoni")


figure
plot(tempo, storia_u + u_ref)
title("Azioni di controllo")
ylabel("Potenza $[W]$" , Interpreter="latex");
xlabel("Tempo $[min]$" , Interpreter="latex");
legend(["Q1"  "Q2"  "Q3"])