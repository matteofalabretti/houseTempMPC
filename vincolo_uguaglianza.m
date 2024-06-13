% In questo script valutiamo MPC con vincolo di disuguaglianza, quindi con
% un controli invariant set (CIS) e con R e Q uguali

%% Richiamiamo lo script di inizzializzazione

clear;
clc;
close all

%Impostiamo il tempo di campionamento
Ts = 60;

inizzializzazione

%% Definizione delle matrici del costo quadratico
Q = 1.e1*eye(6);
R = 1e1*eye(3);

%% N-step controllable set per il  vincolo termiale
G = Hx;
g = [zeros(12,1)];

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


%% simulazione a tempo continuo con il controllo e vincolo terminale

htt=[];
hxx = [];
u_online = [];
x_ini = [284 285 284 0 10 0]';


for i = 1:100


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
figure
subplot(2 , 1 , 1)
plot(htt, hxx(: , 1:3)');
title("Temperature")

subplot(2 , 1 , 2)
plot(htt, hxx(: , 4:6)')
title("Potenza Termosifoni")



figure
plot(u_online+[100 100 100])
title("Ingressi")
