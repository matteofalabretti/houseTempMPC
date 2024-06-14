% In questo script valutiamo MPC con vincolo di disuguaglianza, quindi con
% un controli invariant set (CIS) e con R e Q uguali


clear;
clc;
close all

%% Impostazioni dell script
%Impostiamo il tempo di campionamento
Ts = 30; % [secondi]
% Definizione delle matrici del costo quadratico
Q = 1.e2*eye(6);
R = 1e1*eye(3);


%% Richiamiamo lo script di inizzializzazione
inizzializzazione

%% Verifica dell'esistenza del Controllable Invariant Set
[G, g]= CIS(sys_discretizzato.A, sys_discretizzato.B, zeros(6,1), zeros(3,1), Hx, hx, Hu, hu, Q, R);

CIS_G = Polyhedron(G, g);
CIS_G = minHRep(CIS_G);
% disp("Il Control invariant set è un insieme vuoto? " + boolean(CIS_G.isEmptySet));

CIS_G_T = projection(CIS_G , 1:3);
CIS_G_Q = projection(CIS_G , 4:6);

figure
CIS_G_T.plot();
title("Proiezione del CIS delle temperature nelle stanze")
limitiTemp = [X_vinc_lin(7) X_vinc_lin(1)];
xlim(limitiTemp)
ylim(limitiTemp)
zlim(limitiTemp)

figure;
CIS_G_Q.plot();
title("Proiezione del CIS della potenza termica dei termosifoni")
limitiQ = [X_vinc_lin(10) X_vinc_lin(4)];
xlim(limitiQ)
ylim(limitiQ)
zlim(limitiQ)

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

%% simulazione del sistema discretizzato

[A_cal , A_cal_n , B_cal , B_cal_n,  Q_cal , R_cal] = Calligrafica(sys_discretizzato.A , sys_discretizzato.B , Q , R , Q , Np);

n_sim = 100;
x0_new = x0_centrato;

%calcoliamo H
H = 2 * (B_cal' * Q_cal * B_cal + R_cal);

% impostiamo i vincoli
A_qp = [B_cal; %vincolo di massimo dello stato
        -B_cal; %vincolo di minimo dello stato
        eye(width(B_cal)); % vincolo di massimo dell'ingresso
        -eye(width(B_cal)); % vincolo di minimo dell'ingresso 
        G * B_cal_n]; % Spero che sia questo il vincolo

% definizione vincoli du ingresso e stati centrati
X_max = [];
X_min = [];
U_max = [];
U_min = [];

for i = 1:Np
   X_max = [X_max; X_vinc_lin(1:6)];
   X_min = [X_min; X_vinc_lin(7:end)];
   U_max = [U_max; U_vinc_lin(1:3)];
   U_min = [U_min; U_vinc_lin(4:end)];
end

storia_x = [];
storia_u = [];

for i = 1:n_sim
    % salvo la precedente configurazione degli stati
    storia_x(1:6 , i) = x0_new;

    %calcoliamo f e b_qp
    f = 2* x0_new' * A_cal' * Q_cal * B_cal;
    b_qp = [X_max - A_cal * x0_new;
            -X_min + A_cal * x0_new;
            U_max;
            -U_min;
            g - G * A_cal_n * x0_new]; % vincolo terminale
    
    % % plot dei vincoli
    % figure
    % Vinc_U = Polyhedron('A' , A_qp , 'b' , b_qp);
    % Vinc_U_primo = Vinc_U.projection(1:3);
    % Vinc_U_primo.plot();

    % troviamo il minimo
    [u , ~ , flag] = quadprog(H , f , A_qp , b_qp);
    u_0 = u(1:3);
    
    storia_u(1:3 , i) = u_0; 

    % applichiamo il primo passo
    x0_new = sys_discretizzato.A * x0_new + sys_discretizzato.B * u_0; % applichiamo solo il primo passo

end

%% plot della simulazione

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

%% simulazione a tempo continuo con il controllo

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

    controlAction = MPC(x_run, sys_discretizzato, Q, R, Np, G,g, X_vinc_lin, U_vinc_lin);
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
