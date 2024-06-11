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
u_ref = [100 100  100 ]';


%% ODE del sistema
dxdt = @(t,x) tempCasa(t, x, k, C, tau, T_ext, k_ext, u_ref);

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
B_lin = double(B(u_ref(1) , u_ref(2) , u_ref(3)));

C_lin = eye(6);
D_lin = zeros(6,3);

sys_lineare = ss(A_lin, B_lin, C_lin, D_lin);
% 
% % Verifica della Stabilità del sistema lineare
disp("Autovalori di della matrice A linearizzata:")
disp(eig(A_lin));

% disp("Stati con ingresso costante a: " + u(1));
% disp(-sys_lineare.A^-1 * sys_lineare.B * zeros(3,1));

%% Simulazione Sistema
% 
% tt_sim = 0:1:4999;
% uu_sim = ones(3,5000) * 100;
% xx_0sim = [284 285 284 0 10 0]';
% 
% 
% lsim(sys_lineare , uu_sim , tt_sim , xx_0sim)

%% Discretizziamo
Ts = 60;
sys_discretizzato = c2d(sys_lineare, Ts);
figure
pzmap(sys_discretizzato)

% Verifica della Stabilità del sistema discretizzato
disp("Modulo degli autovalori di della matrice A discretizzata:")
disp(abs(eig(sys_discretizzato.A)));

%punto di equilibrio
% x_eq = (eye(width(sys_discretizzato.A)) - sys_discretizzato.A)^-1 * sys_discretizzato.B * [100 100 100]'; 
% 
% disp("Punto di equilibiro con ingresso pari a 100:")
% disp(x_eq);

%% simulazione sistema discretizzato
% 
% tt_sim = 0:Ts:4999;
% uu_sim = ones(3,length(tt_sim)) * 100;
% xx_0sim = [284 285 284 0 10 0]';
% 
% lsim(sys_discretizzato , uu_sim , tt_sim , xx_0sim)

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

T_vinc = [300 , 282.5];
Q_vinc = [150 , 0];
U_vinc = [150 , 0];

X_vinc = [T_vinc(1)*ones(3,1);
    Q_vinc(1)*ones(3,1);
    T_vinc(2)*ones(3,1);
    Q_vinc(2)*ones(3,1)]; %creazionie vincoli di massimo (prime 6 righe) e di minimo (restanti 6 righe)

U_vinc = [U_vinc(1) * ones(3,1);
    U_vinc(2)*ones(3,1)];

X_vinc_lin = X_vinc - [x_ref ; x_ref]; %calcoliamo i vincoli centrati nel punto di equilibrio
U_vinc_lin = U_vinc - [u_ref ; u_ref];



Hx = [eye(6); -eye(6)];
hx = [ones(6,1) ; -ones(6,1)] .* X_vinc_lin;
Hu = [eye(3); -eye(3)];
hu = [ones(3,1) ; -ones(3,1)] .* U_vinc_lin;

%% definizione delle matrici del costo quadratico
Q = 1.e2*eye(6);
R = 1e1;

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

%% Verifica della fattibilità del n-step controllable invariant set

Np = 5;
[Np_steps_H, Np_steps_h] = controllable_set(Hx, hx, Hu, hu, G, g, sys_discretizzato.A, sys_discretizzato.B, Np);

%%
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

figure;
Np_steps_Q.plot();
title("Proiezione del dominio di attrazione della potenza termica dei termosifoni")
limitiQ = [X_vinc_lin(10) X_vinc_lin(4)];
xlim(limitiQ)
ylim(limitiQ)
zlim(limitiQ)


%% simulazione del sistema

[A_cal , B_cal,  Q_cal , R_cal] = Calligrafica(sys_discretizzato.A , sys_discretizzato.B , Q , R , Q , Np);

n_sim = 83;
x0_new =[284 285 284 0 10 0]';

%calcoliamo H
H = 2 * (B_cal' * Q_cal * B_cal + R_cal);

% impostiamo i vincoli
A_qp = [B_cal; %vincolo di massimo dello stato
        -B_cal; %vincolo di minimo dello stato
        eye(width(B_cal)); % vincolo di massimo dell'ingresso
        -eye(width(B_cal))]; % vincolo di minimo dell'ingresso %vincolo terminale dello stato

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
    
    storia_x(1:6 , i) = x0_new;
    %calcoliamo f e b_qp
    f = 2* x0_new' * A_cal' * Q_cal * B_cal;
    b_qp = [X_max - A_cal * x0_new;
            -X_min + A_cal * x0_new;
            U_max;
            -U_min];%manca il vincolo terminale
    
    % %plot dei vincoli
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




plot(1:n_sim , storia_x(1:3 , :))
title("Evoluzione della temperatura")

figure
plot(1:n_sim , storia_x(4:end , :))
title("Evoluzione della potenza dei termosifoni")

figure
plot(1:n_sim , storia_u)
title("Evoluzione degli ingressi")


