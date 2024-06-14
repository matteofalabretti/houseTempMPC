% Questo script esegue tutte le operzioni preliminari per la progettazione
% del MPC
% In particolare contiene:
% -Definizione delle costanti
% -Simulazione del sistema
% -Linearizzazione
% -Discretizzazione
% -Definizione dei vincoli di X e U

%%

set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
% set(0,'DefaultFigureWindowStyle', 'normal')
set(0,'defaulttextInterpreter','latex')
rng('default');

addpath('funzioni\')
addpath('funzioniAggiuntive\')

%% Definizioni costanti
T_ext = 278;
k_ext = 9;

C = [6300;4600;4200];

k = [16;18;19];
tau = [580;520;540];

% Definizione Obbiettivi di Controllo
x_ref = [289 289 289 100 100 100]';
u_ref = [100 100  100 ]';
x_start = [284 285 284 0 10 0]';


%% ODE del sistema
dxdt = @(t,x) tempCasa(t, x, k, C, tau, T_ext, k_ext, u_ref);

simulazione = 5000;

% simulazione del comportamento del sistema
[tt, xx] = ode45(dxdt, linspace(0, simulazione, simulazione+1), x_start);
figure
hold on

sgtitle("Simulazione del sistema con ingresso constante pari a 100W") 
subplot(2,1,1)
plot(tt, xx(: , 1:3));
title("Temperature delle stanze")
ylabel("Temperatura $[ ^{\circ}C]$" , Interpreter="latex")
xlabel("Tempo $[s]$",  Interpreter="latex")
legend(["T1" "T2" "T3"])

subplot(2, 1,2)
plot(tt , xx(: , 4:6));
title("Potenza termica dei termosifoni")
ylabel("Potenza $[W]$" , Interpreter="latex")
xlabel("Tempo $[s]$",  Interpreter="latex")
legend(["Q1" "Q2" "Q3"])
hold off

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

x0_centrato = x_start - x_ref;

% Verifica della Stabilità del sistema lineare
disp("Autovalori di della matrice A linearizzata:")
disp(eig(A_lin));

%% Discretizziamo
sys_discretizzato = c2d(sys_lineare, Ts);

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

T_vinc = [300 , 282.5];
Q_vinc = [150 , 0];
U_vinc = [150 , 0];

X_vinc = [T_vinc(1)*ones(3,1);
    Q_vinc(1)*ones(3,1);
    T_vinc(2)*ones(3,1);
    Q_vinc(2)*ones(3,1)]; %creazione vincoli di massimo (prime 6 righe) e di minimo (restanti 6 righe)

U_vinc = [U_vinc(1) * ones(3,1);
    U_vinc(2)*ones(3,1)];

X_vinc_lin = X_vinc - [x_ref ; x_ref]; %calcoliamo i vincoli centrati nel punto di equilibrio
U_vinc_lin = U_vinc - [u_ref ; u_ref];



Hx = [eye(6); -eye(6)];
hx = [ones(6,1) ; -ones(6,1)] .* X_vinc_lin;
Hu = [eye(3); -eye(3)];
hu = [ones(3,1) ; -ones(3,1)] .* U_vinc_lin;
