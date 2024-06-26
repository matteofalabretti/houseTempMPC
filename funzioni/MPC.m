function [controlAction , flag] = MPC(x_attuale, A_cal , A_cal_n , B_cal , B_cal_n,  Q_cal , R_cal, Window, G, g, Vinc_X, Vinc_U)


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

for i = 1:Window+1
   X_max = [X_max; Vinc_X(1:6)];
   X_min = [X_min; Vinc_X(7:end)];
   if i ~= Window+1
       U_max = [U_max; Vinc_U(1:3)];
       U_min = [U_min; Vinc_U(4:end)];
   end
end

%calcoliamo f e b_qp
f = 2* x_attuale' * A_cal' * Q_cal * B_cal;
b_qp = [X_max - A_cal * x_attuale;
        -X_min + A_cal * x_attuale;
        U_max;
        -U_min;
        g - G * A_cal_n * x_attuale];

% troviamo il minimo
[controlAction , ~ , flag] = quadprog(H , f , A_qp , b_qp);

end