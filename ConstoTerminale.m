clc

%% LQR
[~ , P] = dlqr(sys_discretizzato.A , sys_discretizzato.B , Q , R);


%%
[A_cal , A_cal_n , B_cal , B_cal_n,  Q_cal , R_cal] = Calligrafica(sys_discretizzato.A , sys_discretizzato.B , Q , R , P , Np);

syms x u

func = (A_cal * x + B_cal* u)' * Q_cal * (A_cal * x + B_cal* u) + u' *R_cal * u + (A_cal_n * x + B_cal_n * u);