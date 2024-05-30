function [G, g] = CIS(A, B, x_ref, u_ref, Hx, hx, Hu, hu, Q, R)

%CIS Calcolo del control invariant set (CIS) di un sistema lineare
%   Questo metodo assume che un controllore LQR venga utilizzato
%   all'interno del CIS
%   Input:
%       - A, B: matrici del sistema
%       - x_ref u_ref: equilibrio attorno al quale costruire il CIS
%       - Hx*x <= hx: vincoli sullo stato
%       - Hu*u <= hu: vincoli sull'ingresso
%       - Q,R: matrici per LQR

%   Controllore LQR
K = -dlqr(A, B, Q, R);

%   Matrice A del sistema controllato con LQR
A_lqr = A + B * K;

%   Vincoli sullo stato e sull'ingresso (H*x <= h)
H = [Hx;Hu * K];
h = [hx;hu + Hu *(K*x_ref - u_ref)];

%   Calcolo del CIS (G*x<=g)
CIS_poly_curr = Polyhedron(H, h);
continua = 1;
i = 0;

tic

while continua
    i = i+1;

    %   Memorizza vecchio candidato
    CIS_poly_prev = CIS_poly_curr;
    
    %   Calcola nuovo candidato (G_hat * x <= g_hat)
    G_hat = [CIS_poly_curr.A * A_lqr;H];
    g_hat = [CIS_poly_curr.b + CIS_poly_curr.A*B*(K*x_ref - u_ref);h];

    CIS_poly_curr= Polyhedron(G_hat, g_hat);
    CIS_poly_curr =  minHRep(CIS_poly_curr); %riduce alla rappresentazione minima del
    %Poliedro
    
    %plot di verifica
    if mod(i , 100) == 0
        disp("Iterazione: " + i)
    end

    continua = CIS_poly_prev ~= CIS_poly_curr; 
end

disp("Numero di iterazioni: " + i)
disp("Tempo di esecuzione: " + toc/1000+  "s")

%   Disequazioni che descrivono il CIS
G = CIS_poly_curr.A;
g = CIS_poly_curr.b;

end