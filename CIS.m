function [G, g] = CIS(A, B, x_ref, u_ref, Hx, hx, Hu, hu, Q, R)
max_iteration = 10000;
%CIS Calcolo del control invariant set (CIS) di un sistema lineare
%   Questo metodo assume che un controllore LQR venga utilizzato
%   all'interno del CIS
%   Input:
%       - A, B: matrici del sistema
%       - x_ref: equilibrio attorno al quale costruire il CIS
%       - Fx*x<=fx: vincoli sullo stato
%       - Fu*u<=fu: vincoli sull'ingresso
%       - Q,R: matrici per LQR

%   Controllore LQR
K = -dlqr(A, B, Q, R);

%   Matrice A del sistema controllato con LQR
A_lrq = A + B*K;

%   Vincoli sullo stato e sull'ingresso (F*x <= f)
H = [Hx;Hu*K];
h = [hx;hu + Hu*(K*x_ref - u_ref)];

%   Calcolo del CIS (G*x<=g)
CIS_poly_prev = Polyhedron();
CIS_poly_curr = Polyhedron(H, h);
i = 0;

while CIS_poly_prev.isEmptySet || CIS_poly_prev ~= CIS_poly_curr || i > max_iteration
    i = i+1;
    
    %   Memorizza vecchio candidato
    CIS_poly_prev = CIS_poly_curr;
    
    %   Calcola nuovo candidato (G_hat * x <= g_hat)
    G_hat = [CIS_poly_curr.A * A_lrq;H];
    g_hat = [CIS_poly_curr.B + CIS_poly_curr.A*B*(K*x_ref - u_ref);h];
    CIS_poly_curr= Polyhedron(G_hat, g_hat);
    %minHRep(CIS_poly_curr); riduce alla rappresentazioneminima del
    %Poliedro
    
end

%   Disequazioni che descrivono il CIS
G = CIS_poly_curr.A;
g = CIS_poly_curr.B;

end