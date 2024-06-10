function [H_nsteps, h_nsteps] = controllable_set(Hx, hx, Hu, hu, H_target, h_target, A, B, N)

    %controllable_set Calcolo del n-step controllable_set di un sistema lineare
%   Input:
%       - Hx*x <= hx: vincoli sullo stato
%       - Hu*u <= hu: vincoli sull'ingresso
%       - H_target <= h_target: la regione obbiettivo
%       - A, B: matrici del sistema
%       - N: il numero di passi di previsione



    n = size(A,2);
    m = size(B, 2);
    
    % candidato iniziale
    H_ii_steps = H_target;
    h_ii_steps = h_target;
    tic;
    for ii = 1:N
        % Calcoliamo il set ad un passo rispetto a quello precedente
        temp = Polyhedron('A', [H_ii_steps*A, H_ii_steps*B; zeros(size(Hu, 1), n), Hu], 'b', [h_ii_steps;hu]);
    
        % Proiezioni in R^n
        temp = projection(temp, 1:n);
        temp = temp.minHRep();
    
        % Intersezione con X := {x | Hx*x <= hx}
        H_ii_steps = [temp.A;Hx];
        h_ii_steps = [temp.b;hx];
    
    end
    disp("Tempo impiegato: " + toc);

    H_nsteps = H_ii_steps;
    h_nsteps = h_ii_steps;
    
end