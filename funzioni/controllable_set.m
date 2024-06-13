function [H_nsteps, h_nsteps , Np] = controllable_set(Hx, hx, Hu, hu, H_target, h_target, A, B, point)

    %controllable_set Calcolo del n-step controllable_set di un sistema lineare
%   Input:
%       - Hx*x <= hx: vincoli sullo stato
%       - Hu*u <= hu: vincoli sull'ingresso
%       - H_target <= h_target: la regione obbiettivo
%       - A, B: matrici del sistema
%       - N: il numero di passi di previsione



    n = size(A,2);
    
    % candidato iniziale
    H_ii_steps = H_target;
    h_ii_steps = h_target;
    tic;
    for ii = 1:1000 

        % Calcoliamo il set ad un passo rispetto a quello precedente
        temp = Polyhedron('A', [H_ii_steps*A, H_ii_steps*B; zeros(size(Hu, 1), n), Hu], 'b', [h_ii_steps;hu]);
    
        % Proiezioni in R^n
        % temp = temp.minHRep();
        temp = projection(temp, 1:n);
        

        % Intersezione con X := {x | Hx*x <= hx}
        H_ii_steps = [temp.A;Hx];
        h_ii_steps = [temp.b;hx];

        
        disp("Fine Calcolo Passo " + ii )

        temp = Polyhedron( 'A' , H_ii_steps , 'b' ,h_ii_steps);

        if temp.contains(point)
            break
        end
    
    end
    disp("Tempo impiegato: " + toc);

    temp = temp.minHRep();
    Np = ii;
    H_nsteps = temp.A;
    h_nsteps = temp.b;
    
    
end